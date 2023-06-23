from collections import defaultdict
import os.path
import sys

from networkx import DiGraph
from pysam import FastaFile
from vcf import VCFReader


class Allele:
    def __init__(self, start, stop, sequence):
        self.start = int(start)
        self.stop = int(stop)
        self.sequence = sequence
        self.samples = set()
        self.is_left_flank = False
        self.is_right_flank = False

    def add_sample(self, s):
        self.samples.add(s)

    def __str__(self):
        return "[%s,%s]\t%s" % (self.start, self.stop, self.sequence)

    def as_comma_separated_str(self):
        return ','.join([str(self.start), str(self.stop), str(len(self.sequence)), ' '.join(list(self.samples)), str(self.is_left_flank), str(self.is_right_flank)])

    def __hash__(self):
        return hash((self.start, self.stop, self.sequence))

    def __eq__(self, other):
        return (self.start, self.stop, self.sequence) == (other.start, other.stop, other.sequence)

    def hash(self):
        return self.__hash__()


"""
Given any number of VCFs, build a set of Allele objects which are merged by their hash(sequence,start,stop), retaining
sample names for all merged alleles
"""
def get_alleles_from_vcfs(ref_path, vcf_paths, chromosome, ref_start=None, ref_stop=None, debug=False):
    region_string = chromosome

    if ref_start is not None and ref_stop is not None:
        region_string += ":" + str(ref_start) + "-" + str(ref_stop)

    ref_sequence = FastaFile(ref_path).fetch(chromosome)

    alleles = dict()

    for vcf_path in vcf_paths:
        sample_name = os.path.basename(vcf_path).split("_")[0]

        with open(vcf_path, 'rb') as file:
            if debug:
                print(vcf_path)

            vcf = VCFReader(file)

            records = vcf.fetch(region_string)

            for record in records:
                if debug:
                    print()
                    print("var_subtype:\t%s" % record.var_subtype)
                    print("start:\t\t%d" % record.start)

                # One VCF per sample, no iterating of samples needed
                call = record.samples[0]
                gt = [int(call.data.GT[0]) if call.data.GT[0] != '.' else 0, int(call.data.GT[-1]) if call.data.GT[-1] != '.' else 0]

                if debug:
                    print("gt:\t\t%d/%d" % (gt[0], gt[1]))
                    print("ref_length:\t%d" % len(record.alleles[0]))
                    print("var_type:\t%s" % record.var_type)
                    print("is_sv_precise:\t%d" % record.is_sv_precise)

                b_length = None
                if record.var_subtype == "DUP":
                    b_length = record.INFO["SVLEN"]
                elif record.var_subtype == "INV":
                    b_length = record.INFO["SVLEN"]
                elif record.var_subtype == "complex":
                    sys.stderr.write("WARNING: skipping break end operation: " + chromosome + ":" + str(record.start) + '\n')
                    continue
                else:
                    try:
                        b_length = len(record.alleles[gt[1]])
                    except Exception as e:
                        sys.stderr.write("WARNING: Skipping unparsable variant:\n")
                        sys.stderr.write(vcf_path)
                        sys.stderr.write('\n')
                        sys.stderr.write(str(record))
                        sys.stderr.write('\n')
                        continue

                # print("a_length:\t%d" % len(record.alleles[gt[0]]))
                if debug:
                    print("b_length:\t%d" % b_length)

                # Iterate unique, non-ref alleles only
                for allele_index in set(gt):
                    if debug:
                        print(record.alleles[allele_index])

                    if allele_index != 0:
                        l = len(record.alleles[0]) if record.alleles[0] != 'N' else 0
                        not_insert = (l >= b_length)

                        # We are sorting by ref coordinates, so we use the ref allele start/stop to keep track of
                        # where the alt allele will be substituted
                        start = int(record.start)
                        stop = start + l + int(not_insert)

                        sequence = str(record.alleles[allele_index])

                        # Occasionally the region can be cut wrong, e.g. directly through a deletion variant,
                        # which means that the graph will contain edges in the flanking regions, which will
                        # be deleted at the end to generate haplotypes... causing issues
                        if ref_start is not None and ref_stop is not None:
                            if ref_stop < start < ref_start or ref_stop < stop < ref_start:
                                raise Exception("ERROR: VCF allele with coords %d-%d extends out of region %s:%d-%d" % (start,stop,chromosome,ref_start,ref_stop))

                        if sequence.strip() == "<DUP>":
                            sequence = ref_sequence[start:start+b_length+1]
                        elif sequence.strip() == "<INV>":
                            sequence = ref_sequence[start:start+b_length+1][::-1]
                        elif sequence == 'N':
                            sequence = ''
                        else:
                            pass

                        # Collapse identical alleles by hashing them as a fn of start,stop,sequence
                        # But keep track of which samples are collapsed together
                        a = Allele(start, stop, sequence)
                        h = a.hash()

                        if h not in alleles:
                            alleles[h] = a

                        alleles[h].add_sample(sample_name)

    # Throw away hashes and keep unique alleles as a list
    return list(alleles.values())


def vcf_to_graph(ref_path, vcf_paths, chromosome, ref_start, ref_stop, ref_sample_name, flank_length):
    ref_sequence = FastaFile(ref_path).fetch(chromosome)

    alleles = get_alleles_from_vcfs(
        ref_path=ref_path,
        vcf_paths=vcf_paths,
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop
    )

    ref_start -= flank_length
    ref_stop += flank_length

    # Construct a list of coordinates along the reference path which contain edges to VCF alleles
    ref_edges = defaultdict(lambda: [[],[]])

    for a,allele in enumerate(alleles):
        ref_edges[allele.start][1].append(a)
        ref_edges[allele.stop][0].append(a)

    # Append dummy item at end of list to make one-pass iteration easier
    ref_edges[ref_stop] = [[],[]]

    # Sort the list by coord so that it can be iterated from left to right
    ref_edges = list(sorted(ref_edges.items(), key=lambda x: x[0]))

    graph = DiGraph()

    ref_id_offset = len(alleles)

    # -- Construct graph and ref alleles --

    # Initialize vars that will be iteration dependent
    prev_coord = ref_start
    in_edges = []

    # First generate nodes for all the known VCF alleles
    for allele_index,allele in enumerate(alleles):
        id = allele_index
        graph.add_node(id)

    # Construct ref backbone nodes with sufficient breakpoints to capture all in/out allele edges
    r = ref_id_offset
    for i,[coord,edges] in enumerate(ref_edges):
        id = r

        sequence = ref_sequence[prev_coord:coord]

        # Create Allele object for this ref node and mimic the allele data structure
        a = Allele(start=prev_coord, stop=coord, sequence=sequence)
        a.add_sample(ref_sample_name)
        alleles.append(a)

        # Create the node in the graph data structure
        graph.add_node(id)

        prev_coord = coord
        r += 1

    # Annotate the ref alleles with flanking information
    alleles[ref_id_offset].is_left_flank = True
    alleles[-1].is_right_flank = True

    # Construct edges from reference backbone to existing alleles and other backbone nodes
    r = ref_id_offset
    for i,[coord,edges] in enumerate(ref_edges):
        id = r

        for allele_index in in_edges:
            other_id = allele_index
            graph.add_edge(other_id,id, weight=0)

        for allele_index in edges[1]:
            other_id = allele_index
            graph.add_edge(id,other_id, weight=0)

        # Add edge to next ref sequence
        if i < len(ref_edges) - 1:
            next_id = r+1
            graph.add_edge(id,next_id, weight=0)

            r += 1

        in_edges = edges[0]

    return graph, alleles

