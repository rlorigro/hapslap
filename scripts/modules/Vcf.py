from collections import defaultdict
from copy import copy,deepcopy
import subprocess
import os.path
import sys

import vcf

if __name__ != "__main__":
    from .Paths import enumerate_paths
else:
    from Paths import enumerate_paths

from networkx import DiGraph
from pysam import FastaFile
from vcf import VCFReader


def decompress_vcf(vcf_path, timeout=60*60, use_cache=False):
    args = ["bgzip", "-d", vcf_path]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    output_vcf_path = vcf_path.replace(".gz","")
    return output_vcf_path


def compress_and_index_vcf(vcf_path, timeout=60*60, use_cache=False):
    output_vcf_path = vcf_path + ".gz"
    output_tbi_path = output_vcf_path + ".tbi"

    # Enable caching by path name
    if use_cache and os.path.exists(output_vcf_path) and os.path.exists(output_tbi_path):
        return output_vcf_path

    args = ["bgzip", "-f", vcf_path]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    args = ["tabix", "-p", "vcf", output_vcf_path]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return output_vcf_path


def merge_vcfs(vcf_paths, output_path, force_samples=False, timeout=60*60):

    if not force_samples:
        args = ["bcftools", "merge", "-m", "none", "-0"] + vcf_paths
    else:
        args = ["bcftools", "merge", "-m", "none", "-0", "--force-samples"] + vcf_paths


    sys.stderr.write(" ".join(args)+'\n')

    with open(output_path, 'w') as file:
        try:
            p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, stdout=file, timeout=timeout)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return None

        except subprocess.TimeoutExpired as e:
            sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return None

    return output_path


# bcftools view [OPTIONS] file.vcf.gz [REGION […]]
def bcftools_view(vcf_path, region_string, output_path, timeout=60*60):
    args = ["bcftools", "view", vcf_path, region_string]

    sys.stderr.write(" ".join(args)+'\n')

    with open(output_path, 'w') as file:
        try:
            p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, stdout=file, timeout=timeout)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return None

        except subprocess.TimeoutExpired as e:
            sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return None

    return output_path


def merge_vcfs_in_directory(vcf_directory, output_path, compress_and_index=True):
    all_vcf_paths = list()

    for filename in os.listdir(vcf_directory):
        if filename.endswith(".vcf.gz"):
            all_vcf_paths.append(os.path.join(vcf_directory, filename))

    output_path = merge_vcfs(vcf_paths=all_vcf_paths, output_path=output_path)

    if compress_and_index:
        compress_and_index_vcf(vcf_path=output_path)

    return output_path


class Allele:
    def __init__(self, start, stop, sequence):
        self.start = int(start)
        self.stop = int(stop)
        self.sequence = sequence
        self.samples = set()
        self.is_left_flank = False
        self.is_right_flank = False
        self.vcf_data = None

    def add_sample(self, s):
        self.samples.add(s)

    def construct_vcf_data(self, vcf_record, alt_index):
        self.vcf_data = {
            "CHROM":str(vcf_record.CHROM),
            "POS":str(vcf_record.POS),
            "ID":str(vcf_record.ID),
            "REF":str(vcf_record.REF),
            "ALT":str(vcf_record.ALT[alt_index-1]),
            "QUAL":str(vcf_record.QUAL),
            "FILTER":"PASS",
            "INFO":'.',
            "FORMAT":"GT",
            "GT":"0/1"
        }

    def get_vcf_data(self):
        return (self.vcf_data[x] for x in ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"])

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
def get_alleles_from_vcfs(ref_path, data_per_sample, chromosome, ref_start=None, ref_stop=None, debug=False, skip_incompatible=False):
    region_string = chromosome

    if ref_start is not None and ref_stop is not None:
        region_string += ":" + str(ref_start) + "-" + str(ref_stop)

    ref_sequence = FastaFile(ref_path).fetch(chromosome)

    alleles = dict()

    for sample_name,data in data_per_sample.items():
        vcf_path = data["vcf"]

        with open(vcf_path, 'rb') as file:
            if debug:
                print(vcf_path)

            vcf = VCFReader(file)

            try:
                records = vcf.fetch(region_string)
            except ValueError:
                print("WARNING: Skipping sample because has no VCF entries: " + sample_name)
                continue

            for record in records:
                if debug:
                    print()
                    print("var_subtype:\t%s" % record.var_subtype)
                    print("start:\t\t%d" % record.start)

                for call in record.samples:
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
                                if start < ref_start - 1 or stop > ref_stop:

                                    if not skip_incompatible:
                                        sys.stderr.write(vcf_path)
                                        sys.stderr.write('\n')
                                        raise Exception("ERROR: VCF allele in sample %s with coords %d-%d extends out of region %s:%d-%d" % (sample_name, start,stop,chromosome,ref_start,ref_stop))
                                    else:
                                        sys.stderr.write(vcf_path)
                                        sys.stderr.write('\n')
                                        sys.stderr.write("WARNING: skipping VCF allele in sample %s with coords %d-%d extending out of region %s:%d-%d\n" % (sample_name, start,stop,chromosome,ref_start,ref_stop))

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
                            a.construct_vcf_data(record,allele_index)
                            h = a.hash()

                            if h not in alleles:
                                alleles[h] = a

                            alleles[h].add_sample(sample_name)

    # Throw away hashes and keep unique alleles as a list
    return list(alleles.values())


def vcf_to_graph(ref_path, data_per_sample, chromosome, ref_start, ref_stop, ref_sample_name, flank_length, skip_incompatible=False):
    ref_sequence = FastaFile(ref_path).fetch(chromosome)

    alleles = get_alleles_from_vcfs(
        ref_path=ref_path,
        data_per_sample=data_per_sample,
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        skip_incompatible=skip_incompatible
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


def remove_empty_nodes_from_variant_graph(graph, alleles):
    empty_nodes = list()
    for n in graph.nodes:
        if len(alleles[n].sequence) == 0:
            empty_nodes.append(n)

    edge_to_deletion_index = dict()
    for n in empty_nodes:
        a_nodes = [e[0] for e in graph.in_edges(n)]
        b_nodes = [e[1] for e in graph.out_edges(n)]

        for a in a_nodes:
            for b in b_nodes:
                graph.add_edge(a,b, weight=0)
                edge_to_deletion_index[(a,b)] = n

        graph.remove_node(n)

    return graph, edge_to_deletion_index


# def get_data(record, n, alleles):
#
#     # CHROM  POS    ID          REF  ALT  QUAL  FILTER  INFO                     FORMAT       NA00001
#     # 20     14370  rs6054257   G    A    29    PASS    NS=3;DP=14;AF=0.5;DB;H2  GT:GQ:DP:HQ  0|0:48:1:51,51
#     alt_index = alleles[n].alt_index
#
#     # Summarize the data and hash by pos + data
#     data = (
#         str(record.CHROM),
#         str(record.POS),
#         str(record.ID),
#         str(record.REF),
#         str(record.ALT[alt_index]),
#         str(record.QUAL),
#         "PASS",
#         '',
#         "GT",
#         "0/1")
#
#     return data


def write_paths_to_vcf(alleles:list, paths:list, output_path:str, sample_name, edge_to_deletion_index=None, compress_and_index=True):
    records_per_position = defaultdict(set)

    with open(output_path, 'w') as file:
        file.write("##fileformat=VCFv4.2\n")
        file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"The genotype of the variant\">\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)

        for path in paths:
            for i,n in enumerate(path):
                # Deletions are implicit in the graph (empty nodes are removed) and therefore need to be searched by
                # pairs of nodes (edges) instead of by single nodes
                if edge_to_deletion_index is not None and i > 0:
                    edge = (path[i-1],n)

                    if edge in edge_to_deletion_index:
                        n_del = edge_to_deletion_index[edge]

                        # If the edge is found in the mapping, it is guaranteed to be in the alleles list
                        pos = int(alleles[n_del].vcf_data["POS"])
                        data = tuple(alleles[n_del].get_vcf_data())
                        records_per_position[pos].add(data)

                # Skip ref nodes which have no record
                if alleles[n].vcf_data is None:
                    continue

                pos = int(alleles[n].vcf_data["POS"])
                data = tuple(alleles[n].get_vcf_data())
                records_per_position[pos].add(data)

        for pos,items in sorted(records_per_position.items(), key=lambda x: x[0]):
            for data in items:
                file.write('\t'.join(data))
                file.write('\n')

    if compress_and_index:
        compress_and_index_vcf(output_path)


def test_vcf():
    project_directory = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_directory = os.path.join(project_directory,"data/test/")

    sample_1_path = os.path.join(data_directory,"n2/HG00438_chr20_sniffles.vcf.gz")
    sample_2_path = os.path.join(data_directory,"n2/HG00621_chr20_sniffles.vcf.gz")
    ref_path = os.path.join(data_directory,"chm13v2.0_chr20.fa")

    data_per_sample = defaultdict(dict)
    data_per_sample["HG00438"]["vcf"] = sample_1_path
    data_per_sample["HG00621"]["vcf"] = sample_2_path

    chromosome = "chr20"
    ref_start = 47474225
    ref_stop = 47476612
    ref_sample_name = "ref"
    flank_length = 200

    # First read the VCF and convert to a graph
    graph, alleles = vcf_to_graph(
        ref_path,
        data_per_sample,
        chromosome,
        ref_start,
        ref_stop,
        ref_sample_name,
        flank_length)

    sample_name = "test"
    output_path = os.path.join(data_directory, "test_output.vcf")
    paths = enumerate_paths(alleles=alleles, graph=graph)

    # Write all alleles to a VCF
    write_paths_to_vcf(
        alleles,
        paths,
        output_path,
        sample_name,
    )

    data_per_sample_2 = defaultdict(dict)
    data_per_sample_2["HG00438"]["vcf"] = output_path + ".gz"

    # Load the newly generated VCF
    graph_2, alleles_2 = vcf_to_graph(
        ref_path,
        data_per_sample_2,
        chromosome,
        ref_start,
        ref_stop,
        ref_sample_name,
        flank_length
    )

    # Delete empty nodes this time
    graph_2, edge_to_deletion_index = remove_empty_nodes_from_variant_graph(graph, alleles)

    for item in edge_to_deletion_index.items():
        print(item)

    sample_name_2 = "test"
    output_path_2 = os.path.join(data_directory, "test_output_2.vcf")
    paths_2 = enumerate_paths(alleles=alleles_2, graph=graph_2)

    # Write the alleles to a VCF again, this time using the edge deletion index for missing nodes
    write_paths_to_vcf(
        alleles_2,
        paths_2,
        output_path_2,
        sample_name_2,
        edge_to_deletion_index=edge_to_deletion_index,
        # edge_to_deletion_index=None,
    )

    decompress_vcf(output_path + ".gz")
    decompress_vcf(output_path_2 + ".gz")

    with open(output_path, 'r') as file_a, open(output_path_2, 'r') as file_b:
        lines_a = [l for l in file_a]
        lines_b = [l for l in file_b]

        for i in range(max(len(lines_a), len(lines_b))):
            if i > len(lines_a) or i > len(lines_b):
                raise Exception("FAIL: lines in a and b VCFs not equal: \n\t" + lines_a[i] + '\n\t' + lines_b[i])

            if lines_a[i] != lines_b[i]:
                raise Exception("FAIL: lines in a and b VCFs not equal: \n\t" + lines_a[i] + '\n\t' + lines_b[i])


if __name__ == "__main__":
    test_vcf()
