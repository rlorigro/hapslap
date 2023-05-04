from collections import defaultdict
from matplotlib import pyplot
from networkx import DiGraph
import networkx
import os.path

from vcf import VCFReader
from pysam import FastaFile
import argparse


class Allele:
    def __init__(self, start, stop, sequence):
        self.start = int(start)
        self.stop = int(stop)
        self.sequence = sequence
        self.samples = set()

    def add_sample(self, s):
        self.samples.add(s)

    def __str__(self):
        return "[%s,%s]\t%s" % (self.start, self.stop, self.sequence)

    def __hash__(self):
        return hash((self.start, self.stop, self.sequence))

    def __eq__(self, other):
        return (self.start, self.stop, self.sequence) == (other.start, other.stop, other.sequence)

    def hash(self):
        return self.__hash__()


def path_recursion(graph, alleles, id, path_sequence=""):
    path_sequence += alleles[id].sequence
    print(id, len(path_sequence))

    out_edges = graph.out_edges(id)
    print(out_edges)

    if len(out_edges) == 0:
        yield path_sequence
        pass
    else:
        for edge in out_edges:
            print(edge[1])
            yield from path_recursion(graph=graph, alleles=alleles, id=edge[1], path_sequence=path_sequence)


def enumerate_paths(alleles, graph, output_directory):
    start_id = next(networkx.topological_sort(graph))

    print("Starting path recursion from %d" % start_id)

    paths = [p for p in path_recursion(graph=graph, alleles=alleles, id=start_id)]

    output_path = os.path.join(output_directory, "paths.fasta")
    with open(output_path, 'w') as file:
        for p,path in enumerate(paths):
            file.write(">%d\n" % p)
            file.write(path)
            file.write('\n')


def main():
    output_directory = "/home/ryan/code/hapslap/data/test/region_1/"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    ref_path = "/home/ryan/data/human/reference/chm13v2.0.fa"

    vcf_paths = [
        # "/home/ryan/code/hapslap/data/test/hprc/HG002_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG00438_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG005_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG00621_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00673_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00733_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00735_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00741_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01071_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01106_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01109_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01123_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01175_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01243_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01258_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01358_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01361_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01891_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01928_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01952_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01978_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02055_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02080_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02109_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02145_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02148_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02257_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02486_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02559_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02572_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02622_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02630_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02717_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02723_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02818_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02886_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03098_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03453_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03486_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03492_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03516_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03540_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03579_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA18906_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA19240_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA20129_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA21309_chr20_sniffles.vcf.gz"
    ]

    chromosome = "chr20"
    region_start = 47474020
    region_stop = 47477018

    region_string = chromosome + ":" + str(region_start) + "-" + str(region_stop)

    gfa_path = os.path.join(output_directory, region_string + ".gfa")
    csv_path = os.path.join(output_directory, region_string + ".csv")

    ref_sequence = FastaFile(ref_path).fetch("chr20")

    alleles = dict()

    for vcf_path in vcf_paths:
        sample_name = os.path.basename(vcf_path).split("_")[0]

        with open(vcf_path, 'rb') as file:
            vcf = VCFReader(file)

            records = vcf.fetch(region_string)

            for record in records:
                print()
                print("var_subtype:\t%s" % record.var_subtype)
                print("start:\t\t%d" % record.start)

                # One VCF per sample, no iterating of samples needed
                call = record.samples[0]
                gt = [int(call.data.GT[0]), int(call.data.GT[-1])]
                print("gt:\t\t%d/%d" % (gt[0], gt[1]))
                print("ref_length:\t%d" % len(record.alleles[0]))
                print("a_length:\t%d" % len(record.alleles[gt[0]]))
                print("b_length:\t%d" % len(record.alleles[gt[1]]))
                print("is_sv_precise:\t%d" % record.is_sv_precise)

                # Iterate unique, non-ref alleles only
                for allele_index in set(gt):
                    print(record.alleles[allele_index])

                    if allele_index != 0:
                        l = len(record.alleles[0]) if record.alleles[0] != 'N' else 0
                        not_insert = (l >= len(record.alleles[allele_index]))

                        # We are sorting by ref coordinates, so we use the ref allele start/stop to keep track of
                        # where the alt allele will be substituted
                        start = int(record.start)
                        stop = start + l + int(not_insert)
                        sequence = str(record.alleles[allele_index])
                        sequence = sequence if sequence != 'N' else ''

                        # Collapse identical alleles by hashing them as a fn of start,stop,sequence
                        # But keep track of which samples are collapsed together
                        a = Allele(start, stop, sequence)
                        h = a.hash()

                        if h not in alleles:
                            alleles[h] = a

                        alleles[h].add_sample(sample_name)

    print()

    # Throw away hashes and keep unique alleles as a list
    alleles = list(alleles.values())

    # sample_to_allele_index = dict()
    #
    # # Build reverse mapping of sample -> allele_index
    # for a,[allele,samples] in enumerate(alleles):
    #     for sample in samples:
    #         sample_to_allele_index[sample] = a
    #
    # # Generate pointers to the allele data and sort them by allele start coord
    # allele_indexes = list(range(len(alleles)))
    # allele_indexes = sorted(allele_indexes, key=lambda x: alleles[x][0].start)
    #
    # overlapping_allele_indexes = list()
    #
    # # Do a sweep over the sorted allele intervals to find overlapping ones
    # stop = -1
    # for i in allele_indexes:
    #     if alleles[i][0].start <= stop:
    #         overlapping_allele_indexes[-1].add(i)
    #
    #         if alleles[i][0].stop > stop:
    #             stop = alleles[i][0].stop
    #     else:
    #         overlapping_allele_indexes.append({i})
    #         stop = alleles[i][0].stop
    #
    # # Print results
    # for i,item in enumerate(overlapping_allele_indexes):
    #     print(i)
    #     for allele_index in item:
    #         allele = alleles[allele_index][0]
    #         sample = alleles[allele_index][1]
    #         print(allele, sample)

    # Construct a list of coordinates along the reference path which contain edges to VCF alleles
    ref_edges = defaultdict(lambda: [[],[]])

    for a,allele in enumerate(alleles):
        print(a)
        print(allele)
        ref_edges[allele.start][1].append(a)
        ref_edges[allele.stop][0].append(a)

    # Append dummy item at end of list to make one-pass iteration easier
    ref_edges[region_stop] = [[],[]]

    # Sort the list by coord so that it can be iterated from left to right
    ref_edges = list(sorted(ref_edges.items(), key=lambda x: x[0]))

    graph = DiGraph()

    ref_id_offset = len(alleles)

    # TODO: replace this with the actual reference name?
    ref_sample_name = "ref"

    # -- Construct graph and ref alleles --

    # Initialize vars that will be iteration dependent
    prev_coord = region_start
    in_edges = []

    # First generate nodes for all the known VCF alleles
    for allele_index,allele in enumerate(alleles):
        id = allele_index
        graph.add_node(id)

    # Construct ref backbone nodes with sufficient breakpoints to capture all in/out allele edges
    r = ref_id_offset
    for i,[coord,edges] in enumerate(ref_edges):
        id = r

        print(id, coord, edges)

        sequence = ref_sequence[prev_coord:coord]

        # Create Allele object for this ref node and mimic the allele data structure
        a = Allele(start=prev_coord, stop=coord, sequence=sequence)
        a.add_sample(ref_sample_name)
        alleles.append(a)

        # Create the node in the graph data structure
        graph.add_node(id)

        prev_coord = coord
        r += 1

    # Construct edges from reference backbone to existing alleles and other backbone nodes
    r = ref_id_offset
    for i,[coord,edges] in enumerate(ref_edges):
        id = r

        for allele_index in in_edges:
            other_id = allele_index
            graph.add_edge(other_id,id)

        for allele_index in edges[1]:
            other_id = allele_index
            graph.add_edge(id,other_id)

        # Add edge to next ref sequence
        if i < len(ref_edges) - 1:
            next_id = r+1
            graph.add_edge(id,next_id)

            r += 1

        in_edges = edges[0]

    # Enumerate paths
    enumerate_paths(alleles=alleles, graph=graph, output_directory=output_directory)

    sample_color = "#007cbe"
    ref_color = "#bebebe"

    # Write a csv that colors the nodes in Bandage
    with open(csv_path, 'w') as csv_file:
        csv_file.write("Name,Color\n")

        for allele_index,allele in enumerate(alleles):
            color = sample_color

            if allele_index >= ref_id_offset:
                color = ref_color

            csv_file.write("%s,%s\n" % (allele_index,color))

    # Write the GFA
    with open(gfa_path, 'w') as gfa_file:
        for allele_index,allele in enumerate(alleles):
            gfa_file.write("S\t%s\t%s\n" % (str(allele_index),allele.sequence))

        for e in graph.edges:
            gfa_file.write("L\t%s\t+\t%s\t+\t0M\n" % (str(e[0]), str(e[1])))

    # -- Plot the graph --
    for layer, nodes in enumerate(networkx.topological_generations(graph)):
        # `multipartite_layout` expects the layer as a node attribute, so add the
        # numeric layer value as a node attribute
        for node in nodes:
            graph.nodes[node]["layer"] = layer

    color_map = []
    for a in range(len(alleles)):
        if a >= ref_id_offset:
            color_map.append(ref_color)
        else:
            color_map.append(sample_color)

    pos = networkx.multipartite_layout(graph, subset_key="layer")
    networkx.draw(graph, pos, connectionstyle="arc3,rad=-0.1", node_color=color_map, font_size=16, with_labels=True)
    pyplot.show()
    pyplot.close()

    return


if __name__ == "__main__":
    main()

