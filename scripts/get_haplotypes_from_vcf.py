from modules.Authenticator import *
from modules.GsUri import *

from collections import defaultdict
from multiprocessing import Pool
from copy import deepcopy,copy
import subprocess
import argparse
import os.path
import re

from matplotlib import pyplot
from networkx import DiGraph
from pysam import FastaFile
from vcf import VCFReader
import networkx


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


def write_graph_to_gfa(output_path, graph, alleles):
    with open(output_path, 'w') as gfa_file:
        for allele_index in graph.nodes:
            gfa_file.write("S\t%s\t%s\n" % (str(allele_index),alleles[allele_index].sequence))

        for e in graph.edges:
            gfa_file.write("L\t%s\t+\t%s\t+\t0M\n" % (str(e[0]), str(e[1])))


def plot_graph(graph, ref_id_offset, ref_color, sample_color, output_path, line_style='-', draw_edge_weight_overlay=False):
    for layer, nodes in enumerate(networkx.topological_generations(graph)):
        # `multipartite_layout` expects the layer as a node attribute, so add the
        # numeric layer value as a node attribute
        for node in nodes:
            graph.nodes[node]["layer"] = layer

    color_map = []
    for n in graph.nodes:
        if n >= ref_id_offset:
            color_map.append(ref_color)
        else:
            color_map.append(sample_color)

    pos = networkx.multipartite_layout(graph, subset_key="layer")

    n_ref = graph.number_of_nodes() - ref_id_offset
    width = 1 + (n_ref*0.7)
    height = 1 + (n_ref*0.1)
    f = pyplot.figure(figsize=(width, height))
    a = pyplot.axes()

    node_size = 30

    networkx.draw(
        graph,
        pos,
        connectionstyle="arc3,rad=-0.35",
        node_color=color_map,
        node_size=node_size,
        style=line_style,
        font_size=4,
        width=0.6,
        arrowsize=3,
        with_labels=True)

    if draw_edge_weight_overlay:
        for edge in graph.edges(data='weight'):
            print(edge)

            networkx.draw_networkx_edges(
                graph,
                pos,
                edgelist=[edge],
                width=edge[2],
                connectionstyle="arc3,rad=-0.35",
                arrowstyle='-',
                node_size=node_size,
                alpha=0.6,
                edge_color="#000000"
            )

    ylim = list(a.get_ylim())
    ylim[0] -= 0.5
    a.set_ylim(ylim)

    pyplot.savefig(output_path, dpi=300)


def get_region_from_bam(output_directory, bam_path, region_string, tokenator, timeout=60*20):
    prefix = os.path.basename(bam_path).split('.')[0]
    filename = prefix + "_" + region_string.replace(":","_") + ".bam"
    local_bam_path = os.path.join(output_directory, filename)

    # Enable caching by path name
    if os.path.exists(local_bam_path):
        return local_bam_path

    tokenator.update_environment()

    args = [
        "samtools",
        "view",
        "-bh",
        "-o", local_bam_path,
        bam_path,
        region_string
    ]

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

    return local_bam_path


def get_reads_from_bam(output_path, bam_path, token):
    samtools_args = ["samtools", "fasta", bam_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(samtools_args)+'\n')

        token.update_environment()
        try:
            p1 = subprocess.run(samtools_args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return output_path


def run_minigraph(output_directory, gfa_path, fasta_path):
    output_path = os.path.join(output_directory, "reads_vs_graph.gaf")

    # minigraph \
    # -cx lr \
    # -o reads_vs_graph.gaf \
    # graph.gfa \
    # reads.fasta \
    args = ["minigraph", "-c", "-x", "lr", "-o", output_path, gfa_path, fasta_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(args)+'\n')

        try:
            p1 = subprocess.run(args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return output_path


def path_recursion(graph, alleles, id, path_sequence=None):
    if path_sequence is None:
        path_sequence = list()

    path_sequence.append(id)
    out_edges = graph.out_edges(id)

    if len(out_edges) == 0:
        yield path_sequence
    else:
        for edge in out_edges:
            yield from path_recursion(graph=graph, alleles=alleles, id=edge[1], path_sequence=copy(path_sequence))


def enumerate_paths(alleles, graph, output_directory):
    start_id = next(networkx.topological_sort(graph))

    print("Starting path recursion from %d" % start_id)

    paths = [p for p in path_recursion(graph=graph, alleles=alleles, id=start_id)]

    output_path = os.path.join(output_directory, "paths.fasta")
    with open(output_path, 'w') as file:
        for path in paths:
            name = '_'.join([str(i) for i in path])
            file.write(">%s\n" % name)
            for i in path:
                file.write("%s" % alleles[i].sequence)

            file.write('\n')


def update_edge_weights_using_alignments(gaf_path, graph, max_width=5.0):
    max_weight = 0.0

    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            path = re.split('[><]', tokens[5][1:])

            for i in range(len(path) - 1):
                a = int(path[i])
                b = int(path[i+1])

                print(a,b)

                if graph.has_edge(a,b):
                    print(graph[a][b]["weight"])
                    graph[a][b]["weight"] += 1
                    print(graph[a][b]["weight"])

                    if graph[a][b]["weight"] > max_weight:
                        max_weight = graph[a][b]["weight"]

                if graph.has_edge(b,a):
                    graph[b][a]["weight"] += 1
                    if graph[b][a]["weight"] > max_weight:
                        max_weight = graph[b][a]["weight"]

    for edge in graph.edges:
        if graph[edge[0]][edge[1]]["weight"] > 0:
            graph[edge[0]][edge[1]]["weight"] /= max_weight/max_width

    return graph


def vcf_to_graph(ref_path, vcf_paths, chromosome, ref_start, ref_stop, ref_sample_name, padding=20000):
    region_string = chromosome + ":" + str(ref_start) + "-" + str(ref_stop)

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

    ref_start -= padding
    ref_stop += padding

    print()

    # Throw away hashes and keep unique alleles as a list
    alleles = list(alleles.values())

    # Construct a list of coordinates along the reference path which contain edges to VCF alleles
    ref_edges = defaultdict(lambda: [[],[]])

    for a,allele in enumerate(alleles):
        print(a)
        print(allele)
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


def main():
    output_directory = "/home/ryan/data/test_hapslap/output2/"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        exit("ERROR: output directory exists already: %s" % output_directory)

    ref_path = "/home/ryan/data/human/reference/chm13v2.0.fa"

    input_directory = "/home/ryan/code/hapslap/data/test/hprc/"

    sample_names = [
        "HG002",
        "HG00438",
        "HG005",
        "HG00621",
        "HG00673",
        "HG00733",
        "HG00735",
        "HG00741",
        "HG01071",
        "HG01106",
        "HG01109",
        "HG01123",
        "HG01175",
        "HG01243",
        "HG01258",
        "HG01358",
        "HG01361",
        "HG01891",
        "HG01928",
        "HG01952",
        "HG01978",
        "HG02055",
        "HG02080",
        "HG02109",
        "HG02145",
        "HG02148",
        "HG02257",
        "HG02486",
        "HG02559",
        "HG02572",
        "HG02622",
        "HG02630",
        "HG02717",
        "HG02723",
        "HG02818",
        "HG02886",
        "HG03098",
        "HG03453",
        "HG03486",
        "HG03492",
        "HG03516",
        "HG03540",
        "HG03579",
        "NA18906",
        "NA19240",
        "NA20129",
        "NA21309"
    ]

    data_per_sample = defaultdict(dict)

    for filename in os.listdir(input_directory):
        path = os.path.join(input_directory, filename)

        for name in sample_names:
            if name in filename:
                if filename.endswith(".bam"):
                    data_per_sample[name]["bam"] = path
                if filename.endswith(".vcf.gz"):
                    data_per_sample[name]["vcf"] = path

    vcf_paths = list()
    bam_paths = list()

    for name in data_per_sample:
        has_vcf = False
        has_bam = False

        if "vcf" in data_per_sample[name]:
            has_vcf = True
            vcf_paths.append(data_per_sample[name]["vcf"])
        if "bam" in data_per_sample[name]:
            has_bam = True
            bam_paths.append(data_per_sample[name]["bam"])

        if (not has_vcf) or (not has_bam):
            exit("ERROR: sample does not have both a BAM and a VCF")

        print(name)
        print("vcf", data_per_sample[name]["vcf"])
        print("bam", data_per_sample[name]["bam"])

    chromosome = "chr20"
    ref_start = 47474020
    ref_stop = 47477018

    # ref_start = 54974920
    # ref_stop = 54977307

    region_string = chromosome + ":" + str(ref_start) + "-" + str(ref_stop)

    tokenator = GoogleToken()

    csv_path = os.path.join(output_directory, "nodes.csv")

    # TODO: make this a function of the actual reference version used
    ref_sample_name = "ref"

    graph,alleles = vcf_to_graph(
        ref_path=ref_path,
        vcf_paths=vcf_paths,
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        ref_sample_name=ref_sample_name)

    ref_id_offset = None
    for i in range(len(alleles)):
        if alleles[i].samples == {ref_sample_name}:
            ref_id_offset = i
            break

    print("ref_id_offset:", ref_id_offset)

    output_gfa_path = os.path.join(output_directory, "graph.gfa")

    sample_color = "#bebebe"
    ref_color = "#007cbe"

    # Write a csv that keeps track of ref coords, is_ref
    with open(csv_path, 'w') as csv_file:
        csv_file.write("id,ref_start,ref_stop,is_ref,color\n")

        for allele_index,allele in enumerate(alleles):
            color = sample_color

            is_ref = False
            if allele_index >= ref_id_offset:
                color = ref_color
                is_ref = True

            csv_file.write(','.join(list(map(str,[allele_index,allele.start,allele.stop,int(is_ref),color]))))
            csv_file.write('\n')

    # Write the GFA
    write_graph_to_gfa(output_path=output_gfa_path, graph=graph, alleles=alleles)

    # Plot the graph
    plot_path = os.path.join(output_directory, "dag.png")
    plot_graph(
        graph=graph,
        ref_id_offset=ref_id_offset,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path)

    # Remove empty nodes
    empty_nodes = list()
    for n in graph.nodes:
        if len(alleles[n].sequence) == 0:
            empty_nodes.append(n)

    for n in empty_nodes:
        a_nodes = [e[0] for e in graph.in_edges(n)]
        b_nodes = [e[1] for e in graph.out_edges(n)]

        for a in a_nodes:
            for b in b_nodes:
                graph.add_edge(a,b, weight=0)

        graph.remove_node(n)

    # Write the GFA
    output_gfa_path = output_gfa_path.replace(".gfa", "_no_empty.gfa")
    write_graph_to_gfa(output_path=output_gfa_path, graph=graph, alleles=alleles)

    # Plot the graph
    plot_path = os.path.join(output_directory, "dag_no_empty.png")
    plot_graph(
        graph=graph,
        ref_id_offset=ref_id_offset,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path)

    # Enumerate paths
    # enumerate_paths(alleles=alleles, graph=graph, output_directory=output_directory)

    output_fasta_path = os.path.join(output_directory, "reads.fasta")

    for path in bam_paths:
        print("bam_path:", path)

        region_bam_path = get_region_from_bam(
            output_directory=output_directory,
            bam_path=path,
            region_string=region_string,
            tokenator=tokenator)

        # This repeatedly appends one FASTA file
        get_reads_from_bam(
            output_path=output_fasta_path,
            bam_path=region_bam_path,
            token=tokenator)

        os.remove(region_bam_path)

    output_gaf_path = run_minigraph(
        output_directory=output_directory,
        gfa_path=output_gfa_path,
        fasta_path=output_fasta_path)

    graph = update_edge_weights_using_alignments(gaf_path=output_gaf_path, graph=graph)

    plot_path = os.path.join(output_directory, "aligned_dag.png")
    plot_graph(
        graph=graph,
        ref_id_offset=ref_id_offset,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path,
        line_style=':',
        draw_edge_weight_overlay=True)


if __name__ == "__main__":
    main()
