import sys

from modules.IterativeHistogram import IterativeHistogram
from modules.IncrementalIdMap import IncrementalIdMap
from modules.Authenticator import *
from modules.Paths import Paths
from modules.GsUri import *

from collections import defaultdict
from multiprocessing import Pool
from copy import deepcopy,copy
import subprocess
import argparse
import os.path
import re

from ortools.graph.python import min_cost_flow
from matplotlib import pyplot
from networkx import DiGraph
from pysam import FastaFile
from vcf import VCFReader
import networkx
import numpy
import pysam


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


def plot_graph(
        graph,
        ref_id_offset,
        ref_color,
        sample_color,
        output_path,
        line_style='-',
        draw_edge_weight_overlay=False,
        figure_width=None,
        figure_height=None,
        connection_style="arc3,rad=-0.35"):

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

    width = None
    height = None

    if figure_width is None:
        width = 1 + (n_ref*0.7)
    else:
        width = figure_width

    if figure_height is None:
        height = 1 + (n_ref*0.1)
    else:
        height = figure_height

    f = pyplot.figure(figsize=(width, height))
    a = pyplot.axes()

    node_size = 30

    networkx.draw(
        graph,
        pos,
        connectionstyle=connection_style,
        node_color=color_map,
        node_size=node_size,
        style=line_style,
        font_size=4,
        width=0.2,
        arrowsize=3,
        with_labels=True)

    if draw_edge_weight_overlay:
        for edge in graph.edges(data='weight'):
            networkx.draw_networkx_edges(
                graph,
                pos,
                edgelist=[edge],
                width=edge[2],
                connectionstyle=connection_style,
                arrowstyle='-',
                node_size=node_size,
                alpha=0.6,
                edge_color="#007cbe")

    ylim = list(a.get_ylim())
    ylim[0] -= 0.5
    a.set_ylim(ylim)

    pyplot.savefig(output_path, dpi=400)


def optimize_with_flow(
        paths,
        read_id_map,
        path_to_read_costs,
        output_directory):

    solver = min_cost_flow.SimpleMinCostFlow()
    source_to_read_arcs = list()
    read_to_path_arcs = list()
    path_to_sink_arcs = list()

    source_id = len(read_id_map) + len(paths)
    sink_id = source_id + 1

    #
    #           [r0]---[h2]
    #          /    \ /    \
    #  [SOURCE]      x      [SINK]
    #          \    / \    /
    #           [r1]---[h1]
    #

    flow_graph = DiGraph()
    flow_graph.add_node(source_id)
    flow_graph.add_node(sink_id)

    # source --> read
    for id,name in read_id_map:
        a = solver.add_arc_with_capacity_and_unit_cost(tail=source_id,head=id,unit_cost=0,capacity=1)
        source_to_read_arcs.append(a)
        solver.set_node_supply(node=id,supply=0)

        flow_graph.add_edge(source_id,id,weight=0)
        flow_graph.add_node(id)

    # path --> sink
    for id,name,_ in paths:
        a = solver.add_arc_with_capacity_and_unit_cost(tail=id,head=sink_id,unit_cost=0,capacity=len(read_id_map))
        path_to_sink_arcs.append(a)
        solver.set_node_supply(node=id,supply=0)

        flow_graph.add_node(id)
        flow_graph.add_edge(id,sink_id,weight=0)

    # read --> path (haplotype)
    for path_id, read_costs in enumerate(path_to_read_costs):
        for read_id,cost in read_costs.items():
            a = solver.add_arc_with_capacity_and_unit_cost(tail=read_id,head=path_id,unit_cost=cost,capacity=1)
            read_to_path_arcs.append(a)

            flow_graph.add_edge(read_id,path_id,weight=cost)

    solver.set_node_supply(node=source_id,supply=len(read_id_map))
    solver.set_node_supply(node=sink_id,supply=-len(read_id_map))

    status = solver.solve_max_flow_with_min_cost()

    if status != solver.OPTIMAL:
        print('There was an issue with the min cost flow input.')
        print(f'Status: {status}')
        exit(1)

    print("Minimum cost: %d" % solver.optimal_cost())
    print("source_to_read_arcs:")
    for arc in source_to_read_arcs:
        a = solver.tail(arc)
        b = solver.head(arc)
        flow = solver.flow(arc)
        print(a, b, flow)
        flow_graph[a][b]["weight"] = flow

    print("read_to_path_arcs:")
    for arc in read_to_path_arcs:
        a = solver.tail(arc)
        b = solver.head(arc)
        flow = solver.flow(arc)
        print(a, b, flow)
        flow_graph[a][b]["weight"] = flow

    print("path_to_sink_arcs:")
    for arc in path_to_sink_arcs:
        a = solver.tail(arc)
        b = solver.head(arc)
        flow = solver.flow(arc)
        print(a, b, flow)
        flow_graph[a][b]["weight"] = flow

    sample_color = "#bebebe"
    ref_color = "#6e6e6e"

    plot_path = os.path.join(output_directory, "flow_graph.png")
    plot_graph(
        graph=flow_graph,
        ref_id_offset=-1,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path,
        line_style=':',
        figure_width=4,
        figure_height=4,
        connection_style="arc3",
        draw_edge_weight_overlay=True,
    )


def enumerate_paths_using_alignments(paths: Paths, graph: DiGraph, alleles, gaf_path, ref_sample_name, output_directory):
    start_node = None
    stop_node = None

    # For this implementation we enforce that reads must span the full graph, so start and end ref nodes are needed.
    # TODO: A smarter implementation might allow partial walks to be extended using all walks from a sample.
    # But it would need to be capable of aborting intractable traversals as a result of too many walks
    for n in graph.nodes():
        if ref_sample_name in alleles[n].samples:
            if len(graph.in_edges(n)) == 0:
                start_node = n
            if len(graph.out_edges(n)) == 0:
                stop_node = n

    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            # Column 5 (0-based) is the path column in a GAF file
            # It can be a forward or reverse alignment, so we will reinterpret them all as forward alignments
            path = None
            if tokens[5][0] == '>':
                path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            if tokens[5][0] == '<':
                path = tuple(map(int,reversed(re.findall(r'\d+', tokens[5]))))

            if len(path) > 1 and path[0] == start_node and path[-1] == stop_node:
                paths.increment_weight(path,1)

    # Write one sequence per FASTA so that all reads can be aligned to each independently
    output_paths = list()
    csv_path = os.path.join(output_directory, "paths.csv")
    with open(csv_path, 'w') as csv_file:
        csv_file.write("index,path,frequency\n")

        for p,path,frequency in paths:
            name = "_".join(map(str,path))

            csv_file.write("%d,%s,%d" % (p,name,frequency) + '\n')

            output_path = os.path.join(output_directory, str(p) + ".fasta")
            output_paths.append(output_path)

            with open(output_path, 'w') as file:
                file.write(">%s\n" % name)
                for i in path:
                    file.write("%s" % alleles[i].sequence)
                file.write("\n")

    return output_paths


def update_edge_weights_using_alignments(gaf_path, graph: DiGraph, min_width=0.5, max_width=5.0):
    max_weight = 0.0

    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            print(path)

            for i in range(len(path) - 1):
                a = int(path[i])
                b = int(path[i+1])

                if graph.has_edge(a,b):
                    print(a,b)
                    graph[a][b]["weight"] += 1
                    if graph[a][b]["weight"] > max_weight:
                        max_weight = graph[a][b]["weight"]

                if graph.has_edge(b,a):
                    print(b,a)
                    graph[b][a]["weight"] += 1
                    if graph[b][a]["weight"] > max_weight:
                        max_weight = graph[b][a]["weight"]

    for edge in graph.edges:
        w = graph[edge[0]][edge[1]]["weight"]
        if w > 0:
            graph[edge[0]][edge[1]]["weight"] = min_width + (w/max_weight)*(max_width - min_width)


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


# Requires samtools installed!
def run_minimap2(ref_fasta_path, reads_fasta_path, preset, n_threads, n_sort_threads, output_directory, filename_prefix="reads_vs_ref"):
    output_filename = os.path.join(output_directory, filename_prefix + ".bam")
    output_path = os.path.join(output_directory,output_filename)

    minimap_args = ["minimap2", "-a", "-x", preset, "--eqx", "-t", str(n_threads), ref_fasta_path, reads_fasta_path]
    sort_args = ["samtools", "sort", "-", "-@", str(n_sort_threads), "-o", output_filename]

    index_args = ["samtools", "index", output_filename]

    with open(output_path, 'w') as file:
        sys.stderr.write(" ".join(minimap_args)+'\n')
        sys.stderr.write(" ".join(sort_args)+'\n')

        p1 = subprocess.Popen(minimap_args, stdout=subprocess.PIPE, cwd=output_directory)
        p2 = subprocess.Popen(sort_args, stdin=p1.stdout, stdout=file, cwd=output_directory)
        p2.communicate()

    success = (p2.returncode == 0)

    if not success:
        sys.stderr.write("ERROR: failed to align: %s to %s\n" % (reads_fasta_path, ref_fasta_path))
        sys.stderr.flush()
        return None

    sys.stderr.write(" ".join(index_args)+'\n')

    try:
        p1 = subprocess.run(index_args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return output_path


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
    n_threads = 30
    n_sort_threads = 4

    chromosome = "chr20"
    ref_start = 47474020
    ref_stop = 47477018

    # ref_start = 54974920
    # ref_stop = 54977307

    region_string = chromosome + ":" + str(ref_start) + "-" + str(ref_stop)

    output_directory = os.path.join("/home/ryan/data/test_hapslap/", region_string.replace(':',"_"))

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
    ref_color = "#6e6e6e"

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

    update_edge_weights_using_alignments(gaf_path=output_gaf_path, graph=graph)

    plot_path = os.path.join(output_directory, "dag_aligned.png")
    plot_graph(
        graph=graph,
        ref_id_offset=ref_id_offset,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path,
        line_style=':',
        draw_edge_weight_overlay=True)

    path_subdirectory = os.path.join(output_directory, "paths")
    os.makedirs(path_subdirectory)

    paths = Paths()

    fasta_paths = enumerate_paths_using_alignments(
        paths=paths,
        graph=graph,
        alleles=alleles,
        gaf_path=output_gaf_path,
        ref_sample_name=ref_sample_name,
        output_directory=path_subdirectory)

    bam_paths = list()
    for p,path in enumerate(fasta_paths):
        bam_path = run_minimap2(
            ref_fasta_path=path,
            reads_fasta_path=output_fasta_path,
            preset="map-hifi",
            n_threads=n_threads,
            n_sort_threads=n_sort_threads,
            output_directory=path_subdirectory,
            filename_prefix=str(p))

        bam_paths.append(bam_path)

    pyplot.close("all")

    fig = pyplot.figure()
    axes = pyplot.axes()

    histogram = IterativeHistogram(start=0, stop=1000, n_bins=200)

    path_to_read_costs = [defaultdict(int) for x in range(len(paths))]
    read_id_map = IncrementalIdMap(offset=len(paths))

    # Because of the potential for supplementary alignments, need to first aggregate costs per read
    for bam_path in bam_paths:
        bam = pysam.AlignmentFile(bam_path, 'rb')

        print(bam_path)

        total = 0.0
        n = 0.0
        for alignment in bam:
            # TODO: verify that the NM tag is not double-counting softclips/hardclips in a supplementary alignment?
            edit_distance = alignment.get_tag("NM")
            path_name = bam.header.get_reference_name(alignment.reference_id)
            path_id = paths.get_path_id(path_name)
            read_name = alignment.query_name
            read_id = read_id_map.add(read_name)

            if not alignment.is_secondary:
                path_to_read_costs[path_id][read_id] += edit_distance
                print(path_id, read_id)

                total += edit_distance
                n += 1

        avg_edit_distance = total/n
        histogram.update(avg_edit_distance)

    axes.plot(histogram.get_bin_centers(), histogram.get_histogram(), color="#007cbe")

    reads_csv_path = os.path.join(output_directory, "reads.csv")
    read_id_map.write_to_file(reads_csv_path)

    # c: cost (edit distance)
    # r: read
    # s: sample
    # h: haplotype

    # cost = sum(c_r_h)
    #
    # boolean constraint: variables are boolean, and indicate read assignment to haplotypes
    # for all r,h: 0 <= r_h <= 1
    #
    # completeness constraint: all reads must be assigned to exactly one haplotype
    # for all r: sum(r_h,h) = 1
    #
    # ploidy_constraint: each sample must be aligned to at most 2 haps
    # Must first group reads r by samples s, then for each sample group:
    # is_active_h = sum(r_h,r) > 0
    # x = sum(r_h,r)
    #
    #


    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
