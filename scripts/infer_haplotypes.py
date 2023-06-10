from threading import Timer
import time
import math
import sys

from modules.IterativeHistogram import IterativeHistogram
from modules.IncrementalIdMap import IncrementalIdMap
from modules.Bam import get_region_from_bam
from modules.Align import run_minimap2
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
from ortools.sat.python import cp_model
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
        self.is_left_flank = False
        self.is_right_flank = False

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


"""
Taken from:
https://stackoverflow.com/questions/73996059/ortools-cp-sat-how-to-fetch-the-feasible-solution-if-the-best-incumbent-solu
"""
class ObjectiveEarlyStopping(cp_model.CpSolverSolutionCallback):
    def __init__(self, timer_limit: int):
        super(ObjectiveEarlyStopping, self).__init__()
        self._timer_limit = timer_limit
        self._timer = None

    def on_solution_callback(self):
        self._reset_timer()

    def _reset_timer(self):
        if self._timer:
            self._timer.cancel()

        self._timer = Timer(self._timer_limit, self.StopSearch)
        self._timer.start()

    def StopSearch(self):
        print(f"{self._timer_limit} seconds without improvement")
        super().StopSearch()
        self._timer.cancel()

    def cancel_timer_thread(self):
        self._timer.cancel()

class Variables:
    def __init__(self):
        # key = (path_id,read_id)
        self.path_to_read = dict()

        # key = (path_id,sample_id)
        self.path_to_sample = dict()

        # key = path_id
        self.path = dict()

        # key = integer of n >= 1
        self.n = dict()

        # Objective a, explicitly minimized
        self.cost_a = None

        # Objective b, marginalized
        self.cost_b = None

        # Return code of optimization
        self.status = None

        # Optimizer stats
        self.response_stats = None

        # Cache contains integer values instead of the cp_model var objects, and can therefore be used downstream
        self.is_cache = False

    '''
    Construct a copy of the variables object, where the values are stored instead of the cp_model variable objects. 
    '''
    def get_cache(self, status, solver: cp_model.CpSolver):
        cache = Variables()

        if not (status == cp_model.OPTIMAL or status == cp_model.FEASIBLE):
            sys.stderr.write("WARNING: non-optimal result from solver")
            cache.cost_a = None
            cache.cost_b = None

        else:
            for key in self.path_to_read:
                cache.path_to_read[key] = solver.Value(self.path_to_read[key])

            for key in self.path_to_sample:
                cache.path_to_sample[key] = solver.Value(self.path_to_sample[key])

            for key in self.path:
                cache.path[key] = solver.Value(self.path[key])

            for key in self.n:
                cache.n[key] = solver.Value(self.n[key])

            cache.cost_a = solver.Value(self.cost_a)
            cache.cost_b = solver.Value(self.cost_b)

        cache.response_stats = solver.ResponseStats()
        cache.status = status
        cache.is_cache = True

        return cache

    def write_results_to_file(self, output_dir: str):
        if not self.is_cache:
            raise Exception("ERROR: cannot write live instance of variables, must get cache instead")

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        output_path = os.path.join(output_dir, "path_to_read.csv")
        with open(output_path, 'w') as file:
            file.write("path_id,read_id,value\n")
            for key,value in self.path_to_read.items():
                file.write("%d,%d,%d" % (key[0],key[1],value))
                file.write('\n')

        output_path = os.path.join(output_dir, "path_to_sample.csv")
        with open(output_path, 'w') as file:
            file.write("path_id,sample_id,value\n")
            for key,value in self.path_to_sample.items():
                file.write("%d,%d,%d" % (key[0],key[1],value))
                file.write('\n')

        output_path = os.path.join(output_dir, "path.csv")
        with open(output_path, 'w') as file:
            file.write("path_id,value\n")
            for item in self.path.items():
                file.write(','.join(list(map(str,item))))
                file.write('\n')

        output_path = os.path.join(output_dir, "n.csv")
        with open(output_path, 'w') as file:
            file.write("n,value\n")
            for item in self.n.items():
                file.write(','.join(list(map(str,item))))
                file.write('\n')

        output_path = os.path.join(output_dir, "stats.csv")
        with open(output_path, 'w') as file:
            file.write(self.response_stats)


def optimize_with_cpsat(
        path_to_read_costs,
        reads: IncrementalIdMap,
        paths: Paths,
        sample_to_reads: dict,
        output_dir: str,
        n_threads: int = 1):

    output_subdir = os.path.join(output_dir, "optimizer")

    vars = Variables()

    model = cp_model.CpModel()

    # Define read assignment variables
    for edge in path_to_read_costs.keys():
        # 'edge' is a tuple with path_id,read_id
        vars.path_to_read[edge] = model.NewIntVar(0, 1, "p%dr%d" % edge)

    # Constraint: each read must map to only one haplotype/path
    for read_id in reads.ids():
        model.Add(sum([vars.path_to_read[(path_id,read_id)] for path_id in paths.ids()]) == 1)

    # Constraint: each sample's reads must map to at most two haplotypes/paths
    # Use a boolean indicator to tell whether any of a sample's reads are assigned to each haplotype/path
    for sample_id,read_group in sample_to_reads.items():
        if not type(sample_id) == int:
            raise Exception("ERROR: non integer ids used for sample: %s" % str(sample_id))

        for id in read_group:
            if not type(id) == int:
                raise Exception("ERROR: non integer ids used for read: %s in sample: %s" % (str(id),str(sample_id)))

        for path_id in paths.ids():
            edge = (path_id,sample_id)
            vars.path_to_sample[edge] = model.NewBoolVar("p%ds%d" % edge)

            s = sum([vars.path_to_read[(path_id,read_id)] for read_id in read_group])
            model.Add(s >= 1).OnlyEnforceIf(vars.path_to_sample[edge])
            model.Add(s == 0).OnlyEnforceIf(vars.path_to_sample[edge].Not())

    # Now that the boolean indicators have been defined, use them to add a constraint on ploidy per sample
    for sample_id,read_group in sample_to_reads.items():
        model.Add(sum([vars.path_to_sample[(path_id,sample_id)] for path_id in paths.ids()]) <= 2)

    for path_id in paths.ids():
        vars.path[path_id] = model.NewBoolVar("p" + str(path_id))

        # Accumulate all possible assignments of this path to any reads
        s = sum([vars.path_to_read[(path_id,read_id)] for read_id in reads.ids()])
        model.Add(s >= 1).OnlyEnforceIf(vars.path[path_id])
        model.Add(s == 0).OnlyEnforceIf(vars.path[path_id].Not())

    # Cost term a: sum of edit distances for all reads assigned to haplotypes
    vars.cost_a = model.NewIntVar(0, 100_000_000, "cost_a")
    model.Add(vars.cost_a == sum([c*vars.path_to_read[e] for e,c in path_to_read_costs.items()]))

    # Cost term b: sum of unique haplotypes used
    vars.cost_b = model.NewIntVar(0, 100_000_000, "cost_b")
    model.Add(vars.cost_b == sum([x for x in vars.path.values()]))

    # n diploid samples can at most fill n*2 haplotypes
    # Sometimes there are may be fewer candidate paths than that
    max_feasible_haplotypes = min(len(paths), len(sample_to_reads)*2) + 1

    for i in range(1,max_feasible_haplotypes):
        vars.n[i] = model.NewBoolVar("n" + str(i))
        model.Add(vars.cost_b == i).OnlyEnforceIf(vars.n[i])
        model.Add(vars.cost_b != i).OnlyEnforceIf(vars.n[i].Not())

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = n_threads
    # solver.parameters.max_time_in_seconds = 300.0
    # solver.parameters.log_search_progress = True
    # solver.log_callback = print

    status = None
    results = dict()

    for i in range(1,max_feasible_haplotypes):
        print(i)
        model.ClearAssumptions()
        model.AddAssumption(vars.n[i])
        model.Minimize(vars.cost_a)

        if status == cp_model.OPTIMAL:
            model.ClearHints()
            for var in vars.path_to_read.values():
                model.AddHint(var, solver.Value(var))

        # It's critical to STOP THE TIMER after this finishes because it relies on a thread which will run on
        # past the optimal solution for as long as the timer is set, potentially starving future iterations of threads.
        o = ObjectiveEarlyStopping(60)
        status = solver.SolveWithSolutionCallback(model, o)
        o.cancel_timer_thread()

        print("=====Stats:======")
        print(solver.SolutionInfo())
        print(solver.ResponseStats())

        results[i] = vars.get_cache(status=status, solver=solver)

        if len(results) > 1 and results[i].cost_a > results[i-1].cost_a:
            sys.stderr.write("Iteration stopped at n=%d because score worsened\n" % i)
            break

    for i,cache in results.items():
        print("%d,%d" % (i, cache.cost_a))
        cache.write_results_to_file(output_dir=os.path.join(output_subdir,str(i)))

    return results


def get_spanning_reads(alleles, gaf_path):
    spanning_reads = set()

    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')
            read_name = tokens[0]

            # Column 5 (0-based) is the path column in a GAF file
            # It can be a forward or reverse alignment, so we will reinterpret them all as forward alignments
            path = None
            if tokens[5][0] == '>':
                path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            if tokens[5][0] == '<':
                path = tuple(map(int,reversed(re.findall(r'\d+', tokens[5]))))

            if len(path) > 1 and alleles[path[0]].is_left_flank and alleles[path[-1]].is_right_flank:
                spanning_reads.add(read_name)

    return spanning_reads


def enumerate_paths_using_alignments(alleles, gaf_path: str, output_directory: str, min_coverage: int, read_subset: set=None):
    paths = Paths()

    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')
            read_name = tokens[0]

            if read_subset is not None and read_name not in read_subset:
                print("skipping read: " + read_name)
                continue

            # Column 5 (0-based) is the path column in a GAF file
            # It can be a forward or reverse alignment, so we will reinterpret them all as forward alignments
            path = None
            if tokens[5][0] == '>':
                path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            elif tokens[5][0] == '<':
                path = tuple(map(int,reversed(re.findall(r'\d+', tokens[5]))))

            else:
                print(tokens[5])
                raise Exception("ERROR: Non GFA character found in path column of GFA")

            if len(path) > 1 and alleles[path[0]].is_left_flank and alleles[path[-1]].is_right_flank:
                paths.increment_weight(path,1)

    filtered_paths = Paths()

    for p,path,frequency in paths:
        if frequency >= min_coverage:
            filtered_paths.add_path(path,frequency)
            print(path)

    paths = filtered_paths

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

    if len(paths) == 0:
        raise Exception("ERROR: no passing paths in alignment graph")

    return paths, output_paths


def update_edge_weights_using_alignments(alleles, gaf_path, graph: DiGraph, min_width=0.5, max_width=5.0):
    max_weight = 0.0

    with open(gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            # Only count spanning reads
            if not (len(path) > 1 and alleles[path[0]].is_left_flank and alleles[path[-1]].is_right_flank):
                continue

            for i in range(len(path) - 1):
                a = int(path[i])
                b = int(path[i+1])

                if graph.has_edge(a,b):
                    graph[a][b]["weight"] += 1
                    if graph[a][b]["weight"] > max_weight:
                        max_weight = graph[a][b]["weight"]

                if graph.has_edge(b,a):
                    graph[b][a]["weight"] += 1
                    if graph[b][a]["weight"] > max_weight:
                        max_weight = graph[b][a]["weight"]

    for edge in graph.edges:
        w = graph[edge[0]][edge[1]]["weight"]
        if w > 0:
            graph[edge[0]][edge[1]]["weight"] = min_width + (w/max_weight)*(max_width - min_width)


def get_read_names_from_bam(bam_path):
    read_names = list()

    bam = pysam.AlignmentFile(bam_path, 'rb')

    print(bam_path)

    for alignment in bam:
        read_names.append(alignment.query_name)

    return read_names


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
    # Get start node
    start_id = None

    for i in range(len(alleles)):
        if alleles[i].is_left_flank:
            start_id = i

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


def vcf_to_graph(ref_path, vcf_paths, chromosome, ref_start, ref_stop, ref_sample_name, flank_length):
    region_string = chromosome + ":" + str(ref_start) + "-" + str(ref_stop)

    ref_sequence = FastaFile(ref_path).fetch("chr20")

    alleles = dict()

    for vcf_path in vcf_paths:
        sample_name = os.path.basename(vcf_path).split("_")[0]

        with open(vcf_path, 'rb') as file:
            print(vcf_path)
            vcf = VCFReader(file)

            records = vcf.fetch(region_string)

            for record in records:
                print()
                print("var_subtype:\t%s" % record.var_subtype)
                print("start:\t\t%d" % record.start)

                # One VCF per sample, no iterating of samples needed
                call = record.samples[0]
                gt = [int(call.data.GT[0]) if call.data.GT[0] != '.' else 0, int(call.data.GT[-1]) if call.data.GT[-1] != '.' else 0]

                print("gt:\t\t%d/%d" % (gt[0], gt[1]))
                print("ref_length:\t%d" % len(record.alleles[0]))

                b_length = None
                if record.var_subtype == "DUP":
                    print(record.INFO)
                    b_length = record.INFO["SVLEN"]
                else:
                    b_length = len(record.alleles[gt[1]])

                print("a_length:\t%d" % len(record.alleles[gt[0]]))
                print("b_length:\t%d" % b_length)
                print("is_sv_precise:\t%d" % record.is_sv_precise)

                # Iterate unique, non-ref alleles only
                for allele_index in set(gt):
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
                        if ref_stop < start < ref_start or ref_stop < stop < ref_start:
                            raise Exception("ERROR: VCF allele with coords %d-%d extends out of region %s:%d-%d" % (start,stop,chromosome,ref_start,ref_stop))

                        if sequence.strip() == "<DUP>":
                            sequence = ref_sequence[start:start+b_length+1]
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

    ref_start -= flank_length
    ref_stop += flank_length

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


def infer_haplotypes(chromosome, ref_start, ref_stop, sample_names):
    n_threads = 30
    n_sort_threads = 4
    flank_length = 20000

    # Used by initial minigraph alignment to select paths for the optimizer to cover
    min_coverage = 2

    region_string = chromosome + ":" + str(ref_start) + "-" + str(ref_stop)

    output_directory = os.path.join("/home/ryan/data/test_hapslap/results", region_string.replace(':',"_"))

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        exit("ERROR: output directory exists already: %s" % output_directory)

    ref_path = "/home/ryan/data/human/reference/chm13v2.0.fa"

    input_directory = "/home/ryan/code/hapslap/data/test/hprc/"

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
        ref_sample_name=ref_sample_name,
        flank_length=flank_length)

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

    sample_to_reads = defaultdict(list)
    for path in bam_paths:
        print("bam_path:", path)

        region_bam_path = get_region_from_bam(
            output_directory=output_directory,
            bam_path=path,
            region_string=region_string,
            tokenator=tokenator)

        read_names = get_read_names_from_bam(region_bam_path)

        # TODO: catalog samples in a more robust way
        sample_name = os.path.basename(path).replace(".bam","").split('_')[0]
        sample_to_reads[sample_name] = read_names

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

    update_edge_weights_using_alignments(
        alleles=alleles,
        gaf_path=output_gaf_path,
        graph=graph)

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

    spanning_reads = get_spanning_reads(
        alleles=alleles,
        gaf_path=output_gaf_path,
    )

    paths, fasta_paths = enumerate_paths_using_alignments(
        alleles=alleles,
        gaf_path=output_gaf_path,
        output_directory=path_subdirectory,
        min_coverage=min_coverage,
        read_subset=spanning_reads
    )

    sys.stderr.write("Optimizing using %d paths with min coverage %d\n" % (len(paths), min_coverage))

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

    read_id_map = IncrementalIdMap(offset=len(paths))
    path_to_read_costs = defaultdict(int)

    # Because of the potential for supplementary alignments, need to first aggregate costs per read
    for bam_path in bam_paths:
        sys.stderr.write(bam_path)
        sys.stderr.write('\n')

        bam = pysam.AlignmentFile(bam_path, 'rb')

        total = 0.0
        n = 0.0
        for alignment in bam:
            read_name = alignment.query_name

            # TODO: [optimize] eventually alignment should be done in memory, and this filtering step could be performed
            # before aligning the reads with minimap2, which would save (potentially considerable) time
            if read_name not in spanning_reads:
                continue

            read_id = read_id_map.add(read_name)

            # TODO: [fix] verify that the NM tag is not double-counting softclips/hardclips in a supplementary alignment?
            path_name = bam.header.get_reference_name(alignment.reference_id)
            path_id = paths.get_path_id(path_name)

            # TODO: more filtering? possible to have unmapped reads?
            if not alignment.is_secondary:
                edit_distance = alignment.get_tag("NM")
                path_to_read_costs[(path_id,read_id)] += edit_distance

                total += edit_distance
                n += 1

        avg_edit_distance = total/n
        histogram.update(avg_edit_distance)

    axes.plot(histogram.get_bin_centers(), histogram.get_histogram(), color="#007cbe")

    reads_csv_path = os.path.join(output_directory, "reads.csv")
    read_id_map.write_to_file(reads_csv_path)

    sample_id_map = IncrementalIdMap()
    for name in sample_names:
        sample_id_map.add(name)

    sample_id_path = os.path.join(output_directory, "samples.csv")
    with open(sample_id_path, 'w') as file:
        for id,name in sample_id_map:
            file.write("%d,%s\n" % (id,name))

    # Make a duplicate of the sample->read mapping, using ids instead of names
    sample_id_to_read_ids = dict()
    for sample_name,read_names in sample_to_reads.items():
        # Sometimes a read may not be in the id_map, because the id_map is generated only using reads that
        # were aligned in the BAM of reads-to-path. We skip any read that doesn't appear in the mappings.
        read_ids = list(filter(None,[read_id_map.get_id(x) for x in read_names]))

        sample_id = sample_id_map.get_id(sample_name)
        sample_id_to_read_ids[sample_id] = read_ids

    results = optimize_with_cpsat(
        path_to_read_costs=path_to_read_costs,
        reads=read_id_map,
        paths=paths,
        sample_to_reads=sample_id_to_read_ids,
        n_threads=30,
        output_dir=output_directory
    )

    # TODO: replace this with average error rate per read?
    c_target = 0

    # TODO: replace this with regional heterozygosity estimate
    n_target = 0

    histogram_path = os.path.join(output_directory, "path_global_score_histogram.png")

    pyplot.savefig(histogram_path, dpi=200)
    pyplot.close('all')

    fig,axes = pyplot.subplots(nrows=1, ncols=1)

    c_start = sys.maxsize
    c_stop = 0
    n_start = 1
    n_stop = len(sample_names)

    for n,cache in results.items():
        if cache.cost_a < c_start:
            c_start = cache.cost_a

        if cache.cost_a > c_stop:
            c_stop = cache.cost_a

    print("n",n_start,n_stop)
    print("c",c_start,c_stop)

    d_min = sys.maxsize
    c_min = None
    n_min = None
    c_per_read_min = None
    n_norm_min = None
    c_norm_min = None

    for n,cache in results.items():
        c_norm = float(cache.cost_a - c_start) / float(c_stop - c_start) + 1
        n_norm = (float(n) - n_start) / float(n_stop - n_start) + 0.001

        distance = math.sqrt((c_norm-c_target)**2 + (n_norm-n_target)**2)
        if distance < d_min:
            d_min = distance
            n_min = n
            c_per_read_min = float(cache.cost_a)/float(len(read_id_map))
            c_min = cache.cost_a

            # For plotting
            n_norm_min = n_norm
            c_norm_min = c_norm

        # axes[0].scatter(n, cache.cost_a, color="C0")
        axes.scatter(n_norm, c_norm, color="C0")
        axes.plot([n_target,n_norm], [c_target,c_norm], color="gray")
        axes.text(n_norm, c_norm, "%.4f" % distance)

    best_result = results[n_min]

    # axes[0].plot(n_min, c_min, color="C1", marker='o', markerfacecolor='none', markersize=6)
    # axes[0].text(n_min, c_min, "n=%d" % n_min, color="red", ha='right', va='top')

    axes.plot(n_norm_min, c_norm_min, color="C1", marker='o', markerfacecolor='none', markersize=6)
    axes.text(n_norm_min, c_norm_min, "n=%d" % n_min, color="red", ha='right', va='top')

    # axes.set_aspect('equal', 'box')
    axes.set_title("Results for %s" % region_string)
    axes.set_xlabel("n_norm")
    axes.set_ylabel("c_norm")

    summary_string = "n,%d\ncost_per_read,%.2f\ndistance,%.4f\n" % (n_min,c_per_read_min,d_min)
    sys.stderr.write(summary_string)

    summary_path = os.path.join(output_directory, "summary.txt")
    with open(summary_path, 'w') as file:
        file.write(summary_string)

    counter = defaultdict(int)

    # Get start and stop nodes
    start_node = None
    stop_node = None

    for i in range(len(alleles)):
        if alleles[i].is_left_flank:
            start_node = i

        if alleles[i].is_right_flank:
            stop_node = i

    print("Trimming start and stop nodes: %d %d " % (start_node,stop_node))
    # Trim flanks off the allele sequences
    s_start = alleles[start_node].sequence
    l_start = len(s_start)
    s_start = s_start[flank_length:]

    s_stop = alleles[stop_node].sequence
    l_stop = len(s_stop)
    s_stop = s_stop[:(l_stop - flank_length)]

    print("before trimming: ", l_start)
    print("after trimming: ", len(s_start))
    print("before trimming: ", l_stop)
    print("after trimming: ", len(s_stop))

    alleles[start_node].sequence = s_start
    alleles[stop_node].sequence = s_stop

    assigned_haplotypes_path = os.path.join(output_directory, "assigned_haplotypes.fasta")
    with open(assigned_haplotypes_path, 'w') as file:
        for [path_id,sample_id],value in best_result.path_to_sample.items():
            # Only consider combinations of path/sample that were actually assigned by the optimizer
            if value == 0:
                continue

            counter[sample_id] += 1
            sample_name = sample_id_map.get_name(sample_id)
            sequence_name = sample_name + "_" + str(counter[sample_id])

            file.write(">" + sequence_name)
            file.write("\n")
            for i in paths.get_path(path_id):
                print(i, len(alleles[i].sequence))
                file.write(alleles[i].sequence)

            file.write("\n")

    optimizer_results_path = os.path.join(output_directory, "optimizer_results.png")
    pyplot.savefig(optimizer_results_path, dpi=200)

    # pyplot.show()
    pyplot.close()

    return summary_string


def main():

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

    chromosome = "chr20"

    coords = [
        # [10437604,10440525],
        # [18259924,18261835],
        # [54975152,54976857], # Nightmare region (tandem)
        # [18689217,18689256],
        # [18828383,18828733],
        # [47475093,47475817],
        # [49404497,49404943],
        # [55000754,55000852],
        [55486867,55492722],
        # [7901318,7901522]
    ]

    summary_strings = list()
    for ref_start, ref_stop in coords:
        summary_string = infer_haplotypes(
            chromosome=chromosome,
            ref_start=ref_start,
            ref_stop=ref_stop,
            sample_names=sample_names
        )

        summary_strings.append(summary_string)

    for s in range(len(summary_strings)):
        sys.stderr.write(str(coords[s]))
        sys.stderr.write('\n')
        sys.stderr.write(summary_strings[s])


if __name__ == "__main__":
    main()
