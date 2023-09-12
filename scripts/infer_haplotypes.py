import ortools.sat.python.cp_model

from modules.Vcf import vcf_to_graph,write_paths_to_vcf,remove_empty_nodes_from_variant_graph,merge_vcfs_in_directory,write_graph_to_gfa,write_node_csv,compress_and_index_vcf
from modules.Align import run_minimap2,run_minigraph,run_mashmap,run_minimap2_on_read_subset,run_panaligner
from modules.Cigar import iterate_cigar,cigar_index_to_char,is_ref_move
from modules.IterativeHistogram import IterativeHistogram
from modules.IncrementalIdMap import IncrementalIdMap
from modules.Bam import get_region_from_bam
from modules.Sequence import iterate_fasta
from modules.Bed import parse_bed_regions
from modules.Authenticator import *
from modules.Paths import Paths
from modules.GsUri import *

from collections import defaultdict,Counter
from itertools import combinations
from multiprocessing import Pool
from copy import deepcopy,copy
from threading import Timer
from uuid import uuid4
import subprocess
import argparse
import os.path
import time
import json
import time
import math
import sys
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

import matplotlib
matplotlib.use('Agg')


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

    dpi = 400

    width = min(width, 4000/dpi)
    height = min(height, 4000/dpi)

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
        with_labels=True
    )

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
                edge_color="#007cbe"
            )

    ylim = list(a.get_ylim())
    ylim[0] -= 0.5
    a.set_ylim(ylim)

    pyplot.savefig(output_path, dpi=dpi)


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
        if self._timer is not None:
            self._timer.cancel()


class Variables:
    def __init__(self):
        # key = (path_id,read_id)
        self.path_to_read = dict()

        # key = (path_id,sample_id)
        self.path_to_sample = dict()

        # key = path_id
        self.path = dict()

        # Used for marginalization/iteration over n (to be deprecated)
        # key = integer of n >= 1
        self.n = dict()

        # Objective a, explicitly minimized
        self.cost_d = None

        # Objective b, marginalized
        self.cost_n = None

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
            cache.cost_d = None
            cache.cost_n = None

        else:
            for key in self.path_to_read:
                cache.path_to_read[key] = solver.Value(self.path_to_read[key])

            for key in self.path_to_sample:
                cache.path_to_sample[key] = solver.Value(self.path_to_sample[key])

            for key in self.path:
                cache.path[key] = solver.Value(self.path[key])

            for key in self.n:
                cache.n[key] = solver.Value(self.n[key])

            cache.cost_d = solver.Value(self.cost_d) if self.cost_d is not None else 0
            cache.cost_n = solver.Value(self.cost_n) if self.cost_n is not None else 0

        cache.response_stats = solver.ResponseStats()
        cache.status = status
        cache.is_cache = True

        return cache

    def infer_n(self, solver, path_to_read_costs, paths, reads):
        n = 0
        for path_id in paths.ids():
            v = [self.path_to_read[(path_id,read_id)] for read_id in reads.ids() if (path_id,read_id) in path_to_read_costs]
            n += sum(v) > 0

        if type(n) == ortools.sat.python.cp_model.IntVar:
            return solver.Value(n)

        elif type(n) == int:
            return n

        else:
            return None

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

    """
    There may be issues with the model which lead it to be infeasible, and CPSAT does not indicate why, so this method
    is a pre-optimizer check
    """
    def validate(
            self,
            sample_id_to_read_ids: dict,
            sample_id_map: IncrementalIdMap,
            path_to_read_costs: dict):

        # First get the reverse map so we can find the sample of each read
        read_id_to_sample_id = dict()
        for sample_id,read_ids in sample_id_to_read_ids.items():
            for read_id in read_ids:
                read_id_to_sample_id[read_id] = sample_id

        passing_samples = set()
        paths_covered_per_sample = defaultdict(set)

        # Iterate all path-read pairs and update the sample-level filter
        for (path_id,read_id) in path_to_read_costs:
            sample_id = read_id_to_sample_id[read_id]
            passing_samples.add(sample_id)
            paths_covered_per_sample[sample_id].add(path_id)

        all_pass = True
        for id,sample in sample_id_map:
            if id not in passing_samples:
                sys.stderr.write("ERROR: sample has no valid path assignments: " + sample)
                sys.stderr.write("\n")
                all_pass = False

        if not all_pass:
            sys.stderr.write("WARNING: Not all samples have valid reads-to-path alignments (spanning)")
            sys.stderr.write("\n")

    def get_genotype_support(
            self,
            paths: Paths,
            sample_id_to_read_ids: dict,
            sample_id_map: IncrementalIdMap,
            read_id_map: IncrementalIdMap):

        # Nested map of sample: path: reads
        # There will only ever be 2 paths per sample, but any number of reads can be assigned to each path,
        # so this datastructure shows the reads that were assigned to each path for each sample, for the purpose of
        # evaluating genotypes
        genotype_support = defaultdict(lambda: defaultdict(list))

        # First get the reverse map so we can find the sample of each read
        read_id_to_sample_id = dict()
        for sample_id,read_ids in sample_id_to_read_ids.items():
            for read_id in read_ids:
                read_id_to_sample_id[read_id] = sample_id

        # Iterate all path,read pairs and build the coverage_per_sample nested map
        for (path_id,read_id),value in self.path_to_read.items():
            if value == 1:
                sample_id = read_id_to_sample_id[read_id]

                sample_name = sample_id_map.get_name(sample_id)
                read_name = read_id_map.get_name(read_id)
                path_name = paths.get_path_name(path_id)

                genotype_support[sample_name][path_name].append(read_name)

        return genotype_support

    def write_genotype_support_to_file(
            self,
            output_path: str,
            paths: Paths,
            sample_id_to_read_ids: dict,
            sample_id_map: IncrementalIdMap,
            read_id_map: IncrementalIdMap):

        support = self.get_genotype_support(paths=paths, sample_id_to_read_ids=sample_id_to_read_ids, sample_id_map=sample_id_map, read_id_map=read_id_map)

        with open(output_path, 'w') as file:
            for sample in support:
                sample_id = sample_id_map.get_id(sample)
                for path in support[sample]:
                    path_id = paths.get_path_id(path)
                    for read in support[sample][path]:
                        read_id = read_id_map.get_id(read)
                        file.write(','.join(list(map(str,[sample_id,sample,path_id,path,read_id,read]))))
                        file.write('\n')

        return


def construct_base_model(
        vars,
        model,
        reads,
        paths,
        sample_to_reads,
        path_to_read_costs):

    # Define read assignment variables
    for edge in path_to_read_costs.keys():
        # 'edge' is a tuple with path_id,read_id
        vars.path_to_read[edge] = model.NewBoolVar("p%dr%d" % edge)

    # Constraint: each read must map to only one haplotype/path
    # Only consider pairs that have a cost assigned to them
    for read_id in reads.ids():
        v = [vars.path_to_read[(path_id,read_id)] for path_id in paths.ids() if (path_id,read_id) in path_to_read_costs]
        if len(v) == 0:
            # print("WARNING: skipping read id %d because has no viable path assignments" % read_id)
            continue

        model.Add(sum(v) == 1)

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

            v = [vars.path_to_read[(path_id,read_id)] for read_id in read_group if (path_id,read_id) in path_to_read_costs]
            s = sum(v)
            if len(v) == 0:
                # print("WARNING: skipping path id %d because has no viable path assignments" % path_id)
                continue

            model.Add(s >= 1).OnlyEnforceIf(vars.path_to_sample[edge])
            model.Add(s == 0).OnlyEnforceIf(vars.path_to_sample[edge].Not())

    # Now that the boolean indicators have been defined, use them to add a constraint on ploidy per sample
    for sample_id,read_group in sample_to_reads.items():
        model.Add(sum([vars.path_to_sample[(path_id,sample_id)] for path_id in paths.ids()]) <= 2)

    return vars, model


def find_all_pareto_solutions(
        path_to_read_costs,
        reads: IncrementalIdMap,
        samples: IncrementalIdMap,
        paths: Paths,
        sample_to_reads: dict,
        output_dir: str,
        n_threads: int = 1):

    output_subdir = os.path.join(output_dir, "optimizer")

    vars = Variables()
    model = cp_model.CpModel()

    vars, model = construct_base_model(
        vars=vars,
        model=model,
        reads=reads,
        paths=paths,
        sample_to_reads=sample_to_reads,
        path_to_read_costs=path_to_read_costs
    )

    # Add boolean indicators which imply that a path has been used (assigned any read)
    for path_id in paths.ids():
        vars.path[path_id] = model.NewBoolVar("p" + str(path_id))

        # Accumulate all possible assignments of this path to any reads
        v = [vars.path_to_read[(path_id,read_id)] for read_id in reads.ids() if (path_id,read_id) in path_to_read_costs]
        s = sum(v)

        if len(v) == 0:
            print("WARNING: skipping path id %d because has no viable path assignments" % path_id)
            continue

        model.Add(s >= 1).OnlyEnforceIf(vars.path[path_id])
        model.Add(s == 0).OnlyEnforceIf(vars.path[path_id].Not())

    # Cost term a: sum of edit distances for all reads assigned to haplotypes
    vars.cost_d = model.NewIntVar(0, 10_000_000_000, "cost_d")
    model.Add(vars.cost_d == sum([c * vars.path_to_read[e] for e,c in path_to_read_costs.items()]))

    # Cost term b: sum of unique haplotypes used
    vars.cost_n = model.NewIntVar(0, 10_000_000_000, "cost_n")
    model.Add(vars.cost_n == sum([x for x in vars.path.values()]))

    # n diploid samples can at most fill n*2 haplotypes
    # Sometimes there are may be fewer candidate paths than that
    max_feasible_haplotypes = min(len(paths), len(sample_to_reads)*2) + 1

    vars.validate(sample_id_to_read_ids=sample_to_reads, sample_id_map=samples, path_to_read_costs=path_to_read_costs)

    for i in range(1,max_feasible_haplotypes):
        vars.n[i] = model.NewBoolVar("n" + str(i))
        model.Add(vars.cost_n == i).OnlyEnforceIf(vars.n[i])
        model.Add(vars.cost_n != i).OnlyEnforceIf(vars.n[i].Not())

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = n_threads
    # solver.parameters.max_time_in_seconds = 300.0
    # solver.parameters.log_search_progress = True
    # solver.log_callback = print

    status = None
    results = dict()
    i_prev = None

    for i in range(1,max_feasible_haplotypes):
        print(i)

        model.ClearAssumptions()
        model.AddAssumption(vars.n[i])
        model.Minimize(vars.cost_d)

        if status == cp_model.OPTIMAL:
            model.ClearHints()
            for var in vars.path_to_read.values():
                model.AddHint(var, solver.Value(var))

        # It's critical to STOP THE TIMER after this finishes because it relies on a thread which will run on
        # past the solution for as long as the timer is set, potentially starving future iterations of threads.
        o = ObjectiveEarlyStopping(60)
        status = solver.SolveWithSolutionCallback(model, o)

        print("----")
        print(solver.SolutionInfo())
        print(solver.ResponseStats())
        sys.stdout.flush()

        if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
            results[i] = vars.get_cache(status=status, solver=solver)

        o.cancel_timer_thread()

        if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
            if i_prev is not None and results[i].cost_d > results[i_prev].cost_d:
                sys.stderr.write("Iteration stopped at n=%d because score worsened\n" % i)
                break

            i_prev = i

    for i,cache in results.items():
        print("%d,%d" % (i, cache.cost_d))
        cache.write_results_to_file(output_dir=os.path.join(output_subdir,str(i)))

    return results


def optimize_n(
        path_to_read_costs,
        reads: IncrementalIdMap,
        samples: IncrementalIdMap,
        paths: Paths,
        sample_to_reads: dict,
        output_dir: str,
        n_threads: int = 1):

    output_subdir = os.path.join(output_dir, "optimizer")

    # First optimize for unconstrained cost first, and fetch the corresponding n value
    vars = Variables()
    model = cp_model.CpModel()

    vars, model = construct_base_model(
        vars=vars,
        model=model,
        reads=reads,
        paths=paths,
        sample_to_reads=sample_to_reads,
        path_to_read_costs=path_to_read_costs
    )

    # Add boolean indicators which imply that a path has been used (assigned any read)
    for path_id in paths.ids():
        vars.path[path_id] = model.NewBoolVar("p" + str(path_id))

        # Accumulate all possible assignments of this path to any reads
        v = [vars.path_to_read[(path_id,read_id)] for read_id in reads.ids() if (path_id,read_id) in path_to_read_costs]
        s = sum(v)

        if len(v) == 0:
            print("WARNING: skipping path id %d because has no viable path assignments" % path_id)
            continue

        model.Add(s >= 1).OnlyEnforceIf(vars.path[path_id])
        model.Add(s == 0).OnlyEnforceIf(vars.path[path_id].Not())

    # Cost term n: sum of unique haplotypes used
    vars.cost_n = model.NewIntVar(0, 10_000_000_000, "cost_n")
    model.Add(vars.cost_n == sum([x for x in vars.path.values()]))

    vars.validate(sample_id_to_read_ids=sample_to_reads, sample_id_map=samples, path_to_read_costs=path_to_read_costs)

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = n_threads

    status = None

    model.Minimize(vars.cost_n)

    # It's critical to STOP THE TIMER after this finishes because it relies on a thread which will run on
    # past the solution for as long as the timer is set, potentially starving future iterations of threads.
    o = ObjectiveEarlyStopping(60)
    status = solver.SolveWithSolutionCallback(model, o)

    print("----")
    print(solver.SolutionInfo())
    print(solver.ResponseStats())
    sys.stdout.flush()

    solution = None
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        solution = vars.get_cache(status=status, solver=solver)
        solution.write_results_to_file(output_dir=os.path.join(output_subdir,"solution"))

    o.cancel_timer_thread()

    return solution


def optimize_d(
        path_to_read_costs,
        reads: IncrementalIdMap,
        samples: IncrementalIdMap,
        paths: Paths,
        sample_to_reads: dict,
        output_dir: str,
        n_threads: int = 1):

    output_subdir = os.path.join(output_dir, "optimizer")

    # First optimize for unconstrained cost first, and fetch the corresponding n value
    vars = Variables()
    model = cp_model.CpModel()

    vars, model = construct_base_model(
        vars=vars,
        model=model,
        reads=reads,
        paths=paths,
        sample_to_reads=sample_to_reads,
        path_to_read_costs=path_to_read_costs
    )

    # # Add boolean indicators which imply that a path has been used (assigned any read)
    # # ONLY used for tracking purposes later..TODO remove this
    # for path_id in paths.ids():
    #     vars.path[path_id] = model.NewBoolVar("p" + str(path_id))
    #
    #     # Accumulate all possible assignments of this path to any reads
    #     v = [vars.path_to_read[(path_id,read_id)] for read_id in reads.ids() if (path_id,read_id) in path_to_read_costs]
    #     s = sum(v)
    #
    #     if len(v) == 0:
    #         print("WARNING: skipping path id %d because has no viable path assignments" % path_id)
    #         continue
    #
    #     model.Add(s >= 1).OnlyEnforceIf(vars.path[path_id])
    #     model.Add(s == 0).OnlyEnforceIf(vars.path[path_id].Not())

    # Cost term d: sum of edit distances for all reads assigned to haplotypes
    vars.cost_d = model.NewIntVar(0, 10_000_000_000, "cost_d")
    model.Add(vars.cost_d == sum([c * vars.path_to_read[e] for e,c in path_to_read_costs.items()]))

    vars.validate(sample_id_to_read_ids=sample_to_reads, sample_id_map=samples, path_to_read_costs=path_to_read_costs)

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = n_threads

    status = None

    model.Minimize(vars.cost_d)

    # It's critical to STOP THE TIMER after this finishes because it relies on a thread which will run on
    # past the solution for as long as the timer is set, potentially starving future iterations of threads.
    o = ObjectiveEarlyStopping(60)
    status = solver.SolveWithSolutionCallback(model, o)

    print("----")
    print(solver.SolutionInfo())
    print(solver.ResponseStats())
    sys.stdout.flush()

    solution = None
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        solution = vars.get_cache(status=status, solver=solver)
        solution.write_results_to_file(output_dir=os.path.join(output_subdir,"solution"))

    o.cancel_timer_thread()

    # Instead of explicitly modeling n, infer it at the end to save complexity during optimization
    solution.cost_n = solution.infer_n(solver=solver, path_to_read_costs=path_to_read_costs, paths=paths, reads=reads)

    return solution


def optimize_d_plus_n(
        path_to_read_costs,
        reads: IncrementalIdMap,
        samples: IncrementalIdMap,
        d_coeff: float,
        n_coeff: float,
        paths: Paths,
        sample_to_reads: dict,
        output_dir: str,
        n_threads: int = 1):

    output_subdir = os.path.join(output_dir, "optimizer")

    # First optimize for unconstrained cost first, and fetch the corresponding n and d value
    vars = Variables()
    model = cp_model.CpModel()

    vars, model = construct_base_model(
        vars=vars,
        model=model,
        reads=reads,
        paths=paths,
        sample_to_reads=sample_to_reads,
        path_to_read_costs=path_to_read_costs
    )

    # Add boolean indicators which imply that a path has been used (assigned any read)
    for path_id in paths.ids():
        vars.path[path_id] = model.NewBoolVar("p" + str(path_id))

        # Accumulate all possible assignments of this path to any reads
        v = [vars.path_to_read[(path_id,read_id)] for read_id in reads.ids() if (path_id,read_id) in path_to_read_costs]
        s = sum(v)

        if len(v) == 0:
            print("WARNING: skipping path id %d because has no viable path assignments" % path_id)
            continue

        model.Add(s >= 1).OnlyEnforceIf(vars.path[path_id])
        model.Add(s == 0).OnlyEnforceIf(vars.path[path_id].Not())

    # Cost term a: sum of edit distances for all reads assigned to haplotypes
    vars.cost_d = model.NewIntVar(0, 10_000_000_000, "cost_d")
    model.Add(vars.cost_d == sum([c * vars.path_to_read[e] for e,c in path_to_read_costs.items()]))

    # Cost term b: sum of unique haplotypes used
    vars.cost_n = model.NewIntVar(0, 10_000_000_000, "cost_n")
    model.Add(vars.cost_n == sum([x for x in vars.path.values()]))

    vars.validate(sample_id_to_read_ids=sample_to_reads, sample_id_map=samples, path_to_read_costs=path_to_read_costs)

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = n_threads

    status = None

    model.Minimize(d_coeff*vars.cost_d + n_coeff*vars.cost_n)

    # It's critical to STOP THE TIMER after this finishes because it relies on a thread which will run on
    # past the solution for as long as the timer is set, potentially starving future iterations of threads.
    o = ObjectiveEarlyStopping(60)
    status = solver.SolveWithSolutionCallback(model, o)

    print("----")
    print(solver.SolutionInfo())
    print(solver.ResponseStats())
    sys.stdout.flush()

    solution = None
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        solution = vars.get_cache(status=status, solver=solver)
        solution.write_results_to_file(output_dir=os.path.join(output_subdir,"solution"))

    o.cancel_timer_thread()

    return solution


def optimize_with_d_norm(
        path_to_read_costs,
        reads: IncrementalIdMap,
        samples: IncrementalIdMap,
        paths: Paths,
        sample_to_reads: dict,
        output_dir: str,
        n_threads: int = 1):

    # Initially we only want n to be a tie-breaker between equally optimal min(d) solutions
    n_coeff = 1
    d_coeff = len(paths)*10

    min_d_solution = optimize_d_plus_n(
        path_to_read_costs=path_to_read_costs,
        reads=reads,
        samples=samples,
        paths=paths,
        d_coeff=d_coeff,
        n_coeff=n_coeff,
        sample_to_reads=sample_to_reads,
        output_dir=output_dir,
        n_threads=n_threads
    )

    # Multiply each factor in the cost function by the other factor's normalization constant to attempt to equalize them
    # Then add an arbitrary additional term to weight the solution
    n_coeff = min_d_solution.cost_d
    d_coeff = min_d_solution.cost_n*6

    print("\n---")
    print("Using n=%d and d=%d" % (min_d_solution.cost_n, min_d_solution.cost_d))
    print("---\n")

    solution = optimize_d_plus_n(
        path_to_read_costs=path_to_read_costs,
        reads=reads,
        samples=samples,
        paths=paths,
        d_coeff=d_coeff,
        n_coeff=n_coeff,
        sample_to_reads=sample_to_reads,
        output_dir=output_dir,
        n_threads=n_threads
    )

    return solution


def get_multiobjective_best_solution(results, n_stop, region_string, output_directory):
    if len(results) == 0:
        sys.stderr.write("WARNING: Region failed due to infeasibility of optimization: %s\n" % region_string)
        return "Region failed due to infeasibility of optimization\n"

    # TODO: replace this with average error rate per read?
    c_target = 0

    # TODO: replace this with regional heterozygosity estimate
    n_target = 0

    fig,axes = pyplot.subplots(nrows=1, ncols=1)

    c_start = sys.maxsize
    c_stop = 0
    n_start = 1
    # n_stop = len(data_per_sample.keys())

    for n,cache in results.items():
        if cache.cost_d < c_start:
            c_start = cache.cost_d

        if cache.cost_d > c_stop:
            c_stop = cache.cost_d

    print("----")
    print("n",n_start,n_stop)
    print("c",c_start,c_stop)
    print()

    d_min = sys.maxsize
    c_min = None
    n_min = None
    c_per_read_min = None
    n_norm_min = None
    c_norm_min = None

    for n,cache in results.items():
        # Normalize the range of outputs for n and c if a range exists
        # A constant of 1 is added to make the c_norm more impactful on total distance from utopia point
        if c_stop == c_start:
            c_norm = c_start
        else:
            c_norm = float(cache.cost_d - c_start) / float(c_stop - c_start) + 1

        if n_stop == n_start:
            n_norm = n_start
        else:
            n_norm = (float(n) - n_start) / float(n_stop - n_start) + 0.001

        distance = math.sqrt((c_norm-c_target)**2 + (n_norm-n_target)**2)
        if distance < d_min:
            d_min = distance
            n_min = n
            c_min = cache.cost_d

            # For plotting
            n_norm_min = n_norm
            c_norm_min = c_norm

        axes.scatter(n_norm, c_norm, color="C0")
        axes.plot([n_target,n_norm], [c_target,c_norm], color="gray")
        axes.text(n_norm, c_norm, "%.4f" % distance)

    best_result = results[n_min]

    axes.plot(n_norm_min, c_norm_min, color="C1", marker='o', markerfacecolor='none', markersize=6)
    axes.text(n_norm_min, c_norm_min, "n=%d" % n_min, color="red", ha='right', va='top')

    axes.set_title("Results for %s" % region_string)
    axes.set_xlabel("n_norm")
    axes.set_ylabel("c_norm")

    optimizer_results_path = os.path.join(output_directory, "optimizer_results.png")
    pyplot.savefig(optimizer_results_path, dpi=200)
    pyplot.close()

    return n_min, best_result


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
                continue

            # Column 5 (0-based) is the path column in a GAF file
            # It can be a forward or reverse alignment, so we will reinterpret them all as forward alignments
            path = None
            if tokens[5][0] == '>':
                path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            elif tokens[5][0] == '<':
                path = tuple(map(int,reversed(re.findall(r'\d+', tokens[5]))))

            else:
                raise Exception("ERROR: Non GFA character found in path column of GFA")

            if len(path) > 1 and alleles[path[0]].is_left_flank and alleles[path[-1]].is_right_flank:
                paths.increment_weight(path,1)

    filtered_paths = Paths()

    for p,path,frequency in paths:
        if frequency >= min_coverage:
            filtered_paths.add_path(path,frequency)

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
        sys.stderr.write("ERROR: no passing paths in alignment graph\n")
        return None, None

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


def trim_flanks_from_ref_alleles(alleles, flank_length):
    # Get start and stop nodes
    start_node = None
    stop_node = None

    for i in range(len(alleles)):
        if alleles[i].is_left_flank:
            start_node = i

        if alleles[i].is_right_flank:
            stop_node = i

    # Trim flanks off the allele sequences
    s_start = alleles[start_node].sequence
    l_start = len(s_start)
    s_start = s_start[flank_length:]

    s_stop = alleles[stop_node].sequence
    l_stop = len(s_stop)
    s_stop = s_stop[:(l_stop - flank_length)]

    alleles[start_node].sequence = s_start
    alleles[stop_node].sequence = s_stop


'''
Very dumb (expensive) way to test a sample's read-to-path edge set for feasibility under the diploid constraint 
'''
def is_feasible(edges):
    all_reads = {read_id for path_id,read_id in edges}
    all_paths = {path_id for path_id,read_id in edges}

    # Start with the easier solution: a single path covers all reads
    feasible = False
    for p in all_paths:
        reads = copy(all_reads)

        for path_id,read_id in edges:
            if path_id == p:
                if read_id in reads:
                    reads.remove(read_id)

        if len(reads) == 0:
            feasible = True
            break

    # If the easy solution works, exit early
    if feasible:
        return True

    # For every possible 2-path solution and check whether all reads are covered
    for a,b in combinations(all_paths,2):
        reads = copy(all_reads)

        for path_id,read_id in edges:
            if path_id == a or path_id == b:
                if read_id in reads:
                    reads.remove(read_id)

        if len(reads) == 0:
            feasible = True
            break

    return feasible


def filter_costs(path_to_read_costs, read_id_map, read_to_sample, max_path_to_read_cost):
    samplewise_costs = defaultdict(list)
    for [path_id,read_id],cost in path_to_read_costs.items():
        read_name = read_id_map.get_name(read_id)
        sample_name = read_to_sample[read_name]
        samplewise_costs[sample_name].append((path_id,read_id,cost))

    for sample_name,costs in samplewise_costs.items():
        edges = {(path_id,read_id) for path_id,read_id,cost in costs}
        edges_to_be_removed = list()

        # Remove edges from this sample in order of greatest until:
        # the cost is below the threshold OR it causes infeasibility under the diploid constraint
        for path_id,read_id,cost in sorted(costs, key=lambda x: x[-1], reverse=True):
            if cost < max_path_to_read_cost:
                break

            edges.remove((path_id,read_id))
            if not is_feasible(edges):
                edges.add((path_id,read_id))
                break

            edges_to_be_removed.append((path_id,read_id))

        for e in edges_to_be_removed:
            del path_to_read_costs[e]

    return path_to_read_costs


def write_costs_to_file(output_directory, path_to_read_costs, paths, read_id_map, read_to_sample, sample_id_map):
    costs_output_path = os.path.join(output_directory, "costs.csv")
    with open(costs_output_path, 'w') as file:
        file.write("path_id,read_id,sample_id,path_name,read_name,sample_name,value")
        file.write('\n')
        for key,value in path_to_read_costs.items():
            path_id = str(key[0])
            path_name = paths.get_path_name(key[0])
            read_id = str(key[1])
            read_name = read_id_map.get_name(key[1])

            sample_name = read_to_sample[read_name]
            sample_id = str(sample_id_map.get_id(sample_name))

            file.write(','.join([path_id,read_id,sample_id,path_name,read_name,sample_name,str(value)]))
            file.write('\n')


def get_read_costs_from_alignments(
        bam_paths,
        read_id_map,
        spanning_reads,
        paths,
        flank_length,
        interval_pad_length,
        path_to_read_costs):

    # For now, double counting softclips is allowed, because we don't really care,
    # the alignment to the true haplotype should contain no clips
    non_match_ops = {'X','I','D','S'}

    min_path_edit_distance = float("inf")

    # Because of the potential for supplementary alignments, need to first aggregate costs per read
    for bam_path in bam_paths:
        sys.stderr.write(bam_path)
        sys.stderr.write('\n')

        bam = pysam.AlignmentFile(bam_path, 'rb')

        total = 0.0
        n = 0.0
        for alignment in bam:
            read_name = alignment.query_name
            read_id = read_id_map.get_id(read_name)

            # TODO: [optimize] eventually alignment should be done in memory, and this filtering step could be performed
            # before aligning the reads with minimap2, which would save (potentially considerable) time
            if read_name not in spanning_reads:
                continue

            # TODO: [fix] verify that the NM tag is not double-counting softclips/hardclips in a supplementary alignment?
            # not sure how much this affects results, if there is a dup, it should be represented in the path already
            path_name = bam.header.get_reference_name(alignment.reference_id)
            path_id = paths.get_path_id(path_name)
            path_length = bam.get_reference_length(path_name)

            window_start = flank_length
            window_stop = path_length - flank_length + 1

            # Sometimes unpredictability in the alignment pushes a large indel just beyond the window...
            # this is added as a margin of safety to catch those. Since the windows are built requiring an empty
            # padding between them, we can use that to decide how much space is safe to look for relevant edits
            window_start -= interval_pad_length
            window_stop += interval_pad_length

            alignment_start = alignment.reference_start
            alignment_stop = alignment.reference_end

            # Check that read is spanning this haplotype (without any softclips)
            # If it's not spanning, add a fixed length penalty to break tie cases between haplotypes
            # TODO: fix this in the region selection process, not by having a constant
            if not alignment_start < window_start < window_stop < alignment_stop:
                path_to_read_costs[(path_id,read_id)] += 10

            if not alignment.is_secondary:
                for ref_start,ref_stop,query_start,query_stop,operation,length in iterate_cigar(alignment):
                    c = cigar_index_to_char[operation]
                    is_non_match = (c in non_match_ops)

                    # This iterator might not give ordered ref start/stop if the read is reversed (TODO: make a ref-oriented version?)
                    if ref_stop < ref_start:
                        r = ref_stop
                        ref_stop = ref_start
                        ref_start = r

                    # Cigar is contained entirely in the window
                    if window_start < ref_start <= ref_stop < window_stop:
                        cost = length * int(is_non_match)
                        path_to_read_costs[(path_id,read_id)] += cost
                        total += cost

                    # Cigar covers entire window
                    elif ref_start <= window_start < window_stop <= ref_stop:
                        cost = (window_stop - window_start) * int(is_non_match)
                        path_to_read_costs[(path_id,read_id)] += cost
                        total += cost

                    # Cigar is overlapping the start of the window
                    elif ref_start <= window_start <= ref_stop:
                        if is_ref_move[operation]:
                            cost = (ref_stop - window_start) * int(is_non_match)
                        else:
                            cost = length * int(is_non_match)

                        path_to_read_costs[(path_id,read_id)] += cost
                        total += cost

                    # Cigar is overlapping the end of the window
                    elif ref_start <= window_stop <= ref_stop:
                        if is_ref_move[operation]:
                            cost = (window_stop - ref_start) * int(is_non_match)
                        else:
                            cost = length * int(is_non_match)

                        path_to_read_costs[(path_id,read_id)] += cost
                        total += cost

                n += 1

        if total < min_path_edit_distance:
            min_path_edit_distance = total

        avg_edit_distance = total/n

    return path_to_read_costs


def write_path_assignments_to_fasta(solution, sample_id_map, paths, alleles, output_directory):
    counter = defaultdict(int)

    assigned_haplotypes_path = os.path.join(output_directory, "assigned_haplotypes.fasta")
    with open(assigned_haplotypes_path, 'w') as file:
        for [path_id,sample_id],value in solution.path_to_sample.items():
            # Only consider combinations of path/sample that were actually assigned by the optimizer
            if value == 0:
                continue

            counter[sample_id] += 1
            sample_name = sample_id_map.get_name(sample_id)
            sequence_name = sample_name + "_" + str(counter[sample_id])

            file.write(">" + sequence_name + " id:" + str(path_id) + " length:"+str(sum(len(alleles[i].sequence) for i in paths.get_path(path_id))))
            file.write("\n")
            for i in paths.get_path(path_id):
                file.write(alleles[i].sequence)

            file.write("\n")


def write_full_haplotypes_to_vcf(chromosome, ref_start, ref_stop, ref_path, solution, sample_id_map, paths, alleles, output_directory):
    output_path = os.path.join(output_directory, "assigned_haplotypes.vcf")

    path_ids_per_sample = defaultdict(list)
    id_to_sequence_index = dict()
    sequences = list()

    chromosome_sequence = FastaFile(ref_path).fetch(chromosome)
    ref_sequence = chromosome_sequence[ref_start:ref_stop]

    for [path_id,sample_id],value in solution.path_to_sample.items():
        # Only consider combinations of path/sample that were actually assigned by the optimizer
        if value == 0:
            continue

        sample_name = sample_id_map.get_name(sample_id)

        # Sample --> path_id
        path_ids_per_sample[sample_name].append(path_id)

        if path_id not in id_to_sequence_index:

            sequence = ""
            for i in paths.get_path(path_id):
                sequence += alleles[i].sequence

            if sequence != ref_sequence:
                # path_id --> sequence_index
                id_to_sequence_index[path_id] = len(sequences)
                sequences.append(sequence)
                print("ADDING path to alt alleles. path=%d alt=%d length=%d" % (path_id, id_to_sequence_index[path_id], len(sequence)))
            else:
                print("SKIPPING sequence identical to reference with path ID: " + str(path_id) + " length " + str(len(sequence)))

    with open(output_path, 'w') as file:
        file.write("##fileformat=VCFv4.2\n")
        file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"The genotype of the variant\">\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(sorted(path_ids_per_sample.keys())))

        gts = list()
        for sample_name,path_ids in sorted(path_ids_per_sample.items(), key=lambda x: x[0]):
            sample_indexes = [id_to_sequence_index[x]+1 for x in path_ids if x in id_to_sequence_index]

            if len(sample_indexes) == 0:
                gts.append("0/0")
            elif len(sample_indexes) == 1:
                gts.append("0/" + str(sample_indexes[0]))
            elif len(sample_indexes) == 2:
                gts.append(str(sample_indexes[0]) + "/" + str(sample_indexes[1]))
            else:
                raise Exception("ERROR: too many alleles specified by sample: " + sample_name)

        prefix = chromosome_sequence[ref_start-1]

        file.write('\n')
        file.write(chromosome)                      # CHROM
        file.write('\t')
        file.write(str(ref_start))                  # POS
        file.write('\t')
        file.write("hapslap_" + uuid4().hex)        # ID
        file.write('\t')
        file.write(ref_sequence)                    # REF
        file.write('\t')
        file.write(','.join([prefix + s for s in sequences]))             # ALT
        file.write('\t')
        file.write("60")                            # QUAL
        file.write('\t')
        file.write("PASS")                          # FILTER
        file.write('\t')
        file.write(".")                             # INFO
        file.write('\t')
        file.write("GT")                            # FORMAT
        file.write('\t')
        file.write('\t'.join(gts))                  # all SAMPLEs

    compress_and_index_vcf(output_path)


def write_full_candidate_haplotypes_to_vcf(chromosome, ref_start, ref_stop, ref_path, paths, alleles, output_directory):
    output_path = os.path.join(output_directory, "candidate_haplotypes.vcf")

    path_ids_per_sample = defaultdict(list)
    id_to_sequence_index = dict()
    sequences = list()

    chromosome_sequence = FastaFile(ref_path).fetch(chromosome)
    ref_sequence = chromosome_sequence[ref_start:ref_stop]

    for path_id in paths.ids():
        sample_name = path_id

        # Sample --> path_id
        path_ids_per_sample[sample_name].append(path_id)

        if path_id not in id_to_sequence_index:
            sequence = ""
            for i in paths.get_path(path_id):
                sequence += alleles[i].sequence

            if sequence != ref_sequence:
                # path_id --> sequence_index
                id_to_sequence_index[path_id] = len(sequences)
                sequences.append(sequence)
            else:
                print("SKIPPING sequence identical to reference with path ID: " + str(path_id))

    with open(output_path, 'w') as file:
        file.write("##fileformat=VCFv4.2\n")
        file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"The genotype of the variant\">\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(list(map(str,sorted(path_ids_per_sample.keys())))))

        gts = list()
        for sample_name,path_ids in sorted(path_ids_per_sample.items(), key=lambda x: x[0]):
            sample_indexes = [id_to_sequence_index[x]+1 for x in path_ids if x in id_to_sequence_index]

            if len(sample_indexes) == 0:
                gts.append("0/0")
            elif len(sample_indexes) == 1:
                gts.append("0/" + str(sample_indexes[0]))
            else:
                raise Exception("ERROR: incorrect number of alleles when iterating paths, ID: " + str(sample_name))

        prefix = chromosome_sequence[ref_start-1]

        file.write('\n')
        file.write(chromosome)                      # CHROM
        file.write('\t')
        file.write(str(ref_start))                  # POS
        file.write('\t')
        file.write("hapslap_" + uuid4().hex)        # ID
        file.write('\t')
        file.write(ref_sequence)                             # REF
        file.write('\t')
        file.write(','.join([prefix + s for s in sequences]))             # ALT
        file.write('\t')
        file.write("60")                            # QUAL
        file.write('\t')
        file.write("PASS")                          # FILTER
        file.write('\t')
        file.write(".")                             # INFO
        file.write('\t')
        file.write("GT")                            # FORMAT
        file.write('\t')
        file.write('\t'.join(gts))                  # all SAMPLEs

    compress_and_index_vcf(output_path)


def write_path_assignments_to_vcf(solution, sample_id_map, paths, alleles, edge_to_allele_index, output_directory, n_threads):
    output_subdirectory = os.path.join(output_directory, "vcf")
    if not os.path.exists(output_subdirectory):
        os.makedirs(output_subdirectory)

    paths_per_sample = defaultdict(set)
    for [path_id,sample_id],value in solution.path_to_sample.items():
        # Only consider combinations of path/sample that were actually assigned by the optimizer
        if value == 0:
            continue

        sample_name = sample_id_map.get_name(sample_id)
        path = paths.get_path(path_id)

        paths_per_sample[sample_name].add(path)

    args = list()

    for sample,paths in paths_per_sample.items():
        output_path = os.path.join(output_subdirectory, sample + ".vcf")
        # alleles, paths, output_path, sample_name, edge_to_deletion_index, compress_and_index
        args.append((alleles, paths, output_path, sample, edge_to_allele_index, True))

    with Pool(n_threads) as pool:
        results = pool.starmap(write_paths_to_vcf, args)

    return output_subdirectory


def infer_haplotypes(
        ref_path,
        data_per_sample,
        chromosome,
        ref_start,
        ref_stop,
        interval_pad_length,
        flank_length,
        min_coverage,
        max_path_to_read_cost,
        output_directory,
        n_threads,
        parameter_size_cutoff):

    time_start = time.time()

    generate_debug_alignments = True
    n_sort_threads = 4

    region_string = chromosome + ":" + str(ref_start) + "-" + str(ref_stop)

    output_directory = os.path.join(output_directory, region_string.replace(':',"_"))

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        exit("ERROR: output directory exists already: %s" % output_directory)

    json_log_path = os.path.join(output_directory, "config.json")
    write_config_to_json(
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        sample_names=list(data_per_sample.keys()),
        ref_path=ref_path,
        flank_length=flank_length,
        min_coverage=min_coverage,
        max_path_to_read_cost=max_path_to_read_cost,
        json_path=json_log_path
    )

    tokenator = GoogleToken()

    # TODO: make this a function of the actual reference version used
    ref_sample_name = "ref"

    graph,alleles = vcf_to_graph(
        ref_path=ref_path,
        data_per_sample=data_per_sample,
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        ref_sample_name=ref_sample_name,
        flank_length=flank_length,
        skip_incompatible=True
    )

    if len(alleles) == 1:
        if alleles[0].is_left_flank and alleles[0].is_right_flank:
            sys.stderr.write("WARNING: no alt alleles found. Skipping region: %s\n" % region_string)
            return "No alt alleles in region\n"
        else:
            raise Exception("ERROR: impossible configuration of alleles, len=1 and not flank")

    alleles_output_path = os.path.join(output_directory,"alleles.csv")
    with open(alleles_output_path, 'w') as file:
        file.write("id,start,stop,length,samples,is_left_flank,is_right_flank\n")
        for i in range(len(alleles)):
            file.write(str(i) + ',' + alleles[i].as_comma_separated_str() + '\n')

    ref_id_offset = None
    for i in range(len(alleles)):
        if alleles[i].samples == {ref_sample_name}:
            ref_id_offset = i
            break

    print("ref_id_offset:", ref_id_offset)

    csv_path = os.path.join(output_directory, "nodes.csv")
    write_node_csv(csv_path, ref_sample_name, alleles)

    output_gfa_path = os.path.join(output_directory, "graph.gfa")

    # Write the GFA
    write_graph_to_gfa(output_path=output_gfa_path, graph=graph, alleles=alleles)

    sample_color = "#0592BA"
    ref_color = "#6e6e6e"

    # Plot the graph
    plot_path = os.path.join(output_directory, "dag.png")
    plot_graph(
        graph=graph,
        ref_id_offset=ref_id_offset,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path
    )

    # Remove empty nodes
    graph, edge_to_allele_index = remove_empty_nodes_from_variant_graph(graph, alleles)

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
        output_path=plot_path
    )

    read_fasta_path = os.path.join(output_directory, "reads.fasta")

    # sample --> [read_name_0,read_name_1,...]
    sample_to_reads = defaultdict(list)
    read_id_map = IncrementalIdMap()

    # Extract all the relevant reads from the ref genome alignment
    args = list()
    sample_list = list()
    for sample_name,data in data_per_sample.items():
        bam_path = data["bam"]
        print("bam_path:", bam_path)

        sample_list.append(sample_name)
        args.append([output_directory,bam_path,region_string,tokenator])

    regional_bams = None
    with Pool(n_threads) as pool:
        regional_bams = pool.starmap(get_region_from_bam, args)

    for region_bam_path, sample_name in zip(regional_bams,sample_list):
        read_names = get_read_names_from_bam(region_bam_path)

        for read_name in read_names:
            read_id_map.add(read_name)

        sample_to_reads[sample_name] = read_names

        # This repeatedly appends one FASTA file
        get_reads_from_bam(
            output_path=read_fasta_path,
            bam_path=region_bam_path,
            token=tokenator
        )

        os.remove(region_bam_path)

    del sample_list

    a = time.time()

    # # Align the relevant reads to the variant graph
    # output_gaf_path = run_minigraph(
    #     output_directory=output_directory,
    #     gfa_path=output_gfa_path,
    #     fasta_path=read_fasta_path,
    #     n_threads=n_threads
    # )

    # args = [
    #     "-c",
    #     "-g", str(20000),
    #     "-k", str(14),
    #     "-r", "10000,20000",
    #     "-n", "3,3",
    #     "-x", "lr",
    # ]

    # Align the relevant haplotypes to the variant graph
    output_gaf_path = run_graphchainer(
        output_directory=output_directory,
        gfa_path=output_gfa_path,
        fasta_path=read_fasta_path,
        n_threads=n_threads,
        args_override=args
    )

    b = time.time()
    time_elapsed_minigraph = b - a

    # Accumulate path weights in the variant graph by parsing the GAF
    update_edge_weights_using_alignments(
        alleles=alleles,
        gaf_path=output_gaf_path,
        graph=graph
    )

    plot_path = os.path.join(output_directory, "dag_aligned.png")
    plot_graph(
        graph=graph,
        ref_id_offset=ref_id_offset,
        ref_color=ref_color,
        sample_color=sample_color,
        output_path=plot_path,
        line_style=':',
        draw_edge_weight_overlay=True
    )

    path_subdirectory = os.path.join(output_directory, "paths")
    os.makedirs(path_subdirectory)

    paths = Paths()

    # Only retain the reads that spanned from source to sink in the graph
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

    if paths is None:
        return "No viable paths in region\n"

    sys.stderr.write("Optimizing using %d paths with min coverage %d\n" % (len(paths), min_coverage))

    a = time.time()

    # mashmap_output_path = run_mashmap(
    #     ref_fasta_paths=fasta_paths,
    #     reads_fasta_path=read_fasta_path,
    #     n_threads=n_threads,
    #     min_identity=0.5,
    #     output_directory=output_directory,
    #     filename_prefix="reads_vs_ref_mashmap"
    # )
    #
    # reads_per_path = defaultdict(set)
    # with open(mashmap_output_path, 'r') as file:
    #     for l,line in enumerate(file):
    #         tokens = line.strip().split()
    #
    #         query_name = tokens[0]
    #         path_name = tokens[5]
    #         path_id = paths.get_path_id(path_name)
    #         reads_per_path[path_id].add(query_name)
    #

    read_subset = set()
    sequences = list()
    for sequence in iterate_fasta(read_fasta_path):
        sequences.append(sequence)
        read_subset.add(sequence.name)

    bam_paths = list()
    # for path_id,read_subset in reads_per_path.items():
    for p,path in enumerate(fasta_paths):
        # path_path = os.path.join(path_subdirectory, str(path_id) + ".fasta")

        bam_path = run_minimap2_on_read_subset(
            ref_fasta_path=path,
            reads=sequences,
            read_names=read_subset,
            preset="map-hifi",
            n_threads=n_threads,
            output_directory=path_subdirectory,
            filename_prefix=str(p)
        )

        bam_paths.append(bam_path)

    del sequences

    b = time.time()

    time_elapsed_minimap = b - a

    # Optionally generate the alignment of paths to the ref, for debug/dev purposes
    if generate_debug_alignments:
        chromosome_fasta_path = os.path.join(output_directory, chromosome + ".fasta")
        with open(chromosome_fasta_path, 'w') as file:
            ref_sequence = FastaFile(ref_path).fetch(chromosome)
            file.write(">")
            file.write(chromosome)
            file.write('\n')
            file.write(ref_sequence)
            file.write('\n')

        combined_fasta_path = os.path.join(output_directory, "all_paths.fasta")

        with open(combined_fasta_path,'w') as out_file:
            for p,path in enumerate(fasta_paths):
                with open(path,'r') as file:
                    for line in file:
                        out_file.write(line)

        run_minimap2(
            ref_fasta_path=chromosome_fasta_path,
            reads_fasta_path=combined_fasta_path,
            preset="asm20",
            n_threads=n_threads,
            n_sort_threads=n_sort_threads,
            output_directory=path_subdirectory,
            filename_prefix="all_paths_vs_ref"
        )

        os.remove(chromosome_fasta_path)

    pyplot.close("all")

    fig = pyplot.figure()
    axes = pyplot.axes()

    histogram = IterativeHistogram(start=0, stop=1000, n_bins=200)

    path_to_read_costs = defaultdict(int)

    # Parse the exhaustive linear alignments to get the cost of assigning each read to each path
    path_to_read_costs = get_read_costs_from_alignments(
        bam_paths=bam_paths,
        read_id_map=read_id_map,
        spanning_reads=spanning_reads,
        paths=paths,
        flank_length=flank_length,
        interval_pad_length=interval_pad_length,
        path_to_read_costs=path_to_read_costs
    )

    reads_csv_path = os.path.join(output_directory, "reads.csv")
    read_id_map.write_to_file(reads_csv_path)

    sample_id_map = IncrementalIdMap()
    for name in data_per_sample.keys():
        sample_id_map.add(name)

    sample_id_path = os.path.join(output_directory, "samples.csv")
    with open(sample_id_path, 'w') as file:
        for id,name in sample_id_map:
            file.write("%d,%s\n" % (id,name))

    # Make a duplicate of the sample->read mapping, using ids instead of names
    sample_id_to_read_ids = dict()
    for sample_name,read_names in sample_to_reads.items():
        read_ids = [read_id_map.get_id(x) for x in read_names]

        sample_id = sample_id_map.get_id(sample_name)
        sample_id_to_read_ids[sample_id] = read_ids

    sample_to_read_output_path = os.path.join(output_directory, "sample_to_read.csv")
    with open(sample_to_read_output_path, 'w') as file:
        for sample_id,read_ids in sample_id_to_read_ids.items():
            for read_id in read_ids:
                file.write("%d,%d" % (sample_id,read_id))
                file.write("\n")

    read_to_sample = dict()
    for sample,reads in sample_to_reads.items():
        for read in reads:
            read_to_sample[read] = sample

    # TODO: do a better job of addressing situations where there is no diploid infeasibility, but all reads of a sample
    # are removed as a result of bad read-to-graph alignment. As an SV merger it is fair to exclude those samples
    # but should they be caught somehow and reported as empty?

    # Filter the edges from path to read that will be used in the model (for saving computational complexity)
    filter_costs(
        path_to_read_costs=path_to_read_costs,
        read_id_map=read_id_map,
        read_to_sample=read_to_sample,
        max_path_to_read_cost=max_path_to_read_cost
    )

    if len(path_to_read_costs) > parameter_size_cutoff:
        print("WARNING: skipping very large optimization problem: " + str(region_string))
        return "Region skipped because it has very large number of variables: " + str(len(path_to_read_costs)) + "\n"

    axes.plot(histogram.get_bin_centers(), histogram.get_histogram(), color="#007cbe")

    write_costs_to_file(
        output_directory=output_directory,
        path_to_read_costs=path_to_read_costs,
        paths=paths,
        read_id_map=read_id_map,
        read_to_sample=read_to_sample,
        sample_id_map=sample_id_map
    )

    a = time.time()

    solution = optimize_with_d_norm(
        path_to_read_costs=path_to_read_costs,
        reads=read_id_map,
        samples=sample_id_map,
        paths=paths,
        sample_to_reads=sample_id_to_read_ids,
        output_dir=output_directory,
        n_threads=n_threads
    )

    b = time.time()
    time_elapsed_optimizer = b - a

    vcf_subdirectory = write_path_assignments_to_vcf(
        solution=solution,
        sample_id_map=sample_id_map,
        paths=paths,
        alleles=alleles,
        edge_to_allele_index=edge_to_allele_index,
        output_directory=output_directory,
        n_threads=n_threads
    )

    merged_vcf_path = os.path.join(output_directory, region_string.replace(":","_") + ".vcf")
    merge_vcfs_in_directory(vcf_subdirectory, output_path=merged_vcf_path)

    # Before writing the final output sequences, trim the flanking sequence (to make evaluation simpler)
    trim_flanks_from_ref_alleles(alleles, flank_length)

    print("Writing haplotypes to VCF")
    write_full_haplotypes_to_vcf(
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        ref_path=ref_path,
        solution=solution,
        sample_id_map=sample_id_map,
        paths=paths,
        alleles=alleles,
        output_directory=output_directory
    )

    print("Writing candidate haplotypes to VCF")
    write_full_candidate_haplotypes_to_vcf(
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        ref_path=ref_path,
        paths=paths,
        alleles=alleles,
        output_directory=output_directory
    )

    write_path_assignments_to_fasta(
        solution=solution,
        sample_id_map=sample_id_map,
        paths=paths,
        alleles=alleles,
        output_directory=output_directory
    )

    genotype_support_output_path = os.path.join(output_directory, "genotype_support.csv")

    solution.write_genotype_support_to_file(
        output_path=genotype_support_output_path,
        paths=paths,
        sample_id_to_read_ids=sample_id_to_read_ids,
        sample_id_map=sample_id_map,
        read_id_map=read_id_map
    )

    time_stop = time.time()
    time_elapsed_total = time_stop - time_start

    summary_string = \
        "n,%d\n" % solution.cost_n + \
        "c,%d\n" % solution.cost_d + \
        "n_candidates,%d\n" % len(paths) + \
        "n_spanning_reads,%d\n" % len(spanning_reads) + \
        "time_elapsed_minigraph_s,%d\n" % time_elapsed_minigraph + \
        "time_elapsed_minimap_s,%d\n" % time_elapsed_minimap + \
        "time_elapsed_optimizer_s,%d\n" % time_elapsed_optimizer + \
        "time_elapsed_total_s,%d\n" % time_elapsed_total

    sys.stderr.write(summary_string)
    sys.stderr.write('\n')

    summary_path = os.path.join(output_directory, "summary.txt")
    with open(summary_path, 'w') as file:
        file.write(summary_string)

    return summary_string


def write_config_to_json(
        json_path,
        sample_names,
        ref_path,
        flank_length,
        min_coverage,
        max_path_to_read_cost,
        input_directory=None,
        output_directory=None,
        bed_path=None,
        chromosome=None,
        ref_start=None,
        ref_stop=None):

    config = dict()

    config["flank_length"] = flank_length
    config["min_coverage"] = min_coverage
    config["max_path_to_read_cost"] = max_path_to_read_cost
    config["sample_names"] = sample_names
    config["ref_path"] = ref_path

    if input_directory is not None:
        config["input_directory"] = input_directory

    if output_directory is not None:
        config["output_directory"] = output_directory

    if bed_path is not None:
        config["bed_path"] = bed_path

    if chromosome is not None:
        config["chromosome"] = chromosome

    if ref_start is not None:
        config["ref_start"] = ref_start

    if ref_stop is not None:
        config["ref_stop"] = ref_stop

    with open(json_path, 'w') as file:
        file.write(json.dumps(config, indent=4))
        file.write('\n')


def read_config_from_json(json_path):
    with open(json_path, 'r') as file:
        config = json.load(file)

    flank_length = config["flank_length"]
    min_coverage = config["min_coverage"]
    max_path_to_read_cost = config["max_path_to_read_cost"]
    bed_path = config["bed_path"]
    sample_names = config["sample_names"]
    ref_path = config["ref_path"]

    # TODO: refactor
    input_directory = config["input_directory"]
    output_directory = config["output_directory"]

    print(type(bed_path),bed_path)
    print(type(sample_names),sample_names)
    print(type(ref_path),ref_path)
    print(type(input_directory),input_directory)
    print(type(output_directory),output_directory)

    return flank_length, min_coverage, max_path_to_read_cost, bed_path, sample_names, ref_path, input_directory, output_directory


def main(json_path, max_parameters, n_threads):
    flank_length, min_coverage, max_path_to_read_cost, bed_path, sample_names, ref_path, input_directory, output_directory = read_config_from_json(json_path)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        exit("ERROR: output directory exists already: %s" % output_directory)

    json_log_path = os.path.join(output_directory, "config.json")
    write_config_to_json(
        bed_path=bed_path,
        sample_names=sample_names,
        ref_path=ref_path,
        flank_length=flank_length,
        min_coverage=min_coverage,
        max_path_to_read_cost=max_path_to_read_cost,
        input_directory=input_directory,
        output_directory=output_directory,
        json_path=json_log_path
    )

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

    regions = parse_bed_regions(bed_path)

    summary_strings = list()
    for region in regions:
        summary_string = infer_haplotypes(
            n_threads=n_threads,
            parameter_size_cutoff=max_parameters,
            interval_pad_length=150,
            ref_path=ref_path,
            data_per_sample=data_per_sample,
            chromosome=region.contig_name,
            ref_start=region.start,
            ref_stop=region.stop,
            flank_length=flank_length,
            min_coverage=min_coverage,
            max_path_to_read_cost=max_path_to_read_cost,
            output_directory=output_directory
        )

        summary_strings.append(summary_string)

    # Group all the results into one message at the end of the output
    for s in range(len(summary_strings)):
        sys.stderr.write(str(regions[s]))
        sys.stderr.write('\n')
        sys.stderr.write(summary_strings[s])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--json",
        required=True,
        type=str,
        help="Input json containing relevant variables to be used during run. Please see repository data/template.json for details"
    )

    parser.add_argument(
        "-t","--n_threads",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use (recommended 15-30)"
    )

    parser.add_argument(
        "-p","--max_parameters",
        required=False,
        default=sys.maxsize,
        type=int,
        help="Maximum size of optimization problem to attempt"
    )

    args = parser.parse_args()

    main(
        json_path=args.json,
        max_parameters=args.max_parameters,
        n_threads=args.n_threads
    )
