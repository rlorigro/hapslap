import multiprocessing

from modules.IterativeHistogram import IterativeHistogram
from modules.IncrementalIdMap import IncrementalIdMap
from modules.Cigar import get_haplotypes_of_region
from modules.Bam import download_regions_of_bam
from modules.Sequence import Sequence

from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from copy import copy, deepcopy
from pathlib import Path
from glob import glob
import argparse
import hashlib
import numpy
import json
import sys
import os
import re

from networkx.drawing.nx_pydot import graphviz_layout
from matplotlib import pyplot,colors
from pywfa import WavefrontAligner
import matplotlib
import networkx
import pandas
import pydot

# matplotlib.use('Agg')


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap


def iterate_fasta(path, force_upper_case=True, normalize_name=True):
    name = ""
    sequence = ""

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                if force_upper_case:
                    sequence = sequence.upper()

                # Output previous sequence
                if l > 0:
                    s = Sequence(name, sequence)

                    if len(s) != 0:
                        if normalize_name:
                            s.normalize_name()

                        yield s

                name = line.strip()[1:].split(' ')[0]
                sequence = ""

            else:
                sequence += line.strip()

    if force_upper_case:
        sequence = sequence.upper()

    s = Sequence(name, sequence)

    if len(s) != 0:
        if normalize_name:
            s.normalize_name()

        yield s


def get_indel_distance_from_string(alignment: dict, minimum_indel_length=1):
    total_indels = 0

    length_token = ""
    for c in alignment["cigars"]:
        if c.isnumeric():
            length_token += c
        else:
            l = int(length_token)
            if l >= minimum_indel_length and (c == 'I' or c == 'D'):
                total_indels += l

    return total_indels


def get_indel_distance_from_tuples(cigar_tuples):
    total_indels = 0

    for c,l in cigar_tuples:
        if c == 1 or c == 2:
            total_indels += l

    return total_indels


def align_and_get_indel_distance(a:Sequence,b:Sequence,output_dir=None):
    if len(a) == 0 or len(b) == 0:
        return max(len(a),len(b))

    aligner = WavefrontAligner(heuristic="adaptive")

    aligner(a.sequence,b.sequence)

    if aligner.status != 0:
        return max(len(a),len(b))

    d = get_indel_distance_from_tuples(aligner.cigartuples)

    if output_dir is not None:
        output_path = os.path.join(output_dir, a.name + "_" + b.name + ".txt")
        aligner.cigar_print_pretty(output_path)

    return d


def cross_align_sequences(sequences:dict, n_threads,output_dir=None):
    pairs = list()
    lengths = list()
    args = list()

    for a,b in combinations(sequences.values(),2):
        if a == b:
            continue

        lengths.append((len(a),len(b)))
        pairs.append((a.name,b.name))
        args.append((a,b,output_dir))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align_and_get_indel_distance, args)

    return pairs,lengths,results


def align_sequences_to_other_sequences(a_seqs:dict, b_seqs:dict, n_threads:int, output_dir:str=None):
    pairs = list()
    lengths = list()
    args = list()

    for a_name,a in a_seqs.items():
        for b_name,b in b_seqs.items():
            lengths.append((len(a),len(b)))
            pairs.append((a,b))
            args.append((a,b,output_dir))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align_and_get_indel_distance, args)

    return pairs,lengths,results


def orient_test_haps_by_best_match(ref_sequences, test_sequences, id_map):
    sample_names = set()
    for id,name in id_map:
        sample_name = name.split('_')[0]
        sample_names.add(sample_name)

    for sample_name in sample_names:
        a = sample_name + "_1"
        b = sample_name + "_2"

        if not (a in ref_sequences and b in ref_sequences and a in test_sequences and b in test_sequences):
            continue

        aa = align_and_get_indel_distance(ref_sequences[a],test_sequences[a])
        bb = align_and_get_indel_distance(ref_sequences[b],test_sequences[b])
        ab = align_and_get_indel_distance(ref_sequences[a],test_sequences[b])
        ba = align_and_get_indel_distance(ref_sequences[b],test_sequences[a])

        cis_distance = aa + bb
        trans_distance = ab + ba

        # print(sample_name)
        # print(aa,ab)
        # print(ba,bb)

        # Flip the sequences around if it is a better match
        if cis_distance > trans_distance:
            a_seq = test_sequences[a].sequence
            b_seq = test_sequences[b].sequence

            test_sequences[a] = Sequence(a,b_seq)
            test_sequences[b] = Sequence(b,a_seq)

    return test_sequences


def duplicate_homozygous_haps(sequences):
    sample_names = set()
    for name in sequences.keys():
        sample_name = name.split('_')[0]
        sample_names.add(sample_name)

    duplicated = set()
    for sample_name in sorted(sample_names):
        a = sample_name + "_1"
        b = sample_name + "_2"

        # If only one hap was generated by the test, just duplicate it (assuming the other one exists)
        if b not in sequences and a in sequences:
            sequences[b] = Sequence(name=b, sequence=sequences[a].sequence)
            duplicated.add(a)
            duplicated.add(b)

        if a not in sequences and b in sequences:
            sequences[a] = Sequence(name=a, sequence=sequences[b].sequence)
            duplicated.add(a)
            duplicated.add(b)

    return sequences, duplicated


def evaluate_upper_bound(ref_sequences, haplotype_fasta_path, output_dir, prefix, flank_size, distance_threshold, n_threads):
    test_sequences = dict()

    for item in iterate_fasta(haplotype_fasta_path, normalize_name=False):
        item.sequence = item.sequence[flank_size:len(item.sequence) - flank_size]
        test_sequences[item.name] = item

    ref_id_map = IncrementalIdMap()
    test_id_map = IncrementalIdMap()

    for key in sorted(ref_sequences):
        ref_id_map.add(key)

    # Order the test sequences by length
    for key in sorted(test_sequences, key=lambda x: len(test_sequences[x]), reverse=True):
        test_id_map.add(key)

    pairs,lengths,results = cross_align_sequences(sequences=ref_sequences,n_threads=n_threads) #,output_dir=output_dir)

    graph = networkx.Graph()

    for i in range(len(results)):
        a,b = pairs[i]
        l_a,l_b = lengths[i]
        d = results[i]

        # w = d
        w = 0
        if l_a + l_b > 0:
            w = float(l_a + l_b - d) / float(l_a + l_b)

        id_a = ref_id_map.get_id(a)
        id_b = ref_id_map.get_id(b)

        if id_a not in graph:
            graph.add_node(id_a)

        if id_b not in graph:
            graph.add_node(id_b)

        if d <= distance_threshold and w != 0:
            graph.add_edge(id_a, id_b, weight=w)

    ref_labels = dict()
    for id,name in ref_id_map:
        ref_labels[id] = name

    components = list(networkx.connected_components(graph))
    components = sorted(components, key=lambda x: len(x), reverse=True)

    component_map = dict()
    for c,component in enumerate(components):
        for n in component:
            component_map[n] = c

    # This will order all IDs so that clusters are contiguous, otherwise randomly ordered
    id_to_cluster_map = dict()
    for id in sorted(range(len(ref_id_map)), key=lambda x: component_map[x]):
        id_to_cluster_map[id] = len(id_to_cluster_map)

    pairs,lengths,results = align_sequences_to_other_sequences(
        a_seqs=ref_sequences,
        b_seqs=test_sequences,
        n_threads=n_threads
        # output_dir=output_dir
    )

    matrix = numpy.zeros([len(ref_id_map),len(test_id_map)])

    feasible_samples = set()

    for i in range(len(results)):
        name_ref = pairs[i][0].name
        name_test = pairs[i][1].name

        l_ref,l_test = lengths[i]
        d = results[i]

        id_ref = ref_id_map.get_id(name_ref)
        id_test = test_id_map.get_id(name_test)

        matrix[id_to_cluster_map[id_ref]][id_test] = d

        if d < distance_threshold:
            feasible_samples.add(name_ref)

    fig = pyplot.figure()
    axes = pyplot.axes()
    p = axes.matshow(matrix,cmap='viridis')
    pyplot.colorbar(p, ax=axes)

    # Construct tick labels based on the clustering order
    labels = [None]*len(ref_id_map.id_to_name)
    for id,name in ref_id_map:
        labels[id_to_cluster_map[id]] = name

    axes.xaxis.set_ticks_position('bottom')

    axes.set_yticks(list(range(len(ref_id_map))))
    axes.set_yticklabels(labels)

    for x in range(matrix.shape[0]):
        for y in range(matrix.shape[1]):
            if matrix[x][y] < distance_threshold:
                pyplot.plot(y,x,marker='o',color='red',markersize=0.7,markeredgecolor='none')

    axes.tick_params(axis='both', which='major', labelsize=3)

    pyplot.tight_layout()
    output_path = os.path.join(output_dir,prefix + "_upper_bound.png")
    fig.savefig(output_path, dpi=500, bbox_inches='tight')
    pyplot.close('all')

    return feasible_samples


def get_all_relevant_chromosome_names(test_dirs):
    names = set()

    for test_dir in test_dirs:
        config_path = os.path.join(test_dir,"config.json")

        if not os.path.exists(config_path):
            continue

        config = None
        with open(config_path, 'r') as file:
            config = json.load(file)

        names.add(config["chromosome"])

    return names


def construct_graph_from_alignments(pairs, lengths, results, id_map, graph, distance_threshold):
    for i in range(len(results)):
        a,b = pairs[i]
        l_a,l_b = lengths[i]
        d = results[i]

        w = 0
        if l_a + l_b > 0:
            w = float(l_a + l_b - d) / float(l_a + l_b)

        id_a = id_map.get_id(a)
        id_b = id_map.get_id(b)

        if id_a not in graph:
            graph.add_node(id_a)

        if id_b not in graph:
            graph.add_node(id_b)

        if d <= distance_threshold and w != 0:
            # print(a, b, d, w, l_a, l_b)
            graph.add_edge(id_a, id_b, weight=d)

    return graph


def evaluate_test_haplotypes(
        input_dir,
        cache_dir,
        output_dir,
        n_threads,
        distance_threshold,
        tsv_path,
        column_names
        ):

    test_dirs = [str(f) for f in Path(input_dir).iterdir() if f.is_dir()]

    chromosome_names = get_all_relevant_chromosome_names(test_dirs)

    print("Downloading relevant chromosomes: " + str(chromosome_names))

    aligned_hap_directory = os.path.join(cache_dir, "aligned_haps")
    bam_paths = download_regions_of_bam(
        regions=chromosome_names,
        tsv_path=tsv_path,
        column_names=column_names,
        output_directory=aligned_hap_directory,
        n_threads=n_threads,
    )

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    summary = dict()

    summary_path = os.path.join(output_dir, "summary.csv")

    for test_dir in test_dirs:
        print("Evaluating: " + test_dir)
        config_path = os.path.join(test_dir,"config.json")
        test_path = os.path.join(test_dir, "assigned_haplotypes.fasta")
        all_haplotypes_path = os.path.join(test_dir, "all_paths.fasta")

        if not os.path.exists(config_path):
            print("Skipping because no config_path found in directory")
            continue

        if not os.path.exists(test_path):
            print("Skipping because no test_path found in directory")
            continue

        config = None
        with open(config_path, 'r') as file:
            config = json.load(file)

        print(config)
        flank_length = config["flank_length"]
        chromosome = config["chromosome"]
        start = config["ref_start"]
        stop = config["ref_stop"]

        prefix = chromosome + "_" + str(start) + "-" + str(stop)

        id_map = IncrementalIdMap()

        results = get_haplotypes_of_region(
            bam_paths=bam_paths,
            chromosome=chromosome,
            start=start,
            stop=stop,
            n_threads=n_threads
        )

        n_results = len([x for x in results if x is not None])
        if n_results == 0:
            print("Skipping because no haplotypes found in truthset")
            continue

        ref_sequences = {x.name:x for x in results if x is not None}
        test_sequences = {x.name:x for x in iterate_fasta(test_path)}

        feasible_samples = evaluate_upper_bound(
            ref_sequences=ref_sequences,
            haplotype_fasta_path=all_haplotypes_path,
            output_dir=output_dir,
            prefix=prefix,
            flank_size=flank_length,
            distance_threshold=distance_threshold,
            n_threads=n_threads)

        for key in sorted(list(ref_sequences.keys()) + list(test_sequences.keys())):
            id_map.add(key)

        # Align ref sequences to themselves (to build a tree)
        pairs,lengths,results = cross_align_sequences(ref_sequences, n_threads)

        ref_graph = networkx.Graph()
        ref_graph = construct_graph_from_alignments(pairs, lengths, results, id_map, ref_graph, distance_threshold)

        ref_graph_cache = ref_graph.copy()
        edges = list(ref_graph_cache.edges(data=True))

        edges_by_weight = defaultdict(list)
        for x in edges:
            edges_by_weight[x[2]["weight"]].append((x[0],x[1]))

        del edges

        tree = networkx.Graph()
        tree_id_map = IncrementalIdMap()

        distances = list()

        # Initialize the root component with all nodes
        component_series = [[{x for x in ref_graph.nodes}]]
        prev_component_map = {x:0 for x in ref_graph.nodes}
        root_name = str(0) + "_" + str(0)
        id = tree_id_map.add(root_name)
        tree.add_node(id, sample='')

        for distance,edges in sorted(edges_by_weight.items(), key=lambda x: x[0], reverse=True):
            affected_components = set()
            for a,b in edges:
                ref_graph_cache.remove_edge(a,b)

                affected_components.add(prev_component_map[a])
                affected_components.add(prev_component_map[b])

            components = list(networkx.connected_components(ref_graph_cache))
            component_map = dict()

            # IF this is the first component set, or if this component set is different from the last
            if len(components) != len(component_series[-1]):
                print(distance, len(components))
                distances.append(distance)

                tier = len(component_series)

                for c,nodes in enumerate(components):
                    name = str(tier) + "_" + str(c)
                    print(name)

                    sample = ""
                    if len(nodes) == 1:
                        sample = id_map.get_name(next(iter(nodes)))

                    id = tree_id_map.add(name)
                    tree.add_node(id, sample=sample)

                    for n in nodes:
                        component_map[n] = c
                        c_prev = prev_component_map[n]

                        name = str(tier) + "_" + str(c)
                        name_prev = str(tier-1) + "_" + str(c_prev)

                        print("new edge", name, name_prev)

                        a = tree_id_map.get_id(name)
                        b = tree_id_map.get_id(name_prev)
                        tree.add_edge(a,b)

                component_series.append(components)
                prev_component_map = component_map

        pos = graphviz_layout(tree, prog="dot")

        # pos = networkx.spring_layout(tree, weight='weight', iterations=50)

        labels = dict()
        for n,[id,values] in enumerate(tree.nodes(data=True)):
            labels[n] = (values["sample"])

        networkx.draw(tree,pos, alpha=0.6, node_size=20, labels=labels)

        # # Align test_sequences to the refs (to tag nodes in the tree)
        # test_pairs,test_lengths,test_results = align_sequences_to_other_sequences(ref_sequences, test_sequences, n_threads)
        #

        pyplot.show()
        pyplot.close()

    return


def main(
        input_dir,
        cache_dir,
        output_dir,
        n_threads,
        distance_threshold,
        tsv_path,
        column_names,
):

    # n_threads = 30
    # distance_threshold = 25
    #
    # output_dir = "/home/ryan/data/test_hapslap/evaluation/test_breakend/"
    # tsv_path = "/home/ryan/data/test_hapslap/terra/hprc_sandbox/subsample.tsv"
    # column_names = ["bam_hap1_vs_chm13","bam_hap2_vs_chm13"]
    #
    # input_directory = Path("/home/ryan/data/test_hapslap/results_breakend/")

    evaluate_test_haplotypes(
        input_dir=input_dir,
        cache_dir=cache_dir,
        output_dir=output_dir,
        n_threads=n_threads,
        distance_threshold=distance_threshold,
        tsv_path=tsv_path,
        column_names=column_names,
    )

    # write_formatted_alignments_per_sample()


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i","--input_dir",
        required=True,
        type=str,
        help="Output directory"
    )

    parser.add_argument(
        "-o","--output_dir",
        required=True,
        type=str,
        help="Output directory which will be created (and must not exist)"
    )

    parser.add_argument(
        "--cache_dir",
        required=True,
        type=str,
        help="Directory where remote BAMs will be stored, and potentially reused for future runs of this script"
    )

    parser.add_argument(
        "-t","--threads",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use"
    )

    parser.add_argument(
        "-d","--distance_threshold",
        required=False,
        default=25,
        type=int,
        help="Default clustering threshold for distance between haplotypes"
    )

    # TODO: refactor
    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="Path to a Terra-style haplotype TSV which contains relevant lookups for aligned hap1/hap2 ref assemblies per sample"
    )

    parser.add_argument(
        "-c","--column_names",
        required=False,
        default="'bam_hap1_vs_chm13','bam_hap2_vs_chm13'",
        type=parse_comma_separated_string,
        help="Which are the relevant columns of the Terra-style TSV (see tsv_path help msg) at the moment, column names must contain the substrings 'hap1' and 'hap2' (sorry)"
    )

    args = parser.parse_args()

    main(
        input_dir=args.input_dir,
        cache_dir=args.cache_dir,
        output_dir=args.output_dir,
        n_threads=args.threads,
        distance_threshold=args.distance_threshold,
        tsv_path=args.tsv,
        column_names=args.column_names
    )