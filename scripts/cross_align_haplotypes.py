from pywfa import WavefrontAligner
import edlib

from modules.IncrementalIdMap import IncrementalIdMap
from modules.Align import run_minimap2
from modules.IterativeHistogram import IterativeHistogram

from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from glob import glob
from copy import copy
import numpy
import os
import re

from matplotlib import pyplot,colors
from edlib import align
import networkx

import matplotlib
matplotlib.use('Agg')


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap


class Sequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def normalize_name(self):
        tokens = re.split("#|_", self.name)
        self.name = '_'.join([tokens[0],self.name[-1]])

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.name == other.name and self.sequence == other.sequence
        else:
            return False

    def __len__(self):
        return len(self.sequence)


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
    aligner = WavefrontAligner()

    aligner(a.sequence,b.sequence)
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

        aa = align_and_get_indel_distance(ref_sequences[a],test_sequences[a])
        bb = align_and_get_indel_distance(ref_sequences[b],test_sequences[b])
        ab = align_and_get_indel_distance(ref_sequences[a],test_sequences[b])
        ba = align_and_get_indel_distance(ref_sequences[b],test_sequences[a])

        cis_distance = aa + bb
        trans_distance = ab + ba

        # Flip the sequences around if it is a better match
        if cis_distance > trans_distance:
            a_seq = test_sequences[a].sequence
            b_seq = test_sequences[b].sequence

            test_sequences[a] = Sequence(a,b_seq)
            test_sequences[b] = Sequence(b,a_seq)

    return test_sequences


def write_sample_alignments(ref_sequences, test_sequences, id_map, output_dir):
    sample_names = set()
    for id,name in id_map:
        sample_name = name.split('_')[0]
        sample_names.add(sample_name)

    for sample_name in sample_names:
        a = sample_name + "_1"
        b = sample_name + "_2"

        aligner = WavefrontAligner()

        output_path = os.path.join(output_dir, a + "_" + a + ".txt")
        aligner(ref_sequences[a].sequence,test_sequences[a].sequence)
        aligner.cigar_print_pretty(output_path)

        output_path = os.path.join(output_dir, b + "_" + b + ".txt")
        aligner(ref_sequences[b].sequence,test_sequences[b].sequence)
        aligner.cigar_print_pretty(output_path)

        output_path = os.path.join(output_dir, a + "_" + b + ".txt")
        aligner(ref_sequences[a].sequence,test_sequences[b].sequence)
        aligner.cigar_print_pretty(output_path)

        output_path = os.path.join(output_dir, b + "_" + a + ".txt")
        aligner(ref_sequences[b].sequence,test_sequences[a].sequence)
        aligner.cigar_print_pretty(output_path)


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


def write_formatted_alignments_per_sample():
    output_dir = "/home/ryan/data/test_hapslap/evaluation/sample_alignments/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    paths = [
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_7901318-7901522.fasta","/home/ryan/data/test_hapslap/results/chr20_7901318-7901522/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_10437604-10440525.fasta","/home/ryan/data/test_hapslap/results/chr20_10437604-10440525/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18259924-18261835.fasta","/home/ryan/data/test_hapslap/results/chr20_18259924-18261835/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18689217-18689256.fasta","/home/ryan/data/test_hapslap/results/chr20_18689217-18689256/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18828383-18828733.fasta","/home/ryan/data/test_hapslap/results/chr20_18828383-18828733/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_47475093-47475817.fasta","/home/ryan/data/test_hapslap/results/chr20_47475093-47475817/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_49404497-49404943.fasta","/home/ryan/data/test_hapslap/results/chr20_49404497-49404943/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55000754-55000852.fasta","/home/ryan/data/test_hapslap/results/chr20_55000754-55000852/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55486867-55492722.fasta","/home/ryan/data/test_hapslap/results/chr20_55486867-55492722/assigned_haplotypes.fasta"]
    ]

    for ref_path,test_path in paths:
        print(ref_path)
        print(test_path)

        id_map = IncrementalIdMap()

        ref_sequences = {x.name:x for x in iterate_fasta(ref_path)}
        test_sequences = {x.name:x for x in iterate_fasta(test_path)}

        duplicated_homozygous_haps = set()
        test_sequences,test_duplicated = duplicate_homozygous_haps(test_sequences)
        ref_sequences,ref_duplicated = duplicate_homozygous_haps(ref_sequences)
        duplicated_homozygous_haps = duplicated_homozygous_haps.union(test_duplicated)
        duplicated_homozygous_haps = duplicated_homozygous_haps.union(ref_duplicated)

        print(duplicated_homozygous_haps)

        for key in list(ref_sequences.keys()):
            if key not in test_sequences:
                print("WARNING: haplotype not in test_sequences: " + key)
                del ref_sequences[key]

        for key in sorted(ref_sequences):
            id_map.add(key)

        write_sample_alignments(ref_sequences, test_sequences, id_map, output_dir)


def evaluate_upper_bound():
    output_dir = "/home/ryan/data/test_hapslap/evaluation/upper_bound"
    n_threads = 30
    flank_size = 20000
    distance_threshold = 25

    paths = [
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_7901318-7901522.fasta","/home/ryan/data/test_hapslap/results/chr20_7901318-7901522/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_10437604-10440525.fasta","/home/ryan/data/test_hapslap/results/chr20_10437604-10440525/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18259924-18261835.fasta","/home/ryan/data/test_hapslap/results/chr20_18259924-18261835/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18689217-18689256.fasta","/home/ryan/data/test_hapslap/results/chr20_18689217-18689256/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18828383-18828733.fasta","/home/ryan/data/test_hapslap/results/chr20_18828383-18828733/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_47475093-47475817.fasta","/home/ryan/data/test_hapslap/results/chr20_47475093-47475817/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_49404497-49404943.fasta","/home/ryan/data/test_hapslap/results/chr20_49404497-49404943/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55000754-55000852.fasta","/home/ryan/data/test_hapslap/results/chr20_55000754-55000852/paths/"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55486867-55492722.fasta","/home/ryan/data/test_hapslap/results/chr20_55486867-55492722/paths/"]
    ]

    for ref_path,test_path in paths:
        test_fasta_paths = [os.path.join(test_path,x) for x in glob("*.fasta", root_dir=test_path)]

        ref_sequences = {x.name:x for x in iterate_fasta(ref_path)}
        test_sequences = dict()

        for path in test_fasta_paths:
            for item in iterate_fasta(path, normalize_name=False):
                item.sequence = item.sequence[flank_size:len(item.sequence) - flank_size]
                test_sequences[item.name] = item

        ref_id_map = IncrementalIdMap()
        test_id_map = IncrementalIdMap()

        for key in sorted(ref_sequences):
            ref_id_map.add(key)

        for key in sorted(test_sequences):
            test_id_map.add(key)

        pairs,lengths,results = cross_align_sequences(sequences=ref_sequences,n_threads=n_threads) #,output_dir=output_dir)

        graph = networkx.Graph()

        for i in range(len(results)):
            a,b = pairs[i]
            l_a,l_b = lengths[i]
            d = results[i]

            # w = d
            w = float(l_a + l_b - d) / float(l_a + l_b)

            id_a = ref_id_map.get_id(a)
            id_b = ref_id_map.get_id(b)

            if id_a not in graph:
                graph.add_node(id_a)

            if id_b not in graph:
                graph.add_node(id_b)

            if d <= distance_threshold:
                graph.add_edge(id_a, id_b, weight=w)

        ref_labels = dict()
        for id,name in ref_id_map:
            ref_labels[id] = name

        components = list(networkx.connected_components(graph))

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
            n_threads=n_threads,
            # output_dir=output_dir
        )

        matrix = numpy.zeros([len(ref_id_map),len(test_id_map)])

        for i in range(len(results)):
            name_ref = pairs[i][0].name
            name_test = pairs[i][1].name

            l_ref,l_test = lengths[i]
            d = results[i]

            id_ref = ref_id_map.get_id(name_ref)
            id_test = test_id_map.get_id(name_test)

            matrix[id_to_cluster_map[id_ref]][id_test] = d

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
        output_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_confusion.png"))
        fig.savefig(output_path, dpi=500, bbox_inches='tight')
        pyplot.close('all')


def evaluate_test_haplotypes():
    n_threads = 30
    distance_threshold = 25

    output_dir = "/home/ryan/data/test_hapslap/evaluation/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    paths = [
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_7901318-7901522.fasta","/home/ryan/data/test_hapslap/results/chr20_7901318-7901522/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_10437604-10440525.fasta","/home/ryan/data/test_hapslap/results/chr20_10437604-10440525/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18259924-18261835.fasta","/home/ryan/data/test_hapslap/results/chr20_18259924-18261835/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18689217-18689256.fasta","/home/ryan/data/test_hapslap/results/chr20_18689217-18689256/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18828383-18828733.fasta","/home/ryan/data/test_hapslap/results/chr20_18828383-18828733/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_47475093-47475817.fasta","/home/ryan/data/test_hapslap/results/chr20_47475093-47475817/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_49404497-49404943.fasta","/home/ryan/data/test_hapslap/results/chr20_49404497-49404943/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55000754-55000852.fasta","/home/ryan/data/test_hapslap/results/chr20_55000754-55000852/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55486867-55492722.fasta","/home/ryan/data/test_hapslap/results/chr20_55486867-55492722/assigned_haplotypes.fasta"]
    ]

    for ref_path,test_path in paths:
        print(ref_path)
        print(test_path)

        id_map = IncrementalIdMap()

        ref_sequences = {x.name:x for x in iterate_fasta(ref_path)}
        test_sequences = {x.name:x for x in iterate_fasta(test_path)}

        duplicated_homozygous_haps = set()
        test_sequences,test_duplicated = duplicate_homozygous_haps(test_sequences)
        ref_sequences,ref_duplicated = duplicate_homozygous_haps(ref_sequences)
        duplicated_homozygous_haps = duplicated_homozygous_haps.union(test_duplicated)
        duplicated_homozygous_haps = duplicated_homozygous_haps.union(ref_duplicated)

        print(duplicated_homozygous_haps)

        for key in list(ref_sequences.keys()):
            if key not in test_sequences:
                print("WARNING: haplotype not in test_sequences: " + key)
                del ref_sequences[key]

        for key in sorted(ref_sequences):
            id_map.add(key)

        test_sequences = orient_test_haps_by_best_match(ref_sequences, test_sequences, id_map)

        pairs,lengths,results = cross_align_sequences(ref_sequences,n_threads)

        graph = networkx.Graph()

        histogram = IterativeHistogram(start=0, stop=2000, n_bins=2000)

        for i in range(len(results)):
            a,b = pairs[i]
            l_a,l_b = lengths[i]
            d = results[i]

            histogram.update(d)

            # w = d
            w = float(l_a + l_b - d) / float(l_a + l_b)

            id_a = id_map.get_id(a)
            id_b = id_map.get_id(b)

            if id_a not in graph:
                graph.add_node(id_a)

            if id_b not in graph:
                graph.add_node(id_b)

            if d <= distance_threshold:
                graph.add_edge(id_a, id_b, weight=w)

        labels = dict()
        for id,name in id_map:
            labels[id] = name

        fig = pyplot.figure()
        axes = pyplot.axes()

        components = list(networkx.connected_components(graph))

        colormap = pyplot.get_cmap("rainbow")
        colormap = truncate_colormap(colormap, 0, 0.7, 100)

        component_map = dict()
        colors = dict()
        for c,component in enumerate(components):
            print(float(c),float(len(components)))
            color_index = float(c+0.001)/float(len(components))
            print(color_index)

            color = colormap(color_index)
            print(c,[id_map.get_name(id) for id in component])
            print('\t',color)

            for n in component:
                component_map[n] = c
                colors[n] = color

        colors_list = list()

        # For fucks sake what is the point of assigning a node ID if Networkx doesn't even use it as an index
        for i,n in enumerate(graph.nodes):
            colors_list.append(colors[n])

        pos = networkx.spring_layout(graph, weight='weight')
        # pos = networkx.spectral_layout(graph, weight='weight')

        networkx.draw(graph,pos,alpha=0.6,node_color=colors_list,node_size=20)

        output_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_graph.png"))
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
        pyplot.close('all')

        fig = pyplot.figure()
        axes = pyplot.axes()

        x = histogram.get_bin_centers()
        y = histogram.get_histogram()

        end = len(y)
        for i in reversed(y):
            if i > 0:
                break

            end -= 1

        x = x[:min(len(y),end+1)]
        y = y[:min(len(y),end+1)]

        axes.plot(x,y)

        axes.set_xlabel("Edit distance")
        axes.set_ylabel("Frequency (#)")
        axes.set_title("Distances in local all-vs-all ref HPRC alignment")

        histogram_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_histogram.png"))
        fig.savefig(histogram_path, dpi=200, bbox_inches='tight')

        pyplot.close('all')

        pairs,lengths,results = align_sequences_to_other_sequences(
            a_seqs=ref_sequences,
            b_seqs=test_sequences,
            n_threads=n_threads
        )

        matrix = numpy.zeros([len(id_map),len(id_map)])

        # This will order all IDs so that clusters are contiguous, otherwise randomly ordered
        id_to_cluster_map = dict()
        for id in sorted(range(len(id_map)), key=lambda x: component_map[x]):
            id_to_cluster_map[id] = len(id_to_cluster_map)

        output_path = os.path.join(output_dir,"alignments.txt")
        with open(output_path, 'w') as file:
            for i in range(len(results)):
                name_ref = pairs[i][0].name
                name_test = pairs[i][1].name

                l_ref,l_test = lengths[i]
                d = results[i]

                id_ref = id_map.get_id(name_ref)
                id_test = id_map.get_id(name_test)

                w = float(l_ref + l_test - d) / float(l_ref + l_test)
                matrix[id_to_cluster_map[id_ref]][id_to_cluster_map[id_test]] = d

                id_test += len(id_map)

                if id_test not in graph:
                    graph.add_node(id_test)

                if d <= distance_threshold:
                    graph.add_edge(id_ref, id_test, weight=w)

        fig = pyplot.figure()
        axes = pyplot.axes()

        for i in range(len(graph.nodes)- len(id_map)) :
            colors_list.append("red")

        pos = networkx.spring_layout(graph, weight='weight')

        networkx.draw(graph,pos,alpha=0.6,node_color=colors_list,linewidths=0,width=0.5,node_size=20)

        output_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_graph_with_test_sequence.png"))
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
        pyplot.close('all')

        fig = pyplot.figure()
        axes = pyplot.axes()
        p = axes.matshow(matrix,cmap='viridis')
        pyplot.colorbar(p, ax=axes)

        cumulative_count = 0
        for c in components:
            cumulative_count += len(c)
            axes.axhline(cumulative_count - 0.5,color="red",linewidth=0.5)
            axes.axvline(cumulative_count - 0.5,color="red",linewidth=0.5)

        axes.tick_params(axis='both', which='major', labelsize=3)

        # Construct tick labels based on the clustering order
        labels = [None]*len(id_map.id_to_name)
        for id,name in id_map:
            labels[id_to_cluster_map[id]] = name

        axes.xaxis.set_ticks_position('bottom')

        axes.set_yticks(list(range(len(id_map))))
        axes.set_yticklabels(labels)
        axes.set_xticks(list(range(len(id_map))))
        axes.set_xticklabels(labels,rotation='vertical')

        axes.set_ylabel("Ref Haplotypes")
        axes.set_xlabel("Test Haplotypes")
        axes.set_title("Edit distance of clustered haplotypes")

        # For debugging purposes, color the samples which were duplicated
        for label in axes.get_xticklabels():
            if str(label.get_text()) in duplicated_homozygous_haps:
                label.set_color("blue")

        pyplot.tight_layout()
        output_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_confusion.png"))
        fig.savefig(output_path, dpi=500, bbox_inches='tight')
        pyplot.close('all')


def main():
    evaluate_test_haplotypes()
    evaluate_upper_bound()
    # write_formatted_alignments_per_sample()


if __name__ == "__main__":
    main()
