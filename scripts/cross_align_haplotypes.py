import edlib

from modules.IncrementalIdMap import IncrementalIdMap
from modules.Align import run_minimap2
from modules.IterativeHistogram import IterativeHistogram

from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from copy import copy
import numpy
import glob
import os
import re

from matplotlib import pyplot
from edlib import align
import networkx

import matplotlib
matplotlib.use('Agg')


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


def cross_align_sequences(sequences:dict, n_threads):
    pairs = list()
    lengths = list()
    args = list()

    for a,b in combinations(sequences.values(),2):
        if a == b:
            continue

        lengths.append((len(a),len(b)))
        pairs.append((a.name,b.name))
        args.append((a.sequence,b.sequence,"NW","distance"))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align, args)

    return pairs,lengths,results


def align_sequences_to_other_sequences(a_seqs:dict, b_seqs:dict, n_threads):
    pairs = list()
    lengths = list()
    args = list()

    for a_name,a in a_seqs.items():
        for b_name,b in b_seqs.items():
            lengths.append((len(a),len(b)))
            pairs.append((a,b))
            args.append((a.sequence,b.sequence,"NW","path"))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align, args)

    return pairs,lengths,results


def orient_test_haps_by_best_match(ref_sequences, test_sequences, id_map):
    sample_names = set()
    for id,name in id_map:
        sample_name = name.split('_')[0]
        sample_names.add(sample_name)

    for sample_name in sample_names:
        a = sample_name + "_1"
        b = sample_name + "_2"

        aa = align(ref_sequences[a].sequence,test_sequences[a].sequence,"NW","distance")["editDistance"]
        bb = align(ref_sequences[b].sequence,test_sequences[b].sequence,"NW","distance")["editDistance"]
        ab = align(ref_sequences[a].sequence,test_sequences[b].sequence,"NW","distance")["editDistance"]
        ba = align(ref_sequences[b].sequence,test_sequences[a].sequence,"NW","distance")["editDistance"]

        cis_distance = aa + bb
        trans_distance = ab + ba

        # print("%.3f %.3f\n%.3f %.3f\n%s" % (float(aa),float(ab),float(ba),float(bb),str(cis_distance > trans_distance)))

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


def main():
    # TODO: do phase assignment step?

    output_dir = "/home/ryan/data/test_hapslap/evaluation/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    paths = [
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_7901318-7901522.fasta","/home/ryan/data/test_hapslap/results/chr20_7901318-7901522/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_10437604-10440525.fasta","/home/ryan/data/test_hapslap/results/chr20_10437604-10440525/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18259924-18261835.fasta","/home/ryan/data/test_hapslap/results/chr20_18259924-18261835/assigned_haplotypes.fasta"],
        ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18689217-18689256.fasta","/home/ryan/data/test_hapslap/results/chr20_18689217-18689256/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_18828383-18828733.fasta","/home/ryan/data/test_hapslap/results/chr20_18828383-18828733/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_47475093-47475817.fasta","/home/ryan/data/test_hapslap/results/chr20_47475093-47475817/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_49404497-49404943.fasta","/home/ryan/data/test_hapslap/results/chr20_49404497-49404943/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55000754-55000852.fasta","/home/ryan/data/test_hapslap/results/chr20_55000754-55000852/assigned_haplotypes.fasta"],
        # ["/home/ryan/data/test_hapslap/haplotypes/haplotypes_chr20_55486867-55492722.fasta","/home/ryan/data/test_hapslap/results/chr20_55486867-55492722/assigned_haplotypes.fasta"]
    ]

    n_threads = 30

    # paths = glob.glob(path + "*.fasta")

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
            d = results[i]["editDistance"]

            histogram.update(d)

            # w = d
            w = float(l_a + l_b - d) / float(l_a + l_b)

            id_a = id_map.get_id(a)
            id_b = id_map.get_id(b)

            if id_a not in graph:
                graph.add_node(id_a)

            if id_b not in graph:
                graph.add_node(id_b)

            if d <= 25:
                graph.add_edge(id_a, id_b, weight=w)

        labels = dict()
        for id,name in id_map:
            labels[id] = name

        fig = pyplot.figure()
        axes = pyplot.axes()

        components = list(networkx.connected_components(graph))

        colormap = pyplot.get_cmap("rainbow",len(components))

        component_map = dict()
        colors = dict()
        for c,component in enumerate(components):
            color = colormap(float(c)/float(len(components)))
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
                d = results[i]["editDistance"]

                id_ref = id_map.get_id(name_ref)
                id_test = id_map.get_id(name_test)

                w = float(l_ref + l_test - d) / float(l_ref + l_test)
                matrix[id_to_cluster_map[id_ref]][id_to_cluster_map[id_test]] = d

                # file.write(','.join([name_ref,str(id_ref),str(l_ref),name_test,str(id_test),str(l_test),str(d)]))
                # file.write('\n')
                # file.write(pairs[i][0].sequence)
                # file.write('\n')
                # file.write(pairs[i][1].sequence)
                # file.write('\n')
                # alignment = edlib.getNiceAlignment(results[i],pairs[i][0].sequence,pairs[i][1].sequence)
                # file.write(alignment['query_aligned'])
                # file.write('\n')
                # file.write(alignment['matched_aligned'])
                # file.write('\n')
                # file.write(alignment['target_aligned'])
                # file.write('\n')
                # file.write('\n')

        # sums = numpy.sum(matrix, axis=1)
        # matrix /= sums

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


if __name__ == "__main__":
    main()

