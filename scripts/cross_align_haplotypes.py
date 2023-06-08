import edlib

from modules.IncrementalIdMap import IncrementalIdMap
from modules.Align import run_minimap2
from modules.IterativeHistogram import IterativeHistogram

from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from copy import deepcopy
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


def iterate_fasta(path, force_upper_case=True):
    name = ""
    sequence = ""

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                if force_upper_case:
                    sequence = sequence.upper()

                # Output previous sequence
                if l > 0:
                    yield Sequence(name, sequence)

                name = line.strip()[1:].split(' ')[0]
                sequence = ""

            else:
                sequence += line.strip()

    if force_upper_case:
        sequence = sequence.upper()

    yield Sequence(name, sequence)


def cross_align_sequences(sequences, n_threads):
    pairs = list()
    lengths = list()
    args = list()

    for a,b in combinations(sequences,2):
        if a == b:
            continue

        lengths.append((len(a),len(b)))
        pairs.append((a.name,b.name))
        args.append((a.sequence,b.sequence,"NW","distance"))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align, args)

    return pairs,lengths,results


def align_sequences_to_other_sequences(a_seqs, b_seqs, n_threads):
    pairs = list()
    lengths = list()
    args = list()

    for a in a_seqs:
        for b in b_seqs:
            lengths.append((len(a),len(b)))
            pairs.append((a,b))
            args.append((a.sequence,b.sequence,"NW","path"))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align, args)

    return pairs,lengths,results


def main():
    # TODO: do phase assignment step?

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

    n_threads = 30

    # paths = glob.glob(path + "*.fasta")

    for ref_path,test_path in paths:
        print(ref_path)
        print(test_path)

        id_map = IncrementalIdMap()

        ref_sequences = [x for x in iterate_fasta(ref_path)]
        test_sequences = [x for x in iterate_fasta(test_path)]

        for i in range(len(ref_sequences)):
            ref_sequences[i].normalize_name()

        for i in range(len(test_sequences)):
            test_sequences[i].normalize_name()

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

            id_a = id_map.add(a)
            id_b = id_map.add(b)

            if id_a not in graph:
                graph.add_node(id_a)

            if id_b not in graph:
                graph.add_node(id_b)

            if d < 25:
                graph.add_edge(id_a, id_b, weight=w)

        labels = dict()
        for id,name in id_map:
            labels[id] = name

        fig2 = pyplot.figure()
        axes2 = pyplot.axes()

        components = list(networkx.connected_components(graph))

        colormap = pyplot.get_cmap("rainbow",len(id_map))

        component_map = dict()
        colors = dict()
        for c,component in enumerate(components):
            color = colormap(float(c)/float(len(components)))
            print(c,component)

            for n in component:
                component_map[n] = c
                colors[n] = color

        colors = [colors[x] for x in sorted(colors)]

        pos = networkx.spring_layout(graph, weight='weight')
        # pos = networkx.spectral_layout(graph, weight='weight')

        networkx.draw(graph,pos,alpha=0.6,node_color=colors,node_size=20)

        output_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_graph.png"))
        pyplot.savefig(output_path, dpi=200)
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

        pyplot.savefig(os.path.basename(ref_path).replace(".fasta", "_histogram.png"), dpi=200)

        pyplot.close('all')

        pairs,lengths,results = align_sequences_to_other_sequences(
            a_seqs=ref_sequences,
            b_seqs=test_sequences,
            n_threads=n_threads
        )

        matrix = numpy.zeros([len(id_map),len(id_map)])

        new_id_map = dict()
        for id in sorted(range(len(id_map)), key=lambda x: component_map[x]):
            new_id_map[id] = len(new_id_map)
            print(id_map.get_name(id), component_map[id])

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
                matrix[new_id_map[id_ref]][new_id_map[id_test]] = d

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

        sums = numpy.sum(matrix, axis=0)
        matrix /= sums

        pyplot.figure()
        axes = pyplot.axes()
        axes.matshow(matrix)

        print(','.join(sorted(id_map.id_to_name)))
        print(','.join(sorted(id_map.id_to_name)))

        axes.tick_params(axis='both', which='major', labelsize=4)

        axes.set_yticks(list(range(len(id_map))))
        axes.set_yticklabels(id_map.id_to_name)
        axes.set_xticks(list(range(len(id_map))))
        axes.set_xticklabels(id_map.id_to_name,rotation='vertical')

        pyplot.tight_layout()
        output_path = os.path.join(output_dir,os.path.basename(ref_path).replace(".fasta", "_confusion.png"))
        pyplot.savefig(output_path, dpi=200)


if __name__ == "__main__":
    main()

