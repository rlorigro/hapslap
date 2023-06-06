from modules.IncrementalIdMap import IncrementalIdMap
from modules.Align import run_minimap2

from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from copy import deepcopy
import glob
import os

from matplotlib import pyplot
from edlib import align
import networkx

import matplotlib
matplotlib.use('Agg')


def iterate_fasta(path):
    name = ""
    sequence = ""

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                # Output previous sequence
                if l > 0:
                    yield name, sequence

                name = line.strip()[1:].split(' ')[0]
                sequence = ""

            else:
                sequence += line.strip()

    yield name, sequence


def cross_align_fasta(path, n_threads):
    pairs = list()
    lengths = list()
    args = list()

    seqs = [x for x in iterate_fasta(path)]

    for a,b in combinations(seqs,2):
        if a == b:
            continue

        lengths.append((len(a[1]),len(b[1])))
        pairs.append((a[0],b[0]))
        args.append((a[1],b[1],"NW","distance"))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align, args)

    return pairs,lengths,results


def main():
    path = "/home/ryan/data/test_hapslap/haplotypes/"
    n_threads = 30

    paths = glob.glob(path + "*.fasta")

    for path in paths:
        name_id_map = IncrementalIdMap()

        pairs,lengths,results = cross_align_fasta(path,n_threads)

        graph = networkx.Graph()

        for i in range(len(results)):
            a,b = pairs[i]
            l_a,l_b = lengths[i]
            d = results[i]["editDistance"]

            w = float(l_a + l_b - d) / float(l_a + l_b)
            print(a,b,w)

            id_a = name_id_map.add(a)
            id_b = name_id_map.add(b)

            if id_a not in graph:
                graph.add_node(id_a)

            if id_b not in graph:
                graph.add_node(id_b)

            graph.add_edge(id_a, id_b, weight=w)

        # pos = networkx.spring_layout(graph, weight='weight')
        pos = networkx.spectral_layout(graph, weight='weight')

        labels = dict()
        for id,name in name_id_map:
            labels[id] = name

        fig = pyplot.figure()

        networkx.draw(graph,pos,
                      alpha=0.6
                      )

        pyplot.savefig(os.path.basename(path).replace(".fasta", ".png"), dpi=200)

        pyplot.close('all')


if __name__ == "__main__":
    main()

