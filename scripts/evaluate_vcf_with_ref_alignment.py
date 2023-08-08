from collections import defaultdict
import argparse
import numpy
import sys
import re
import os

from matplotlib import pyplot

from modules.Vcf import vcf_to_graph,compress_and_index_vcf,remove_empty_nodes_from_variant_graph
from modules.Cigar import get_haplotypes_of_region,char_is_query_move,char_is_ref_move
from modules.IterativeHistogram import IterativeHistogram
from modules.Bam import download_regions_of_bam
from modules.Align import run_minigraph
from modules.Sequence import Sequence


class Alignment:
    def __init__(self, ref_name, query_name, ref_edit_distance, query_edit_distance, ref_start, ref_stop, query_start, query_stop):
        self.ref_name = ref_name
        self.query_name = query_name
        self.ref_edit_distance = ref_edit_distance
        self.query_edit_distance = query_edit_distance
        self.ref_start = ref_start
        self.ref_stop = ref_stop
        self.query_start = query_start
        self.query_stop = query_stop

    def __str__(self):
        s = ""
        s += "ref_name:\t" + self.ref_name + '\n'
        s += "query_name:\t" + self.query_name + '\n'
        s += "ref_edit_distance:\t" + str(self.ref_edit_distance) + '\n'
        s += "query_edit_distance:\t" + str(self.query_edit_distance) + '\n'
        s += "ref_start:\t" + str(self.ref_start) + '\n'
        s += "ref_stop:\t" + str(self.ref_stop) + '\n'
        s += "query_start:\t" + str(self.query_start) + '\n'
        s += "query_stop:\t" + str(self.query_stop) + '\n'

        return s


def write_graph_to_gfa(output_path, graph, alleles):
    with open(output_path, 'w') as gfa_file:
        for allele_index in graph.nodes:
            gfa_file.write("S\t%s\t%s\n" % (str(allele_index),alleles[allele_index].sequence))

        for e in graph.edges:
            gfa_file.write("L\t%s\t+\t%s\t+\t0M\n" % (str(e[0]), str(e[1])))


def parse_region_string(s):
    tokens = re.split(r'[{_:\-}]+', s.strip())
    chromosome = tokens[0]
    start = int(tokens[1])
    stop = int(tokens[2])

    return chromosome, start, stop


def iter_tuples_of_cigar_string(s):
    l_token = ""

    for c in s:
        if c.isnumeric():
            l_token += c
        else:
            l = int(l_token)
            l_token = ""

            yield l,c


def get_bounded_cigar_data(c, l, window_start, window_stop, cigar_start, cigar_stop, is_query):
    cost = 0
    min_index = None
    max_index = None

    # Some calculations need to know if the cigar operation is a move operation for ref or query
    is_move = char_is_ref_move
    if is_query:
        is_move = char_is_query_move

    # print("A", cigar_start, cigar_stop, window_start, window_stop)

    # Cigar is contained entirely in window
    if window_start <= cigar_start <= cigar_stop <= window_stop:
        # print("CONTAINED")
        min_index = cigar_start
        max_index = cigar_stop
        cost += l

    # Cigar covers entire window
    elif cigar_start <= window_start <= window_stop <= cigar_stop:
        # print("COVERS")
        min_index = window_start
        max_index = window_stop
        cost += (window_stop - window_start)

    # Cigar is overlapping the start of the window
    elif cigar_start <= window_start <= cigar_stop:
        min_index = window_start
        max_index = cigar_stop
        # print("LEFT")

        if is_move[c]:
            cost += (cigar_stop - window_start)
        else:
            cost += l

    # Cigar is overlapping the end of the window
    elif cigar_start <= window_stop <= cigar_stop:
        min_index = cigar_start
        max_index = window_stop
        # print("RIGHT")

        if is_move[c]:
            cost += (window_stop - cigar_start)
        else:
            cost += l

    return cost, min_index, max_index


class Results:
    def __init__(self, query_coverage, query_lengths, n_alignments, n_nodes, n_edges, nodes_covered, edges_covered):
        self.query_coverage = query_coverage
        self.query_lengths = query_lengths
        self.n_alignments = n_alignments
        self.n_nodes = n_nodes
        self.n_edges = n_edges
        self.nodes_covered = nodes_covered
        self.edges_covered = edges_covered

    def get_percent_query_coverage(self, query_name):
        return float(self.query_coverage[query_name]) / float(self.query_lengths[query_name])

    def iter_query_names(self):
        yield from self.query_lengths.keys()

    def __str__(self):
        s = ""
        s += "query_coverage:\t" + str(len(self.query_coverage)) + '\n'
        s += "query_lengths:\t" + str(len(self.query_lengths)) + '\n'
        s += "n_alignments:\t" + str(self.n_alignments) + '\n'
        s += "n_nodes:\t" + str(self.n_nodes) + '\n'
        s += "n_edges:\t" + str(self.n_edges) + '\n'
        s += "nodes_covered:\t" + str(self.nodes_covered) + '\n'
        s += "edges_covered:\t" + str(self.edges_covered) + '\n'

        return s


def evaluate_vcf_with_ref_alignment(vcf_path, ref_path, ref_sequences: dict, flank_length, buffer_length, output_directory, n_threads):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if vcf_path.endswith(".vcf"):
        vcf_path = compress_and_index_vcf(vcf_path, use_cache=True)

    sample_name = "test"

    data_per_sample = defaultdict(dict)
    data_per_sample[sample_name]["vcf"] = vcf_path

    region_string = os.path.basename(vcf_path).replace(".vcf.gz","")
    chromosome, ref_start, ref_stop = parse_region_string(region_string)

    # First read the VCF and convert to a graph
    graph, alleles = vcf_to_graph(
        ref_path,
        data_per_sample,
        chromosome=chromosome,
        ref_start=ref_start,
        ref_stop=ref_stop,
        ref_sample_name="ref",
        flank_length=flank_length,
        skip_incompatible=True
    )

    graph, edge_to_deletion_index = remove_empty_nodes_from_variant_graph(graph, alleles)

    output_gfa_path = os.path.join(output_directory, "variant_graph.gfa")
    write_graph_to_gfa(output_gfa_path, graph, alleles)

    fasta_path = os.path.join(output_directory, "ref_haps.fasta")

    with open(fasta_path, 'w') as file:
        for s in ref_sequences.values():
            file.write(s.to_fasta_string())
            file.write('\n')

    args = [
        "-c",
        "-g", str(20000),
        "-k", str(14),
        # "-f", "0.5",
        "-r", "10000,20000",
        "-n", "3,3",
        # "-m", "30,20",
        "-x", "lr",
    ]

    # Align the relevant haplotypes to the variant graph
    output_gaf_path = run_minigraph(
        output_directory=output_directory,
        gfa_path=output_gfa_path,
        fasta_path=fasta_path,
        n_threads=n_threads,
        args_override=args
    )

    alignments = defaultdict(list)
    ref_distance = 0
    n_alignments = 0
    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()
    nodes_covered = set()
    edges_covered = set()

    # Iterate the alignments and gather some stats
    with open(output_gaf_path, 'r') as file:
        for l,line in enumerate(file):
            n_alignments += 1
            tokens = line.strip().split('\t')
            query_name = tokens[0]

            # Column 5 (0-based) is the path column in a GAF file
            # It can be a forward or reverse alignment, so we will temporarily interpret them all as forward alignments
            path = None
            is_reverse = False
            if tokens[5][0] == '>':
                path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            elif tokens[5][0] == '<':
                is_reverse = True
                path = tuple(map(int,reversed(re.findall(r'\d+', tokens[5]))))

            else:
                raise Exception("ERROR: Non GFA character found in path column of GFA")

            for i in range(len(path)):
                # Give the flanking coverage for free
                if (alleles[i].is_left_flank or alleles[i].is_right_flank):
                    nodes_covered.add(path[i])

                # Assume that a node is fully covered if a path fully traverses it
                elif i > 0 and i < len(path) - 1:
                    nodes_covered.add(path[i])

                if i > 0:
                    edge = (path[i-1], path[i])
                    edges_covered.add(edge)

            is_spanning = False
            if len(path) > 1 and alleles[path[0]].is_left_flank and alleles[path[-1]].is_right_flank:
                is_spanning = True

            if not is_spanning:
                sys.stderr.write("WARNING: non-spanning alignment: " + query_name + '\n')

            cigar_string = None
            for token in tokens[11:]:
                if token.startswith("cg:Z:"):
                    cigar_string = token[5:]

            cigar_tuples = [x for x in iter_tuples_of_cigar_string(cigar_string)]

            ref_length = int(tokens[6])
            query_length = int(tokens[1])

            ref_start = 0
            if alleles[path[0]].is_left_flank:
                ref_start = (flank_length - buffer_length)

            ref_stop = ref_length
            if alleles[path[-1]].is_right_flank:
                ref_stop = ref_length - (flank_length - buffer_length)

            query_start = (flank_length - buffer_length)
            query_stop = query_length - (flank_length - buffer_length)

            # Don't have to deal with clips (thank the lord)
            # If the alignment is in reverse, ref coordinates decrease
            ref_cigar_start = int(tokens[7]) if not is_reverse else int(tokens[8])
            ref_cigar_stop = ref_cigar_start

            # Query coordinates always increase
            query_cigar_start = int(tokens[2])
            query_cigar_stop = 0

            ref_cost = 0
            query_cost = 0

            min_query_index = sys.maxsize
            max_query_index = -1

            for l,c in cigar_tuples if not is_reverse else reversed(cigar_tuples):
                ref_cigar_stop = ref_cigar_start + l*char_is_ref_move[c] * (-1 if is_reverse else 1)
                query_cigar_stop = query_cigar_start + l*char_is_query_move[c]

                # print(l,c,is_reverse,query_cigar_start,query_cigar_stop,query_start,query_stop)

                # Flip window bounds if it's a reverse alignment
                if is_reverse:
                    t = ref_cigar_stop
                    ref_cigar_stop = ref_cigar_start
                    ref_cigar_start = t

                is_match = True
                # if c == 'X' or c == 'D' or c == 'I':
                if c == 'D' or c == 'I':
                    is_match = False

                c, min_index, max_index = get_bounded_cigar_data(
                    c=c,
                    l=l,
                    window_start=query_start,
                    window_stop=query_stop,
                    cigar_start=query_cigar_start,
                    cigar_stop=query_cigar_stop,
                    is_query=True
                )

                # print("IN WINDOW" if min_index is not None and max_index is not None else "", c, min_index, max_index)

                if min_index is not None and min_index < min_query_index:
                    min_query_index = min_index

                if max_index is not None and max_index > max_query_index:
                    max_query_index = max_index

                if not is_match:
                    query_cost += c

                # Unflip window bounds
                if is_reverse:
                    t = ref_cigar_stop
                    ref_cigar_stop = ref_cigar_start
                    ref_cigar_start = t

                ref_cigar_start = ref_cigar_stop
                query_cigar_start = query_cigar_stop

            # Only count alignments that actually overlap the window to some extent
            if max_query_index != -1 and min_query_index != sys.maxsize:
                a = Alignment(
                    ref_name="ref",
                    query_name=query_name,
                    ref_edit_distance=ref_cost,
                    query_edit_distance=query_cost,
                    ref_start=None,
                    ref_stop=None,
                    query_start=min_query_index,
                    query_stop=max_query_index
                )

                # print()
                # print(a)

                alignments[query_name].append(a)
                ref_distance += ref_cost

    # Once the alignment blocks are grouped by sample, they need to be sorted and overlap needs to be subtracted.
    # Edit distance is not considered for overlap subtraction, which is crude, but easier to deal with for now.
    # Worst case, a near-boundary overlapping non-match gets double counted, reducing the score
    query_coverage = defaultdict(int)
    for query_name in alignments:
        prev_stop = None
        for a in sorted(alignments[query_name], key=lambda x: x.query_start):
            query_coverage[query_name] += a.query_stop - a.query_start + 1 - a.query_edit_distance
            if prev_stop is not None:
                overlap = prev_stop - a.query_start
                if overlap > 0:
                    query_coverage[query_name] -= overlap

            prev_stop = a.query_stop

    print(vcf_path)
    for name,coverage in query_coverage.items():
        l = len(ref_sequences[name]) - (flank_length - buffer_length)*2 + 1
        print(name,coverage,l,float(coverage)/float(l))

    print(ref_distance)
    print()

    query_lengths = dict()
    for s in ref_sequences.values():
        query_lengths[s.name] = len(s) - (flank_length - buffer_length)*2 + 1

    results = Results(
        query_coverage=query_coverage,
        query_lengths=query_lengths,
        n_alignments=n_alignments,
        n_nodes=n_nodes,
        n_edges=n_edges,
        nodes_covered=len(nodes_covered),
        edges_covered=len(edges_covered)
    )

    print(str(results))

    return results


def get_region_haps(region_string, cache_dir, tsv_path, column_names, n_threads):
    chromosome, ref_start, ref_stop = parse_region_string(region_string)
    flank_length = 800
    buffer_length = 10

    aligned_hap_directory = os.path.join(cache_dir, "aligned_haps")
    bam_paths = download_regions_of_bam(
        regions={chromosome},
        tsv_path=tsv_path,
        column_names=column_names,
        output_directory=aligned_hap_directory,
        n_threads=n_threads,
    )

    results = get_haplotypes_of_region(
        bam_paths=bam_paths,
        chromosome=chromosome,
        start=ref_start - flank_length,
        stop=ref_stop + flank_length,
        n_threads=n_threads
    )

    n_results = len([x for x in results if x is not None])
    if n_results == 0:
        print("Skipping because no haplotypes found in truthset")
        return

    ref_sequences = {x.name:x for x in results if x is not None}

    return ref_sequences


def evaluate_directories(input_directories: list, ref_path, tsv_path, column_names, cache_dir, output_directory, n_threads):
    flank_length = 800
    buffer_length = 10

    fig,axes = pyplot.subplots(nrows=2,ncols=3)

    # First gather all the regions to be evaluated
    region_strings = defaultdict(set)
    for i, input_directory in enumerate(input_directories):
        for name in sorted(os.listdir(input_directory)):
            path = os.path.join(input_directory, name)

            # if "16275838-16277739" not in path and "1893503-1907430" not in path:
            #     continue

            if path.endswith(".vcf.gz") or path.endswith(".vcf"):
                region_string = os.path.basename(path).split('.')[0]
                region_strings[region_string].add(i)

    identities_per_dir = [IterativeHistogram(start=0.0,stop=1.0,n_bins=100) for x in range(len(input_directories))]
    node_coverage_per_dir = [IterativeHistogram(start=0.0,stop=1.0,n_bins=100,unbounded_lower_bin=True) for x in range(len(input_directories))]
    edge_coverage_per_dir = [IterativeHistogram(start=0.0,stop=1.0,n_bins=100,unbounded_lower_bin=True) for x in range(len(input_directories))]
    n_nodes_per_dir = [IterativeHistogram(start=0.0,stop=800,n_bins=100,unbounded_upper_bin=True) for x in range(len(input_directories))]
    n_alignments_per_dir = [IterativeHistogram(start=0.0,stop=400,n_bins=100,unbounded_upper_bin=True) for x in range(len(input_directories))]

    # Iterate the regions, extract the ref haplotypes only once for each, and evaluate all input dirs for each
    for region_string,input_set in region_strings.items():
        print(region_string, input_set)

        if len(input_set) != len(input_directories):
            sys.stderr.write("WARNING: skipping region not found in all directories: " + region_string + '\n')
            continue

        ref_sequences = get_region_haps(region_string, cache_dir, tsv_path, column_names, n_threads)

        for i, input_directory in enumerate(input_directories):
            vcf_path = os.path.join(input_directory, region_string + ".vcf.gz")
            if not os.path.exists(vcf_path):
                vcf_path = os.path.join(input_directory, region_string + ".vcf")

            output_subdirectory = os.path.join(output_directory, str(i))
            results = evaluate_vcf_with_ref_alignment(
                vcf_path=vcf_path,
                ref_path=ref_path,
                ref_sequences=ref_sequences,
                flank_length=flank_length,
                buffer_length=buffer_length,
                output_directory=output_subdirectory,
                n_threads=n_threads)

            for name in results.iter_query_names():
                identity = results.get_percent_query_coverage(name)
                node_coverage = float(results.nodes_covered) / float(results.n_nodes) if results.n_nodes > 0 else 1.0
                edge_coverage = float(results.edges_covered) / float(results.n_edges) if results.n_edges > 0 else 1.0

                identities_per_dir[i].update(identity)
                node_coverage_per_dir[i].update(node_coverage)
                edge_coverage_per_dir[i].update(edge_coverage)
                n_nodes_per_dir[i].update(results.n_nodes)
                n_alignments_per_dir[i].update(results.n_alignments)

    for h,histogram in enumerate(identities_per_dir):
        label = os.path.basename(input_directories[h])

        x = histogram.get_bin_centers()
        y = 1 - numpy.cumsum(histogram.get_normalized_histogram())
        axes[0][0].plot(x,y,alpha=0.6, color="C" + str(h), label=label)
        axes[0][0].set_xlabel("identities")
        axes[0][0].set_ylabel("Cumulative density")
        axes[0][0].set_xlim([-0.05,1.05])
        axes[0][0].invert_xaxis()

    for h,histogram in enumerate(n_nodes_per_dir):
        label = os.path.basename(input_directories[h])

        x = histogram.get_bin_centers()
        y = numpy.cumsum(histogram.get_normalized_histogram())
        axes[0][2].plot(x,y,alpha=0.6, color="C" + str(h), label=label)
        axes[0][2].set_xlabel("n_nodes")
        axes[0][2].set_ylabel("Cumulative density")

    for h,histogram in enumerate(node_coverage_per_dir):
        label = os.path.basename(input_directories[h])

        x = histogram.get_bin_centers()
        y = 1 - numpy.cumsum(histogram.get_normalized_histogram())
        axes[1][0].plot(x,y,alpha=0.6, color="C" + str(h), label=label)
        axes[1][0].set_xlabel("node_coverage")
        axes[1][0].set_ylabel("Cumulative density")
        axes[1][0].set_xlim([-0.05,1.05])
        axes[1][0].invert_xaxis()

    for h,histogram in enumerate(edge_coverage_per_dir):
        label = os.path.basename(input_directories[h])

        x = histogram.get_bin_centers()
        y = 1 - numpy.cumsum(histogram.get_normalized_histogram())
        axes[1][1].plot(x,y,alpha=0.6, color="C" + str(h), label=label)
        axes[1][1].set_xlabel("edge_coverage")
        axes[1][1].set_ylabel("Cumulative density")
        axes[1][1].set_xlim([-0.05,1.05])
        axes[1][1].invert_xaxis()

    for h,histogram in enumerate(n_alignments_per_dir):
        label = os.path.basename(input_directories[h])

        x = histogram.get_bin_centers()
        y = numpy.cumsum(histogram.get_normalized_histogram())
        axes[0][1].plot(x,y,alpha=0.6, color="C" + str(h), label=label)
        axes[0][1].set_xlabel("n_alignments")
        axes[0][1].set_ylabel("Cumulative density")

    axes[1][2].axis("off")

    fig.set_size_inches(14,8)
    axes[0][0].legend()

    pyplot.savefig("comparison.png", dpi=300)
    pyplot.show()
    pyplot.close()


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i","--input_directories",
        required=True,
        type=parse_comma_separated_string,
        help="Multiple directories, each containing a set of VCFs to be evaluated"
    )

    parser.add_argument(
        "-r","--ref_path",
        required=True,
        type=str,
        help="Path to ref sequence fasta which was used for alignment in VCF and truth assemblies"
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
        help="Directory where remote asm-to-ref BAMs will be stored, and potentially reused for future runs of this script"
    )

    parser.add_argument(
        "-t","--n_threads",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use"
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

    evaluate_directories(
        input_directories=args.input_directories,
        ref_path=args.ref_path,
        tsv_path=args.tsv,
        column_names=args.column_names,
        cache_dir=args.cache_dir,
        output_directory=args.output_dir,
        n_threads=args.n_threads
    )

