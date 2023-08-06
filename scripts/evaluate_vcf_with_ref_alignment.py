from collections import defaultdict
import argparse
import sys
import re
import os

from modules.Vcf import vcf_to_graph,compress_and_index_vcf,remove_empty_nodes_from_variant_graph
from modules.Cigar import get_haplotypes_of_region,char_is_query_move,char_is_ref_move
from modules.Bam import download_regions_of_bam
from modules.Align import run_minigraph
from modules.Sequence import Sequence


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


def evaluate_vcf_with_ref_alignment(vcf_path, ref_path, tsv_path, column_names, cache_dir, output_directory, n_threads):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        exit("ERROR: output directory exists already: %s" % output_directory)

    if vcf_path.endswith(".vcf"):
        vcf_path = compress_and_index_vcf(vcf_path, use_cache=True)

    sample_name = "test"

    data_per_sample = defaultdict(dict)
    data_per_sample[sample_name]["vcf"] = vcf_path

    region_string = os.path.basename(vcf_path).replace(".vcf.gz","")

    chromosome, ref_start, ref_stop = parse_region_string(region_string)
    ref_sample_name = "ref"
    flank_length = 500

    # First read the VCF and convert to a graph
    graph, alleles = vcf_to_graph(
        ref_path,
        data_per_sample,
        chromosome,
        ref_start,
        ref_stop,
        ref_sample_name,
        flank_length)

    graph, edge_to_deletion_index = remove_empty_nodes_from_variant_graph(graph, alleles)

    output_gfa_path = os.path.join(output_directory, "variant_graph.gfa")
    write_graph_to_gfa(output_gfa_path, graph, alleles)

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

    fasta_path = os.path.join(output_directory, "ref_haps.fasta")

    with open(fasta_path, 'w') as file:
        for s in ref_sequences.values():
            file.write(s.to_fasta_string())
            file.write('\n')

    # Align the relevant haplotypes to the variant graph
    output_gaf_path = run_minigraph(
        output_directory=output_directory,
        gfa_path=output_gfa_path,
        fasta_path=fasta_path
    )

    distance_per_query = defaultdict(int)
    query_coords = defaultdict(list)
    ref_distance = 0
    covered_nodes = set()

    with open(output_gaf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')
            query_name = tokens[0]

            # Column 5 (0-based) is the path column in a GAF file
            # It can be a forward or reverse alignment, so we will reinterpret them all as forward alignments
            path = None
            is_reverse = False
            if tokens[5][0] == '>':
                path = tuple(map(int,re.findall(r'\d+', tokens[5])))

            elif tokens[5][0] == '<':
                is_reverse = True
                path = tuple(map(int,reversed(re.findall(r'\d+', tokens[5]))))

            else:
                raise Exception("ERROR: Non GFA character found in path column of GFA")

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
                ref_start = flank_length

            ref_stop = ref_length
            if alleles[path[-1]].is_right_flank:
                ref_stop = ref_length - flank_length

            query_start = flank_length
            query_stop = query_length - flank_length

            # Don't have to deal with clips (thank the lord)
            # If the alignment is in reverse, ref coordinates decrease
            ref_cigar_start = int(tokens[7]) if not is_reverse else int(tokens[8])
            ref_cigar_stop = ref_cigar_start

            # Query coordinates always increase
            query_cigar_start = int(tokens[2])
            query_cigar_stop = 0

            ref_cost = 0
            query_cost = 0

            print(query_name, query_length)
            for l,c in cigar_tuples if not is_reverse else reversed(cigar_tuples):
                ref_cigar_stop = ref_cigar_start + l*char_is_ref_move[c] * (-1 if is_reverse else 1)
                query_cigar_stop = query_cigar_start + l*char_is_query_move[c]

                print(l,c,is_reverse,ref_cigar_start,ref_cigar_stop,ref_start,ref_stop)

                # Flip window bounds if it's a reverse alignment
                if is_reverse:
                    t = ref_cigar_stop
                    ref_cigar_stop = ref_cigar_start
                    ref_cigar_start = t

                is_match = True
                if c == 'X' or c == 'D' or c == 'I':
                    is_match = False

                if not is_match:
                    # Cigar is contained entirely in the window
                    if ref_start < ref_cigar_start <= ref_cigar_stop < ref_stop:
                        print("in window")
                        ref_cost += l

                    # Cigar covers entire window
                    elif ref_cigar_start <= ref_start < ref_stop <= ref_cigar_stop:
                        print("in window")
                        ref_cost += (ref_stop - ref_start)

                    # Cigar is overlapping the start of the window
                    elif ref_cigar_start <= ref_start <= ref_cigar_stop:
                        print("in window")
                        if char_is_ref_move[c]:
                            ref_cost += (ref_stop - ref_start)
                        else:
                            ref_cost += l

                    # Cigar is overlapping the end of the window
                    elif ref_cigar_start <= ref_stop <= ref_cigar_stop:
                        print("in window")
                        if char_is_ref_move[c]:
                            ref_cost += (ref_stop - ref_start)
                        else:
                            ref_cost += l

                if not is_match:
                    # Cigar is contained entirely in the window
                    if query_start < query_cigar_start <= query_cigar_stop < query_stop:
                        print("in window")
                        query_cost += l

                    # Cigar covers entire window
                    elif query_cigar_start <= query_start < query_stop <= query_cigar_stop:
                        print("in window")
                        query_cost += (query_stop - query_start)

                    # Cigar is overlapping the start of the window
                    elif query_cigar_start <= query_start <= query_cigar_stop:
                        print("in window")
                        if char_is_query_move[c]:
                            query_cost += (query_stop - query_start)
                        else:
                            query_cost += l

                    # Cigar is overlapping the end of the window
                    elif query_cigar_start <= query_stop <= query_cigar_stop:
                        print("in window")
                        if char_is_query_move[c]:
                            query_cost += (query_stop - query_start)
                        else:
                            query_cost += l

                # Unflip window bounds
                if is_reverse:
                    t = ref_cigar_stop
                    ref_cigar_stop = ref_cigar_start
                    ref_cigar_start = t

                ref_cigar_start = ref_cigar_stop
                query_cigar_start = query_cigar_stop

            distance_per_query[query_name] += query_cost
            ref_distance += ref_cost

    print(vcf_path)
    total_query_cost = 0
    for name,distance in distance_per_query.items():
        print(name,distance)
        total_query_cost += distance

    print(total_query_cost)
    print(ref_distance)
    print()


def evaluate_directory(input_directory, ref_path, tsv_path, column_names, cache_dir, output_directory, n_threads):
    for name in os.listdir(input_directory):
        path = os.path.join(input_directory, name)

        if path.endswith(".vcf.gz") or path.endswith(".vcf"):
            evaluate_vcf_with_ref_alignment(
                vcf_path=path,
                ref_path=ref_path,
                tsv_path=tsv_path,
                column_names=column_names,
                cache_dir=cache_dir,
                output_directory=output_directory,
                n_threads=n_threads)

            exit()


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i","--input_directory",
        required=True,
        type=str,
        help="Path to vcf to be evaluated"
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

    evaluate_directory(
        input_directory=args.input_directory,
        ref_path=args.ref_path,
        tsv_path=args.tsv,
        column_names=args.column_names,
        cache_dir=args.cache_dir,
        output_directory=args.output_dir,
        n_threads=args.n_threads
    )

