from modules.Vcf import vcf_to_graph,write_graph_to_gfa,remove_empty_nodes_from_variant_graph
import argparse
import re
import os


def parse_region_string(s):
    tokens = re.split(r'[{:\-}]+', s.strip())
    chromosome = tokens[0]

    start = None
    stop = None
    if len(tokens) > 1:
        start = int(tokens[1])
        stop = int(tokens[2])

    return chromosome, start, stop


def main(vcf_path, ref_path, region_string, output_directory, keep_empty):
    data_per_sample = {"sample":{"vcf":vcf_path}}
    ref_sample_name = "ref"
    flank_length = 500

    chromosome, ref_start, ref_stop = parse_region_string(region_string)

    graph, alleles = vcf_to_graph(
        ref_path,
        data_per_sample,
        chromosome,
        ref_start,
        ref_stop,
        ref_sample_name,
        flank_length)

    if not keep_empty:
        graph, edge_to_allele_index = remove_empty_nodes_from_variant_graph(graph, alleles)

    name = os.path.basename(vcf_path)
    name = name.split('.')[0] + ".gfa"

    gfa_path = os.path.join(output_directory, name)

    print(gfa_path)
    write_graph_to_gfa(gfa_path, graph, alleles)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--vcf",
        required=True,
        type=str,
        help="Input vcf to be converted"
    )

    parser.add_argument(
        "--ref",
        required=True,
        type=str,
        help="FASTA file path of reference sequenced used to create VCF"
    )

    parser.add_argument(
        "--region",
        required=True,
        type=str,
        help="any valid htslib region string, e.g. 'chr20:69000-4200000'"
    )

    parser.add_argument(
        "--output_dir","-o",
        required=True,
        type=str,
        help="where to write the GFA output"
    )

    parser.add_argument('--keep_empty', action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    main(
        vcf_path=args.vcf,
        ref_path=args.ref,
        output_directory=args.output_dir,
        region_string=args.region,
        keep_empty=args.keep_empty,
    )
