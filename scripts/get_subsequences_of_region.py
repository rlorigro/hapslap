from modules.Cigar import iter_query_sequences_of_region
from modules.Bam import get_region_from_bam,index_bam
from modules.Authenticator import *

from multiprocessing import Pool
import argparse
import pandas
import os
import re


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


def parse_region_string(s):
    tokens = re.split(r'[{:\-}]+', s.strip())
    chromosome = tokens[0]
    start = int(tokens[1])
    stop = int(tokens[2])

    return chromosome, start, stop


def get_region(tsv_path, n_threads, samples, column_names, region_string, output_directory):
    chromosome, start, stop = parse_region_string(region_string)

    token = GoogleToken()
    index_threads = 1

    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    df = pandas.read_table(tsv_path, sep='\t', header=0)

    n_rows, n_cols = df.shape

    print(n_rows)

    args = list()

    samples = set(samples)
    for column_name in column_names:
        for i in range(n_rows):
            sample_name = df.iloc[i][0]

            if sample_name not in samples:
                continue

            print(sample_name)

            gs_uri = df.iloc[i][column_name]

            # Each tool downloads its regions to its own subdirectory to prevent overwriting (filenames are by region)
            output_subdirectory = os.path.join(output_directory, sample_name)
            if not os.path.exists(output_subdirectory):
                os.makedirs(output_subdirectory)

            filename = sample_name + "_" + column_name + "_" + region_string.replace(":","_") + ".bam"
            args.append([output_subdirectory,gs_uri,region_string,token,600,filename,True])

    with Pool(n_threads) as pool:
        download_results = pool.starmap(get_region_from_bam, args)

    output_fasta = os.path.join(output_directory, "haplotypes" + "_" + region_string.replace(":","_") + ".fasta")
    sequences = list()

    with open(output_fasta, 'w') as file:
        for path in download_results:
            suffix = None

            if "hap1" in path:
                suffix = "hap1"

            if "hap2" in path:
                suffix = "hap2"

            index_path = path + ".bai"

            if not os.path.exists(index_path):
                index_bam(path, index_threads)

            iter = iter_query_sequences_of_region(bam_path=path, chromosome=chromosome, ref_start=start, ref_stop=stop)

            for item in iter:
                sequences.append(item)

        for query_name, is_reverse, query_start, query_stop, seq in sorted(sequences, key=lambda x: len(x[-1]), reverse=True):
            file.write(">" + query_name + "_" + suffix)
            file.write('\n')
            file.write(seq)
            file.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="Input tsv containing URIs, linking to tarballs to be parsed"
    )

    parser.add_argument(
        "-t",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use"
    )

    parser.add_argument(
        "-s",
        required=True,
        type=parse_comma_separated_string,
        help="Which samples to download, as a comma separated list"
    )

    parser.add_argument(
        "-c",
        required=True,
        type=parse_comma_separated_string,
        help="Column names in tsv which pertain to relevant data to download. Each column will be given its own directory"
    )

    parser.add_argument(
        "-r",
        required=True,
        type=str,
        help="SAMtools formatted region, e.g. chr1:2000-4000 "
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()
    get_region(
        tsv_path=args.tsv,
        n_threads=args.t,
        samples=args.s,
        column_names=args.c,
        region_string=args.r,
        output_directory=args.o
    )
