from modules.GsUri import download_gs_uri

from multiprocessing import Pool
import argparse
import pandas
import os
import re


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


def download_sample_data(tsv_path, n_threads, samples, column_names, output_directory):
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

            # This is going to get mangled if anyone uses '.' in their filename
            suffix = '.'.join(gs_uri.split('.')[1:])
            output_filename = column_name + '.' + suffix
            args.append([gs_uri,output_subdirectory,output_filename])

    for item in args:
        print(item)

    with Pool(n_threads) as pool:
        download_results = pool.starmap(download_gs_uri, args)


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
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()
    download_sample_data(
        tsv_path=args.tsv,
        n_threads=args.t,
        samples=args.s,
        column_names=args.c,
        output_directory=args.o
    )
