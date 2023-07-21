from modules.IterativeHistogram import IterativeHistogram
from statistics import median,mean,mode
from matplotlib import pyplot
from pathlib import Path
import argparse
import numpy
import os


def main(input_path):
    y = dict()

    with open(input_path, 'r') as file:
        for line in file:
            tokens = line.strip().split(',')
            y[tokens[0]] = float(tokens[1])

    y = y.values()

    fig = pyplot.figure(layout='constrained',figsize=[5,5])

    ax = pyplot.axes()

    y_bins = numpy.arange(0, 1 + 0.01, 0.01, dtype=float)

    ax.hist(y,bins=y_bins)

    ax.set_xlabel("Portion of feasible haplotypes correctly assigned by optimizer")
    ax.set_ylabel("Frequency")

    pyplot.savefig("sample_support.png",dpi=200)

    pyplot.show()
    pyplot.close()

    print(len(y))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="input summary.csv"
    )

    args = parser.parse_args()

    main(
        input_path=args.i
    )
