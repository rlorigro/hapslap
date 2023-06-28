from modules.IterativeHistogram import IterativeHistogram
from statistics import median,mean,mode
from matplotlib import pyplot
from pathlib import Path
import argparse
import numpy
import os


def main(input_dir):
    test_dirs = [str(f) for f in Path(input_dir).iterdir() if f.is_dir()]

    names = list()
    n = list()
    n_candidates = list()
    n_spanning_reads = list()
    n_valid_pairs = list()
    time_elapsed_minigraph_s = list()
    time_elapsed_minimap_s = list()
    time_elapsed_optimizer_s = list()
    time_elapsed_other_s = list()
    time_elapsed_total_s = list()

    for test_dir in test_dirs:
        summary_path = os.path.join(test_dir, "summary.txt")
        costs_path = os.path.join(test_dir, "costs.csv")

        names.append(test_dir.split('/')[-1])

        if not os.path.exists(summary_path):
            continue

        p = 0
        with open(costs_path, 'r') as file:
            for line in file:
                p += 1

        n_valid_pairs.append(p)

        with open(summary_path, 'r') as file:
            for line in file:
                tokens = line.split(',')

                if len(tokens) != 2:
                    raise Exception("ERROR: Unparsable line: " + line)

                key,value = tokens

                # n,6
                # n_candidates,4
                # n_spanning_reads,452
                # cost_per_read,1.01
                # distance,28.0000
                # time_elapsed_minigraph_s,0
                # time_elapsed_minimap_s,9
                # time_elapsed_optimizer_s,34
                # time_elapsed_total_s,55
                if key == "n":
                    n.append(int(value))

                elif key == "n_candidates":
                    n_candidates.append(int(value))

                elif key == "n_spanning_reads":
                    n_spanning_reads.append(int(value))

                elif key == "time_elapsed_minigraph_s":
                    time_elapsed_minigraph_s.append(int(value))

                elif key == "time_elapsed_minimap_s":
                    time_elapsed_minimap_s.append(int(value))

                elif key == "time_elapsed_optimizer_s":
                    time_elapsed_optimizer_s.append(int(value))

                elif key == "time_elapsed_total_s":
                    total = int(value)
                    time_elapsed_total_s.append(total)
                    t_other = total - time_elapsed_minigraph_s[-1] - time_elapsed_minimap_s[-1] - time_elapsed_optimizer_s[-1]
                    time_elapsed_other_s.append(int(t_other))

    n_parameters = [x*y for x,y in zip(n_candidates,n_spanning_reads)]

    print(names[:20])
    print(n_candidates[:20])
    print(n_spanning_reads[:20])
    print(n_parameters[:20])
    print(time_elapsed_optimizer_s[:20])

    y = time_elapsed_optimizer_s
    x = n_valid_pairs

    fig = pyplot.figure(layout='constrained',figsize=[10,10])

    ax = fig.add_gridspec(top=0.75, right=0.75).subplots()

    ax_histx = ax.inset_axes([0, 1.15, 1, 0.25],sharex=ax)
    ax_histy = ax.inset_axes([1.15, 0, 0.25, 1],sharey=ax)

    ax.scatter(x=x,y=y,alpha=0.6,edgecolors="none")

    x_max = int(ax.get_xlim()[-1])
    x_width = x_max/100

    y_max = int(ax.get_ylim()[-1])
    y_width = y_max/100

    x_bins = numpy.arange(0, x_max + x_width, x_width, dtype=float)
    y_bins = numpy.arange(0, y_max + y_width, y_width, dtype=float)

    ax_histx.hist(x,bins=x_bins)
    ax_histy.hist(y,bins=y_bins,orientation='horizontal')

    ax.set_xlabel("# of valid read/path pairs")
    ax.set_ylabel("Time (s)")

    median_time_total = median(y)
    mean_time_total = mean(y)
    mode_time_total = mode(y)

    print("Total samples:", len(y))
    print("Average time:", mean_time_total)
    print("Median time:", median_time_total)
    print("Mode time:", mode_time_total)

    pyplot.savefig("run_time.png",dpi=200)

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i","--input_dir",
        required=True,
        type=str,
        help="Output directory"
    )

    args = parser.parse_args()

    main(
        input_dir=args.input_dir
    )
