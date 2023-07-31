from collections import defaultdict
import sys
import os
import re

from vcf import VCFReader

from modules.IntervalGraph import IntervalGraph
from modules.Bed import iter_bed_items


def parse_region_string(s):
    tokens = re.split(r'[{:\-}]+', s.strip())
    chromosome = tokens[0]
    start = int(tokens[1])
    stop = int(tokens[2])

    return chromosome, start, stop


def get_overlap_sites(output_dir, region_string, vcfs_per_sample, bed_path, bed_name, padding, max_interval_length):
    # The IntervalGraph data structure has no concept of chromosomes, only coordinates, so chr must be factored out
    intervals_per_chromosome = defaultdict(list)

    ref_start = None
    ref_stop = None

    if ref_start is not None and ref_stop is not None:
        region_string += ":" + str(ref_start) + "-" + str(ref_stop)

    for sample_name in vcfs_per_sample:
        vcf_path = vcfs_per_sample[sample_name]

        with open(vcf_path, 'rb') as file:
            vcf = VCFReader(file)

            try:
                records = vcf.fetch(region_string)
            except ValueError:
                print("WARNING: Skipping sample because has no VCF entries: " + sample_name)
                continue

            for record in records:
                chromosome = record.CHROM
                start = record.affected_start
                stop = record.affected_end

                if 'END' in record.INFO:
                    stop = record.sv_end

                if stop - start > max_interval_length:
                    continue

                intervals_per_chromosome[chromosome].append((start, stop+padding, sample_name))

    if bed_path is not None:
        for chr,start,stop in iter_bed_items(bed_path):
            intervals_per_chromosome[chr].append((start, stop+padding, bed_name))

    output_path = os.path.join(output_dir,"sites.bed")

    with open(output_path, 'w') as file:
        for chromosome,intervals in intervals_per_chromosome.items():
            graph = IntervalGraph(intervals)

            components = graph.get_connected_components()

            sizes = list()
            n_samples = list()

            for c,component in enumerate(components):
                s = len(component)

                samples = set()

                left_coord = sys.maxsize
                right_coord = -1

                for x in component:
                    samples = samples.union(graph.graph[x].values)

                    if x[0] < left_coord:
                        left_coord = x[0]

                    if x[1] > right_coord:
                        right_coord = x[1]

                # Remove the bed name from the set of samples which have overlapping intervals in this component,
                if bed_name in samples:
                    samples.remove(bed_name)
                    s -= 1

                right_coord -= padding

                n = len(samples)

                sizes.append(s)
                n_samples.append(n)

                # Only include regions that have at least one real interval (not from the BED file)
                if len(samples) >= 1 and right_coord - left_coord < max_interval_length:
                    file.write('\t'.join([chromosome,str(left_coord),str(right_coord)]))
                    file.write('\n')

        # pyplot.figure()
        # axes = pyplot.axes()
        #
        # n_max = max(n_samples)
        # size_max = max(sizes)
        #
        # axes.hist2d(x=n_samples,y=sizes,bins=(n_max,size_max),norm=colors.LogNorm())
        # axes.set_xlabel("# samples")
        # axes.set_ylabel("Component size (# overlapping intervals)")
        #
        # pyplot.show()
        # pyplot.close()

    return output_path


def test_locally():
    region_string = "chr20"

    input_directory = "/home/ryan/code/hapslap/data/test/hprc/"
    bed_path = "/home/ryan/data/human/reference/human_chm13v2.0_maskedY_rCRS.trf.bed"
    bed_name = "tandems"

    padding = 150
    max_interval_length = 20000

    sample_names = [
        "HG002",
        "HG00438",
        "HG005",
        "HG00621",
        "HG00673",
        "HG00733",
        "HG00735",
        "HG00741",
        "HG01071",
        "HG01106",
        "HG01109",
        "HG01123",
        "HG01175",
        "HG01243",
        "HG01258",
        "HG01358",
        "HG01361",
        "HG01891",
        "HG01928",
        "HG01952",
        "HG01978",
        "HG02055",
        "HG02080",
        "HG02109",
        "HG02145",
        "HG02148",
        "HG02257",
        "HG02486",
        "HG02559",
        "HG02572",
        "HG02622",
        "HG02630",
        "HG02717",
        "HG02723",
        "HG02818",
        "HG02886",
        "HG03098",
        "HG03453",
        "HG03486",
        "HG03492",
        "HG03516",
        "HG03540",
        "HG03579",
        "NA18906",
        "NA19240",
        "NA20129",
        "NA21309"
    ]

    data_per_sample = dict()

    for filename in os.listdir(input_directory):
        path = os.path.join(input_directory, filename)

        for name in sample_names:
            if name in filename:
                if filename.endswith(".vcf.gz"):
                    data_per_sample[name] = path

    output_dir = "/home/ryan/data/test_hapslap/test_regions/test_refactor"

    get_overlap_sites(
        output_dir,
        region_string,
        data_per_sample,
        bed_path,
        bed_name,
        padding,
        max_interval_length
    )


if __name__ == "__main__":
    test_locally()
