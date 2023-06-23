from collections import defaultdict
import sys
import os

from matplotlib import pyplot,colors
from vcf import VCFReader

from modules.IntervalGraph import IntervalGraph
from modules.Vcf import get_alleles_from_vcfs


def iter_bed_items(path):
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split()
            contig = tokens[0]
            start = int(tokens[1])
            stop = int(tokens[2])

            yield contig, start, stop


def main():
    output_dir = "/home/ryan/data/test_hapslap/test_regions/automated"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        exit("ERROR: output directory exists already: %s" % output_dir)

    chromosome = "chr20"

    input_directory = "/home/ryan/code/hapslap/data/test/hprc/"
    ref_path = "/home/ryan/data/human/reference/chm13v2.0.fa"
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

    data_per_sample = defaultdict(dict)

    for filename in os.listdir(input_directory):
        path = os.path.join(input_directory, filename)

        for name in sample_names:
            if name in filename:
                if filename.endswith(".bam"):
                    data_per_sample[name]["bam"] = path
                if filename.endswith(".vcf.gz"):
                    data_per_sample[name]["vcf"] = path

    vcf_paths = list()
    bam_paths = list()

    for name in data_per_sample:
        has_vcf = False
        has_bam = False

        if "vcf" in data_per_sample[name]:
            has_vcf = True
            vcf_paths.append(data_per_sample[name]["vcf"])
        if "bam" in data_per_sample[name]:
            has_bam = True
            bam_paths.append(data_per_sample[name]["bam"])

        if (not has_vcf) or (not has_bam):
            exit("ERROR: sample does not have both a BAM and a VCF")

        # print(name)
        # print("vcf", data_per_sample[name]["vcf"])
        # print("bam", data_per_sample[name]["bam"])

    intervals = list()

    ref_start = None
    ref_stop = None

    # FOR DEBUGGING PURPOSES ONLY, WILL NOT PROPERLY INCLUDE PARTIALLY CONTAINED INTERVALS
    # ref_start = 10351130
    # ref_stop = 10351874

    region_string = chromosome

    if ref_start is not None and ref_stop is not None:
        region_string += ":" + str(ref_start) + "-" + str(ref_stop)

    for vcf_path in vcf_paths:
        sample_name = os.path.basename(vcf_path).split("_")[0]

        with open(vcf_path, 'rb') as file:
            vcf = VCFReader(file)

            records = vcf.fetch(region_string)

            for record in records:
                # print()
                # print("var_subtype:\t%s" % record.var_subtype)
                # print("start:\t\t%d" % record.start)
                # print("affected_start:\t%d" % record.affected_start)
                # print("affected_end:\t%d" % record.affected_end)
                # # print("sv_end:\t\t%d" % record.sv_end)
                #
                # # One VCF per sample, no iterating of samples needed
                # call = record.samples[0]
                # gt = [int(call.data.GT[0]) if call.data.GT[0] != '.' else 0, int(call.data.GT[-1]) if call.data.GT[-1] != '.' else 0]
                #
                # print("gt:\t\t%d/%d" % (gt[0], gt[1]))
                # print("ref_length:\t%d" % len(record.alleles[0]))
                # print("var_type:\t%s" % record.var_type)
                # print("is_sv_precise:\t%d" % record.is_sv_precise)

                start = record.affected_start
                stop = record.affected_end

                if 'END' in record.INFO:
                    stop = record.sv_end

                if stop - start > max_interval_length:
                    continue

                intervals.append((start, stop + padding, sample_name))

    # for interval in intervals:
    #     print(interval)

    bed_intervals = [(start,stop+padding,bed_name) for chr,start,stop in iter_bed_items(bed_path) if chr == chromosome] # and start >= ref_start and stop <= ref_stop]

    intervals += bed_intervals

    graph = IntervalGraph(intervals)

    components = graph.get_connected_components()

    sizes = list()
    n_samples = list()

    output_path = os.path.join(output_dir,"sites.bed")

    with open(output_path, 'w') as file:
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

            if bed_name in samples:
                samples.remove(bed_name)
                s -= 1

            right_coord -= padding

            n = len(samples)

            sizes.append(s)
            n_samples.append(n)

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


if __name__ == "__main__":
    main()
