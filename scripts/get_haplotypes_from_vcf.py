from collections import defaultdict
import os.path

from vcf import VCFReader
from pysam import FastaFile
import argparse


class Allele:
    def __init__(self, start, stop, sequence):
        self.start = int(start)
        self.stop = int(stop)
        self.sequence = sequence

    def __str__(self):
        return "[%s,%s]\t%s" % (self.start, self.stop, self.sequence)

    def __hash__(self):
        return hash((self.start, self.stop, self.sequence))

    def __eq__(self, other):
        return (self.start, self.stop, self.sequence) == (other.start, other.stop, other.sequence)


def main():
    ref_path = "/home/ryan/data/human/reference/chm13v2.0.fa"

    vcf_paths = [
        # "/home/ryan/code/hapslap/data/test/hprc/HG002_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG00438_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG005_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG00621_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00673_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00733_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00735_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG00741_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01071_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01106_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01109_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01123_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01175_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01243_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01258_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01358_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01361_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01891_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01928_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01952_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG01978_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02055_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02080_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02109_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02145_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02148_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02257_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02486_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02559_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02572_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02622_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02630_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02717_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02723_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02818_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG02886_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03098_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03453_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03486_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03492_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03516_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03540_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/HG03579_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA18906_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA19240_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA20129_chr20_sniffles.vcf.gz",
        # "/home/ryan/code/hapslap/data/test/hprc/NA21309_chr20_sniffles.vcf.gz"
    ]

    chromosome = "chr20"
    region_start = 47474020
    region_stop = 47477018

    region_string = chromosome + ":" + str(region_start) + "-" + str(region_stop)

    gfa_path = region_string + ".gfa"

    ref_sequence = FastaFile(ref_path).fetch("chr20")

    alleles = defaultdict(set)

    for vcf_path in vcf_paths:
        sample_name = os.path.basename(vcf_path).split("_")[0]

        with open(vcf_path, 'rb') as file:
            vcf = VCFReader(file)

            records = vcf.fetch(region_string)

            for record in records:
                print()
                print("var_subtype:\t%s" % record.var_subtype)
                print("start:\t\t%d" % record.start)

                # One VCF per sample, no iterating of samples needed
                call = record.samples[0]
                gt = [int(call.data.GT[0]), int(call.data.GT[-1])]
                print("gt:\t\t%d/%d" % (gt[0], gt[1]))
                print("ref_length:\t%d" % len(record.alleles[0]))
                print("a_length:\t%d" % len(record.alleles[gt[0]]))
                print("b_length:\t%d" % len(record.alleles[gt[1]]))
                print("is_sv_precise:\t%d" % record.is_sv_precise)

                # Iterate unique, non-ref alleles only
                for allele_index in set(gt):
                    print(record.alleles[allele_index])

                    if allele_index != 0:
                        l = len(record.alleles[0]) if record.alleles[0] != 'N' else 0

                        # We are sorting by ref coordinates, so we use the ref allele start/stop to keep track of
                        # where the alt allele will be substituted
                        start = int(record.start)
                        stop = start + l
                        sequence = str(record.alleles[allele_index])
                        sequence = sequence if sequence != 'N' else ''

                        # Collapse identical alleles by hashing them as a fn of start,stop,sequence
                        # But keep track of which samples are collapsed together
                        a = Allele(start, stop, sequence)
                        alleles[a].add(sample_name)

    print()

    alleles = list(alleles.items())

    # sample_to_allele_index = dict()
    #
    # # Build reverse mapping of sample -> allele_index
    # for a,[allele,samples] in enumerate(alleles):
    #     for sample in samples:
    #         sample_to_allele_index[sample] = a
    #
    # # Generate pointers to the allele data and sort them by allele start coord
    # allele_indexes = list(range(len(alleles)))
    # allele_indexes = sorted(allele_indexes, key=lambda x: alleles[x][0].start)
    #
    # overlapping_allele_indexes = list()
    #
    # # Do a sweep over the sorted allele intervals to find overlapping ones
    # stop = -1
    # for i in allele_indexes:
    #     if alleles[i][0].start <= stop:
    #         overlapping_allele_indexes[-1].add(i)
    #
    #         if alleles[i][0].stop > stop:
    #             stop = alleles[i][0].stop
    #     else:
    #         overlapping_allele_indexes.append({i})
    #         stop = alleles[i][0].stop
    #
    # # Print results
    # for i,item in enumerate(overlapping_allele_indexes):
    #     print(i)
    #     for allele_index in item:
    #         allele = alleles[allele_index][0]
    #         sample = alleles[allele_index][1]
    #         print(allele, sample)

    # Construct a list of coordinates along the reference path which contain edges to VCF alleles
    ref_edges = defaultdict(lambda: [[],[]])

    for a,[allele,_] in enumerate(alleles):
        print(a)
        print(allele)
        ref_edges[allele.start][1].append(a)
        ref_edges[allele.stop][0].append(a)

    # Append dummy item at end of list to make one-pass iteration easier
    ref_edges[region_stop] = [[],[]]

    # Sort the list by coord so that it can be iterated from left to right
    ref_edges = list(sorted(ref_edges.items(), key=lambda x: x[0]))

    # Save all the edges so they can be written separately at the end
    gfa_edge_lines = list()

    # Construct GFA
    with open(gfa_path, 'w') as file:
        # Initialize vars that will be iteration dependent
        prev_coord = region_start
        in_edges = []

        # First generate nodes for all the known VCF alleles
        for allele_index,[allele,samples] in enumerate(alleles):
            name = str(allele_index)
            file.write("S\t%s\t%s\n" % (name,allele.sequence))

        r = 0
        for i,[coord,edges] in enumerate(ref_edges):
            id = "ref_" + str(r)

            for allele_index in in_edges:
                name = str(allele_index)
                gfa_edge_lines.append("L\t%s\t+\t%s\t+\t0M\n" % (name,id))

            for allele_index in edges[1]:
                name = str(allele_index)
                gfa_edge_lines.append("L\t%s\t+\t%s\t+\t0M\n" % (id,name))

            print(id, coord, edges)

            if i < len(ref_edges) - 1:
                # Write ref sequence to GFA
                sequence = ref_sequence[prev_coord:coord]

                # Add edge to next ref sequence
                next_id = "ref_" + str(r+1)
                file.write("S\t%s\t%s\n" % (id,sequence))
                gfa_edge_lines.append("L\t%s\t+\t%s\t+\t0M\n" % (id,next_id))

                prev_coord = coord
                r += 1
            else:
                sequence = ref_sequence[prev_coord:coord]
                file.write("S\t%s\t%s\n" % (id,sequence))

            in_edges = edges[0]

        for line in gfa_edge_lines:
            file.write(line)

    return


if __name__ == "__main__":
    main()

