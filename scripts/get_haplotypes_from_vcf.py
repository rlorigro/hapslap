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
        "/home/ryan/code/hapslap/data/test/hprc/HG002_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG005_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG00438_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG00621_chr20_sniffles.vcf.gz",
        "/home/ryan/code/hapslap/data/test/hprc/HG01243_chr20_sniffles.vcf.gz"
    ]

    chromosome = "chr20"
    region_start = 47474020
    region_stop = 47477018

    region_string = chromosome + ":" + str(region_start) + "-" + str(region_stop)

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

                # Iterate unique, non-ref alleles only
                for allele_index in set(gt):
                    print(record.alleles[allele_index])

                    if allele_index != 0:
                        l = len(record.alleles[allele_index])

                        # We are sorting by ref coordinates, so we use the ref allele start/stop to keep track of
                        # where the alt allele will be substituted
                        start = int(record.start)
                        stop = start + len(record.alleles[0])
                        sequence = str(record.alleles[allele_index])

                        # Collapse identical alleles by hashing them as a fn of start,stop,sequence
                        a = Allele(start, stop, sequence)
                        alleles[a].add(sample_name)

    print()

    alleles = list(alleles.items())
    overlapping_allele_indexes = list()

    # Generate pointers to the allele data and sort them by allele start coord
    allele_indexes = list(range(len(alleles)))
    allele_indexes = sorted(allele_indexes, key=lambda x: alleles[x][0].start)

    # Do a sweep over the sorted allele intervals to find overlapping ones
    stop = -1
    for i in allele_indexes:
        if alleles[i][0].start <= stop:
            overlapping_allele_indexes[-1].add(i)

            if alleles[i][0].stop > stop:
                stop = alleles[i][0].stop
        else:
            overlapping_allele_indexes.append({i})
            stop = alleles[i][0].stop

    # # Print results
    # for i,item in enumerate(overlapping_allele_indexes):
    #     print(i)
    #     for allele_index in item:
    #         allele = alleles[allele_index][0]
    #         sample = alleles[allele_index][1]
    #         print(allele, sample)

    ref_edges = defaultdict(lambda: [[],[]])

    for a,[allele,_] in enumerate(alleles):
        print(a)
        print(allele)

        ref_edges[allele.start][1].append(a)
        ref_edges[allele.stop][0].append(a)

    ref_edges[region_stop] = [[],[]]

    gfa_path = region_string + ".gfa"

    gfa_edge_lines = list()

    with open(gfa_path, 'w') as file:
        # Construct GFA
        prev_coord = region_start

        for a,[allele,samples] in enumerate(alleles):
            file.write("S\t%s\t%s\n" % ('_'.join([str(a)] + list(samples)),allele.sequence))

        ref_edges = list(sorted(ref_edges.items(), key=lambda x: x[0]))
        in_edges = []

        r = 0
        for i,[coord,edges] in enumerate(ref_edges):
            id = "ref_" + str(r)

            for allele_index in in_edges:
                samples = alleles[allele_index][1]
                name = '_'.join([str(allele_index)] + list(samples))
                gfa_edge_lines.append("L\t%s\t+\t%s\t+\t0M\n" % (name,id))

            for allele_index in edges[1]:
                samples = alleles[allele_index][1]
                name = '_'.join([str(allele_index)] + list(samples))
                gfa_edge_lines.append("L\t%s\t+\t%s\t+\t0M\n" % (id,name))

            print(id, coord, edges)

            if i < len(ref_edges) - 1:
                next_edges = ref_edges[i+1][1]

                # if len(edges[1]) > 0 or len(next_edges[0]) > 0:
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

