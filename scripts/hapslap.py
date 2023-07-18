from get_overlap_sites import *
from generate_sv_calls import *
from infer_haplotypes import *

import argparse
import re


def infer_haplotypes_for_all_regions(data_per_sample, bed_path, flank_length, min_coverage, max_path_to_read_cost, ref_path, output_directory):
    regions = parse_bed_regions(bed_path)

    summary_strings = list()
    for region in regions:
        summary_string = infer_haplotypes(
            ref_path=ref_path,
            data_per_sample=data_per_sample,
            chromosome=region.contig_name,
            ref_start=region.start,
            ref_stop=region.stop,
            flank_length=flank_length,
            min_coverage=min_coverage,
            max_path_to_read_cost=max_path_to_read_cost,
            output_directory=output_directory
        )

        summary_strings.append(summary_string)

    # Group all the results into one message at the end of the output
    for s in range(len(summary_strings)):
        sys.stderr.write(str(regions[s]))
        sys.stderr.write('\n')
        sys.stderr.write(summary_strings[s])


def localize_bams_per_sample(cache_directory, regions, tsv_path, column_names, n_threads):
    tokenator = GoogleToken()
    tokenator.update_environment()

    # Terminal directory names "should" be guaranteed to be sample names in the output of this operation
    bam_paths_per_sample = download_regions_of_bam(
        regions=regions,
        tsv_path=tsv_path,
        column_names=column_names,
        output_directory=cache_directory,
        n_threads=n_threads,
        as_dict=True
    )

    for path in bam_paths_per_sample.values():
        # index BAMs but don't over-parallelize it (anecdotally worsens performance?)
        index_bam(path,n_threads=min(16,n_threads))

    return bam_paths_per_sample


def generate_sv_calls(
        bams_per_sample,
        ref_path,
        n_threads,
        output_directory):

    output_subdirectory = os.path.join(output_directory, "vcfs")

    vcf_paths_per_sample = dict()
    for sample,path in bams_per_sample.items():
        # TODO: add tandem bed path argument for sniffles
        vcf_path = run_sniffles(ref_path=ref_path, output_dir=output_subdirectory, bam_path=path, n_threads=n_threads)
        indexed_vcf_path = compress_and_index_vcf(vcf_path=vcf_path)
        vcf_paths_per_sample[sample] = indexed_vcf_path

    return vcf_paths_per_sample


def run_hapslap(n_threads, tsv_path, column_name, ref_path, region_string, interval_bed_path, output_directory, cache_directory):
    interval_padding = 150
    max_interval_length = 15000
    flank_length = 5000
    min_coverage = 2
    max_path_to_read_cost = 250

    bams_per_sample = localize_bams_per_sample(
        cache_directory,
        region_string,
        tsv_path,
        column_name,
        n_threads
    )

    vcfs_per_sample = generate_sv_calls(
        bams_per_sample=bams_per_sample,
        ref_path=ref_path,
        n_threads=n_threads,
        output_directory=output_directory
    )

    data_per_sample = defaultdict(dict)
    for sample in bams_per_sample.keys():
        data_per_sample[sample]["vcf"] = vcfs_per_sample[sample]
        data_per_sample[sample]["bam"] = bams_per_sample[sample]

    bed_path = get_overlap_sites(
        output_dir=output_directory,
        region_string=region_string,
        vcfs_per_sample=vcfs_per_sample,
        bed_path=interval_bed_path,
        bed_name="tandems",
        padding=interval_padding,
        max_interval_length=max_interval_length)

    infer_haplotypes_for_all_regions(
        data_per_sample=data_per_sample,
        bed_path=bed_path,
        flank_length=flank_length,
        min_coverage=min_coverage,
        max_path_to_read_cost=max_path_to_read_cost,
        ref_path=ref_path,
        output_directory=output_directory)


def parse_comma_separated_string(s):
    return re.split(r'[{\'\",}]+', s.strip("\"\'[]{}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ref",
        required=True,
        type=str,
        help="Path to reference which reads are aligned to (e.g. the chm13 v2.0 reference)"
    )

    parser.add_argument(
        "--bed",
        required=False,
        default=None,
        type=str,
        help="Path to BED file which should be used to merge overlapping VCF intervals (e.g. tandem regions)"
    )

    parser.add_argument(
        "--region",
        required=True,
        type=str,
        help="Any valid samtools region string (e.g. 'chr20' or 'chr1 chr2' or 'chr3:5000-10000')"
    )

    parser.add_argument(
        "-o","--output_dir",
        required=True,
        type=str,
        help="Output directory which will be created (and must not exist)"
    )

    parser.add_argument(
        "-t","--threads",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use (recommended 15-30)"
    )

    parser.add_argument(
        "--cache_dir",
        required=True,
        type=str,
        help="Directory where remote BAMs will be stored, and potentially reused for future runs of this script"
    )

    parser.add_argument(
        "--tsv",
        required=True,
        type=str,
        help="Path to a Terra-style haplotype TSV which contains relevant lookups for aligned hap1/hap2 ref assemblies per sample"
    )

    parser.add_argument(
        "-c","--column_name",
        required=False,
        default="'bam_vs_chm13'",
        type=parse_comma_separated_string,
        help="Which is the relevant column of the Terra-style TSV which contains the BAM of aligned reads to the reference"
    )

    args = parser.parse_args()

    run_hapslap(
        n_threads=args.threads,
        tsv_path=args.tsv,
        column_name=args.column_name,
        ref_path=args.ref,
        output_directory=args.output_dir,
        cache_directory=args.cache_dir,
        region_string=args.region,
        interval_bed_path=args.bed,
    )
