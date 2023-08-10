from modules.Vcf import merge_vcfs,bcftools_view,compress_and_index_vcf
from modules.Bed import parse_bed_regions
import os


def main():
    bed_path = "/home/ryan/data/test_hapslap/test_regions/automated/sites_100.bed"

    vcf_paths = [
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02818_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02559_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG00733_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02257_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG005_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG002_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02486_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG00438_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/NA20129_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01123_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01175_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG00673_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG00621_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01071_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02055_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01258_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG00735_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01928_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01952_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01891_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01978_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01358_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02145_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG00741_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02080_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03492_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02886_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02148_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01361_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/NA21309_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/NA19240_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03486_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02109_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02622_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01106_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/NA18906_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02717_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02630_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02572_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01243_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03516_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG01109_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03579_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03098_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03453_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG02723_bam_vs_chm13_chr20_sniffles.vcf.gz",
        "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/HG03540_bam_vs_chm13_chr20_sniffles.vcf.gz"
    ]

    output_directory = "/home/ryan/data/test_hapslap/results/competitors/sniffles_default/merged"
    output_filename = "hprc_47.vcf.gz"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    merged_vcf_path = os.path.join(output_directory, output_filename)
    merged_vcf_path = merge_vcfs(vcf_paths=vcf_paths, output_path=merged_vcf_path, force_samples=True)
    merged_vcf_path = compress_and_index_vcf(vcf_path=merged_vcf_path)

    regions = parse_bed_regions(bed_path)

    output_subdirectory = os.path.join(output_directory, "regions")

    if not os.path.exists(output_subdirectory):
        os.makedirs(output_subdirectory)

    for region in regions:
        output_path = os.path.join(output_subdirectory, str(region) + ".vcf.gz")
        bcftools_view(vcf_path=merged_vcf_path, region_string=region.to_samtools_string(), output_path=output_path)


if __name__ == "__main__":
    main()
