from modules.Bam import download_regions_of_bam, index_bam
from modules.Authenticator import *
from modules.GsUri import *

from multiprocessing import Pool
import subprocess
import argparse
import re


def run_sniffles(ref_path, output_dir, bam_path, n_threads, timeout=60*3):
    output_filename = os.path.basename(bam_path).replace(".bam", "_sniffles.vcf")
    output_path = os.path.join(output_dir, output_filename)

    # Enable caching by path name
    if os.path.exists(output_path):
        return output_path

    # sniffles \
    # --input /home/ryan/data/test_hapslap/hg00733_1fc_chr20.bam \
    # --vcf /home/ryan/data/test_hapslap/hg00733_1fc_chr20_sniffles.vcf \
    # --reference /home/ryan/data/human/reference/chm13v2.0.fa \
    # --threads 30 \
    # --output-rnames
    args = [
        "sniffles",
        "--input", bam_path,
        "--vcf", output_path,
        "--reference", ref_path,
        "--threads", str(n_threads),
        "--output-rnames"
    ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return output_path


def compress_and_index_vcf(vcf_path, timeout=60*3):
    output_vcf_path = vcf_path + ".gz"
    output_tbi_path = output_vcf_path + ".tbi"

    # Enable caching by path name
    if os.path.exists(output_vcf_path) and os.path.exists(output_tbi_path):
        return output_vcf_path

    args = ["bgzip", vcf_path]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    args = ["tabix", "-p", "vcf", output_vcf_path]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return output_vcf_path


def main(
        tsv_path,
        column_names,
        ref_path,
        regions,
        n_threads,
        output_directory):

    output_directory = os.path.abspath(output_directory)
    cache_directory = os.path.join(output_directory, "input_bams")

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    tokenator = GoogleToken()
    tokenator.update_environment()

    bam_paths = download_regions_of_bam(
        regions=regions,
        tsv_path=tsv_path,
        column_names=column_names,
        output_directory=cache_directory,
        n_threads=n_threads,
    )

    for p in bam_paths:
        index_bam(p,n_threads)

    output_subdirectory = os.path.join(output_directory, "vcfs")
    for path in bam_paths:
            # TODO: add tandem bed path argument for sniffles
            vcf_path = run_sniffles(ref_path=ref_path, output_dir=output_subdirectory, bam_path=path, n_threads=n_threads)
            indexed_vcf_path = compress_and_index_vcf(vcf_path=vcf_path)

    return


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
        "-c","--column_names",
        required=False,
        default="'bam_vs_chm13'",
        type=parse_comma_separated_string,
        help="Which are the relevant columns of the Terra-style TSV (see tsv_path help msg) at the moment"
    )

    args = parser.parse_args()

    main(
        n_threads=args.threads,
        tsv_path=args.tsv,
        column_names=args.column_names,
        ref_path=args.ref,
        output_directory=args.output_dir,
        regions=args.region
    )
