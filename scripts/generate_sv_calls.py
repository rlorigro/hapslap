from modules.Authenticator import *
from modules.GsUri import *
from multiprocessing import Pool
import subprocess


def download_bam(output_directory, bam_path, region_string, tokenator, index=True, timeout=60*20):
    sample_name = os.path.basename(bam_path).split('.')[0]
    local_bam_filename = sample_name + "_" + region_string.replace(":","_") + ".bam"
    local_bam_path = os.path.join(output_directory,local_bam_filename)

    # Don't re-download
    if os.path.exists(local_bam_path):
        return local_bam_path

    tokenator.update_environment()

    # sniffles \
    # --input /home/ryan/data/test_hapslap/hg00733_1fc_chr20.bam \
    # --vcf /home/ryan/data/test_hapslap/hg00733_1fc_chr20_sniffles.vcf \
    # --reference /home/ryan/data/human/reference/chm13v2.0.fa \
    # --threads 30 \
    # --output-rnames
    args = [
        "samtools",
        "view",
        "-bh",
        "-o", local_bam_path,
        bam_path,
        region_string
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

    if index:
        samtools_index_args = ["samtools", "index", local_bam_path]
        sys.stderr.write(" ".join(samtools_index_args)+'\n')

        try:
            p1 = subprocess.run(samtools_index_args, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return local_bam_path


def run_sniffles(ref_path, bam_path, n_threads, timeout=60*3):
    output_path = bam_path.replace(".bam", "_sniffles.vcf")

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


def main():
    data = [
        ("HG002","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/9ee7caa4-3d4c-492f-8ce5-fd1f36a513ae/call-alignAndSortBAM/cacheCopy/HG002.bam"),
        ("HG00438","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/18f8defc-6bc9-4aaf-8844-6c26349b54e9/call-alignAndSortBAM/cacheCopy/HG00438.bam"),
        ("HG005","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/61ead59f-0665-40e9-ac01-baf6fc4d2878/call-alignAndSortBAM/cacheCopy/HG005.bam"),
        ("HG00621","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/696d4743-1c6e-43ca-86d0-ce03027ab67d/call-alignAndSortBAM/cacheCopy/HG00621.bam"),
        ("HG00673","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/b7c96642-49fe-4ab4-a92c-d5ade67373a6/call-alignAndSortBAM/cacheCopy/HG00673.bam"),
        ("HG00733","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/b08b607d-aaf9-4f9e-8d5b-b228afb4fc05/call-alignAndSortBAM/cacheCopy/HG00733.bam"),
        ("HG00735","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/97dc767b-9f61-45f9-9ae2-7449f7dd4f32/call-alignAndSortBAM/cacheCopy/HG00735.bam"),
        ("HG00741","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/d0c46f48-9ae3-4aa6-8e86-f8de7c9e5293/call-alignAndSortBAM/cacheCopy/HG00741.bam"),
        ("HG01071","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/18affbe5-d1f3-43bd-b993-69d167e974b2/call-alignAndSortBAM/cacheCopy/HG01071.bam"),
        ("HG01106","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/a67e5ae8-7db5-4c97-aa89-4f86796d37be/call-alignAndSortBAM/cacheCopy/HG01106.bam"),
        ("HG01109","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/88abf65a-4fe8-4633-a328-c996a4747a66/call-alignAndSortBAM/cacheCopy/HG01109.bam"),
        ("HG01123","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/5fe9d177-975e-4969-88b2-1d666ab85461/call-alignAndSortBAM/cacheCopy/HG01123.bam"),
        ("HG01175","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/368d0ee4-787e-4dd2-bb91-ef47991d4faa/call-alignAndSortBAM/cacheCopy/HG01175.bam"),
        ("HG01243","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/8f2495ee-e407-4cfa-bbdb-3d072273f300/call-alignAndSortBAM/cacheCopy/HG01243.bam"),
        ("HG01258","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/1383d647-7a9f-4d74-a490-cee8c07419b4/call-alignAndSortBAM/cacheCopy/HG01258.bam"),
        ("HG01358","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/1b622d1c-ae90-416a-878f-111cc546cfca/call-alignAndSortBAM/cacheCopy/HG01358.bam"),
        ("HG01361","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/eed50abf-25c8-41c1-92ef-b3aa56eb7a7b/call-alignAndSortBAM/cacheCopy/HG01361.bam"),
        ("HG01891","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/0c2c93a1-15a7-408a-967b-8534459c1282/call-alignAndSortBAM/cacheCopy/HG01891.bam"),
        ("HG01928","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/ae4488e8-80dc-45e6-8105-b602c817bbdb/call-alignAndSortBAM/cacheCopy/HG01928.bam"),
        ("HG01952","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/ad5908d6-f6cb-4152-8ce8-972d232b709a/call-alignAndSortBAM/cacheCopy/HG01952.bam"),
        ("HG01978","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/348b80f7-b98d-4872-b506-f1afb17b953f/call-alignAndSortBAM/cacheCopy/HG01978.bam"),
        ("HG02055","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/ae7377f3-c332-4f5b-83dc-c509ea135c4f/call-alignAndSortBAM/cacheCopy/HG02055.bam"),
        ("HG02080","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/1431ea16-ed73-41f7-9033-68648994433a/call-alignAndSortBAM/cacheCopy/HG02080.bam"),
        ("HG02109","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/64d92792-96b5-4034-8ce9-967a357a6fd3/call-alignAndSortBAM/cacheCopy/HG02109.bam"),
        ("HG02145","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/6edec7eb-34c3-417b-bb0e-650ddf0bb7ca/call-alignAndSortBAM/cacheCopy/HG02145.bam"),
        ("HG02148","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/cba38e18-36dc-4e3b-8d5c-9fc605f727f6/call-alignAndSortBAM/cacheCopy/HG02148.bam"),
        ("HG02257","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/8f8a9a90-2c16-46d1-a642-5421e972be2c/call-alignAndSortBAM/cacheCopy/HG02257.bam"),
        ("HG02486","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/8e59944d-4df3-451e-9c60-ea14943d4b19/call-alignAndSortBAM/cacheCopy/HG02486.bam"),
        ("HG02559","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/6c6f3783-8cfc-4d3e-bb42-4dae8f021329/call-alignAndSortBAM/cacheCopy/HG02559.bam"),
        ("HG02572","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/0a6d8d3b-06d5-472d-8752-24567363c572/call-alignAndSortBAM/cacheCopy/HG02572.bam"),
        ("HG02622","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/71123442-ff82-45e9-a131-545be2fedfaa/call-alignAndSortBAM/cacheCopy/HG02622.bam"),
        ("HG02630","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/aac31522-0cf0-421d-9b0a-7c722450b16d/call-alignAndSortBAM/cacheCopy/HG02630.bam"),
        ("HG02717","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/71ed6298-285d-4f3b-97da-c066383d7767/call-alignAndSortBAM/cacheCopy/HG02717.bam"),
        ("HG02723","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/1d3b20b6-d5c5-4bfa-9c6f-55948484c7e6/call-alignAndSortBAM/cacheCopy/HG02723.bam"),
        ("HG02818","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/756fbbf4-329a-4ea7-b9b2-c7a927297143/call-alignAndSortBAM/cacheCopy/HG02818.bam"),
        ("HG02886","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/28afcd4f-dcaa-445d-9441-c60e936c42de/call-alignAndSortBAM/cacheCopy/HG02886.bam"),
        ("HG03098","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/b12faf1f-54c2-46b2-aef9-c2c9114b3c35/call-alignAndSortBAM/cacheCopy/HG03098.bam"),
        ("HG03453","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/3c7d8151-be23-45d0-a4ac-dc12dbd2f3e8/call-alignAndSortBAM/cacheCopy/HG03453.bam"),
        ("HG03486","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/cbfcdaa6-d980-4af4-ab25-8fc99ea6b976/call-alignAndSortBAM/cacheCopy/HG03486.bam"),
        ("HG03492","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/f1958600-b7b1-4c88-819c-213a526e1f6f/call-alignAndSortBAM/cacheCopy/HG03492.bam"),
        ("HG03516","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/6e6ef333-14cf-47d6-b9e8-02e80375bdb3/call-alignAndSortBAM/cacheCopy/HG03516.bam"),
        ("HG03540","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/98050da4-446d-4b11-9dd3-ec5c1e2eff7b/call-alignAndSortBAM/cacheCopy/HG03540.bam"),
        ("HG03579","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/455c9671-dad5-453c-81e1-2da37568b8fd/call-alignAndSortBAM/cacheCopy/HG03579.bam"),
        ("NA18906","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/34c1d44f-7af0-4858-8f4f-e1469a7c5c7f/call-alignAndSortBAM/cacheCopy/NA18906.bam"),
        ("NA19240","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/137782f4-67a9-479d-893a-9f2d020b376a/call-alignAndSortBAM/cacheCopy/NA19240.bam"),
        ("NA20129","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/d13f295a-f8e6-46a6-938c-d381666363a3/call-alignAndSortBAM/cacheCopy/NA20129.bam"),
        ("NA21309","gs://fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/9a6c205f-381e-4a27-b445-f27627b7ad44/minimap2/35ca7981-e1fb-4180-8c11-8a7dd0a88a6a/call-alignAndSortBAM/cacheCopy/NA21309.bam"),
    ]

    n_threads = 30

    output_directory = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data/test/hprc/")
    print(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    ref_path = "/home/ryan/data/human/reference/chm13v2.0.fa"

    tokenator = GoogleToken()
    tokenator.update_environment()

    region_string = "chr20"

    args = list()
    bam_results = None

    for name,path in data:
        args.append([output_directory, path, region_string, tokenator])

    with Pool(processes=n_threads) as pool:
        bam_results = pool.starmap(download_bam, args)

    for path in bam_results:
        vcf_path = run_sniffles(ref_path=ref_path, bam_path=path, n_threads=n_threads)
        indexed_vcf_path = compress_and_index_vcf(vcf_path=vcf_path)

    return


if __name__ == "__main__":
    main()
