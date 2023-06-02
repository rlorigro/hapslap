import subprocess
import sys
import os


def index_bam(bam_path, n_threads):
    args = [
        "samtools",
        "index",
        "-@", str(n_threads),
        bam_path,
    ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return


def get_region_from_bam(output_directory, bam_path, region_string, tokenator, timeout=60*20, output_filename=None, remove_local_index=False):
    if output_filename is None:
        prefix = os.path.basename(bam_path).split('.')[0]
        output_filename = prefix + "_" + region_string.replace(":","_") + ".bam"

    local_bam_path = os.path.join(output_directory, output_filename)

    # Enable caching by path name
    if os.path.exists(local_bam_path):
        return local_bam_path

    tokenator.update_environment()

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

        # Always remove the local index if there is an error
        os.remove(os.path.basename(bam_path) + ".bai")

        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()

        # Always remove the local index if there is an error
        os.remove(os.path.basename(bam_path) + ".bai")

        return None

    if remove_local_index:
        os.remove(os.path.basename(bam_path) + ".bai")

    return local_bam_path

