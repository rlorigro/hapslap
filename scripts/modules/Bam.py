import subprocess
import sys
import os
import re


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


def get_region_from_bam(output_directory, bam_path, region_string, tokenator, timeout=60*20, output_filename=None):
    # samtools creates a local copy of the index file wherever the CWD is. But it is named not based on the path of
    # the online BAM, only on the filename. This leads to conflicts, where 2 different oneline BAMs have the same index,
    # and therefore causes a crash when fetching regions. The solution is to make a temp workdir that is unique to the
    # full online path of the BAM... As long as absolute paths are specified for input and output, the workdir does
    # not matter. This also means that future fetch operations can correctly locate and reuse the index if it has
    # been downloaded on the machine previously (since last clearing of the /tmp/ dir)

    # Replace all non-alphanumeric chars in the path with a '_'
    name = re.sub('[^0-9a-zA-Z]+', '_', bam_path)
    temp_working_dir = os.path.join("/tmp", name)

    if not os.path.exists(temp_working_dir):
        os.makedirs(temp_working_dir)

    if output_filename is None:
        prefix = os.path.basename(bam_path).split('.')[0]
        output_filename = prefix + "_" + region_string.replace(":","_") + ".bam"

    local_bam_path = os.path.join(output_directory, output_filename)
    index_path = bam_path + ".bai"

    # Enable caching by path name
    if os.path.exists(local_bam_path):
        print("exists! using cached file")
        return local_bam_path

    tokenator.update_environment()

    args = [
        "samtools",
        "view",
        "-bh",
        "-o", local_bam_path,
        "-X",
        bam_path,
        index_path,
        region_string
    ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, cwd=temp_working_dir, check=True, stderr=subprocess.PIPE, timeout=timeout)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()

        return None

    except subprocess.TimeoutExpired as e:
        sys.stderr.write("Status: FAIL due to timeout " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()

        return None

    return local_bam_path

