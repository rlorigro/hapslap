from .Authenticator import GoogleToken
import pandas

from multiprocessing import Pool
import subprocess
import hashlib
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
    name = bam_path
    name = re.sub('[^0-9a-zA-Z]+', '_', name)
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


def download_regions_of_bam(regions, tsv_path, column_names, output_directory, n_threads, samples=None, as_dict=False):
    token = GoogleToken()

    output_directory = os.path.abspath(output_directory)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    df = pandas.read_table(tsv_path, sep='\t', header=0)

    n_rows, n_cols = df.shape

    if type(regions) == str:
        regions = regions.split(' ')

    elif type(regions) == set:
        regions = list(regions)

    # Samtools compatible formatting of regions
    region_string = ' '.join(regions)

    file_tag = region_string
    file_tag = file_tag.replace(' ','_').replace(':','_')

    # If too many regions were specified to make a sensible file path, use a deterministic hash identifier instead
    if len(file_tag) > 64:
        sha = hashlib.sha256()
        sha.update(region_string.encode())
        file_tag = sha.hexdigest()

        sys.stderr.write("Using substitute hash for regions:%s\n\tregions:%s\n" % (file_tag, region_string))

    args = list()
    result_samples = list()

    if samples is not None:
        samples = set(samples)

    for column_name in column_names:
        for i in range(n_rows):
            sample_name = df.iloc[i][0]

            if samples is not None and sample_name not in samples:
                continue

            print(sample_name)

            gs_uri = df.iloc[i][column_name]

            # Each sample downloads its regions to its own subdirectory to prevent overwriting (filenames are by region)
            output_subdirectory = os.path.join(output_directory, sample_name)

            if not os.path.exists(output_subdirectory):
                os.makedirs(output_subdirectory)

            filename = sample_name + "_" + column_name + "_" + file_tag + ".bam"

            print(gs_uri)

            args.append([output_subdirectory,gs_uri,region_string,token,600,filename])
            result_samples.append(sample_name)

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(get_region_from_bam, args)

    if as_dict:
        results = {a:b for a,b in zip(result_samples, results)}

    return results
