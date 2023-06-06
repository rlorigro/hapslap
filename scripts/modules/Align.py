import time
import sys

from modules.GsUri import *

import subprocess
import os.path


# Requires samtools installed!
def run_minimap2(ref_fasta_path, reads_fasta_path, preset, n_threads, n_sort_threads, output_directory, filename_prefix="reads_vs_ref"):
    output_filename = os.path.join(output_directory, filename_prefix + ".bam")
    output_path = os.path.join(output_directory,output_filename)

    minimap_args = ["minimap2", "-a", "-x", preset, "--eqx", "-t", str(n_threads), ref_fasta_path, reads_fasta_path]
    sort_args = ["samtools", "sort", "-", "-@", str(n_sort_threads), "-o", output_filename]

    index_args = ["samtools", "index", output_filename]

    with open(output_path, 'w') as file:
        sys.stderr.write(" ".join(minimap_args)+'\n')
        sys.stderr.write(" ".join(sort_args)+'\n')

        p1 = subprocess.Popen(minimap_args, stdout=subprocess.PIPE, cwd=output_directory)
        p2 = subprocess.Popen(sort_args, stdin=p1.stdout, stdout=file, cwd=output_directory)
        p2.communicate()

    success = (p2.returncode == 0)

    if not success:
        sys.stderr.write("ERROR: failed to align: %s to %s\n" % (reads_fasta_path, ref_fasta_path))
        sys.stderr.flush()
        return None

    sys.stderr.write(" ".join(index_args)+'\n')

    try:
        p1 = subprocess.run(index_args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    return output_path
