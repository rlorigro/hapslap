from pywfa import WavefrontAligner
import time
import sys

if __name__ != "__main__":
    from .Sequence import Sequence,iterate_fasta
    from .GsUri import *
else:
    from Sequence import Sequence,iterate_fasta
    from GsUri import *

from itertools import combinations
from multiprocessing import Pool
import tempfile
import subprocess
import os.path


def get_indel_distance_from_string(alignment: dict, minimum_indel_length=1):
    total_indels = 0

    length_token = ""
    for c in alignment["cigars"]:
        if c.isnumeric():
            length_token += c
        else:
            l = int(length_token)
            if l >= minimum_indel_length and (c == 'I' or c == 'D'):
                total_indels += l

    return total_indels


def get_indel_distance_from_tuples(cigar_tuples):
    total_indels = 0

    for c,l in cigar_tuples:
        if c == 1 or c == 2:
            total_indels += l

    return total_indels


def align_and_get_indel_distance(a:Sequence,b:Sequence,output_dir=None):
    if len(a) == 0 or len(b) == 0:
        return max(len(a),len(b))

    aligner = WavefrontAligner(heuristic="adaptive")

    aligner(a.sequence,b.sequence)

    if aligner.status != 0:
        return max(len(a),len(b))

    d = get_indel_distance_from_tuples(aligner.cigartuples)

    if output_dir is not None:
        output_path = os.path.join(output_dir, a.name + "_" + b.name + ".txt")
        aligner.cigar_print_pretty(output_path)

    return d


def cross_align_sequences(sequences:dict, n_threads,output_dir=None):
    pairs = list()
    lengths = list()
    args = list()

    for a,b in combinations(sequences.values(),2):
        if a == b:
            continue

        lengths.append((len(a),len(b)))
        pairs.append((a.name,b.name))
        args.append((a,b,output_dir))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align_and_get_indel_distance, args)

    return pairs,lengths,results


def align_sequences_to_other_sequences(a_seqs:dict, b_seqs:dict, n_threads:int, output_dir:str=None):
    pairs = list()
    lengths = list()
    args = list()

    for a_name,a in a_seqs.items():
        for b_name,b in b_seqs.items():
            lengths.append((len(a),len(b)))
            pairs.append((a,b))
            args.append((a,b,output_dir))

    results = None
    with Pool(n_threads) as pool:
        results = pool.starmap(align_and_get_indel_distance, args)

    return pairs,lengths,results


def run_minigraph(output_directory, gfa_path, fasta_path, n_threads, args_override=None):
    output_path = os.path.join(output_directory, "reads_vs_graph.gaf")

    # minigraph \
    # -cx lr \
    # -o reads_vs_graph.gaf \
    # graph.gfa \
    # reads.fasta \
    # args = ["minigraph", "-c", "-x", "lr", "-o", output_path, gfa_path, fasta_path]

    if args_override is None:
        args = [
            "minigraph",
            "-c",
            "-g", str(10000),
            "-k", str(14),
            "-f", "0.25",
            "-r", "1000,20000",
            "-n", "3,3",
            "-p", str(0.5),
            # "-j", str(0.85),  # <-- this alone causes horrific slowdown in some regions, no idea why
            "-x", "lr",
            "-t", str(n_threads),
            "-o", output_path,
            gfa_path,
            fasta_path]
    else:
        args = \
            ["minigraph"] + \
            args_override + \
            ["-o", output_path,
            gfa_path,
            fasta_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(args)+'\n')

        try:
            p1 = subprocess.run(args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return output_path


def run_graphchainer(output_directory, gfa_path, fasta_path, n_threads, args_override=None):
    output_path = os.path.join(output_directory, "reads_vs_graph.gaf")

    # GraphChainer -t 30 -f ./bad.fasta -g ./bad.gfa -a graphchain.gaf

    if args_override is None:
        args = [
            "GraphChainer",
            "-c",
            "-x", "lr",
            "-t", str(n_threads),
            "-a", output_path,
            "-g", gfa_path,
            "-f", fasta_path]
    else:
        args = \
            ["minichain"] + \
            args_override + \
            [
                "-t", str(n_threads),
                "-a", output_path,
                "-g", gfa_path,
                "-f", fasta_path
            ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return False
    except Exception as e:
        sys.stderr.write(str(e))
        return False

    return output_path


def run_panaligner(output_directory, gfa_path, fasta_path, n_threads, args_override=None):
    output_path = os.path.join(output_directory, "reads_vs_graph.gaf")

    # minigraph \
    # -cx lr \
    # -o reads_vs_graph.gaf \
    # graph.gfa \
    # reads.fasta \
    # args = ["minigraph", "-c", "-x", "lr", "-o", output_path, gfa_path, fasta_path]

    if args_override is None:
        args = [
            "PanAligner",
            "-c",
            "-g", str(10000),
            "-k", str(14),
            "-f", "0.25",
            "-r", "1000,20000",
            "-n", "3,3",
            "-p", str(0.5),
            # "-j", str(0.85),  # <-- this alone causes horrific slowdown in some regions, no idea why
            "-x", "lr",
            "-t", str(n_threads),
            "-o", output_path,
            gfa_path,
            fasta_path]
    else:
        args = \
            ["PanAligner"] + \
            args_override + \
            ["-o", output_path,
            gfa_path,
            fasta_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(args)+'\n')

        try:
            p1 = subprocess.run(args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    return output_path


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


def run_minimap2_on_read_subset(
        ref_fasta_path,
        reads,
        read_names,
        preset,
        n_threads,
        output_directory,
        filename_prefix="reads_vs_ref"):

    output_filename = os.path.join(output_directory, filename_prefix + ".bam")
    output_path = os.path.join(output_directory,output_filename)

    minimap_args = ["minimap2", "-a", "-x", preset, "--eqx", "-t", str(n_threads), ref_fasta_path, "-"]

    fasta_input_string = list()

    for sequence in reads:
        if sequence.name in read_names:
            fasta_input_string.append(">%s\n%s\n" % (sequence.name, sequence.sequence))

    fasta_input_string = ''.join(fasta_input_string)

    sys.stderr.write("WRITING: %s\n" % output_path)
    with open(output_path, 'w') as file:
        sys.stderr.write(" ".join(minimap_args)+'\n')

        p1 = subprocess.Popen(minimap_args, stdin=subprocess.PIPE, stdout=file, cwd=output_directory)
        p1.communicate(input=fasta_input_string.encode("utf-8"))

    success = (p1.returncode == 0)

    if not success:
        sys.stderr.write("ERROR: failed to align: %s\n" % ref_fasta_path)
        sys.stderr.flush()
        return None

    return output_path


def run_mashmap(ref_fasta_paths: list, reads_fasta_path, min_identity, n_threads, output_directory, filename_prefix="reads_vs_ref_mash"):
    output_filename = os.path.join(output_directory, filename_prefix + ".paf")
    output_path = os.path.join(output_directory,output_filename)

    temp_file_path = os.path.join(output_directory, "refs.txt")
    with open(temp_file_path, 'w') as temp_file:
        for path in ref_fasta_paths:
            print(path)
            temp_file.write(path)
            temp_file.write('\n')

    # mashmap --rl paths/ref_paths.txt -q reads.fasta --dense -s 500 --pi 70 -t 30
    mashmap_args = [
        "mashmap",
        "--rl", temp_file_path,
        "-q", reads_fasta_path,
        "--dense",
        "-s", str(500),
        "--pi", str(int(round(min_identity*100))),
        "-t", str(n_threads),
        "-o", output_path
    ]

    sys.stderr.write(" ".join(mashmap_args)+'\n')

    try:
        p1 = subprocess.run(mashmap_args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        raise e
    except Exception as e:
        sys.stderr.write(str(e))
        raise e

    return output_path


def test():
    project_directory = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    modules_directory = os.path.join(project_directory,"scripts/modules/")
    data_directory = os.path.join(project_directory,"data/test/")

    ref_path = os.path.join(data_directory,"test_ref.fasta")
    query_path = os.path.join(data_directory,"test_query.fasta")

    sequences = list()
    for sequence in iterate_fasta(query_path):
        sequences.append(sequence)

    print(len(sequences))

    read_names = {
        "exact_match_reverse",
        "softclip_2000_reverse",
        "del_5_at_33_reverse",
        "ins_5_at_38_reverse"
    }

    run_minimap2_on_read_subset(
        ref_fasta_path=ref_path,
        reads=sequences,
        read_names=read_names,
        preset="map-hifi",
        n_threads=4,
        output_directory=data_directory,
        filename_prefix="reads_vs_ref_subset")

    return


if __name__ == "__main__":
    test()
