from .Sequence import Sequence
from .Bam import index_bam

from collections import defaultdict
from multiprocessing import Pool
import sys
import os

import pysam
from pysam import AlignedSegment


cigar_index_to_char = [
    'M',    # 0  M  MATCH
    'I',    # 1  I  INS
    'D',    # 2  D  DEL
    'N',    # 3  N  REF_SKIP
    'S',    # 4  S  SOFT_CLIP
    'H',    # 5  H  HARD_CLIP
    'P',    # 6  P  PAD
    '=',    # 7  =  EQUAL
    'X',    # 8  X  DIFF
]

cigar_index_to_formatted_alignment_char = [
    '|',    # 0  M  MATCH
    '-',    # 1  I  INS
    '-',    # 2  D  DEL
    'N',    # 3  N  REF_SKIP
    '#',    # 4  S  SOFT_CLIP
    '#',    # 5  H  HARD_CLIP
    'P',    # 6  P  PAD
    '|',    # 7  =  EQUAL
    '*',    # 8  X  DIFF
]


is_ref_move = [
    True,   # 'M'
    False,  # 'I'
    True,   # 'D'
    True,   # 'N'
    False,  # 'S'
    False,  # 'H'
    False,  # 'P'
    True,   # '='
    True,   # 'X'
]


is_query_move = [
    True,   # 'M'
    True,   # 'I'
    False,  # 'D'
    False,  # 'N'
    True,   # 'S'
    False,  # 'H'
    False,  # 'P'
    True,   # '='
    True,   # 'X'
]


complement = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A'
}


def get_reverse_complement(seq):
    rc = list()
    for s in reversed(seq):
        rc.append(complement[s])

    return ''.join(rc)


"""
From the perspective of the read (query) in its forward direction, iterate the alignment, keeping track of ref and 
query coords for each step (cigar operation)

As a result of using the query as the "ref" in this iterator, ref coords can occasionally increment in reverse, i.e.
ref start > ref_stop, in which case any downstream method would need to account for the reverse complement of the ref
"""
def iterate_cigar(alignment: AlignedSegment):
    # The output is a list of (operation, length) tuples, such as ``[(0, 30)]``.
    cigar_operations = alignment.cigartuples

    ref_start = alignment.reference_start
    query_start = 0

    if alignment.is_reverse:
        cigar_operations = reversed(cigar_operations)
        ref_start = alignment.reference_end
        query_start = 0

    for operation,length in cigar_operations:
        ref_stop = ref_start + length*(1-2*int(alignment.is_reverse))*is_ref_move[operation]
        query_stop = query_start + length*is_query_move[operation]

        yield ref_start,ref_stop,query_start,query_stop,operation,length

        if is_ref_move[operation]:
            ref_start = ref_stop

        if is_query_move[operation]:
            query_start = query_stop


"""
This will get the coordinates of the query sequence in the coordinate space and orientation of the original query 
sequence, NOT the ref-forward-reoriented query sequence which is stored in the BAM/SAM.
"""
def get_query_coord_of_ref_coord(alignment, ref_start, ref_stop, offset_by_hardclip=False):
    query_start = None
    query_stop = None
    first_valid_query_coord = None
    last_valid_query_coord = None

    hardclip_length = 0
    has_hardclip = False
    spans_ref_start = False
    spans_ref_stop = False

    # print("---", alignment.query_name, "is supp: ", alignment.is_supplementary)

    # first_operation,first_length = alignment.cigartuples[0]
    # last_operation,last_length = alignment.cigartuples[-1]
    #
    # if first_operation == 5:
    #     has_hardclip = True
    #     hardclip_length = first_length

    # print(first_operation,first_length)
    # print(last_operation,last_length)


    for i,[cigar_ref_start,cigar_ref_stop,cigar_query_start,cigar_query_stop,operation,length] in enumerate(iterate_cigar(alignment)):
        # TODO: debug HG00741#1#JAHALY010000028.1

        if i == 0 and operation == 5:
            hardclip_length = length
            # print("has hardclip: " + str(hardclip_length))
            has_hardclip = True

        if operation == 5 or operation == 4:
            continue

        in_window = False

        # Ref decreases, query increases
        if alignment.is_reverse:
            # Keep track of any valid coordinates for the case when the read is not spanning
            if ref_start <= cigar_ref_start <= ref_stop:
                if first_valid_query_coord is None:
                    first_valid_query_coord = cigar_query_start

                # Use this indicator to skip unnecessary operations
                in_window = True

            # Keep track of any valid coordinates for the case when the read is not spanning
            if ref_start <= cigar_ref_stop <= ref_stop:
                last_valid_query_coord = cigar_query_stop

                # Use this indicator to skip unnecessary operations
                in_window = True

            # Cigar contains entire window
            #
            #          stop     start
            # ref          |-----|
            # cigar-ref  |-+-----+--|
            #         stop        start
            # query   start        stop
            if cigar_ref_stop <= ref_start <= ref_stop <= cigar_ref_start and is_query_move:
                query_stop = cigar_query_start + (cigar_ref_start - ref_start)
                query_start = cigar_query_stop - (ref_stop - cigar_ref_stop)
                spans_ref_start = True
                spans_ref_stop = True
                break

            if not in_window:
                if ((last_valid_query_coord is not None) or (first_valid_query_coord is not None)):
                    break
                continue
            # else:
            #     if alignment.query_name == "HG00741#1#JAHALY010000028.1":
            #         print(' '.join(list(map(str,['w',ref_start,ref_stop,'r',cigar_ref_start,cigar_ref_stop,'q',cigar_query_start,cigar_query_stop,alignment.is_reverse,cigar_index_to_char[operation],length]))))

            # Cigar overlaps window bound
            #
            #            start   stop
            # ref          |--+--|
            # cigar-ref       |--+--|
            #              stop    start
            # query        start    stop
            if cigar_ref_stop <= ref_stop <= cigar_ref_start:
                query_start = cigar_query_stop - int(is_query_move[operation])*(ref_stop - cigar_ref_stop)
                spans_ref_stop = True

            # Cigar overlaps window bound
            #
            #               start   stop
            # ref             |--+--|
            # cigar-ref    |--+--|
            #            stop    start
            # query      start    stop
            if cigar_ref_stop <= ref_start <= cigar_ref_start:
                query_stop = cigar_query_stop - int(is_query_move[operation])*(ref_start - cigar_ref_stop)
                spans_ref_start = True

        # Ref increases, query increases
        else:
            # Keep track of any valid coordinates for the case when the read is not spanning
            if ref_start <= cigar_ref_start <= ref_stop:
                if first_valid_query_coord is None:
                    first_valid_query_coord = cigar_query_start

                # Use this indicator to skip unnecessary operations
                in_window = True

            # Keep track of any valid coordinates for the case when the read is not spanning
            if ref_start <= cigar_ref_stop <= ref_stop:
                last_valid_query_coord = cigar_query_stop

                # Use this indicator to skip unnecessary operations
                in_window = True

            # Cigar contains entire window
            #
            #          start     stop
            # ref          |-----|
            # cigar-ref  |-+-----+--|
            #         start        stop
            if cigar_ref_start <= ref_start <= ref_stop <= cigar_ref_stop and is_query_move:
                query_start = cigar_query_start + (ref_start - cigar_ref_start)
                query_stop = cigar_query_stop - (cigar_ref_stop - ref_stop)
                spans_ref_start = True
                spans_ref_stop = True
                break

            if not in_window:
                if ((last_valid_query_coord is not None) or (first_valid_query_coord is not None)):
                    break
                continue
            # else:
                # if alignment.query_name == "HG00741#1#JAHALY010000028.1":
                #     print(' '.join(list(map(str,['w',ref_start,ref_stop,'r',cigar_ref_start,cigar_ref_stop,'q',cigar_query_start,cigar_query_stop,alignment.is_reverse,cigar_index_to_char[operation],length]))))

            # Cigar overlaps window bound
            #
            #               start   stop
            # ref          |  |--+--|
            # cigar-ref    |--+--|
            #           start   stop
            if cigar_ref_start <= ref_start <= cigar_ref_stop:
                query_start = cigar_query_start + int(is_query_move[operation])*(ref_start - cigar_ref_start)
                spans_ref_start = True

            # Cigar overlaps window bound
            #
            #          start   stop
            # ref        |--+--|  |
            # cigar-ref     |--+--|
            #             start    stop
            #
            if cigar_ref_start <= ref_stop <= cigar_ref_stop:
                query_stop = cigar_query_start + int(is_query_move[operation])*(ref_stop - cigar_ref_start)
                spans_ref_stop = True

        # If the alignment fully covers the window, stop as soon as both bounds are reached
        if query_start is not None and query_stop is not None:
            break

    # This line may be reached even if the window is not fully covered, in which case we return the last visited coords
    # for whichever side was not spanned
    if query_start is None and query_stop is not None:
        query_start = first_valid_query_coord

    if query_stop is None and query_start is not None:
        query_stop = last_valid_query_coord

    if query_start is None and query_stop is None:
        query_start = first_valid_query_coord
        query_stop = last_valid_query_coord

    # print("before hardclip offset: ", query_start, query_stop, "hardclip length:", hardclip_length)

    if offset_by_hardclip:
        query_start += hardclip_length
        query_stop += hardclip_length

    # print("after hardclip offset: ", query_start, query_stop, "hardclip length:", hardclip_length)
    # sys.stdout.flush()

    return query_start, query_stop, spans_ref_start, spans_ref_stop


def iter_query_sequences_of_region(bam_path, chromosome, ref_start, ref_stop, skip_non_spanning=False):
    sam = pysam.AlignmentFile(bam_path)

    # Need to track multiple alignment records in cases of supplementaries:
    sequence_coords = defaultdict(lambda: [sys.maxsize,0])
    is_spanning = defaultdict(lambda: [False,False])
    is_reverse = defaultdict(lambda: False)
    sequences = dict()

    for alignment in sam.fetch(contig=chromosome, start=ref_start, stop=ref_stop):
        if alignment.is_secondary or alignment.is_unmapped or alignment.mapping_quality < 5:
            continue

        is_hardclipped = False
        first_operation,first_length = alignment.cigartuples[0]
        last_operation,last_length = alignment.cigartuples[-1]

        if first_operation == 5 or last_operation == 5:
            is_hardclipped = True

        query_name = alignment.query_name

        # Only store the sequence if there are no hardclips that will mess up the sequence
        # (otherwise it will be found in an exhaustive search later)
        if not is_hardclipped:
            sequences[query_name] = alignment.query_sequence

        query_start, query_stop, spans_ref_start, spans_ref_stop = get_query_coord_of_ref_coord(
            alignment=alignment,
            ref_start=ref_start,
            ref_stop=ref_stop,
            offset_by_hardclip=True
        )

        is_spanning[query_name][0] += spans_ref_start
        is_spanning[query_name][1] += spans_ref_stop

        # print(query_start,query_stop,alignment.is_reverse,alignment.query_name,alignment.is_supplementary)

        # If we want to use the sequence field of the BAM, we annoyingly have to reconvert back to forward
        # reference orientation because the BAM will store only the forward reference oriented sequence even if the
        # query is a reverse strand read
        if alignment.is_reverse:
            l = 0
            if not is_hardclipped:
                l = len(alignment.query_sequence)
            else:
                l = alignment.infer_read_length()

            start = l - query_stop
            stop = l - query_start

            # Extract the minimum and maximum observed query coordinates that are relevant to this ref region
            # print(query_name, start, stop, stop - start, spans_ref_start, spans_ref_stop)
            sequence_coords[query_name][0] = min(start,sequence_coords[query_name][0])
            sequence_coords[query_name][1] = max(stop,sequence_coords[query_name][1])
            is_reverse[query_name] += True
        else:
            # Extract the minimum and maximum observed query coordinates that are relevant to this ref region
            # print(query_name, query_start, query_stop, query_stop - query_start, spans_ref_start, spans_ref_stop)
            sequence_coords[query_name][0] = min(query_start,sequence_coords[query_name][0])
            sequence_coords[query_name][1] = max(query_stop,sequence_coords[query_name][1])
            is_reverse[query_name] += False

    # Only look for the spanning reads/contigs in the alignment (if specified)
    if skip_non_spanning:
        for name in list(sequences.keys()):
            if not (is_spanning[name][0] and is_spanning[name][1]):
                print("WARNING: skipping non-spanning sequence: " + name)
                del is_reverse[name]
                del sequences[name]
                del sequence_coords[name]

    # It may be that we never iterated the primary sequence, and therefore don't have the unclipped version of it...
    remaining_sequences = set(sequence_coords.keys()) - set(sequences.keys())

    # Initially search just this chromosome
    if len(remaining_sequences) > 0:
        sam = pysam.AlignmentFile(bam_path)

        # Iterate every alignment (in this chromosome) to find it... :(
        for alignment in sam.fetch(contig=chromosome):
            if len(remaining_sequences) == 0:
                break

            if not alignment.is_secondary and not alignment.is_supplementary:
                if alignment.query_name in remaining_sequences:
                    sequences[alignment.query_name] = alignment.query_sequence
                    remaining_sequences.remove(alignment.query_name)

    # Search entire BAM (please no)
    if len(remaining_sequences) > 0:
        sam = pysam.AlignmentFile(bam_path)

        # Iterate every alignment to find it... :(
        for alignment in sam:
            if len(remaining_sequences) == 0:
                break

            if not alignment.is_secondary and not alignment.is_supplementary:
                if alignment.query_name in remaining_sequences:
                    sequences[alignment.query_name] = alignment.query_sequence
                    remaining_sequences.remove(alignment.query_name)

    # Finally iterate the coordinates that were covered by this region and extract them
    for name,[start,stop] in sequence_coords.items():
        if name in sequences:
            s = sequences[name][start:stop]
            r = is_reverse[name]
        else:
            print("WARNING: skipping sequence with no primary alignment in chromosome: " + name)
            continue

        if stop < start:
            raise Exception("ERROR: negative sequence coords obtained from %s:%d-%d in %s" % (chromosome,start,stop,bam_path))

        print("final sequence:", name, stop-start)

        yield name, r, start, stop, s


def print_formatted_alignment_of_query(ref_path, query_path, sam_path, output_directory):
    ref_fasta = pysam.FastaFile(ref_path)
    query_fasta = pysam.FastaFile(query_path)
    sam = pysam.AlignmentFile(sam_path)

    ref_sequence = ref_fasta.fetch('a')

    output_path = os.path.join(output_directory, "test.txt")
    with open(output_path, 'w') as file:
        for alignment in sam:
            file.write('\n')
            file.write('\n')
            file.write(alignment.query_name)
            file.write('\n')

            formatted_query = ""
            formatted_ref = ""
            formatted_cigar = ""

            query_sequence = query_fasta.fetch(alignment.query_name)

            for ref_start,ref_stop,query_start,query_stop,operation,length in iterate_cigar(alignment):
                file.write("operation:\t\t\t%s" % cigar_index_to_char[operation])
                file.write('\n')
                file.write("length:\t\t\t\t%d" % length)
                file.write('\n')
                file.write("ref_coord:\t\t\t%d" % ref_start)
                file.write('\n')
                file.write("ref_coord_stop:\t\t%d" % ref_stop)
                file.write('\n')
                file.write("query_coord:\t\t%d" % query_start)
                file.write('\n')
                file.write("query_coord_stop:\t%d" % query_stop)
                file.write('\n')
                file.write('\n')

                c = cigar_index_to_formatted_alignment_char[operation]

                if is_ref_move[operation]:
                    if alignment.is_reverse:
                        formatted_ref += get_reverse_complement(ref_sequence[ref_stop:ref_start])
                    else:
                        formatted_ref += (ref_sequence[ref_start:ref_stop])
                else:
                    formatted_ref += c*length

                if is_query_move[operation]:
                    formatted_query += (query_sequence[query_start:query_stop])
                else:
                    formatted_query += c*length

                if is_query_move or is_ref_move:
                    formatted_cigar += c*length

            file.write(formatted_ref)
            file.write('\n')
            file.write(formatted_cigar)
            file.write('\n')
            file.write(formatted_query)
            file.write('\n')


def get_sample_haplotypes_of_region(bam_path, chromosome, start, stop, index_threads):
    suffix = None

    if "hap1" in bam_path:
        suffix = "hap1"

    elif "hap2" in bam_path:
        suffix = "hap2"

    else:
        raise Exception("ERROR: 'hap1' or 'hap2' not in file name, unparsable: " + bam_path)

    index_path = bam_path + ".bai"

    if not os.path.exists(index_path):
        index_bam(bam_path, index_threads)

    iter = iter_query_sequences_of_region(
        bam_path=bam_path,
        chromosome=chromosome,
        ref_start=start,
        ref_stop=stop,
        skip_non_spanning=True
    )

    sequences = list()
    for query_name, is_reverse, query_start, query_stop, seq in iter:
        s = Sequence(query_name + "_" + suffix, seq.upper())
        s.normalize_name()
        sequences.append(s)

    if len(sequences) > 1:
        print("WARNING: multiple spanning query sequences in region: %s %s" % (str([x.name for x in sequences]), bam_path))
        print("arbitrarily selecting first sequence...")
    elif len(sequences) == 0:
        print("WARNING: no spanning query sequences found in region: %s" % bam_path)
        return None

    return sequences[0]


def get_haplotypes_of_region(bam_paths, chromosome, start, stop, n_threads):
    args = list()

    for path in bam_paths:
        print(path)
        args.append([path, chromosome, start, stop, True])

    with Pool(n_threads) as pool:
        results = pool.starmap(get_sample_haplotypes_of_region, args)

    return results


"""
This is where regression tests live for the window-specific haplotype fetching problem
"""
def test_haplotype_fetching(data_directory):

    # Multiple spanning sequences...
    bam_path = os.path.join(data_directory,"HG01175_bam_hap1_vs_chm13_chr20.bam")
    chromosome="chr20"
    start=64108857
    stop=64113364
    get_haplotypes_of_region(bam_paths=[bam_path],chromosome=chromosome, start=start, stop=stop, n_threads=1)
    # no solution to this problem yet

    # Oversize sequence
    bam_path = os.path.join(data_directory,"HG00741_bam_hap1_vs_chm13_chr20.bam")
    chromosome="chr20"
    start=64108857
    stop=64113364
    results = get_haplotypes_of_region(bam_paths=[bam_path],chromosome=chromosome, start=start, stop=stop, n_threads=1)

    for item in results:
        print(item.name, len(item))
        if len(item) > 8000:
            raise Exception("FAIL: oversize sequence in HG00741_bam_hap1_vs_chm13_chr20 64108857-64113364")

    exit()




def test_region(test_sam_path, chromosome, start, stop):
    data = defaultdict(list)

    prefixes = set()
    for name,is_reverse,start,stop,sequence in iter_query_sequences_of_region(bam_path=test_sam_path, chromosome=chromosome, ref_start=start, ref_stop=stop):
        data[name].append((is_reverse, start, stop, sequence))
        prefixes.add(name.replace("_reverse",""))

    for prefix in sorted(prefixes):
        a = prefix
        b = prefix + "_reverse"

        a_items = data[a]
        b_items = data[b]

        seqs = set()

        print(prefix)
        for item in a_items:
            is_reverse,start,stop,sequence = item
            print("%d\t%d\t%d\t%s" % (int(is_reverse),start,stop,sequence))
            seqs.add(sequence)

        for item in b_items:
            is_reverse,start,stop,sequence = item
            print("%d\t%d\t%d\t%s" % (int(is_reverse),start,stop,sequence))
            seqs.add(sequence)

        if len(seqs) != 1:
            raise Exception("ERROR: F and R complement not equal")


def test():
    project_directory = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_directory = os.path.join(project_directory,"data/test/")
    test_sam_path = os.path.join(data_directory,"test_alignment_hardclip_sorted.bam")

    ref_path = os.path.join(data_directory,"test_ref.fasta")
    query_path = os.path.join(data_directory,"test_query.fasta")

    test_haplotype_fetching(data_directory)

    print_formatted_alignment_of_query(
        ref_path=ref_path,
        query_path=query_path,
        sam_path=test_sam_path,
        output_directory=data_directory)

    chromosome='a'
    start=2010
    stop=2060
    print("\nTESTING: ", chromosome, start, stop)
    test_region(test_sam_path=test_sam_path, chromosome=chromosome, start=start, stop=stop)

    chromosome='a'
    start=2038
    stop=2040
    print("\nTESTING: ", chromosome, start, stop)
    test_region(test_sam_path=test_sam_path, chromosome=chromosome, start=start, stop=stop)

    chromosome='a'
    start=2048
    stop=2050
    print("\nTESTING: ", chromosome, start, stop)
    test_region(test_sam_path=test_sam_path, chromosome=chromosome, start=start, stop=stop)



if __name__ == "__main__":
    test()
