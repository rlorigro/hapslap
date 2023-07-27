import sys
import re


class Sequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def normalize_name(self):
        try:
            tokens = re.split("#|_", self.name)
            self.name = '_'.join([tokens[0],self.name[-1]])
        except Exception as e:
            sys.stderr.write("ERROR: cannot normalize sequence name: '%s'\n" % self.name)
            raise e

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.name == other.name and self.sequence == other.sequence
        else:
            return False

    def __len__(self):
        return len(self.sequence)


def iterate_fasta(path, force_upper_case=False, normalize_name=False):
    name = ""
    sequence = ""

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                if force_upper_case:
                    sequence = sequence.upper()

                # Output previous sequence
                if l > 0:
                    s = Sequence(name, sequence)

                    if len(s) != 0:
                        if normalize_name:
                            s.normalize_name()

                        yield s

                name = line.strip()[1:].split(' ')[0]
                sequence = ""

            else:
                sequence += line.strip()

    if force_upper_case:
        sequence = sequence.upper()

    s = Sequence(name, sequence)

    if len(s) != 0:
        if normalize_name:
            s.normalize_name()

        yield s
