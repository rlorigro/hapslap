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

