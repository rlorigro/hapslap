import re


class Sequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def normalize_name(self):
        tokens = re.split("#|_", self.name)
        self.name = '_'.join([tokens[0],self.name[-1]])

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.name == other.name and self.sequence == other.sequence
        else:
            return False

    def __len__(self):
        return len(self.sequence)

