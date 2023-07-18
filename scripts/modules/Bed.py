class Region:
    def __init__(self, chromosome_name, start, stop):
        self.contig_name = chromosome_name
        self.start = start
        self.stop = stop

    def __str__(self):
        return "%s_%d-%d" % (self.contig_name, self.start, self.stop)


def iter_bed_items(path):
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split()
            contig = tokens[0]
            start = int(tokens[1])
            stop = int(tokens[2])

            yield contig, start, stop


def parse_bed_regions(bed_path):
    regions = list()
    for chromosome, start, stop in iter_bed_items(bed_path):
            regions.append(Region(chromosome, start, stop))

    return regions

