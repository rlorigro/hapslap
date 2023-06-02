from collections import defaultdict


def main():
    path = "/home/ryan/code/hapslap/data/test/hprc_hifiasm/regional_haplotpyes/abpoa_chr20_47474337-47476723_sorted_hi-lo.gfa"

    output_path = path.replace(".gfa", "_with_coverage.gfa")

    coverages = defaultdict(int)
    s_lines = dict()
    l_lines = list()
    p_lines = list()
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            if tokens[0] == 'S':
                name = tokens[1]
                s_lines[name] = line

            if tokens[0] == 'L':
                l_lines.append(line)

            if tokens[0] == 'P':
                p_lines.append(line)

                data = tokens[2].split(',')
                for item in data:
                    name = item.strip('-').strip('+')
                    coverages[name] += 1

    with open(output_path, 'w') as file:
        for name,line in s_lines.items():
            line = line.strip('\n')

            file.write(line)
            file.write('\t')
            file.write("RC:i:")
            file.write(str(coverages[name]))
            file.write('\n')

        for item in l_lines:
            file.write(item)

        for item in p_lines:
            file.write(item)




if __name__ == "__main__":
    main()
