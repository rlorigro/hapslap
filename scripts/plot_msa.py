from matplotlib import pyplot
import numpy


def main():
    path = "/home/ryan/code/hapslap/data/test/hprc_hifiasm/regional_haplotpyes/abpoa_chr20_47474337-47476723_sorted_hi-lo.pir"

    colormap = {
        'A': [1,0,0],
        'C': [0.9,0.8,0],
        'G': [0,1,0],
        'T': [0,0,1],
        '-': [0,0,0]
    }

    n = 0
    names = list()
    rows = list()
    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if line.startswith('>'):
                name = line.strip()[1:]
                names.append(name)
                rows.append([])
            else:
                colors = [colormap[x] for x in line.strip()]
                rows[-1].extend(colors)

                print(len(rows[-1]))
                n += 1

    rows = sorted(rows, key=lambda x: x.count(colormap['-']))

    pyplot.imshow(rows)
    pyplot.tight_layout()

    output_path = '.'.join(path.split('.')[:-1]) + '.pdf'
    pyplot.savefig(output_path, dpi=1200)

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()