from modules.Vcf import decompress_vcf
import os


def run():
    input_dir = "/home/ryan/data/test_hapslap/results/competitors/svimmer"
    output_dir = "/home/ryan/data/test_hapslap/results/competitors/svimmer_gt"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for name in os.listdir(input_dir):
        subpath = os.path.join(input_dir, name)

        if subpath.endswith(".vcf.gz"):
            subpath = decompress_vcf(subpath)

        print(subpath)

        lines = list()
        if subpath.endswith(".vcf"):
            with open(subpath, 'r') as file:
                for l,line in enumerate(file):

                    if line.startswith("#") and not line.startswith("##"):
                        line = line.strip() + "\tFORMAT\tGT\n"

                    elif not line.startswith("#") and len(line) > 1:
                        line = line.strip() + "\tGT\t0/1\n"

                    lines.append(line)

            out_path = os.path.join(output_dir, name)

            print(subpath)
            print(out_path)

            with open(out_path.replace(".gz",""), 'w') as file:
                for line in lines:
                    file.write(line)





if __name__ == "__main__":
    run()
