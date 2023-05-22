import random


def main():
    bases = ['A','C','G','T']

    l = 2000
    seq = list()
    for i in range(l):
        seq.append(random.choice(bases))

    print(''.join(seq))


if __name__ == "__main__":
    main()