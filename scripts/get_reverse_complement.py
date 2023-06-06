import random


def main():
    complement = {
        'A':'T',
        'C':'G',
        'G':'C',
        'T':'A'
    }

    seq = "CAAAAGATACCATAGCCTATGACAGACTACAAGATATGAACGATGGGGGCAATATGGTTCGCGACAGAATCGCTGTCCTAGTGGTTACGCGCTTACGGAAGTCGGACCACCCTATCGCAAAAACCCCACTCGCCATAGTTCGCCGATGGGAAAAAATTAGACAATCCGCCTAGGATTACTTACGAATGAATTTATGCCATAGCAAATGAAATTGATGACAAACCGGAAGAAAAAGAATCACAAACCCAATGTGTGAAAGCGTCTATTAAGCCATGGTTAGAGGGAAGGAATTGTGTTGATTTTGAAGTAAAAGTGTTTGAAACATTCCAGCCCCGTCACTGTCCTCATCTCCATGGCTATGCAGTCCCGTACACTGCAGCGTTCGCCGACCCAGCGTAAACTCATTGATCTGTTTCCTAGTTAGTTGGGCACTATTCTTGGACTGTTGAGCAATAACCACTACTATAGAGCAGTTGGGTATTCTGGCGGGATGTCTTATCACACCGTGCTCTATTACCTAAAACAAAGCTCCAAGAAGCATAGTTTTCGGCTAAAAACGTGTGTGGGGGCGGGTGCTATGGTTCGGATGACAGAGCTCTATATAAACTACATCGCCTTCTTCGGTTGGCCTATCGGGGGTTCAAGTCCGGGCTCGTTCTGTAACTATAATGCCGACGACTGCTGGGCAGCGGCTAACTGTATTCGTAAGTAATGTGAAGAACTAGTGCCGCTCCTACGCGACAGTGATAGGTCTTCTCTCGTATATTTTGCAACGGCTCTAGTGCGCTACTTGGTAGAGGTTACGCCGTTACGTGAGGTTAGCAGCGCTGCCAGAAAGGGGGCCTGAAGAGAAACTATAACATCCCCTGTGTCTTGCTATTCGGGGGCCCCTAGCCTCGCCCGCTGGCAGAGAGCTTAAATATTGTTATGTAACAGGAGCAGTTGCTACGTGATGGACGATCGAAACAGCGAATTAAACTACAGACCCACGTTGAGTGGGCCCTAGGGTCGTAGCCCATAACTATGATATAACTACGGATGATGGTTTTTCTTCCTAGGATTTTCCATCCGCCGGACCGCTGCCAAGCTAAAGACCCTGATGATAAAATTTCCGTTTGCGTTCATGTCTCGGTCGCAAATGTGGGGTCTTAGGGGTTTTTGCCTTATTGCATCCTGGTCTCATTGCCTTGTGTTTGGTCCATTACGCAGTCAAGTCTAGGGTTGTTACCGAACTCTTCCTCCGATCACACACGTACGGACTTCGCAGGCGATGCATAGTATGATTCCGGTCTCCACAAAGAGCATAACTATAAGGTTTAGGACCAAAATTGGAAGTAACGACTTTCTTGGTGCGACCTTCCACTTAGACCACTTAAGCCAACGATACGACGGTTCGAAAAGAGATATACTTAAGCGTTGAGGTCCTTGTGCGGATACCTCTGCTACAACTTAGGTATTGCACACGGTTACTTATTCAAATCGGGAGGTAACCAGTCCAAAAAGATAGTAGATAGAGATATTGAACAGCACGGCGCTAACATCACCGAAGATCCTCTATTTCCATGCAGCTTCACACTTAAAGAGGTGAGTATACAGTAATACGTTGAATGCCGTTCAGACTACCTCGCTTTCACACACCGTAATTAATCTACTGAAAAACCTGATTGAGTATGCGGGGCTTGGGAGATCCTAGCCGGCGGTTATTTTAAGTGCGTGTGCCTTTTGATTCTTCCCCCGGAGGCTAACTTTACCTCCTACAATACCATGTAAGAGCGCCCATATCCAAGAGGCTGGCCTTAACTGAAGGCAGTTGGTGAATTTTATCTTAGTATTGCCTTGTTCGAGGTATCTCAGTCTACTGAACCGTTAACCTGAGAGGGTGCGAATCGTCGCAAGACGCCAAGGGAACACACTATCAGATTCGGTACCGTACGAACCGTAACTATTTTCTTCAGTAACGGAACCACCAACTGAGAACATGTTCTTCCAGCAAGGACACG"
    rc = list()
    for s in reversed(seq):
        rc.append(complement[s])

    print(''.join(rc))


if __name__ == "__main__":
    main()