#!/usr/bin/env python3
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def calculate_hamming_distance(seq_1, seq_2):
    hamming_dist = 0
    if len(seq_1) != len(seq_2):
        hamming_dist = -1
    else:
        for element1, element2 in zip(seq_1, seq_2):
            if element1 != element2:
                hamming_dist += 1
    return hamming_dist


def main():
    seq_1 = Seq('TTATCGCGCTTTCTTCCAA', generic_dna)
    seq_2 = Seq('ATACCGCGCGTTCGACCAA', generic_dna)

    hamming_dist = calculate_hamming_distance(seq_1, seq_2)
    print(hamming_dist)


if __name__ == '__main__':
    main()