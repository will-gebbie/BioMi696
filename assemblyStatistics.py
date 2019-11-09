import argparse
from Bio import SeqIO
import numpy as np


def calculate_N50(contig_lengths):
    contig_sum = sum(contig_lengths)
    half_sum = int(contig_sum / 2)

    contigs = []
    i = 0
    while sum(contigs) <= half_sum:
        contigs.append(contig_lengths[i])
        i += 1
    n50 = contigs[-2]
    return n50


def get_num_contigs(file):
    count = 0
    contig_lengths = []
    for record in SeqIO.parse(file, "fasta"):
        contig_lengths.append(len(record.seq))
        count += 1

    sorted_lengths = np.sort(contig_lengths)
    return count, sorted_lengths[::-1]


def arg_parser():
    parser = argparse.ArgumentParser(description='Tabulate the number of contigs, and calculate the N50 of '
                                                 'an input FASTA file')
    parser.add_argument('fasta_file', metavar='<file.fasta>', help='A fasta formatted file')
    args = parser.parse_args()

    return args


def main():
    args = arg_parser()

    file = args.fasta_file

    num_contigs, contig_lengths = get_num_contigs(file)
    print('Number of contigs: ', num_contigs)
    print('Contig lengths: ', contig_lengths)

    n50 = calculate_N50(contig_lengths)
    print('N50: ', n50)


if __name__ == '__main__':
    main()
