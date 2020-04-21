#!/usr/bin/env python3
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from numpy import sort


def create_index_matrix(seq):
    str_T = '$' + str(seq)
    list_T = []
    for char in str_T:
        list_T.append(char)

    index_matrix = ['' for i in range(len(list_T))]
    for i in range(len(list_T) - 1, -1, -1):
        index_matrix[i] = str_T
        left_dollar = list_T[:len(list_T) - 1]
        right_dollar = list_T[len(list_T) - 1:]
        list_T = right_dollar + left_dollar
        str_T = ''.join(list_T)

    return index_matrix


def bwm(index_matrix):
    sorted_index_matrix = sort(index_matrix)
    return sorted_index_matrix


def get_first_and_last(bw_matrix):
    first_col = []
    last_col = []
    for row in bw_matrix:
        first_col.append(row[0])
        last_col.append(row[-1])
    return first_col, last_col


def create_index(index_matrix, bw_matrix):
    index_dict = dict(enumerate(index_matrix))
    index_all = []
    index_evens = []

    for seq in bw_matrix:
        for key, value in index_dict.items():
            if seq == value:
                index_all.append(key)
                if key % 2 == 0:
                    index_evens.append(key)
                else:
                    index_evens.append(None)

    return index_all, index_evens


def index_list_to_file(index_list, filename):
    with open(filename, 'w') as filehandle:
        for index in index_list:
            filehandle.write(str(index) + '\n')


def first_last_to_file(first, last, filename):
    with open(filename, 'w') as filehandle:
        for f, l in zip(first, last):
            filehandle.write(f + '    ' + l + '\n')


def file_to_seq(filename):
    with open(filename, 'r') as file:
        seq_string = file.read()

    seq = Seq(seq_string[1:], generic_dna)
    return seq


def main():
    parser = argparse.ArgumentParser(description='Create suffix array index files and '
                                                 'first and last column file.')
    parser.add_argument('T', metavar='T_File',
                        help='A file containing a Sequence T string in the format'
                             ' ${Sequence}.')
    args = parser.parse_args()

    seq_1 = file_to_seq(args.T)

    seq_1_index_matrix = create_index_matrix(seq_1)
    seq_1_bwm = bwm(seq_1_index_matrix)
    all_indexes, even_indexes = create_index(seq_1_index_matrix, seq_1_bwm)
    first_col, last_col = get_first_and_last(seq_1_bwm)

    # Output to files
    index_list_to_file(all_indexes, 'reference-all.sa')
    index_list_to_file(even_indexes, 'reference-even.sa')
    first_last_to_file(first_col, last_col, 'reference.fl')


if __name__ == '__main__':
    main()
