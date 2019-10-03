#!/usr/bin/env python3
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


def bwt(bw_matrix):
    bwt_list = []
    for row in bw_matrix:
        bwt_list.append(row[-1])

    bwt_string = ''.join(bwt_list)
    return bwt_string


def main():
    seq_1 = Seq('AATTGCGCGG', generic_dna)

    seq_1_index_matrix = create_index_matrix(seq_1)
    seq_1_bwm = bwm(seq_1_index_matrix)
    seq_1_bwt = bwt(seq_1_bwm)

    print(seq_1_bwt)


if __name__ == '__main__':
    main()
