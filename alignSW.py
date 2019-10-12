#!/usr/bin/env python3
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from pandas import DataFrame


def create_sw_table(seq_1, seq_2, input_gap, input_match, input_mismatch):
    # Added gap element
    rows = len(seq_2) + 1
    cols = len(seq_1) + 1

    # Initialize 2d array
    scores = [[(0, '') for i in range(cols)] for j in range(rows)]
    gap = input_gap
    match = input_match
    mismatch = input_mismatch

    # Initialize first row and first column of array with gap values (0).
    scores[0][0] = (0, 'H')
    for i in range(1, len(scores[0])):
        scores[0][i] = (0, 'H')
    for i in range(1, len(scores)):
        scores[i][0] = (0, 'V')

    # Calculate Scores of matrix
    for i in range(1, len(scores)):
        for j in range(1, len(scores[0])):
            horizontal_score = scores[i][j - 1][0] + gap
            vertical_score = scores[i - 1][j][0] + gap

            if seq_1[j - 1] == seq_2[i - 1]:
                diagonal_score = scores[i - 1][j - 1][0] + match
            else:
                diagonal_score = scores[i - 1][j - 1][0] + mismatch

            if diagonal_score > horizontal_score and diagonal_score >= vertical_score:
                scores[i][j] = (diagonal_score, 'D')
            elif vertical_score > horizontal_score and vertical_score > diagonal_score:
                scores[i][j] = (vertical_score, 'V')
            elif horizontal_score >= vertical_score and horizontal_score >= diagonal_score:
                scores[i][j] = (horizontal_score, 'H')

            # Added condition for SW math
            d_v_h = scores[i][j][1]
            if scores[i][j][0] < 0:
                scores[i][j] = (0, d_v_h)

    return scores


def find_max_index(sw_matrix):
    max_i = 0
    max_j = 0
    for i in range(len(sw_matrix)):
        for j in range(len(sw_matrix[0])):
            if sw_matrix[max_i][max_j] < sw_matrix[i][j]:
                max_i = i
                max_j = j

    return max_i, max_j


def traceback_length(sw_matrix, max_i, max_j):
    sw_i = max_i
    sw_j = max_j
    length = 0
    while sw_matrix[sw_i][sw_j][0] > 0:
        if sw_matrix[sw_i][sw_j][1] == 'H':
            sw_j -= 1
            length += 1
        elif sw_matrix[sw_i][sw_j][1] == 'V':
            sw_i -= 1
            length += 1
        elif sw_matrix[sw_i][sw_j][1] == 'D':
            sw_i -= 1
            sw_j -= 1
            length += 1
    return length


def alignment_figure(seq1, seq2, sw_matrix, align_length, max_i, max_j):
    figure_cols = align_length
    figure_rows = 3

    figure_matrix = [['' for i in range(figure_cols)] for j in range(figure_rows)]
    sw_i = max_i
    sw_j = max_j
    fig_j = figure_cols - 1

    while sw_matrix[sw_i][sw_j][0] > 0:
        if sw_matrix[sw_i][sw_j][1] == 'H':
            figure_matrix[0][fig_j] = seq1[sw_j - 1]
            figure_matrix[2][fig_j] = '-'
            fig_j -= 1
            sw_j -= 1
        elif sw_matrix[sw_i][sw_j][1] == 'V':
            figure_matrix[0][fig_j] = '-'
            figure_matrix[2][fig_j] = seq2[sw_i - 1]
            fig_j -= 1
            sw_i -= 1
        elif sw_matrix[sw_i][sw_j][1] == 'D':
            figure_matrix[0][fig_j] = seq1[sw_j - 1]
            figure_matrix[2][fig_j] = seq2[sw_i - 1]
            if sw_matrix[sw_i][sw_j][0] > sw_matrix[sw_i - 1][sw_j - 1][0]:
                figure_matrix[1][fig_j] = '|'
            elif sw_matrix[sw_i][sw_j][0] < sw_matrix[sw_i - 1][sw_j - 1][0]:
                figure_matrix[1][fig_j] = ''
            fig_j -= 1
            sw_i -= 1
            sw_j -= 1
    return figure_matrix


def pretty2d(arr):
    pretty_data = DataFrame(arr)
    print(pretty_data.to_string(index=False, header=False))


def file_to_seq(filename):
    with open(filename, 'r') as file:
        seq_string = file.read()

    seq = Seq(seq_string, generic_dna)
    return seq


def main():
    parser = argparse.ArgumentParser(description='Smith Waterman Algorithm implementation '
                                                 'to find largest relative sub-sequence '
                                                 'between two sequences.')
    parser.add_argument('T1', metavar='T file #1',
                        help='A file containing a Sequence T string. (w/o $)')
    parser.add_argument('T2', metavar='T file #2',
                        help='A file containing a Sequence T string. (w/o $)')
    parser.add_argument('--gap', action='store', dest='gap', default=-2, type=int,
                        help='Value for gap penalty. Default = -2')
    parser.add_argument('--match', action='store', dest='match', default=1, type=int,
                        help='Value for match score. Default = 1')
    parser.add_argument('--mismatch', action='store', dest='mismatch', default=-1, type=int,
                        help='Value for mismatch penalty. Default = -1')
    args = parser.parse_args()

    seq_1 = file_to_seq(args.T1)
    seq_2 = file_to_seq(args.T2)
    input_gap = args.gap
    input_match = args.match
    input_mismatch = args.mismatch

    sw_matrix = create_sw_table(seq_1, seq_2, input_gap, input_match, input_mismatch)
    max_row, max_col = find_max_index(sw_matrix)
    sw_align_length = traceback_length(sw_matrix, max_row, max_col)
    sw_align_fig = alignment_figure(seq_1, seq_2, sw_matrix, sw_align_length,
                                    max_row, max_col)

    print('------------------------------------------------')
    pretty2d(sw_align_fig)
    print('------------------------------------------------')


if __name__ == '__main__':
    main()
