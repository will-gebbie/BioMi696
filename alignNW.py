#!/usr/bin/env python3
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from pandas import DataFrame


def create_nw_table(seq_1, seq_2):
    # Added gap element
    rows = len(seq_2) + 1
    cols = len(seq_1) + 1

    # Initialize 2d array
    scores = [[(0, '') for i in range(cols)] for j in range(rows)]
    gap = -2
    match = 1
    mismatch = -1

    # Initialize first row and first column of array with gap values.
    scores[0][0] = (0, 'H')
    for i in range(1, len(scores[0])):
        scores[0][i] = (gap * i, 'H')
    for i in range(1, len(scores)):
        scores[i][0] = (gap * i, 'V')

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

    return scores


def alignment_score(nw_matrix):
    last_row = len(nw_matrix) - 1
    last_col = len(nw_matrix[0]) - 1

    return nw_matrix[last_row][last_col][0]


def alignment_figure(seq1, seq2, nw_matrix):
    cols = len(nw_matrix[0])
    rows = len(nw_matrix)
    figure_cols = max(cols, rows)
    figure_rows = 3

    # Initialize matrix to display alignment with '|'s
    figure_matrix = [['' for i in range(figure_cols)] for j in range(figure_rows)]
    nw_i = rows - 1
    nw_j = cols - 1
    fig_j = figure_cols - 1

    # Fill in figure with proper traceback route for NW algorithm
    while fig_j >= 0 and nw_i > 0:
        if nw_matrix[nw_i][nw_j][1] == 'H':
            figure_matrix[0][fig_j] = seq1[nw_j - 1]
            figure_matrix[2][fig_j] = '-'
            fig_j -= 1
            nw_j -= 1
        elif nw_matrix[nw_i][nw_j][1] == 'V':
            figure_matrix[0][fig_j] = '-'
            figure_matrix[2][fig_j] = seq2[nw_i - 1]
            fig_j -= 1
            nw_i -= 1
        elif nw_matrix[nw_i][nw_j][1] == 'D':
            figure_matrix[0][fig_j] = seq1[nw_j - 1]
            figure_matrix[2][fig_j] = seq2[nw_i - 1]
            if nw_matrix[nw_i][nw_j][0] > nw_matrix[nw_i - 1][nw_j - 1][0]:
                figure_matrix[1][fig_j] = '|'
            elif nw_matrix[nw_i][nw_j][0] < nw_matrix[nw_i - 1][nw_j - 1][0]:
                figure_matrix[1][fig_j] = ''
            fig_j -= 1
            nw_i -= 1
            nw_j -= 1

    return figure_matrix


def pretty2d(arr):
    pretty_data = DataFrame(arr)
    print(pretty_data.to_string(index=False, header=False))


def main():
    seq_1 = Seq('TATCGCGCTTT', generic_dna)
    seq_2 = Seq('ATTACCGCCGTT', generic_dna)

    nw_matrix = create_nw_table(seq_1, seq_2)
    alignment_score_seq1_seq2 = alignment_score(nw_matrix)
    alignment_seq1_seq2 = alignment_figure(seq_1, seq_2, nw_matrix)

    print('Alignment Score =', alignment_score_seq1_seq2)
    print('------------------------------------------------')
    pretty2d(alignment_seq1_seq2)
    print('------------------------------------------------')


if __name__ == '__main__':
    main()
