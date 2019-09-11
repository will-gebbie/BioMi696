#!/usr/bin/env
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


# practicing getting in the habit of writing functions
def transcribe_seq(original_seq):
    transcribed_rna_seq = original_seq.transcribe()
    return transcribed_rna_seq


# practicing getting in the habit of writing functions
def translate_rna(transcribed_seq):
    amino_acid = transcribed_seq.translate()
    return amino_acid


def print_sequences(original_seq):
    rna_seq = transcribe_seq(original_seq)
    single_letter_aa = translate_rna(rna_seq)

    print('Original Sequence: {}'.format(original_seq))
    print('Transcribed RNA: {}'.format(rna_seq))
    print('Single Letter Amino Acid Translation: {}'.format(single_letter_aa))


def main():
    seq_1 = Seq('ATGATTGGCCCGGTTTTTTAA', generic_dna)
    seq_2 = Seq('GTGGTGGGGAAATTCCGCTGA', generic_dna)

    print_sequences(seq_1)
    print()
    print_sequences(seq_2)


if __name__ == '__main__':
    main()
