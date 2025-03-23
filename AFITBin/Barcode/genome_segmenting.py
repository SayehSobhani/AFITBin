#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 19:35:09 2024

@author: sayeh
"""

from Bio import SeqIO

def split_sequences(input_file, output_prefix, size):
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        sequence_length = len(sequence)
        num_splits = sequence_length // size

        for i in range(num_splits):
            start = i * size
            end = start + size
            split_seq = sequence[start:end]
            output_file = f"{output_prefix}_{record.id}_split{i+1}.fasta"
            with open(output_file, "w") as f:
                f.write(f">{record.id}_split{i+1}\n{split_seq}\n")

# Replace 'input.fasta' with the name of your input FASTA file
# Replace 'output_prefix' with the desired prefix for the output files
split_sequences('/home/sayeh/Ecoli/2/Escherichia coli O100:H21 strain Res13-Lact-ER07-18 chromosome, complete genome.fasta', '/home/sayeh/Ecoli/2/precoli2', 1000)
