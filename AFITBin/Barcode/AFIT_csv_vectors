#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def novel_frequency(sequences, K=8, just_count=False):
    """
    Calculate the Novel Frequency for sequences.

    Parameters
    ----------
    sequences : SeqIO.FastaIO.FastaIterator, list of Bio.SeqRecord.SeqRecord,
        str
        path of fasta file or a list of dna strings.
    K : int, optional
        Length of longest blank area in sub-sequences.
        default is 8.
    just_count : bool, optional
        Output will not divide.
        default False.
    Returns
    -------
    out : pandas.DataFrame
        It has 10 * (K+1) columns.
        Indexes are Bio.SeqRecord.SeqRecord.id.
    """
    import Bio.SeqIO as SeqIO
    import numpy as np
    import pandas as pd
    import re

    if type(sequences) is SeqIO.FastaIO.FastaIterator:
        sequences = list(sequences)
    elif type(sequences) is str:
        sequences = list(SeqIO.parse(sequences, "fasta"))

    patterns = []  # (K+1) list of 10
    pd_columns = []
    for i in range(K + 1):
        pats = []
        dots = i * '.'
        pats.append('(?=(A' + dots + 'A|T' + dots + 'T))')
        pats.append('(?=(A' + dots + 'C|G' + dots + 'T))')
        pats.append('(?=(A' + dots + 'G|C' + dots + 'T))')
        pats.append('(?=(A' + dots + 'T))')
        pats.append('(?=(C' + dots + 'A|T' + dots + 'G))')
        pats.append('(?=(C' + dots + 'C|G' + dots + 'G))')
        pats.append('(?=(C' + dots + 'G))')
        pats.append('(?=(G' + dots + 'A|T' + dots + 'C))')
        pats.append('(?=(G' + dots + 'C))')
        pats.append('(?=(T' + dots + 'A))')
        patterns.append(pats)
        pd_columns.append('A' + str(i) + 'A')
        pd_columns.append('A' + str(i) + 'C')
        pd_columns.append('A' + str(i) + 'G')
        pd_columns.append('A' + str(i) + 'T')
        pd_columns.append('C' + str(i) + 'A')
        pd_columns.append('C' + str(i) + 'C')
        pd_columns.append('C' + str(i) + 'G')
        pd_columns.append('G' + str(i) + 'A')
        pd_columns.append('G' + str(i) + 'C')
        pd_columns.append('T' + str(i) + 'A')

    out = pd.DataFrame(columns=pd_columns)
    for rec in sequences:
        k_frequancy = []
        sq = str(rec.seq).upper()
        for i in range(K + 1):
            for pat in patterns[i]:
                k_frequancy.append(len(re.findall(pat, sq)))
        if not just_count:
            for i in range(K + 1):
                sumofall = sum(k_frequancy[i * 10: (i + 1) * 10])
                for j in range(10):
                    k_frequancy[i * 10 + j] = k_frequancy[i * 10 + j] /\
                        sumofall
        out.loc[rec.id] = k_frequancy
    return out.astype('float32')


def nucleotide_frequency(sequences, K=4, just_count=False):
    """
    Calculate the Nucleotide_frequency Frequancy for sequences.

    Parameters
    ----------
    sequences : SeqIO.FastaIO.FastaIterator, list of Bio.SeqRecord.SeqRecord,
        str
        path of fasta file or a list of dna strings.
    K : int, optional
        Size of k-mers.
        default is 4.
    just_count : bool, optional
        Output will not divide.
        default False.
    Returns
    -------
    out : pandas.DataFrame
        It has (4^K - palindromic) / 2 + palindromic columns.
        Indexes are Bio.SeqRecord.SeqRecord.id.
    """
    def revers_dna(dna):
        rev = ''
        for noc in dna.upper():
            if noc == 'A':
                rev += 'T'
            elif noc == 'T':
                rev += 'A'
            elif noc == 'C':
                rev += 'G'
            elif noc == 'G':
                rev += 'C'
        return rev[::-1]

    def next_step(dna):
        ret = ''
        for i in range(len(dna)):
            if dna[i] == 'T':
                ret += 'A'
            else:
                if dna[i] == 'A':
                    ret += 'C'
                elif dna[i] == 'C':
                    ret += 'G'
                else:
                    ret += 'T'
                break
        return ret + dna[len(ret):]

    import Bio.SeqIO as SeqIO
    import numpy as np
    import pandas as pd
    import re

    if type(sequences) is SeqIO.FastaIO.FastaIterator:
        sequences = list(sequences)
    elif type(sequences) is str:
        sequences = list(SeqIO.parse(sequences, "fasta"))

    patterns = []
    seen = []  # patterns that seen
    pd_columns = []
    step = 'A' * K
    while step != 'T' * K:
        if step in seen:
            step = next_step(step)
            continue
        pd_columns.append(step)
        rev_pat = revers_dna(step)
        seen.append(step)
        seen.append(rev_pat)
        if rev_pat == step:
            patterns.append('(?=(' + step + '))')
        else:
            patterns.append('(?=(' + step + '|' + rev_pat + '))')

    out = pd.DataFrame(columns=pd_columns)
    for rec in sequences:
        frequancy = []
        sq = str(rec.seq).upper()
        for pat in patterns:
            frequancy.append(len(re.findall(pat, sq)))
        if not just_count:
            sumofall = sum(frequancy)
            frequancy = [f / sumofall for f in frequancy]
        out.loc[rec.id] = frequancy
    return out.astype('float32')


def combine_frequency(sequences, Knovel=8, Knuc=4, just_count=False):
    """
    """
    import pandas as pd
    if Knovel < 4:
        return None
    novel = novel_frequency(sequences, K=Knovel, just_count=just_count)
    novel = novel[novel.columns[30:]]
    tnf = novel_frequency(sequences, K=Knuc, just_count=just_count)
    if not just_count:
        tnf = tnf * (Knovel - 3)
    return pd.concat([novel, tnf], 1)


from glob import glob
import os 
import pandas as pd

split_files = os.listdir('/home/sayeh/fasta_BMC/fol/')
print(split_files)
# Get a list of all split FASTA files
#split_files = glob('/home/sayeh/fasta_BMC/fol/pref_NC_083852.1*.fasta')

# Initialize an empty DataFrame to store concatenated results
concatenated_results = pd.DataFrame()

# Process each split file using the provided functions
for file in split_files:
    result = novel_frequency('/home/sayeh/fasta_BMC/fol/' + file, K=8, just_count=False)
    # Do something with the result, such as saving it to a file or further analysis
    # For example, to save the result to a CSV file:
    result.to_csv(f"/home/sayeh/results.csv")
    # Append the result to the concatenated DataFrame
    concatenated_results = pd.concat([concatenated_results, result])

# Save the concatenated results to a single CSV file
concatenated_results.to_csv('/home/sayeh/combined_results.csv')


import pandas as pd
from natsort import index_natsorted

# Read the CSV file into a DataFrame
data = pd.read_csv('/home/sayeh/combined_results.csv')

# Get the index (second column) as a list
index_list = list(data.iloc[:, 0])

# Perform natural sorting on the index list
sorted_index = index_natsorted(index_list)

# Reorder the DataFrame based on the sorted index
data_sorted = data.loc[sorted_index]

# Save the sorted DataFrame back to a CSV file
data_sorted.to_csv('/home/sayeh/combined_results_sorted.csv', index=False)

t = pd.read_csv('/home/sayeh/combined_results.csv').head()
j = pd.read_csv('/home/sayeh/combined_results_sorted.csv').head()

print(t)
print
