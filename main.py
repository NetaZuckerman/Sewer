
from collections import Counter
import pandas as pd
import numpy as np
from sys import argv
import glob
import os
"""
for sewer samples.
iterate over mapped and sorted bam files,
calculate mutations frequencies by pileup files.
"""
# TODO: results table to calculate overall mutations frequencies.
# TODO: for each sample: csv file of mutations positions (mutation = position total depth > 10, and mutation(not N)
#       is > 1% frequency
#

def frequency(mut_val, pos, pileup_df, depth_threshold):
    """
    return frequency of mut_val base in specified position.
    :param mut_val: mutation nucleotide (A,C,T,G,-)
    :param pos: position
    :param pileup_df: the pileup dataframe
    :param depth_threshold: minimum depth threshold. below that, frequency is 0.
    :return: frequency of mutation nucleotide in position (mut_depth/sum*100)
    """
    mut_val = 'del' if mut_val == '-' else mut_val
    total = pileup_df.loc[pos-1]['sum']
    if total:
        count = pileup_df.loc[pos-1][mut_val]
        if count > depth_threshold:
            freq = (count / total) * 100
        else:
            freq = 0.0
    else:
        freq = 0.0
    return freq


if __name__ == '__main__':
    muttable = pd.read_csv("novelMutTable.csv")
    uniq_lineages = set()
    for lin in muttable.lineage:
        for x in lin.split(','):
            uniq_lineages.add(x.strip())
    muttable_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}
    for lin, table in muttable_by_lineage.items():
        table.lineage = lin

    final_df = pd.concat([frame for frame in muttable_by_lineage.values()])

    pileup_table = pd.read_csv('temp_pileuptable.csv')
    pileup_table['ref_freq'] = pileup_table.apply(
        lambda row: (row[row['ref']] / row['sum']) * 100 if row['sum'] else 0.0, axis=1)

    final_df['file_name'] = final_df.apply(lambda row: frequency(row['mut'], row['pos'], pileup_table, 5), axis=1)

    for x in final_df:
        print(x)





