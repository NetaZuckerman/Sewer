from collections import Counter
import pandas as pd
import numpy as np
from sys import argv
import glob
import os
import re
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
    total = pileup_df.loc[pos]['sum']
    freq = None
    if total and total >= depth_threshold:
        count = pileup_df.loc[pos][mut_val]  # specific mutation frequency
        if count >= depth_threshold:
            freq = (count / total) * 100
        else:
            freq = 0.0
            # freq = None

    return freq


def keysort(elem):
    try:
        a = re.findall('\d+|\D+', elem)
        return int(a[1])
    except ValueError:
        return 0


def sortAndTranspose(df):
    # 'Unnamed: 0',
    df = df.reindex(columns=[
        'B.1.1.7 - UK avg', 'B.1.351 - SA avg', 'P.1 - Manaus avg', 'P.2 - Rio de jeneiro avg',
        'B.1.429 - California avg',
        'B.1.525 - Global avg', 'B.1.526 - New york avg', 'A.23.1 - Uganda avg',
        '20C/H655Y - Brittany avg', 'B.1.1.7 - UK freq', 'B.1.351 - SA freq', 'P.1 - Manaus freq'
        , 'P.2 - Rio de jeneiro freq', 'B.1.429 - California freq', 'B.1.525 - Global freq', 'B.1.526 - New york freq',
        'A.23.1 - Uganda freq', '20C/H655Y - Brittany freq', 'VOI-18.02 - WHO freq', 'VUI_L452R/L1063F_Israel freq',
        'VUI_N481K_Israel freq', 'VUI_P681H_Israel freq', 'VOI-18.02 - WHO avg', 'VUI_L452R/L1063F_Israel avg',
        'VUI_N481K_Israel avg', 'VUI_P681H_Israel avg'])
    df = df.transpose()
    try:
        df = df[sorted(df.columns, key=keysort)]
    except:
        print("error in columns sorting")
    df = df.transpose()
    return df

if __name__ == '__main__':
    pileup_table = pd.read_csv("Env1.csv")
    min_depth = 5
    refseq_path = "REF_NC_045512.2.fasta"
    muttable = pd.read_csv("novelMutTable.csv")
    muttable = muttable.drop(muttable[muttable['type'] == 'Insertion'].index)
    uniq_lineages = set()
    for lin in muttable.lineage:
        for x in lin.split(','):
            uniq_lineages.add(x.strip())
    muttable_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}
    for lin, table in muttable_by_lineage.items():
        table.lineage = lin
    final_df = pd.concat([frame for frame in muttable_by_lineage.values()])
    all_tables = {}
    final_df = final_df.apply(lambda row: frequency(row['mut'], row['pos']-1, pileup_table, min_depth), axis=1)