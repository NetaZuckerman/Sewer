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

def no_uk_calculate(no_uk_df, other_variants):
    no_uk_df = no_uk_df[(no_uk_df.AA.isin(other_variants))]
    # create another surveillance table
    lineage_avg = no_uk_df.drop('pos', axis=1).groupby('lineage').mean().transpose()
    # calculate frequency
    lineage_num_muts = no_uk_df.groupby('lineage')['lineage'].count().to_frame().rename(columns={'lineage': 'total'})
    lineage_non_zero_count = no_uk_df.drop(columns=['nucleotide', 'AA', 'gene', 'type', 'pos', 'REF', 'mut']) \
        .groupby('lineage').agg(lambda x: x.ne(0).sum())
    lineage_freq = lineage_num_muts.join(lineage_non_zero_count)
    return lineage_freq, lineage_avg


def uk_calculate(uk_df, uk_variant_mutations):
    uk_df = uk_df[(uk_df.AA.isin(uk_variant_mutations))]
    # create another surveillance table
    lineage_avg = uk_df.drop('pos', axis=1).groupby('lineage').mean().transpose()
    # calculate frequency
    lineage_num_muts = uk_df.groupby('lineage')['lineage'].count().to_frame().rename(columns={'lineage': 'total'})
    lineage_non_zero_count = uk_df.drop(columns=['nucleotide', 'AA', 'gene', 'type', 'pos', 'REF', 'mut']) \
        .groupby('lineage').agg(lambda x: x.ne(0).sum())
    lineage_freq = lineage_num_muts.join(lineage_non_zero_count)
    lineage_avg = lineage_avg['B.1.1.7 - UK']
    uk_total= lineage_freq['total']['B.1.1.7 - UK']
    lineage_freq = lineage_freq.loc['B.1.1.7 - UK', :].transpose()
    lineage_freq = lineage_freq.astype(int).astype(str) + '\\' + uk_total.astype(str)
    return lineage_freq, lineage_avg

if __name__ == '__main__':
    # pileup_table = pd.read_csv("Env1.csv")
    min_depth = 5
    # refseq_path = "REF_NC_045512.2.fasta"
    muttable = pd.read_csv("novelMutTable.csv")
    muttable = muttable.drop(muttable[muttable['type'] == 'Insertion'].index)
    uniq_lineages = set()
    for lin in muttable.lineage:
        for x in lin.split(','):
            uniq_lineages.add(x.strip())
    all_mutations = set([x for x in muttable.AA])
    muttable_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}
    for lin, table in muttable_by_lineage.items():
        table.assign(lineage=lin)
    final_df = pd.read_csv('monitored_mutations.csv', index_col=0)
    all_tables = {}
    uk_variant_mutations = set(muttable_by_lineage['B.1.1.7 - UK']['AA'])  # list of mutations of uk variant
    other_variants = all_mutations - uk_variant_mutations

    no_uk_lineage_freq, no_uk_lineage_avg = no_uk_calculate(final_df.copy(), other_variants)
    uk_lineage_freq, uk_lineage_avg = uk_calculate(final_df.copy(), uk_variant_mutations)

    for name in all_tables.keys():
        # lineage_freq[name] /= lineage_freq['total']/100
        no_uk_lineage_freq[name] = no_uk_lineage_freq[name].astype(int).astype(str) + '\\' + no_uk_lineage_freq[
            'total'].astype(str)

    lineage_freq = no_uk_lineage_freq.drop(columns='total').transpose()
    surv_table = lineage_freq.add_suffix(' freq').join(no_uk_lineage_avg.add_suffix(' avg'))
    surv_table = sortAndTranspose(surv_table)
    surv_table['B.1.1.7 - UK avg'] = uk_lineage_avg
    surv_table['B.1.1.7 - UK freq'] = uk_lineage_freq
    surv_table.to_csv('surveillance_table.csv')

