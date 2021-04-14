import pysam
from collections import Counter
import pandas as pd
import numpy as np
from sys import argv
import glob
import os
import errno
import re

"""
for sewer samples.
iterate over mapped and sorted bam files,
calculate mutations frequencies by pileup files.
"""


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
        'B.1.1.7 - UK avg', 'B.1.1.7 - UK freq', 'B.1.351 - SA avg', 'B.1.351 - SA freq', 'P.1 - Manaus avg',
        'P.1 - Manaus freq' 'P.2 - Rio de jeneiro avg', 'P.2 - Rio de jeneiro freq',
        'B.1.429 - California avg', 'B.1.429 - California freq',
        'B.1.525 - Global avg', 'B.1.525 - Global freq', 'B.1.526 - New york avg', 'B.1.526 - New york freq',
        'A.23.1 - Uganda avg', 'A.23.1 - Uganda freq',
        '20C/H655Y - Brittany avg', '20C/H655Y - Brittany freq',
        'VOI-18.02 - WHO freq', 'VUI_L452R/L1063F_Israel freq',
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
    no_uk_df.fillna(-1, inplace=True)
    lineage_non_zero_count = no_uk_df.drop(columns=['nucleotide', 'AA', 'gene', 'type', 'pos', 'REF', 'mut']) \
        .groupby('lineage').agg(lambda x: x.gt(0).sum())
    no_uk_df.replace(-1, None, inplace=True)
    lineage_freq = lineage_num_muts.join(lineage_non_zero_count)
    return lineage_freq, lineage_avg


def uk_calculate(uk_df, uk_variant_mutations):
    uk_df = uk_df[(uk_df.AA.isin(uk_variant_mutations))]
    # create another surveillance table
    lineage_avg = uk_df.drop('pos', axis=1).groupby('lineage').mean().transpose()
    # calculate frequency
    lineage_num_muts = uk_df.groupby('lineage')['lineage'].count().to_frame().rename(columns={'lineage': 'total'})
    uk_df.fillna(-1, inplace=True)
    lineage_non_zero_count = uk_df.drop(columns=['nucleotide', 'AA', 'gene', 'type', 'pos', 'REF', 'mut']) \
        .groupby('lineage').agg(lambda x: x.gt(0).sum())
    uk_df.replace(-1, None, inplace=True)
    lineage_freq = lineage_num_muts.join(lineage_non_zero_count)
    lineage_avg = lineage_avg['B.1.1.7 - UK']
    uk_total = lineage_freq['total']['B.1.1.7 - UK']
    lineage_freq = lineage_freq.loc['B.1.1.7 - UK', :].transpose()
    lineage_freq = lineage_freq.astype(int).astype(str) + '\\' + uk_total.astype(str)
    return lineage_freq, lineage_avg


if __name__ == '__main__':
    # user input
    bam_dir = argv[1]
    min_depth = int(argv[2])
    refseq_path = argv[3]
    # preparations
    refseq_name = os.path.basename(refseq_path).strip('.fasta')
    # index refseq
    pysam.faidx(refseq_path)
    refseq_series = pd.Series([x for x in pysam.Fastafile(refseq_path).fetch(reference=refseq_name)])
    muttable = pd.read_csv("/data/projects/Dana/scripts/covid19/novelMutTable.csv")  # TODO: get from other location!
    # muttable = pd.read_csv("novelMutTable.csv") # TODO change before commit
    muttable = muttable.drop(muttable[muttable['type'] == 'Insertion'].index)
    uniq_lineages = set()
    for lin in muttable.lineage:
        for x in lin.split(','):
            uniq_lineages.add(x.strip())
    muttable_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}
    for lin, table in muttable_by_lineage.items():
        table.lineage = lin
    # table = table.assign(lineage=lin)

    final_df = pd.concat([frame for frame in muttable_by_lineage.values()])
    all_mutations = set([x for x in muttable.AA])
    uk_variant_mutations = set(muttable_by_lineage['B.1.1.7 - UK']['AA'])  # list of mutations of uk variant
    other_variants = all_mutations - uk_variant_mutations
    other_variants.remove('L5F')

    all_tables = {}

    files_list = glob.glob(bam_dir + '/*.mapped.sorted.bam')

    # iterate all bam files:
    for file in files_list:
        pileup_table = pd.DataFrame(np.empty(shape=(29903, 6)) * np.nan, columns=['C', 'A', 'G', 'T', 'N', 'del'],
                                    index=list(range(29903)))
        bam = pysam.AlignmentFile(file, 'rb')
        pileup_iter = bam.pileup(stepper='nofilter')
        # iterate over reads in each position and count nucleotides, Ns and deletions.
        for position in pileup_iter:
            c = Counter({'C': 0, 'A': 0, 'G': 0, 'T': 0, 'N': 0, 'del': 0})
            for pileupread in position.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    c[pileupread.alignment.query_sequence[pileupread.query_position].upper()] += 1
                elif pileupread.is_del:
                    c['del'] += 1
                elif pileupread.is_refskip:  # N?
                    c['N'] += 1
            pileup_table.loc[position.reference_pos] = pd.Series(c)
        # produce pileup table(for each bam): pos,A,C,T,G,N,del,totaldepth,
        pileup_table.index.name = 'pos'
        pileup_table['sum'] = pileup_table['A'] + pileup_table['C'] + pileup_table['T'] + pileup_table['G'] + \
                              pileup_table['del'] + pileup_table['N']

        pileup_table['ref'] = refseq_series
        # pileup_table.to_csv('temp_pileuptable.csv')  # to remove after debug
        pileup_table['ref_freq'] = pileup_table.apply(
            lambda row: (row[row['ref']] / row['sum']) * 100 if row['sum'] else None,
            axis=1)  # if not row['sum'] then no coverage at all.
        pileup_table['C_freq'] = pileup_table.apply(
            lambda row: (row['C'] / row['sum']) * 100 if row['sum'] else None, axis=1)
        pileup_table['A_freq'] = pileup_table.apply(
            lambda row: (row['A'] / row['sum']) * 100 if row['sum'] else None, axis=1)
        pileup_table['G_freq'] = pileup_table.apply(
            lambda row: (row['G'] / row['sum']) * 100 if row['sum'] else None, axis=1)
        pileup_table['T_freq'] = pileup_table.apply(
            lambda row: (row['T'] / row['sum']) * 100 if row['sum'] else None, axis=1)
        pileup_table['N_freq'] = pileup_table.apply(
            lambda row: (row['N'] / row['sum']) * 100 if row['sum'] else None, axis=1)
        pileup_table['del_freq'] = pileup_table.apply(
            lambda row: (row['del'] / row['sum']) * 100 if row['sum'] else None, axis=1)
        # add sample to table
        file_name = file.strip('BAM/').strip('.mapped.sorted.bam')
        all_tables[file_name] = pileup_table
        final_df[file_name] = final_df.apply(lambda row: frequency(row['mut'], row['pos'] - 1, pileup_table, min_depth),
                                             axis=1)

    final_df = final_df.sort_values(["lineage", "gene"], ascending=(True, False))  # sort by:(1)lineage (2)gene(S first)

    sortednames = sorted([x for x in final_df.columns.values if "nv" in x], key=keysort)
    sorted_cols = [c for c in final_df.columns.values if c not in sortednames] + sortednames
    final_df = final_df.reindex(columns=sorted_cols)

    if not os.path.exists('results'):
        os.mkdir('results/')
    monitoredfile = final_df.copy()
    monitoredfile.fillna(-1, inplace=True)
    monitoredfile.replace(-1, "No Coverage", inplace=True)
    monitoredfile.to_csv("results/monitored_mutations.csv")

    try:
        os.makedirs('results/mutationsPileups')
    except OSError as e:
        if e.errno == errno.EEXIST:
            print("directory results/mutationsPileups already exists, continuing.")
        else:
            raise
    try:
        os.makedirs('results/fullPileups')
    except OSError as e:
        if e.errno == errno.EEXIST:
            print("directory results/mutationsPileups already exists, continuing.")
        else:
            raise
    # write pileup files that contain only positions mutations
    for name, table in all_tables.items():
        # keep only lines that: >1% frequency of non refseq mutation AND >=10 depth (line.sum)
        table['N_freq'] = table.apply(lambda row: (row['N'] / row['sum']) * 100 if row['sum'] else 0.0, axis=1)
        table = table.dropna(thresh=3)
        table.to_csv('results/fullPileups/' + name + '.csv')
        indexNames = table[(table['sum'] < 10) | (table['ref_freq'] > 99) | (table['N_freq'] > 99)].index
        table = table.drop(index=indexNames, columns=['N_freq'])
        table.reset_index(level=0, inplace=True)
        table.to_csv('results/mutationsPileups/' + name + '.csv', index=False)

    no_uk_lineage_freq, no_uk_lineage_avg = no_uk_calculate(final_df.copy(), other_variants)
    uk_lineage_freq, uk_lineage_avg = uk_calculate(final_df.copy(), uk_variant_mutations)

    for name in all_tables.keys():
        # lineage_freq[name] /= lineage_freq['total']/100
        no_uk_lineage_freq[name] = no_uk_lineage_freq[name].astype(int).astype(str) + '\\' + no_uk_lineage_freq[
            'total'].astype(str)+"("+(no_uk_lineage_freq[name]/no_uk_lineage_freq['total'] / 100).astype(str)+")"

    lineage_freq = no_uk_lineage_freq.drop(columns='total').transpose()
    surv_table = lineage_freq.add_suffix(' freq').join(no_uk_lineage_avg.add_suffix(' avg'))
    surv_table = sortAndTranspose(surv_table)
    surv_table['B.1.1.7 - UK avg'] = uk_lineage_avg
    surv_table['B.1.1.7 - UK freq'] = uk_lineage_freq
    surv_table.to_csv('results/surveillance_table.csv')
