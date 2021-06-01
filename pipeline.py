import pysam
from collections import Counter
import pandas as pd
import numpy as np
from sys import argv
import glob
import os
import errno
import re
import traceback


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
    """
    this is the sort function for the surveillance_table file
    the function sort the env samples from low to high
    the function called from sortAndTranspose function
    """
    try:
        a = re.findall('\d+|\D+', elem)
        return int(a[1])
    except ValueError:
        return 0


def sortAndTranspose(df):
    """
    This function sort and transpose the surveillance_table file by the keysort function
    """
    # taking the desired columns order from a file in the server
    df = df.reindex(columns=[x[:-1] for x in open("/data/projects/Dana/scripts/Sewer/Lineages_ordered.txt", "r")])
    df = df.transpose()
    try:
        df = df[sorted(df.columns, key=keysort)]
    except:
        print("error in columns sorting")
    df = df.transpose()
    return df


def no_uk_calculate(no_uk_df, other_variants):
    """
    this function calculate freq and average for the non-uk mutations and lineages
    (the non-uk doesn't include b.1.1.7 mutations)
    For the surveillance_table.csv
    :param no_uk_df: the monitored_mutations DataFrame, about to filter out the b.1.1.7 mutations
    :param other_variants:  list of all the mutations that doesn't belong to b.1.1.7 mutations
    :return:
    """
    # filtering mutations
    no_uk_df = no_uk_df[(no_uk_df.variant.isin(other_variants))]
    # calculate avg and std
    # avg and std doesn't include NA's in the calculation
    lineage_avg = no_uk_df.drop('Position', axis=1).groupby('lineage').mean().transpose()
    # round to 2 digits after decimal point
    lineage_avg = round(lineage_avg, 2)
    lineage_std = no_uk_df.drop('Position', axis=1).groupby('lineage').std()
    # calculate frequency
    # number of total mutations (include NA's) per lineage
    lineage_num_muts = no_uk_df.groupby('lineage')['lineage'].count().to_frame().rename(columns={'lineage': 'total'})
    # temporarily Changes NA's to -1 to help the counting
    no_uk_df.fillna(-1, inplace=True)
    # number of non zero mutations (Greater than 0, not including NA's) per lineage
    lineage_non_zero_count = no_uk_df.drop(columns=['varname', 'variant', 'protein', 'Mutation type', 'Position', 'Reference', 'Mutation']) \
        .groupby('lineage').agg(lambda x: x.gt(0).sum())
    # number of zeros per lineage
    lineage_zero_count = no_uk_df.drop(columns=['varname', 'variant', 'protein', 'Mutation type', 'Position', 'Reference', 'Mutation']).groupby(
        'lineage').agg(lambda x: x.eq(0).sum())
    # number of NA's per lineage
    lineage_na_count = no_uk_df.drop(columns=['varname', 'variant', 'protein', 'Mutation type', 'Position', 'Reference', 'Mutation']).groupby(
        'lineage').agg(lambda x: x.eq(-1).sum())
    no_uk_df.replace(-1, None, inplace=True)
    lineage_freq = lineage_num_muts.join(lineage_non_zero_count)
    return lineage_freq, lineage_avg, lineage_std, lineage_zero_count, lineage_na_count


def uk_calculate(uk_df, uk_variant_mutations):
    """
    this function calculate freq and average for the uk lineage mutations
    (only B.1.1.7 mutations)
    For the surveillance_table.csv
    :param uk_df: the monitored_mutations DataFrame, about to remain with only uk mutations
    :param uk_variant_mutations:  the mutations of the uk variant
    :return:
    """
    # filtering mutations
    uk_df = uk_df[(uk_df.variant.isin(uk_variant_mutations))]
    # calculate avg and std
    # avg and std doesn't include NA's in the calculation
    lineage_avg = uk_df.drop('Position', axis=1).groupby('lineage').mean().transpose()
    # round to 2 digits after decimal point
    lineage_avg = round(lineage_avg, 2)
    lineage_std = uk_df.drop('Position', axis=1).groupby('lineage').std()
    # calculate frequency
    # number of total mutations (include NA's) per lineage
    lineage_num_muts = uk_df.groupby('lineage')['lineage'].count().to_frame().rename(columns={'lineage': 'total'})
    # temporarily Changes NA's to -1 to help the counting
    uk_df.fillna(-1, inplace=True)
    # number of non zero mutations (Greater than 0, not including NA's) per lineage
    lineage_non_zero_count = uk_df.drop(columns=['varname', 'variant', 'protein', 'Mutation type', 'Position', 'Reference', 'Mutation']) \
        .groupby('lineage').agg(lambda x: x.gt(0).sum())
    # number of zeros per lineage
    lineage_zero_count = uk_df.drop(columns=['varname', 'variant', 'protein', 'Mutation type', 'Position', 'Reference', 'Mutation']) \
        .groupby('lineage').agg(lambda x: x.eq(0).sum())
    # number of NA's per lineage
    lineage_na_count = uk_df.drop(columns=['varname', 'variant', 'protein', 'Mutation type', 'Position', 'Reference', 'Mutation']) \
        .groupby('lineage').agg(lambda x: x.eq(-1).sum())
    uk_df.replace(-1, None, inplace=True)
    lineage_freq = lineage_num_muts.join(lineage_non_zero_count)
    lineage_avg = lineage_avg['B.1.1.7']
    uk_total = lineage_freq['total']['B.1.1.7']
    lineage_std = lineage_std.loc['B.1.1.7', :].transpose()
    lineage_freq = lineage_freq.loc['B.1.1.7', :].transpose()
    lineage_zero_count = lineage_zero_count.loc['B.1.1.7', :].transpose()
    lineage_na_count = lineage_na_count.loc['B.1.1.7', :].transpose()
    # lineage_freq = lineage_freq.astype(int).astype(str) + '\\' + uk_total.astype(str)
    # Build the value in B.1.1.7 Freq column
    lineage_freq = "all:" + lineage_freq.astype(int).astype(str) + '\\' + uk_total.astype(str) + " ; (" + round(
        (lineage_freq / uk_total * 100), 2).astype(str) + "%," + " sd: " + round(lineage_std, 2).astype(
        str) + "); zero: " + lineage_zero_count.astype(int).astype(str) + '\\' + uk_total.astype(
        str) + "; NA:" + lineage_na_count.astype(int).astype(str) + '\\' + uk_total.astype(str)
    return lineage_freq, lineage_avg


def addVerdict(survTable):
    try:
        survTable.insert(1, 'verdict', "")
    except:
        print("verdict column is already exist")
    try:
        for index, row in survTable.iterrows():
            verList = []
            for (columnName, columnData) in row.iteritems():
                if "freq" in columnName:
                    # check if not nan
                    if columnData == columnData and columnData != 0 and columnData != '0':
                        freq = columnData.split(";")[1].split("%")[0][2:]
                        if float(freq) >= 60:
                            lineageName = str(columnName).split(" ")[0]
                            avgColName = lineageName + " avg"
                            lineageAvg = row[avgColName]
                            verList.append(lineageName + " " + str(lineageAvg) + "%")
                        elif 60 > float(freq) >= 35:
                            numOfZeros = int(columnData.split(";")[2].split(":")[1].split("\\")[0])
                            total = int(columnData.split(";")[2].split(":")[1].split("\\")[1])
                            if numOfZeros / total * 100 < 10:
                                lineageName = str(columnName).split(" ")[0]
                                avgColName = lineageName + " avg"
                                lineageAvg = row[avgColName]
                                verList.append("Suspect: " + lineageName + " " + str(lineageAvg) + "%")
            if len(verList) > 0:
                toSurv = ' '.join(verList)
            else:
                toSurv = "Undetermined"
            survTable["verdict"][index] = toSurv
        return survTable
    except:
        print("Data: "+str(columnData))
        traceback.print_exc()

        
        
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
    excel_mutTable = pd.read_excel("/data/projects/Dana/scripts/covid19/mutationsTable.xlsx", sheet_name=None,engine='openpyxl')
    # muttable = pd.read_csv("novelMutTable.csv") # TODO change before commit
    mutTable_copy = excel_mutTable.copy()
    for name, frame in mutTable_copy.items():
        frame['Mutation type'] = frame['Mutation type'].str.lower()  # make sure lower case to prevent mistakes
        excel_mutTable[name] = frame[frame['Mutation type'] != 'insertion']
        excel_mutTable[name]['lineage'] = name  # add a lineage column to all variant's tables
    muttable = pd.concat(excel_mutTable.values(), ignore_index=True)
    uniq_lineages = set()
    for lin in muttable.lineage:
        lin=lin.rsplit('_',1)
        uniq_lineages.add(lin[0].strip())
    muttable_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}
    for lin, table in muttable_by_lineage.items():
        table.lineage = lin
    # table = table.assign(lineage=lin)

    final_df = pd.concat([frame for frame in muttable_by_lineage.values()])
    # getting list of all mutations
    all_mutations = set([x for x in muttable.variant])
    # only uk mutations
    uk_variant_mutations = set(muttable_by_lineage['B.1.1.7']['variant'])  # list of mutations of uk variant
    # all of the non-uk mutations
    other_variants = all_mutations - uk_variant_mutations
    # L5F doesn't included in both lists
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
        final_df[file_name] = final_df.apply(lambda row: frequency(row['Mutation'], row['Position'] - 1, pileup_table, min_depth),
                                             axis=1)

    # sorting and organizing the monitored_mutations file
    final_df = final_df.sort_values(["lineage", "protein"], ascending=(True, False))  # sort by:(1)lineage (2)gene(S first)
    sortednames = sorted([x for x in final_df.columns.values if "nv" in x], key=keysort)
    sorted_cols = [c for c in final_df.columns.values if c not in sortednames] + sortednames
    final_df = final_df.reindex(columns=sorted_cols)
    # creating folders
    if not os.path.exists('results'):
        os.mkdir('results/')
    final_df = final_df.drop(['Unnamed: 6','% of sequences'], axis=1)
    monitoredfile = final_df.copy()
    # replacing NA's with "No Coverage" Text
    monitoredfile.fillna(-1, inplace=True)
    monitoredfile.replace(-1, "No Coverage", inplace=True)
    monitoredfile.to_csv("results/monitored_mutations.csv")
    # Folders for the pileups
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

    # Start Handling the surveillance_table.csv
    no_uk_lineage_freq, no_uk_lineage_avg, no_uk_lineage_std, no_uk_zero, no_uk_na = no_uk_calculate(final_df.copy(),
                                                                                                     other_variants)
    uk_lineage_freq, uk_lineage_avg = uk_calculate(final_df.copy(), uk_variant_mutations)

    for name in all_tables.keys():
        # lineage_freq[name] /= lineage_freq['total']/100
        try:
            # Build the value in Lineage Freq column
            no_uk_lineage_freq[name] = "all:" + no_uk_lineage_freq[name].astype(int).astype(str) + '\\' + \
                                       no_uk_lineage_freq[
                                           'total'].astype(str) + " ; (" + round(
                (no_uk_lineage_freq[name] / no_uk_lineage_freq['total'] * 100),
                2).astype(str) + "%, sd:" + round(no_uk_lineage_std[name], 2).astype(str) + "); zero: " + no_uk_zero[
                                           name].astype(int).astype(str) + '\\' + no_uk_lineage_freq[
                                           'total'].astype(str) + "; NA:" + no_uk_na[name].astype(int).astype(
                str) + '\\' + no_uk_lineage_freq[
                                           'total'].astype(str)
        except:
            print(name)

    lineage_freq = no_uk_lineage_freq.drop(columns='total').transpose()
    surv_table = lineage_freq.add_suffix(' freq').join(no_uk_lineage_avg.add_suffix(' avg'))
    #surv_table = sortAndTranspose(surv_table)
    surv_table['B.1.1.7 avg'] = uk_lineage_avg
    surv_table['B.1.1.7 freq'] = uk_lineage_freq
    #addVerdict(surv_table)
    surv_table.to_csv('results/surveillance_table.csv')
