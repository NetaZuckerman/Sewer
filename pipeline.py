import pysam
from collections import Counter
import pandas as pd
import numpy as np
from sys import argv
import glob
"""
for sewer samples.
iterate over mapped and sorted bam files,
calculate mutations frequencies by pileup files.
"""


def frequency(mut_val, pos, pileup_df):
    """
    return frequency of mut_val base in specified position.
    :param mut_val: mutation nucleotide (A,C,T,G,-)
    :param pos: position
    :param pileup_df: the pileup dataframe
    :return: frequency of mutation nucleotide in position (mut_depth/sum*100)
    """
    mut_val = 'del' if mut_val == '-' else mut_val
    if pileup_df.loc[pos-1]['sum']:
        val = pileup_df.loc[pos-1][mut_val] / pileup_df.loc[pos-1]['sum']
    else:
        val = 0.0
    return val*100

# TODO: is there a faster way of pileup?


if __name__ == '__main__':
    # inputs
    bam_dir = argv[1]
    out_file = argv[2]
    muttable = pd.read_csv("novelMutTable.csv")  # TODO: get from other location!

    uniq_lineages = set()
    for lin in muttable.lineage:
        for x in lin.split(','):
            uniq_lineages.add(x.strip())
    muttable_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}
    for lin, table in muttable_by_lineage.items():
        table.lineage = lin

    final_df = pd.concat([frame for frame in muttable_by_lineage.values()])

    all_tables = []
    # get BAM location and iterate through it
    files_list = glob.glob(bam_dir + '/*.mapped.sorted.bam')

    for file in files_list:
        pileup_table = pd.DataFrame(np.zeros(shape=(29903, 5)), columns=['C', 'A', 'G', 'T', 'del'],
                                    index=list(range(29903)))
        bam = pysam.AlignmentFile(file, 'rb')
        pileup_iter = bam.pileup(stepper='nofilter')
        for position in pileup_iter:
            c = Counter({'C': 0, 'A': 0, 'G': 0, 'T': 0, 'del': 0})
            for pileupread in position.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    c[pileupread.alignment.query_sequence[pileupread.query_position].upper()] += 1
                if pileupread.is_del:
                    c['del'] += 1
            pileup_table.loc[position.reference_pos + 1] = pd.Series(c)
        pileup_table.index.name = 'pos'
        pileup_table['sum'] = pileup_table['A'] + pileup_table['C'] + pileup_table['T'] + pileup_table['G'] + \
                              pileup_table['del']

        all_tables.append(pileup_table)
        file_name = file.strip('.mapped.sorted.bam')
        final_df[file_name] = final_df.apply(lambda row: frequency(row['mut'], row['pos'], pileup_table), axis=1)

    final_df.to_csv(out_file)

