import pysam
from collections import Counter
import pandas as pd
import numpy as np
from sys import argv
import glob
import os
import errno



if __name__ == '__main__':
    file = 'BAM/Env1.mapped.sorted.bam'
    refseq_series = pd.Series([x for x in pysam.Fastafile('refs/REF_NC_045512.2.fasta').fetch(reference="REF_NC_045512.2")])

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
    pileup_table.to_csv('temp_pileuptable.csv')  # to remove after debug
    pileup_table['ref_freq'] = pileup_table.apply(
        lambda row: (row[row['ref']] / row['sum']) * 100 if row['sum'] else 0.0, axis=1)
    pileup_table = pileup_table.dropna()
    pileup_table.to_csv()