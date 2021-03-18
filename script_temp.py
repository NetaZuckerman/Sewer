from collections import Counter
import pandas as pd
import numpy as np

pileup_df = pd.read_csv("pileup.csv")
muttable = pd.read_csv("novelMutTable.csv")
uniq_lineages = set()
for lin in muttable.lineage:
    for x in lin.split(','):
        uniq_lineages.add(x.strip())

mutations_by_lineage = {x: muttable[muttable.lineage.str.contains(x)] for x in uniq_lineages}

pileup_df['postotal'] = pileup_df['A'] + pileup_df['C'] + pileup_df['T'] + pileup_df['G'] + pileup_df['isindel']

frequency_a = pileup_df