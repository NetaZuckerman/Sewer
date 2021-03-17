import pysam
from collections import Counter
import pandas as pd
import numpy as np

file = pysam.AlignmentFile('s6519_Ashdod_Jan1.mapped.sorted.bam', "rb")

iter = file.pileup()
# positions_dict = {i: dict() for i in range(0, 29903)}
# mystring = ''
# mylist = []
# for pileupcol in file.pileup():
#     mystring = ''
#     for pileupread in pileupcol.pileups:
#         if not pileupread.is_del and not pileupread.is_refskip:
#             pos = pileupread.query_position
#             # print(pos)
#             if pileupread.alignment.query_sequence[pos] not in positions_dict[pos]:
#                 positions_dict[pos][pileupread.alignment.query_sequence[pos]] = 0
#             positions_dict[pos][pileupread.alignment.query_sequence[pos]] += 1
#             # mystring = mystring + pileupread.alignment.query_sequence[pileupread.query_position]
#
#     # mystring = "".join(sorted(list(mystring)))
#     # mylist.append((pileupcol.pos+1, mystring))
#     #
# print(positions_dict)
df_pysam = pd.DataFrame(np.zeros(shape=(29903, 5)), columns=['C', 'A', 'G', 'T', 'isindel'], index=list(range(29903)))

MTpileup = file.pileup(stepper='nofilter')
for position in MTpileup:
    c = Counter({'C': 0, 'A': 0, 'G': 0, 'T': 0, 'isindel': 0})
    for pileupread in position.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            c[pileupread.alignment.query_sequence[pileupread.query_position].upper()] += 1
        if pileupread.indel:
            c['isindel'] += 1
    df_pysam.loc[position.reference_pos+1] = pd.Series(c)

file.close
print(c)