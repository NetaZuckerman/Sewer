import pandas as pd
from sys import argv
"""
merge monitored mutations output tables.
get two tables as input and create a new one called mergedTable.csv
"""

table1 = pd.read_csv('all.csv')
table2 = pd.read_csv('ngs107.csv').set_index(["Position", "Reference", "Mutation", "protein", "variant", "Mutation type",
                                         "annotation", "varname", "lineage"])

merged_table = table1.join(table2, how='outer',
                           on=["Position", "Reference", "Mutation", "protein", "variant", "Mutation type",
                               "annotation", "varname", "lineage"]).fillna('X')

merged_table.to_csv('mergedTable.csv')
