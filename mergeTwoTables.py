import pandas as pd
from sys import argv
"""
merge compressed and not compressed tables.
get two tables as input and create a new one called mergedTable.csv
"""

table1 = pd.read_csv(argv[1])
table2 = pd.read_csv(argv[2]).set_index(["Position", "Reference", "Mutation", "protein", "variant", "Mutation type",
                                         "annotation", "varname", "lineage"])

merged_table = table1.join(table2, how='outer', on=["Position", "Reference", "Mutation", "protein", "variant",
                                                    "Mutation type", "annotation", "varname", "lineage"])
merged_table.to_csv('mergedTable.csv')
