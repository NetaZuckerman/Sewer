from pathlib import Path
import pandas as pd
from argparse import ArgumentParser
"""
merge monitored mutations output tables.
get two tables as input and create a new one called mergedTable.csv
"""

SEWER_PATH = Path('/data3/sewer')
new_runs = ['NGS110_18082021', 'NGS111_20082021', 'NGS113_27082021']



def configure_parser():
    parser = ArgumentParser()
    parser.add_argument('--runs',
                        help="Run NUMBER",
                        nargs='+',
                        type=str
                        )
    
    parser.add_argument(
        '--monitored', '-m',
        action='store_true',
        default=False
        )
    
    parser.add_argument(
        '--surveillance', '-s',
        action='store_true',
        default=False
        )
    
    return parser


def merge_monitored():
    print('Merging monitored variants')
    monitored_path = SEWER_PATH / 'monitored.xlsx'
    table1 = pd.read_excel(monitored_path, engine='openpyxl')
    
    for run in new_runs:
        run_path = SEWER_PATH / run / 'results' / 'monitored_mutations.csv'
        table2 = pd.read_csv(run_path).set_index(["Position", "Reference", "Mutation", "protein", "variant", "Mutation type", "annotation", "varname", "lineage"])
        table1 = table1.join(table2, how='outer',
                               on=["Position", "Reference", "Mutation", "protein", "variant", "Mutation type",
                                   "annotation", "varname", "lineage"]).fillna('X')
    
    table1.to_csv('mergedTable.csv', index=False)
    

def merge_surv():
    print('Merging surveillanced variants')
    env_surv_path = SEWER_PATH / 'EnvSurv_excel.xlsx'
    env_surv_df = pd.read_excel(env_surv_path)
    samples_ind = [
        num[0] for num in env_surv_df['sample number'] 
            .str.strip()
            .str.findall('\D(\d+)$')
            ]
    
    env_surv_df = env_surv_df.set_axis(samples_ind, axis=0)
    
    for run in new_runs:
        run_path = SEWER_PATH / run / 'results' / 'surveillance_table.csv'
        surv = pd.read_csv(run_path, index_col=0)
        samples_ind = [
        num[0] for num in surv.index
            .str.strip()
            .str.findall('\D(\d+)$')
            ]
        
        surv = surv.set_axis(samples_ind, axis=0)
        missing_ind = surv.index[~surv.index.isin(env_surv_df.index)]
        if len(missing_ind) > 0:
            missing_s = ' '.join(missing_ind)
            s = f"{run}: Folloing samples were not found in main table: {missing_s}"
            print(s)
            surv = surv.loc[~missing_ind]
        
        env_surv_df.loc[surv.index, surv.columns] = surv
    output_path = SEWER_PATH / 'updated_envsurv.xlsx'
    env_surv_df.to_excel(output_path)

# if __name__ == "__main__":
#     parser = configure_parser()
#     args = parser.parse_args()
    
merge_monitored()
merge_surv()
    