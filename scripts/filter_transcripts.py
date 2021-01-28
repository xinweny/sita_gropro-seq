#### Packages ####
import pandas as pd
import argparse

#### Functions ####
def calculate_rpm(prom_df):
    # Make counts table
    counts_df = prom_df.drop(prom_df.iloc[:, 0:7], axis=1)

    # Calculate RPM
    sf = counts_df.sum(axis=0) / 1000000
    rpm_df = counts_df.div(sf, axis=1)

    # Replace raw counts with RPM values
    for col in list(rpm_df):
        prom_df[col] = rpm_df[col]

    return prom_df

def filter_transcripts(prom_df, rpm_thresh):
    # Make counts table
    counts_df = prom_df.drop(prom_df.iloc[:, 0:7], axis=1)
    samples = len(counts_df.columns)

def main():
    # Parser
    parser = argparse.ArgumentParser(description='Filter transripts based on promoter activeness and downstreamness')
    parser.add_argument('-i', required=True, help='Path to the promoter counts table.')
    parser.add_argument('-r', required=True, help='RPM threshold for filtering')

    args = parser.parse_args()

    prom_df = pd.read_csv(args.i, header=0, sep='\t')
    print(prom_df)
    rpm_thresh = float(args.r)

    rpm_df = calculate_rpm(prom_df)
    print(rpm_df)


#### Execute code ####
if __name__ == "__main__":
    main()