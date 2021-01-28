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

def filter_transcripts(prom_df, gene_body_df, rpm_thresh):
    # Make counts table
    counts_df = prom_df.drop(prom_df.iloc[:, 0:7], axis=1)
    samples = len(counts_df.columns)

    # Filter on promoter activeness
    prom_df['mean_rpm'] = counts_df.mean(axis=1)
    filt_df = prom_df[prom_df['mean_rpm'] > rpm_thresh]

    # Select the most active, downstream transcript per gene
    gene_ids = filt_df['name'].unique()

    transcripts = []

    for g in gene_ids:
        transcripts_df = filt_df[filt_df['name'] == g]

        if len(transcripts_df) > 1:
            strand = transcripts_df['strand'].unique()[0]

            if strand == '+':
                transcripts_df = transcripts_df.sort_values(['mean_rpm', 'start'], ascending=[False, False])
            elif strand == '-':
                transcripts_df = transcripts_df.sort_values(['mean_rpm', 'end'], ascending=[False, True])

            transcripts.append(transcripts_df.iloc[0, :]['transcript_id'])
        else:
            transcripts.append(transcripts_df['transcript_id'].values[0])

    # Filter out gene body BED file with filtered transcript ID list
    filt_bed_df = gene_body_df[gene_body_df['transcript_id'].isin(transcripts)]
    filt_bed_df = filt_bed_df.drop(['transcript_id'], axis=1)

    return filt_bed_df
    
def main():
    # Parser
    parser = argparse.ArgumentParser(description='Filter transripts based on promoter activeness and downstreamness')
    parser.add_argument('-p', required=True, help='Path to the promoter counts table.')
    parser.add_argument('-b', required=True, help='Path to gene body BED file')
    parser.add_argument('-r', required=True, help='RPM threshold for filtering.')
    parser.add_argument('-o', required=True, help="Path to output file.")

    args = parser.parse_args()

    prom_df = pd.read_csv(args.p, header=0, sep='\t')
    gene_body_df = pd.read_csv(args.b, sep='\t', names=prom_df.columns[0:7])

    rpm_thresh = float(args.r)

    # Calculate RPM across raw counts table
    rpm_df = calculate_rpm(prom_df)

    # Filter transcripts to select the most active and downstream TSS
    filt_df = filter_transcripts(prom_df, gene_body_df, rpm_thresh)

    # Write to BED output
    filt_df.to_csv(args.o, header=False, sep='\t')

#### Execute code ####
if __name__ == "__main__":
    main()