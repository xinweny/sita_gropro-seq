#### Packages ####
import pandas as pd
import pybedtools as pbt
import argparse

#### Functions ####
def calculate_rpm(prom_df):
    # Make counts table
    counts_df = prom_df.drop(prom_df.iloc[:, 0:7], axis=1)

    # Calculate RPM
    print("Calculating RPM...")

    print("RPM scaling factors:")
    sf = counts_df.sum(axis=0) / 1000000

    print(sf)

    rpm_df = counts_df.div(sf, axis=1)

    # Replace raw promoter counts with RPM values
    for col in list(rpm_df.columns):
        prom_df[col] = rpm_df[col]

    return prom_df

def filter_transcripts(prom_df, gb_df, rpm_thresh, dist_thresh, gene_bed, chrom_sizes, sl, sr):
    # Make counts table
    counts_df = prom_df.drop(prom_df.iloc[:, 0:7], axis=1)
    samples = len(counts_df.columns)

    # Filter on promoter RPM in both control and treatmen
    print(f"Filtering transcripts by mean RPM > {rpm_thresh}...")

    prom_df['mean_rpm'] = counts_df.mean(axis=1)
    filt_rpm_df = prom_df[(prom_df['mean_rpm'] > rpm_thresh)]

    print(f"Keeping {len(filt_rpm_df)} out of {len(prom_df)} promoters ({len(filt_rpm_df) / len(prom_df) * 100}%).")

    # Select the most active, downstream transcript per gene
    print("Selecting the most active, downstream transcript per gene...")

    gene_ids = filt_rpm_df['name'].unique()

    transcripts = []

    # Sort transcripts of each gene by mean RPM first, then most downstream depending on strand if > 1 transcript
    for g in gene_ids:
        transcripts_df = filt_rpm_df[filt_rpm_df['name'] == g]

        if len(transcripts_df) > 1:
            strand = transcripts_df['strand'].unique()[0]

            if strand == '+':
                transcripts_df = transcripts_df.sort_values(['mean_rpm', 'start'], ascending=[False, False])
            elif strand == '-':
                transcripts_df = transcripts_df.sort_values(['mean_rpm', 'end'], ascending=[False, True])

            transcripts.append(transcripts_df.iloc[0, :]['transcript_id'])
        else:
            transcripts.append(transcripts_df['transcript_id'].values[0])

    # Filter gene body counts with filtered transcript ID list
    print("Filtering gene body counts file with filtered transcript ID list...")
    filt_gb_df = gb_df[gb_df['transcript_id'].isin(transcripts)]

    print(f"Keeping {len(filt_gb_df)} out of {len(gb_df)} transcripts ({len(filt_gb_df) / len(gb_df) * 100}%).")

    # Make pybedtools objects from BED files and slop to restore original transcript coordinates
    print(f"Keeping gene body transcripts with min distance to closest gene > {dist_thresh}...")
    filt_gb_bt = pbt.BedTool.from_dataframe(filt_gb_df.iloc[:, 0:7])
    filt_gb_bt = filt_gb_bt.slop(g=chrom_sizes, s=True, l=sl, r=sr)
    filt_gb_bt = filt_gb_bt.sort()

    # Load gene BED file as data frame
    gene_df = pd.read_csv(gene_bed, sep='\t', names=['chr', 'start', 'end', 'name', 'score', 'strand'])

    # Remove genes which have transcripts in filt_gb_df
    gene_df = gene_df[~gene_df['name'].isin(filt_gb_df['name'])]

    # Calculate distance to closest transcript using gene BED file 
    gene_bt = pbt.BedTool.from_dataframe(gene_df).sort()
    closest_bt = filt_gb_bt.closest(gene_bt, D='a', N=True)
    closest_df = closest_bt.to_dataframe()

    # Filter for genes with closest distance > threshold
    closest_df = closest_df[(closest_df[closest_df.columns[-1]] > dist_thresh) |
                            (closest_df[closest_df.columns[-1]] > -dist_thresh)]

    transcripts = list(closest_df[closest_df.columns[6]])

    # Filter gene body counts with filtered transcript ID list
    filt_gb_df = filt_gb_df[filt_gb_df['transcript_id'].isin(transcripts)]

    print(f"Keeping {len(filt_gb_df)} out of {len(gb_df)} transcripts ({len(filt_gb_df) / len(gb_df) * 100}%).")

    filt_gb_df = filt_gb_df.drop(['transcript_id'], axis=1)

    return filt_gb_df
    
def main():
    # Parser
    parser = argparse.ArgumentParser(description='Filter gene body transcripts for PRO-seq analysis.')
    parser.add_argument('-p', required=True, help='Path to the promoter counts table.')
    parser.add_argument('-b', required=True, help='Path to gene body counts file')
    parser.add_argument('-r', required=True, type=float, help='RPM threshold for filtering.')
    parser.add_argument('-g', required=True, help='Path to gene BED file for distance filtering.')
    parser.add_argument('-c', required=True, help='Path to genome chrom.sizes file.')
    parser.add_argument('-d', required=True, type=int, help='Min distance (bp) to closest transcript for filtering.')
    parser.add_argument('-sl', required=True, type=int, help='TSS shift (bp).')
    parser.add_argument('-sr', required=True, type=int, help='TES shift (bp).')
    parser.add_argument('-o', required=True, help='Path to output file.')

    args = parser.parse_args()

    prom_df = pd.read_csv(args.p, header=0, sep='\t')
    gb_df = pd.read_csv(args.b, header=0, sep='\t')

    # Calculate RPM across raw counts table
    rpm_df = calculate_rpm(prom_df)

    # Filter transcripts to select the most active and downstream TSS
    filt_df = filter_transcripts(prom_df=prom_df, gb_df=gb_df,
                                 rpm_thresh=args.r, dist_thresh=args.d,
                                 gene_bed=args.g, chrom_sizes=args.c,
                                 sl=args.sl, sr=args.sr)

    # Write filtered count table to output
    print(f"Saving output {args.o}...")

    filt_df.to_csv(args.o, index=False, sep='\t')

#### Execute code ####
if __name__ == "__main__":
    main()