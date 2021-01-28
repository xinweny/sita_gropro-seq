#### Packages ####
import argparse
import os
import re
import pandas as pd
import glob
from requests import get
from bs4 import BeautifulSoup

#### Functions ####
def scrape_sample_names(gse):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    response = get(url)

    html_soup = BeautifulSoup(response.text, 'html.parser')
    tds = [td.text for td in html_soup.findAll("td")]
    tables = html_soup.findAll('table')

    start = 18
    end = 20

    if 'NIH grant(s)' in tds:
        start += 1
        end += 1

    sample_table = tables[start:end]

    gsm_sample = {}

    for table in sample_table:
        for i, row in enumerate(table.findAll("tr")):
            cells = row.findAll("td")
            text = [cell.text for cell in cells]
            name = text[0]

            gsm_sample[name] = re.sub(r'[\/:]', "", text[1])

    return gsm_sample

def aspera_download(srr, ext, paired=False):
    if paired:
        sys_value = 0

        for i in [1, 2]:
            sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/00{srr[-1]}/{srr}/{srr}_{i}.fastq.gz \
                        ./fastq/{srr}_{i}.{ext}")

            if sys_val == 256:
                print("Retrying new ftp...")
                sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                            era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/0{srr[-2:]}/{srr}/{srr}_{i}.fastq.gz \
                            ./fastq/{srr}_{i}.{ext}")

            if sys_val == 256:
                print("Retrying new ftp...")
                sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                            era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/{srr}/{srr}_{i}.fastq.gz \
                            ./fastq/{srr}_{i}.{ext}")

            sys_value += sys_val

    else:
        sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                    era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/00{srr[-1]}/{srr}/{srr}.fastq.gz \
                    ./fastq/{srr}.{ext}")

        if sys_val == 256:
            print("Retrying new ftp...")
            sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/0{srr[-2:]}/{srr}/{srr}.fastq.gz \
                        ./fastq/{srr}.{ext}")

        if sys_val == 256:
            print("Retrying new ftp...")
            sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srr[0:6]}/{srr}/{srr}.fastq.gz \
                        ./fastq/{srr}.{ext}")

        sys_value = sys_val

    return sys_value

def samplename_from_metadata(gsm, metadata_df, columns):
    info = [metadata_df.loc[metadata_df['Sample Name'] == gsm, col].iloc[0] for col in columns]
    name = re.sub(r'[\/:]', "", '_'.join(info))

    return name

def is_paired(metadata, col_name, element):
    return metadata.loc[metadata[col_name] == element, 'LibraryLayout'].iloc[0] == "PAIRED"

def check_exitcodes(sys_values):
    if sys_values == {}:
        return True

    if sum(sys_values.values()) > 0:
        print("WARNING: Some files not fully downloaded:")
        for srr, value in sys_values.items():
            if value != 0:
                print(srr)

        return False
    else:
        print("All files downloaded completely.")

        return True

# Main
def main():
    # Parser
    parser = argparse.ArgumentParser(description='Download fq.gz files from SRA')
    parser.add_argument('-m', required=True, help='Path to SRA run table metadata.')
    parser.add_argument('-c', nargs='?', const='', help='Columns in metadata table for naming, separated by |. If given sample name on GSE page is preferred, leave as an empty string.')
    parser.add_argument('-e', nargs='?', const='fq.gz', default='fq.gz', help='Extension name (default: fq.gz)')

    args = parser.parse_args()
    metadata_path = args.m
    ext = args.e

    columns = (args.c).split('|')

    if args.c == '':
        gsm_samples = scrape_sample_names(os.getcwd().split("/")[-1])

    metadata_df = pd.read_csv(metadata_path, header=0, sep=',')
    gsms = metadata_df['Sample Name'].unique()

    # Set up file architecture
    os.system("mkdir fastq")

    for i, gsm in enumerate(gsms):
        print(f"({i + 1}/{len(gsms)}) Processing {gsm}...")

        if len([file for file in glob.glob(f"{gsm}*.{ext}")]) > 0:
            print(f"{gsm}*.fq.gz already downloaded. Skipping...")
            continue

        gsm_srrs = metadata_df.loc[metadata_df['Sample Name'] == gsm, 'Run']
        name = gsm_samples[gsm] if args.c == '' else samplename_from_metadata(gsm, metadata_df, columns)

        not_all_downloaded = True

        while not_all_downloaded:
            sys_values = {}

            for srr in list(gsm_srrs):
                print(f"Downloading {srr}...")
                sys_value = aspera_download(srr, ext, paired=is_paired(metadata_df, 'Run', srr))

                sys_values[srr] = sys_value

            not_all_downloaded = not check_exitcodes(sys_values)

        if len(gsm_srrs) == 1:
            print(f"Renaming fastq files for {gsm}...")
            srr = gsm_srrs.iloc[0]

            if is_paired(metadata_df, 'Run', srr):
                for i in [1, 2]:
                    os.system(f"mv -v '{srr}_{i}.{ext}' '{gsm}_{name}_{i}.{ext}'")
            else:
                os.system(f"mv -v '{srr}.{ext}' '{gsm}_{name}.{ext}'")
        else:
            print(f"Merging technical runs and renaming fastq files for {gsm}...")
            if is_paired(metadata_df, 'Run', srr):
                for i in [1, 2]:
                    single_fqs = ' '.join(gsm_srrs + f"_{i}.{ext}")
                    os.system(f"cat {single_fqs} > '{gsm}_{name}_{i}.{ext}' && rm {single_fqs}")
            else:
                srr_fqs = ' '.join(f"{gsm_srrs}.{ext}")
                os.system(f"cat {srr_fqs} > '{gsm}_{name}.{ext}' && rm {srr_fqs}")

#### Execute code ####
if __name__ == "__main__":
    main()
