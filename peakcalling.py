import os
import pandas as pd
import numpy as np

# Initial settings
input_path = "results"
output_file = "totalcombined_snps_summary_HN_sorghum.csv"
gene_info_file = "Sbicolor_genes_positions.csv"  # File with gene information

# Total number of SNPs for threshold of significance
total_snps = 4560009

# Statistical significance threshold
p_cutoff = -np.log10(0.05 / total_snps)

# Load gene information
gene_info = pd.read_csv(gene_info_file)
# Extend gene regions by 1 Mb on both sides
gene_info['start_extended'] = gene_info['start'] - 1e6
gene_info['end_extended'] = gene_info['end'] + 1e6

# Ensure chromosome columns are of the same type (string)
gene_info['chr'] = gene_info['chr'].astype(str)
new_chr = []
for achr in gene_info['chr']:
    if not "scaffold" in achr:
        new_chr.append("Chr{0}".format(str(int(achr)).zfill(2)))
    else:
        new_chr.append(achr)
gene_info['chr'] = new_chr

# Peak identification function
def identify_combined_peaks(df):
    df = df.sort_values(by=['CHROM', 'POS'], inplace=False)
    peaks = []
    current_peak_start = None
    current_peak_end = None
    current_chrom = None
    peak_id = 1
    
    for _, row in df.iterrows():
        if current_peak_start is None:
            current_peak_start = row['POS']
            current_peak_end = row['POS']
            current_chrom = row['CHROM']
            peaks.append(peak_id)
        else:
            if row['CHROM'] == current_chrom and row['POS'] - current_peak_end <= 1e6:
                current_peak_end = row['POS']
                peaks.append(peak_id)
            else:
                peak_id += 1
                current_peak_start = row['POS']
                current_peak_end = row['POS']
                current_chrom = row['CHROM']
                peaks.append(peak_id)
    
    df['peak_id'] = peaks
    return df

# Create a list to store the peak summaries
peaks_summary_list = []

# Process each file in the results directory
for filename in os.listdir(input_path):
    if filename.endswith("_2023_HN_333filtered.csv.gz"):
        file_path = os.path.join(input_path, filename)
        print(f"Processing file: {filename}")
        
        # Read compressed CSV file
        df = pd.read_csv(file_path, compression='gzip')
        
        # Check if the expected columns are present
        required_columns = {'SNP', 'CHROM', 'POS', 'REF', 'ALT', 'MAF','Effect', 'SE', 'pval', 'BPcum'}
        if not required_columns.issubset(df.columns):
            print(f"Warning: {filename} does not contain all the necessary columns. Skipping.")
            continue
        
        # Calculate -log10(p-value)
        df['-log10Pvalue'] = -np.log10(df['pval'].astype(float))
        
        # Filter out significant SNPs
        significant_snps = df[df['-log10Pvalue'] > p_cutoff].copy()
        
        if significant_snps.empty:
            print(f"No significant SNPs in {filename}")
            continue
        
        significant_snps.reset_index(drop=True, inplace=True)
        
        # Add a column to identify the gene (based on the file name)
        gene_id = filename.split("_")[0]
        significant_snps['gene_id'] = gene_id
        
        # Identifying peaks
        significant_snps = identify_combined_peaks(significant_snps)
        
        # Select the most significant SNP per peak
        top_snps_per_peak = significant_snps.loc[significant_snps.groupby('peak_id')['pval'].idxmin()].rename(
            columns={'SNP': 'top_SNP', 'pval': 'top_Pvalue'}
        )[['peak_id', 'CHROM', 'POS', 'top_SNP', 'top_Pvalue']]
        
        # Compute peak range, associated traits, and SNP count
        peak_ranges = significant_snps.groupby('peak_id').agg(
            pStart=('POS', 'min'),
            pStop=('POS', 'max'),
            pLength=('POS', lambda x: x.max() - x.min()),
            num_SNPs=('POS', 'count'),
            traits=('gene_id', lambda x: ', '.join(x.unique()))
        ).reset_index()

        # Merge top SNP info with peak ranges
        peaks_summary = top_snps_per_peak.merge(peak_ranges, on='peak_id')

        # Merge gene information based on 'traits'
        peaks_summary = peaks_summary.merge(
            gene_info, left_on='traits', right_on='gene', how='left'
        )

        # Ensure chromosome columns are of the same type (string)
        peaks_summary['CHROM'] = peaks_summary['CHROM'].astype(str)

        # Determine if the peak is in cis or trans
        peaks_summary['cis_trans'] = peaks_summary.apply(
            lambda row: 'cis' if (
                row['CHROM'] == row['chr'] and
                row['pStart'] <= row['end_extended'] and
                row['pStop'] >= row['start_extended']
            ) else 'trans', axis=1
        )

        # Append summary to the list
        peaks_summary_list.append(peaks_summary)

# Combine all peak summaries into a single DataFrame and save to CSV
if peaks_summary_list:
    combined_peaks_summary = pd.concat(peaks_summary_list)
    combined_peaks_summary.to_csv(output_file, index=False)
    print(f"Summary file saved at: {output_file}")
else:
    print("No significant SNPs were found in the processed files.")
