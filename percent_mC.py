ce#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.ndimage import gaussian_filter1d
from matplotlib.backends.backend_pdf import PdfPages

def add_cpg_methylation_plot(input_csv, output_pdf):
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv)
        
        # Calculate percent C (unmethylated) and percent T (methylated)
        total_coverage = df['A'] + df['C'] + df['G'] + df['T']
        df['percent_C'] = df['C'] / total_coverage * 100  # Unmethylated
        df['percent_T'] = df['T'] / total_coverage * 100  # Methylated
        
        # Create a PDF file
        with PdfPages(output_pdf) as pdf:
            # Plot CpG-specific C percentage (smoothed)
            df = df.sort_values('POS')  # Ensure positions are ordered
            df['next_pos'] = df['POS'].shift(-1)
            df['next_ref'] = df['REF'].shift(-1)
            
            # Find CpGs (C followed by G at next position)
            cpg_df = df[(df['REF'] == 'C') & 
                       (df['next_ref'] == 'G') &
                       (df['next_pos'] == df['POS'] + 1)]
            
            if len(cpg_df) > 1:
                fig2, ax2 = plt.subplots(figsize=(15, 5))
                
                # Smooth the C percentage data
                sigma = 1.5  # Smoothing factor
                smoothed_C = gaussian_filter1d(cpg_df['percent_C'], sigma)
                
                ax2.plot(cpg_df['POS'], smoothed_C, 'b-', linewidth=2, label='Smoothed %C')
                ax2.scatter(cpg_df['POS'], cpg_df['percent_C'], 
                           c='blue', alpha=0.3, label='Raw %C')
                
                # Optional: Add T percentage for comparison
                # smoothed_T = gaussian_filter1d(cpg_df['percent_T'], sigma)
                # ax2.plot(cpg_df['POS'], smoothed_T, 'r--', linewidth=1, label='Smoothed %T')
                
                ax2.set_xlabel('Position')
                ax2.set_ylabel('Percent Methylation')
                ax2.set_title('Percent methylation per CpG')
                ax2.set_ylim(0, 100)
                ax2.grid(True, linestyle='--', alpha=0.3)
                ax2.legend()
                
                pdf.savefig(fig2)
                plt.close(fig2)
            else:
                print(f"Found only {len(cpg_df)} CpG sites - not enough for smoothing")
                print("Debug info:")
                print("Total positions:", len(df))
                print("C positions:", len(df[df['REF'] == 'C']))
                print("G positions:", len(df[df['REF'] == 'G']))
    
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        raise

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python3 percent_mC.py <input_csv> <output_pdf>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_pdf = sys.argv[2]
    
    add_cpg_methylation_plot(input_csv, output_pdf)
