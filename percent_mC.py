#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.ndimage import gaussian_filter1d
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import binomtest

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate CpG methylation plots from CSV data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python3 percent_mC.py input.csv output.pdf
  python3 percent_mC.py input.csv output.pdf -snp
        '''
    )
    
    parser.add_argument('input_csv', help='Input CSV file containing methylation data')
    parser.add_argument('output_pdf', help='Output PDF file for the plot')
    parser.add_argument('-snp', action='store_true', 
                       help='Add vertical lines at predefined SNP positions (397, 492) with zygosity analysis')
    
    return parser.parse_args()

def add_cpg_methylation_plot(input_csv, output_pdf, add_snp_lines=False, positions=None):
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv)
        
        # Calculate percent C (unmethylated) and percent T (methylated)
        total_coverage = df['A'] + df['C'] + df['G'] + df['T']
        df['percent_C'] = df['C'] / total_coverage * 100  # Unmethylated
        df['percent_T'] = df['T'] / total_coverage * 100  # Methylated
        
        zygosity_results = []
        
        # Only perform zygosity analysis if SNP lines are requested
        if add_snp_lines and positions:
            # Analyze A, G content using binomial test
            for pos in positions:
                row = df[df['POS'] == pos]
                if not row.empty:
                    a, g = row[['A', 'G']].values[0]
                    total = a + g
                    if total > 0:
                        result = binomtest(a, n=total, p=0.5)
                        p_val = result.pvalue
                        if p_val < 0.05:
                            genotype = "A-G"
                        else:
                            genotype = "A-A" if a > g else "G-G"
                        zygosity_results.append((pos, genotype, p_val))
                    else:
                        zygosity_results.append((pos, "No A/G", None))
                else:
                    zygosity_results.append((pos, "Not found", None))

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
                fig2, ax2 = plt.subplots(figsize=(10, 5))
                fig2.subplots_adjust(bottom=0.4)
                
                # Smooth the C percentage data
                sigma = 1.5  # Smoothing factor
                smoothed_C = gaussian_filter1d(cpg_df['percent_C'], sigma)
                
                ax2.plot(cpg_df['POS'], smoothed_C, 'b-', linewidth=2, label='Smoothed %C')
                ax2.scatter(cpg_df['POS'], cpg_df['percent_C'], 
                           c='blue', alpha=0.3, label='Raw %C')
                
                # Add vertical lines and labels at specified positions (only if requested)
                if add_snp_lines:
                    for pos, genotype, p_val in zygosity_results:
                        if isinstance(p_val, float):
                            if genotype == "A-G":
                                color = "green"
                            elif genotype == "A-A":
                                color = "blue"
                            elif genotype == "G-G":
                                color = "red"
                            else:
                                color = "gray"

                            ax2.axvline(x=pos, color=color, linestyle='--', linewidth=1.5, alpha=0.7)

                            label = f"{genotype}\nPosition: {pos}\np = {p_val:.3g}"
                            ax2.annotate(label,
                                xy=(pos, 0), xycoords='data',
                                xytext=(0, -50), textcoords='offset points',
                                ha='center', va='top',
                                fontsize=9, color=color,
                                bbox=dict(facecolor='white', edgecolor=color, boxstyle='round,pad=0.3'))
                
                # Optional: Add T percentage for comparison
                # smoothed_T = gaussian_filter1d(cpg_df['percent_T'], sigma)
                # ax2.plot(cpg_df['POS'], smoothed_T, 'r--', linewidth=1, label='Smoothed %T')
                
                ax2.set_xlabel('Position')
                ax2.set_ylabel('Percent Methylation')
                ax2.set_title('Percent methylation per CpG')
                ax2.set_ylim(-10, 100)
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

def main():
    args = parse_arguments()
    
    # Use predefined positions if SNP analysis is requested
    predefined_positions = [397, 492]
    add_snp_lines = args.snp
    positions = predefined_positions if add_snp_lines else None
    
    add_cpg_methylation_plot(args.input_csv, args.output_pdf, add_snp_lines, positions)

if __name__ == "__main__":
    main()