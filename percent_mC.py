
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
        description='Generate CpG methylation plots from multiple CSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 percent_mC.py output.pdf --case file1.csv file2.csv --control file3.csv file4.csv
  python3 percent_mC.py output.pdf --case file1.csv --control file2.csv -snp
        """
    )
    
    parser.add_argument('output_pdf', help='Output PDF file for the plot')
    parser.add_argument('--case', nargs='+', required=True,
                       help='Input CSV files for case group')
    parser.add_argument('--control', nargs='+', required=True,
                       help='Input CSV files for control group')
    parser.add_argument('-snp', action='store_true', 
                       help='Add vertical lines at predefined SNP positions (397, 492) with zygosity analysis')
    parser.add_argument('--color', action='store_true',
                       help='Use different colors for each file instead of grouping by case/control')
    
    return parser.parse_args()

def process_csv_file(input_csv):
    """Process a single CSV file and return CpG data"""
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv)
        
        # Calculate percent C (unmethylated) and percent T (methylated)
        total_coverage = df['A'] + df['C'] + df['G'] + df['T']
        df['percent_C'] = df['C'] / total_coverage * 100  # Unmethylated
        df['percent_T'] = df['T'] / total_coverage * 100  # Methylated
        
        # Sort and find CpGs
        df = df.sort_values('POS')
        df['next_pos'] = df['POS'].shift(-1)
        df['next_ref'] = df['REF'].shift(-1)
        
        # Find CpGs (C followed by G at next position)
        cpg_df = df[(df['REF'] == 'C') & 
                   (df['next_ref'] == 'G') &
                   (df['next_pos'] == df['POS'] + 1)]
        
        return cpg_df
    
    except Exception as e:
        print(f"Error processing file {input_csv}: {str(e)}")
        return None

def analyze_snp_zygosity(df, positions):
    """Analyze zygosity at specified SNP positions"""
    zygosity_results = []
    
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
    
    return zygosity_results

def create_methylation_plot(case_files, control_files, output_pdf, add_snp_lines=False, use_colors=False):
    """Create combined methylation plot for all files"""
    predefined_positions = [397, 492]
    
    # Process all files
    case_data = []
    control_data = []
    
    print("Processing case files...")
    for f in case_files:
        cpg_df = process_csv_file(f)
        if cpg_df is not None and len(cpg_df) > 1:
            case_data.append((f, cpg_df))
    
    print("Processing control files...")
    for f in control_files:
        cpg_df = process_csv_file(f)
        if cpg_df is not None and len(cpg_df) > 1:
            control_data.append((f, cpg_df))
    
    if not case_data and not control_data:
        print("No valid CpG data found in any files")
        return
    
    # Create PDF
    with PdfPages(output_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(12, 6))
        fig.subplots_adjust(bottom=0.4)
        
        sigma = 1.5  # Smoothing factor
        
        # Define color palette for --color mode
        color_palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        if use_colors:
            # Use different color for each file
            all_files = [(f, cpg_df, 'Case') for f, cpg_df in case_data] + \
                       [(f, cpg_df, 'Control') for f, cpg_df in control_data]
            
            for i, (filename, cpg_df, group) in enumerate(all_files):
                smoothed_C = gaussian_filter1d(cpg_df['percent_C'], sigma)
                color = color_palette[i % len(color_palette)]
                ax.plot(cpg_df['POS'], smoothed_C, 
                       color=color, linewidth=2,
                       label=f'{group}: {filename.split("/")[-1]}')
        else:
            # Use group colors (blue for case, red for control)
            # Plot case group (blue shades)
            for i, (filename, cpg_df) in enumerate(case_data):
                smoothed_C = gaussian_filter1d(cpg_df['percent_C'], sigma)
                alpha = 0.6 if len(case_data) > 1 else 1.0
                ax.plot(cpg_df['POS'], smoothed_C, 
                       color='blue', linewidth=2, alpha=alpha,
                       label=f'Case: {filename.split("/")[-1]}')
            
            # Plot control group (red shades)
            for i, (filename, cpg_df) in enumerate(control_data):
                smoothed_C = gaussian_filter1d(cpg_df['percent_C'], sigma)
                alpha = 0.6 if len(control_data) > 1 else 1.0
                ax.plot(cpg_df['POS'], smoothed_C, 
                       color='red', linewidth=2, alpha=alpha,
                       label=f'Control: {filename.split("/")[-1]}')
        
        # Add SNP analysis if requested (using first file from each group)
        if add_snp_lines:
            # Use first case file for SNP analysis
            if case_data:
                first_file = case_files[0]
                df = pd.read_csv(first_file)
                total_coverage = df['A'] + df['C'] + df['G'] + df['T']
                df['percent_C'] = df['C'] / total_coverage * 100
                
                zygosity_results = analyze_snp_zygosity(df, predefined_positions)
                
                for pos, genotype, p_val in zygosity_results:
                    if isinstance(p_val, float):
                        if genotype == "A-G":
                            color = "green"
                        elif genotype == "A-A":
                            color = "purple"
                        elif genotype == "G-G":
                            color = "orange"
                        else:
                            color = "gray"
                        
                        ax.axvline(x=pos, color=color, linestyle='--', 
                                 linewidth=1.5, alpha=0.7)
                        
                        label = f"{genotype}\nPosition: {pos}\np = {p_val:.3g}"
                        ax.annotate(label,
                            xy=(pos, 0), xycoords='data',
                            xytext=(0, -50), textcoords='offset points',
                            ha='center', va='top',
                            fontsize=9, color=color,
                            bbox=dict(facecolor='white', edgecolor=color, 
                                    boxstyle='round,pad=0.3'))
        
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Percent Methylation', fontsize=12)
        ax.set_title('CpG Methylation: Case vs Control', fontsize=14, fontweight='bold')
        ax.set_ylim(-10, 100)
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.legend(loc='best', fontsize=9)
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    
    print(f"Plot saved to {output_pdf}")

def main():
    args = parse_arguments()
    
    create_methylation_plot(
        case_files=args.case,
        control_files=args.control,
        output_pdf=args.output_pdf,
        add_snp_lines=args.snp,
        use_colors=args.color
    )

if __name__ == "__main__":
    main()