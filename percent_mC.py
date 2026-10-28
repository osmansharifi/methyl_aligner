#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import warnings
from scipy.ndimage import gaussian_filter1d
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import binomtest, ttest_ind, fisher_exact

# Journal figure style (Arial on most systems, editable text in PDFs)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': False, 'legend.frameon': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.dpi': 300,
})
# Default reference: NHIP_pro_forward_HY (801 bp); SNP positions are coordinates on it only
DEFAULT_REFERENCE_SEQ = (
    "GCAGACTCTGGAGAAGCCCAGGAGCGCGTGTCCCCGCCGGCCCCGCTCAGCCTGACCTGCTTTTCCCCAC"
    "AGAGGTGCCGGGGTTCGGTGCAGTCTGGCCGTGTCCTCCTCTTCTCGCGGTTTCAGATTCTTCTTGCGTG"
    "ACATTTTATCACACAAGACCGTGTTCGGGGATGGGGTGTTCTGAGCACGCCCTTACCGAGCCAGGCTCAC"
    "CCAGTTCCCCATCTAGAGCTCCCTGCTGGCATCTCTCACCTGGGACCCCAGCCTGCAAGACTCCGACACC"
    "CACGCCCTTAGGGTGATGTGGCCGCGCCCAGCCTGCGCTCAGCACCACGGGAGGGAGGGCGGCAGGCGGA"
    "AGCCCCGGCCATCTTGGGAAGGGCCTTTTGGGACCGTGTGTTCCGCTGAGAACCACAGACCCATGCTTCA"
    "GACGCTACAGATTTCTCTCTCAAATGCCTGGATTCAAAATATAAAAGTGTGTTTGAGCAAACAGCAAAAT"
    "ATTTGGGAAATTCTCCCAGCGAGGCTTCTAATGACCCCGTAGGAAAAGCCAGGCTGCTCCCTGCTGCGCC"
    "TCCAGCGAACTTCCCACCGGGCAGGACTCTGGGAAGAAGAGCAAAGGGTTCTATTTTGCGTAGAATTTAA"
    "ATAAAATTGCTTCACAGACGGGAAACGTCAAGCAAGCAGCACTCTGAGCCCATGCGTGCTTTTAATGCCT"
    "GAATGACATCGCAATTCAGTCTGAGTTCCCGCCTGGGCTCAGGGCCTTCTGAACTGCCGGGGGCTGTGAA"
    "ATGGTTTATAATCAGAAGAGAAGCAGTGAGA"
)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate CpG methylation plots from multiple CSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 percent_mC.py output.pdf --case file1.csv file2.csv --control file3.csv file4.csv
  python3 percent_mC.py output.pdf --case file1.csv --control file2.csv -snp
  python3 percent_mC.py output.pdf --sample file1.csv -snp
        """
    )
    
    parser.add_argument('output_pdf', help='Output PDF file for the plot')
    parser.add_argument('--case', nargs='+',
                       help='Input CSV files for case group (differential methylation needs --case and --control)')
    parser.add_argument('--control', nargs='+',
                       help='Input CSV files for control group (differential methylation needs --case and --control)')
    parser.add_argument('--sample',
                       help='A single CSV file for SNP calling with -snp, no case or control needed')
    parser.add_argument('-snp', action='store_true', 
                       help='Add vertical lines at rs4823981 (G>T, 238) and rs4823982 (G>A, 333) with zygosity analysis '
                            '(default NHIP_pro_forward_HY reference only)')
    parser.add_argument('--color', action='store_true',
                       help='Use different colors for each file instead of grouping by case/control')
    
    args = parser.parse_args()
    
    if args.sample:
        if args.case or args.control:
            parser.error('use --sample on its own, or --case and --control for differential methylation')
        if not args.snp:
            parser.error('--sample is for SNP calling, so add -snp')
    elif not (args.case and args.control):
        parser.error('differential methylation needs both --case and --control (or use --sample FILE -snp for one sample)')
    
    return args

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
    """Analyze zygosity at specified SNP positions, given as (rsid, pos, ref, alt)"""
    zygosity_results = []
    
    for rsid, pos, ref, alt in positions:
        row = df[df['POS'] == pos]
        if not row.empty:
            r, a = row[[ref, alt]].values[0]
            total = r + a
            if total >= 10:  # minimum ref + alt reads to make a call
                result = binomtest(a, n=total, p=0.5)
                p_val = result.pvalue
                # Not significantly different from 50:50 means heterozygous
                if p_val >= 0.05:
                    genotype = f"{ref}/{alt}"
                else:
                    genotype = f"{alt}/{alt}" if a > r else f"{ref}/{ref}"
                zygosity_results.append((rsid, pos, ref, alt, r, a, genotype, p_val))
            else:
                zygosity_results.append((rsid, pos, ref, alt, r, a, "Low coverage", None))
        else:
            zygosity_results.append((rsid, pos, ref, alt, 0, 0, "Not found", None))
    
    return zygosity_results

def differential_methylation(case_data, control_data, min_cov=10):
    """Test case vs control methylation at each CpG and across the whole region.

    With 2 or more samples per group, a Welch t-test on per-sample percent methylation
    (samples are the unit, so PCR duplicates do not inflate p values). With 1 sample in
    either group, Fisher's exact test on pooled C and T read counts.
    Per-CpG p values are Benjamini-Hochberg corrected.
    """
    def counts(data):
        # one row per CpG, C and T counts for each sample; CpGs below min_cov set to missing
        C = pd.concat([d.set_index('POS')['C'] for _, d in data], axis=1, keys=range(len(data)))
        T = pd.concat([d.set_index('POS')['T'] for _, d in data], axis=1, keys=range(len(data)))
        low = (C + T) < min_cov
        return C.mask(low), T.mask(low)

    case_C, case_T = counts(case_data)
    ctrl_C, ctrl_T = counts(control_data)
    case_pct = case_C / (case_C + case_T) * 100
    ctrl_pct = ctrl_C / (ctrl_C + ctrl_T) * 100
    replicates = len(case_data) >= 2 and len(control_data) >= 2
    test = 'Welch t-test' if replicates else "Fisher's exact (reads, no replicates)"

    rows = []
    for pos in case_pct.index.intersection(ctrl_pct.index):
        a, b = case_pct.loc[pos].dropna(), ctrl_pct.loc[pos].dropna()
        if replicates:
            if len(a) < 2 or len(b) < 2:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)  # identical values at a CpG
                p_val = ttest_ind(a, b, equal_var=False).pvalue
        else:
            if len(a) == 0 or len(b) == 0:
                continue
            table = [[case_C.loc[pos].sum(), case_T.loc[pos].sum()],
                     [ctrl_C.loc[pos].sum(), ctrl_T.loc[pos].sum()]]
            p_val = fisher_exact(table).pvalue
        p_val = 1.0 if np.isnan(p_val) else p_val  # no variance in either group
        rows.append((pos, a.mean(), b.mean(), a.mean() - b.mean(), p_val))

    stats_df = pd.DataFrame(rows, columns=['POS', 'case_pct', 'control_pct', 'delta', 'p_value'])
    if stats_df.empty:
        return stats_df, None

    # Benjamini-Hochberg FDR
    p = stats_df['p_value'].values
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    fdr = np.empty(len(p))
    fdr[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    stats_df['fdr'] = fdr

    # Whole region: mean methylation across tested CpGs
    tested = stats_df['POS']
    case_region = case_pct.loc[tested].mean()
    ctrl_region = ctrl_pct.loc[tested].mean()
    if replicates:
        region_p = ttest_ind(case_region, ctrl_region, equal_var=False).pvalue
    else:
        region_p = fisher_exact([[case_C.loc[tested].sum().sum(), case_T.loc[tested].sum().sum()],
                                 [ctrl_C.loc[tested].sum().sum(), ctrl_T.loc[tested].sum().sum()]]).pvalue
    region_p = 1.0 if np.isnan(region_p) else region_p
    region = {'test': test, 'n_case': len(case_data), 'n_control': len(control_data),
              'delta': case_region.mean() - ctrl_region.mean(), 'p_value': region_p,
              'n_sig': int((stats_df['fdr'] < 0.05).sum()), 'n_tested': len(stats_df)}
    return stats_df, region

def create_methylation_plot(case_files, control_files, output_pdf, add_snp_lines=False, use_colors=False, sample_mode=False):
    """Create combined methylation plot for all files"""
    predefined_positions = [('rs4823981', 238, 'G', 'T'), ('rs4823982', 333, 'G', 'A')]
    
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
        fig, ax = plt.subplots(figsize=(7.2, 3.2))  # double-column journal width
        
        sigma = 1.5  # Smoothing factor
        notes = []   # figure notes printed under the plot
        
        # Okabe-Ito colorblind-safe palette, used in this fixed order for --color mode
        color_palette = ['#0072B2', '#D55E00', '#009E73', '#CC79A7',
                        '#E69F00', '#56B4E9', '#000000']
        group_colors = {'Case': '#D55E00', 'Control': '#0072B2', 'Sample': '#000000'}
        case_label = 'Sample' if sample_mode else 'Case'  # --sample runs use the case slot
        
        def sample_name(f):
            return f.split("/")[-1].rsplit(".", 1)[0]
        
        if use_colors:
            # Use different color for each file
            all_files = [(f, cpg_df, case_label) for f, cpg_df in case_data] + \
                       [(f, cpg_df, 'Control') for f, cpg_df in control_data]
            
            for i, (filename, cpg_df, group) in enumerate(all_files):
                smoothed_C = gaussian_filter1d(cpg_df['percent_C'], sigma)
                ax.plot(cpg_df['POS'], smoothed_C,
                       color=color_palette[i % len(color_palette)], linewidth=1.2,
                       label=f'{group}: {sample_name(filename)}')
        else:
            # Group colors: one line per sample, one legend entry per group
            for group, data in ((case_label, case_data), ('Control', control_data)):
                for i, (filename, cpg_df) in enumerate(data):
                    if len(data) == 1:
                        label = f'{group}: {sample_name(filename)}'
                    else:
                        label = f'{group} (n = {len(data)})' if i == 0 else None
                    ax.plot(cpg_df['POS'], gaussian_filter1d(cpg_df['percent_C'], sigma),
                           color=group_colors[group], linewidth=1.0, alpha=0.8, label=label)
        
        # CpG sites as tick marks along the x axis
        cpg_positions = pd.concat([d['POS'] for _, d in case_data + control_data]).unique()
        ax.plot(cpg_positions, np.full(len(cpg_positions), 1.5), linestyle='none', marker='|',
               markersize=4, markeredgewidth=0.8, color='0.55', label='CpG site')
        
        # Add SNP analysis if requested (using first file from each group)
        if add_snp_lines:
            # Use first case file for SNP analysis
            if case_data:
                first_file = case_files[0]
                df = pd.read_csv(first_file)
                total_coverage = df['A'] + df['C'] + df['G'] + df['T']
                df['percent_C'] = df['C'] / total_coverage * 100
                
                # SNP positions only apply to the default reference
                default_ref_used = ''.join(df.sort_values('POS')['REF'].astype(str)).upper() == DEFAULT_REFERENCE_SEQ
                if not default_ref_used:
                    print(f"WARNING: {first_file} was not aligned to the default NHIP_pro_forward_HY reference.\n"
                          f"SNP calling does not apply to a custom reference, so -snp is ignored.")
                
                zygosity_results = analyze_snp_zygosity(df, predefined_positions) if default_ref_used else []
                
                snp_notes = []
                for rsid, pos, ref, alt, r, a, genotype, p_val in zygosity_results:
                    print(f"{rsid} ({ref}>{alt}) position {pos}: {genotype}  {ref}={r} {alt}={a}")
                    if isinstance(p_val, float):
                        ax.axvline(x=pos, color='0.45', linestyle='--', linewidth=0.8, zorder=0)
                        ax.text(pos, 1.01, f"{rsid}\n{genotype}", transform=ax.get_xaxis_transform(),
                               ha='center', va='bottom', fontsize=7)
                        snp_notes.append(f"{rsid} {genotype} (P = {p_val:.2g})")
                    else:
                        snp_notes.append(f"{rsid} {genotype.lower()}")
                if snp_notes:
                    notes.append(f"Genotypes ({sample_name(first_file)}): " + ", ".join(snp_notes)
                                 + "; P < 0.05 means homozygous")
        
        # Differential methylation: case vs control
        if case_data and control_data:
            stats_df, region = differential_methylation(case_data, control_data)
            if region:
                stats_csv = output_pdf.rsplit('.', 1)[0] + '_differential_methylation.csv'
                stats_df.to_csv(stats_csv, index=False)
                print(f"Differential methylation ({region['test']}): mean difference "
                      f"{region['delta']:+.1f}%, p = {region['p_value']:.3g}, "
                      f"{region['n_sig']} of {region['n_tested']} CpGs FDR < 0.05")
                print(f"Per-CpG statistics saved to {stats_csv}")
                
                # Mark significant CpGs
                sig = stats_df[stats_df['fdr'] < 0.05]
                for pos in sig['POS']:
                    ax.text(pos, 94, '*', ha='center', va='center', fontsize=12, fontweight='bold')
                
                test = 't-test' if region['test'] == 'Welch t-test' else "Fisher's exact on reads"
                notes.insert(0, f"* Significant CpG (FDR < 0.05): {region['n_sig']} of {region['n_tested']}. "
                                f"Overall, case is {region['delta']:+.1f}% vs control "
                                f"(P = {region['p_value']:.2g}, {test}, n = {region['n_case']} vs {region['n_control']})")
        
        ax.set_xlabel('Position in amplicon (bp)')
        ax.set_ylabel('CpG methylation (%)')
        ax.set_ylim(0, 100)
        ax.set_xlim(left=0)
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        if notes:
            ax.text(0, -0.24, '\n'.join(notes), transform=ax.transAxes, ha='left', va='top',
                   fontsize=6.5, color='0.25', linespacing=1.5)
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    
    print(f"Plot saved to {output_pdf}")

def main():
    args = parse_arguments()
    
    create_methylation_plot(
        case_files=[args.sample] if args.sample else args.case,
        control_files=[] if args.sample else args.control,
        output_pdf=args.output_pdf,
        add_snp_lines=args.snp,
        use_colors=args.color,
        sample_mode=bool(args.sample)
    )

if __name__ == "__main__":
    main()