#!/usr/bin/env python3
"""methyl_aligner: bisulfite amplicon alignment, CpG methylation plots,
SNP genotyping and case vs control differential methylation.

Commands:
  run     start here: align FASTQ files and make the figures in one step
  align   align FASTQ files only (step 1 of 2)
  plot    make the figures from aligned samples (step 2 of 2)

Run "python3 methyl_aligner.py <command> -h" for the options of each command.
"""

import argparse
import datetime
import gzip
import os
import sys
import subprocess
import tempfile
import shlex
import shutil
import threading
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from queue import Queue
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage import gaussian_filter1d
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

# Default reference: NHIP_pro_forward_HY (801 bp), built in so no FASTA is needed
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

# SNPs genotyped for each --sample: (rsID, position on NHIP_pro_forward_HY, REF, ALT)
SNP_POSITIONS = [('rs4823981', 238, 'G', 'T'), ('rs4823982', 333, 'G', 'A')]

#######################
## Argument Parsing  ##
#######################

EXAMPLES = """
Examples:
  # One or more samples: methylation figure and SNP genotypes for each
  python3 methyl_aligner.py run --sample sample1.fastq sample2.fastq

  # Case vs control: figure and differential methylation statistics
  python3 methyl_aligner.py run --case case1.fastq case2.fastq --control ctrl1.fastq ctrl2.fastq

  # Two steps instead: align once, then re-plot any time without realigning
  python3 methyl_aligner.py align case1.fastq case2.fastq ctrl1.fastq ctrl2.fastq
  python3 methyl_aligner.py plot --case results/samples/case1/case1_matrix.csv results/samples/case2/case2_matrix.csv \\
                                 --control results/samples/ctrl1/ctrl1_matrix.csv results/samples/ctrl2/ctrl2_matrix.csv

Outputs:
  results/samples/<sample>/   count matrices and nucleotide count plots
  results/analyses/<name>/    figure, SNP genotypes or statistics, summary.txt
"""

class FriendlyParser(argparse.ArgumentParser):
    """Argument parser whose errors say what went wrong and where to get help."""
    def error(self, message):
        self.print_usage(sys.stderr)
        sys.stderr.write(f"\nError: {message}\nFor help: python3 {self.prog} -h\n")
        sys.exit(2)

def parse_arguments(argv=None):
    """Parse and return all command line arguments."""
    argv = sys.argv[1:] if argv is None else list(argv)
    parser = FriendlyParser(
        prog='methyl_aligner.py',
        description='Bisulfite amplicon long read alignment, CpG methylation figures, '
                    'SNP genotyping and case vs control differential methylation.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EXAMPLES)
    commands = parser.add_subparsers(dest='command', metavar='command', title='commands')
    
    def add_alignment_options(p):
        group = p.add_argument_group('alignment options')
        group.add_argument('-r', '--reference', dest='fasta', metavar='FASTA',
                           help='custom reference (default: built-in NHIP_pro_forward_HY). '
                                'At your own risk: SNP genotyping is skipped for a custom reference.')
        group.add_argument('-p', '--percent', type=int, default=70, metavar='N',
                           help='minimum percent identity to keep a read (default: 70)')
        group.add_argument('-t', '--threads', type=int, default=4, metavar='N',
                           help='number of threads (default: 4)')
    
    def add_sample_options(p, kind):
        group = p.add_argument_group('samples (use --sample, or --case and --control)')
        group.add_argument('--sample', nargs='+', metavar=kind,
                           help='one or more samples: a methylation figure and SNP genotypes '
                                '(rs4823981, rs4823982) for each')
        group.add_argument('--case', nargs='+', metavar=kind,
                           help='case group, compared with --control for differential methylation')
        group.add_argument('--control', nargs='+', metavar=kind,
                           help='control group, compared with --case')
    
    def add_output_options(p, figures=True):
        group = p.add_argument_group('output options')
        if figures:
            group.add_argument('-o', '--name',
                               help='name of the analysis folder (default: the sample name for --sample, '
                                    'case_vs_control for --case/--control)')
            group.add_argument('--color', action='store_true',
                               help='give each case and control sample its own color')
        group.add_argument('--outdir', default='results', metavar='DIR',
                           help='folder for all results (default: results)')
    
    subparsers = {}
    p = commands.add_parser('run', help='start here: align FASTQ files and make the figures in one step',
                            description='Align FASTQ files, then make the figures. Use --sample for SNP '
                                        'genotypes, or --case and --control for differential methylation.',
                            formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EXAMPLES)
    add_sample_options(p, 'FASTQ')
    add_alignment_options(p)
    add_output_options(p)
    subparsers['run'] = p
    
    p = commands.add_parser('align', help='align FASTQ files only (step 1 of 2)',
                            description='Align FASTQ files (.fastq, .fq, or gzipped) to the reference. '
                                        'Writes results/samples/<sample>/ for each file.',
                            formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EXAMPLES)
    p.add_argument('fastq', nargs='+', metavar='FASTQ', help='one or more FASTQ files')
    add_alignment_options(p)
    add_output_options(p, figures=False)
    subparsers['align'] = p
    
    p = commands.add_parser('plot', help='make the figures from aligned samples (step 2 of 2)',
                            description='Make figures from the matrix CSVs that "align" wrote to '
                                        'results/samples/. Use --sample for SNP genotypes, or --case '
                                        'and --control for differential methylation.',
                            formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EXAMPLES)
    add_sample_options(p, 'CSV')
    add_output_options(p)
    subparsers['plot'] = p
    
    if not argv:
        parser.print_help()
        sys.exit(0)
    if '-snp' in argv or '--snp' in argv:
        argv = [a for a in argv if a not in ('-snp', '--snp')]
        sys.stderr.write("Note: -snp is no longer needed; --sample always genotypes the SNPs.\n")
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    p = subparsers[args.command]
    
    # Check the sample groups
    if args.command == 'align':
        inputs = args.fastq
    else:
        if args.sample and (args.case or args.control):
            p.error('use --sample, or --case with --control, but not both')
        if not args.sample and not (args.case and args.control):
            if args.case or args.control:
                p.error('differential methylation needs both --case and --control')
            p.error('tell me which samples to use: --sample FILE [FILE ...], or --case ... --control ...')
        if args.sample and args.name and len(args.sample) > 1:
            p.error('-o names a single analysis; with several --sample files each one is named after its sample')
        inputs = (args.sample or []) + (args.case or []) + (args.control or [])
    
    # Check the files: they must exist and be the right kind for the command
    missing = [f for f in inputs if not os.path.isfile(f)]
    if missing:
        p.error('file not found: ' + ', '.join(missing))
    if args.command == 'plot':
        fastqs = [f for f in inputs if not f.lower().endswith('.csv')]
        if fastqs:
            p.error('plot takes the matrix CSVs made by align (results/samples/<sample>/<sample>_matrix.csv). '
                    'To start from FASTQ files, use run instead.')
    else:
        csvs = [f for f in inputs if f.lower().endswith('.csv')]
        if csvs:
            p.error(f'{args.command} takes FASTQ files. To use matrix CSVs from an earlier run, use plot instead.')
    if getattr(args, 'fasta', None) and not os.path.isfile(args.fasta):
        p.error(f'reference not found: {args.fasta}')
    
    return args

def get_sample_name(fastq_path):
    """Extract sample name from FASTQ filename."""
    base = os.path.basename(fastq_path)
    # Remove common extensions
    for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    return base

#########################
## Alignment Functions ##
#########################

def anti(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[nt] for nt in seq[::-1]])

def run_needle(ref_file, seq, strand, wdir, tid, mat_file):
    """Run needle alignment and return alignment results."""
    fasta = f"{wdir}/{tid}.{strand}.fa"
    out = f"{wdir}/{tid}.{strand}.fa.out"
    
    with open(fasta, 'w') as f:
        f.write(f">s\n{seq}\n")
    
    params = [
        "needle",
        "-asequence", ref_file,
        "-bsequence", fasta,
        "-gapopen", "3",
        "-gapextend", "1",
        "-datafile", mat_file,
        "-endweight",
        "-outfile", out
    ]
    
    try:
        subprocess.run(params, check=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        sys.stderr.write("Error running needle\n")
        sys.exit(1)
    
    ra, sa, score = "", "", 0
    with open(out, 'r') as fh:
        for line in fh:
            if line.startswith('# Score:'):
                score = float(line.split()[2])
            elif line.startswith('r') or line.startswith('s'):
                parts = line.split()
                if parts[0] == 'r':
                    ra += parts[2]
                elif parts[0] == 's':
                    sa += parts[2]
    
    # Percent identity from A and G only: bisulfite conversion changes C to T
    match, total = 0, 0
    for r, s in zip(ra, sa):
        if r == '-' or s == '-':
            continue
        if r in ['A', 'G']:
            if r == s:
                match += 1
            total += 1
    
    return {'ref': ra, 'seq': sa, 'pct': match/total if total > 0 else 0}

def align_worker(q, results, ref_file, wdir, mat_file):
    """Worker function for alignment threads."""
    tid = threading.get_ident()
    alignments = []
    while True:
        try:
            plus = q.get_nowait()
        except:
            break
        
        anti_seq = anti(plus)
        a1 = run_needle(ref_file, plus, f"{tid}+", wdir, tid, mat_file)
        a2 = run_needle(ref_file, anti_seq, f"{tid}-", wdir, tid, mat_file)
        max_align = a1 if a1['pct'] > a2['pct'] else a2
        alignments.append(max_align)
        q.task_done()
    
    for item in alignments:
        results.put(item)

#############################
## Visualization Functions ##
#############################

def filter_cpg_sites(df):
    """Filter dataframe to only CpG sites (where REF = C and next REF = G)."""
    return df[(df['REF'] == 'C') & (df['REF'].shift(-1) == 'G')]

def calculate_percentages(df):
    """Calculate C and T percentages from nucleotide counts."""
    total_counts = df[['A', 'C', 'G', 'T']].sum().sum()
    percent_C = (df['C'].sum() / total_counts) * 100
    percent_T = (df['T'].sum() / total_counts) * 100
    return percent_C, percent_T

def create_plot(df, title, output_path, percentages):
    """Create and save a line plot of nucleotide counts."""
    fig, ax = plt.subplots(figsize=(7.2, 3.0))  # double-column journal width
    colors = {'A': '#009E73', 'C': '#0072B2', 'G': '#000000', 'T': '#D55E00'}  # Okabe-Ito
    for nucleotide in ['A', 'C', 'G', 'T']:
        ax.plot(df['POS'], df[nucleotide], label=nucleotide, color=colors[nucleotide], linewidth=0.7)

    ax.set_xlabel('Position in amplicon (bp)')
    ax.set_ylabel('Read count')
    ax.set_title(title, loc='left')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    percent_C, percent_T = percentages
    ax.legend(title='Base', bbox_to_anchor=(1.01, 1), loc='upper left', alignment='left')

    ax.text(1.02, 0.42, f'C: {percent_C:.1f}%\nT: {percent_T:.1f}%',
            transform=ax.transAxes, va='top', ha='left')

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved as {output_path}")

#################################
## Methylation Plot Functions ##
#################################

def process_csv_file(input_csv):
    """Process a single CSV file and return CpG data"""
    try:
        # Read the CSV file
        df = pd.read_csv(input_csv)
        
        # Calculate percent C (methylated, protected from conversion) and percent T (unmethylated)
        total_coverage = df['A'] + df['C'] + df['G'] + df['T']
        df['percent_C'] = df['C'] / total_coverage * 100  # Methylated
        df['percent_T'] = df['T'] / total_coverage * 100  # Unmethylated
        
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
              'delta': round(case_region.mean() - ctrl_region.mean(), 1) + 0.0, 'p_value': region_p,
              'n_sig': int((stats_df['fdr'] < 0.05).sum()), 'n_tested': len(stats_df)}
    return stats_df, region

def display_name(f):
    """Sample name for figures: file name without folder, extension or the _matrix suffix."""
    name = f.split("/")[-1].rsplit(".", 1)[0]
    return name[:-len('_matrix')] if name.endswith('_matrix') else name

def create_methylation_plot(case_files, control_files, output_pdf, add_snp_lines=False, use_colors=False, sample_mode=False):
    """Create combined methylation plot for all files; returns a report of the results"""
    predefined_positions = SNP_POSITIONS
    report = {'snps': [], 'snp_sample': None, 'region': None, 'sig_positions': [], 'files': [output_pdf]}
    case_label = 'Sample' if sample_mode else 'Case'  # --sample runs use the case slot
    
    # Process all files
    case_data = []
    control_data = []
    
    print(f"Reading {case_label.lower()} files...")
    for f in case_files:
        cpg_df = process_csv_file(f)
        if cpg_df is not None and len(cpg_df) > 1:
            case_data.append((f, cpg_df))
    
    if control_files:
        print("Reading control files...")
    for f in control_files:
        cpg_df = process_csv_file(f)
        if cpg_df is not None and len(cpg_df) > 1:
            control_data.append((f, cpg_df))
    
    if not case_data and not control_data:
        print("No valid CpG data found in any files")
        return None
    
    # Create PDF
    with PdfPages(output_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(7.2, 3.2))  # double-column journal width
        
        sigma = 1.5  # Smoothing factor
        notes = []   # figure notes printed under the plot
        
        # Okabe-Ito colorblind-safe palette, used in this fixed order for --color mode
        color_palette = ['#0072B2', '#D55E00', '#009E73', '#CC79A7',
                        '#E69F00', '#56B4E9', '#000000']
        group_colors = {'Case': '#D55E00', 'Control': '#0072B2', 'Sample': '#000000'}
        
        sample_name = display_name
        
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
                          f"SNP genotyping is skipped for a custom reference.")
                
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
                if zygosity_results:
                    snp_csv = output_pdf.rsplit('.', 1)[0] + '_snp_genotypes.csv'
                    pd.DataFrame([(sample_name(first_file),) + tuple(r) for r in zygosity_results],
                                 columns=['sample', 'rsid', 'position', 'ref', 'alt', 'ref_reads',
                                          'alt_reads', 'genotype', 'p_value']).to_csv(snp_csv, index=False)
                    print(f"SNP genotypes saved to {snp_csv}")
                    report.update(snps=zygosity_results, snp_sample=sample_name(first_file))
                    report['files'].append(snp_csv)
        
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
                
                report.update(region=region)
                report['files'].append(stats_csv)
                
                # Mark significant CpGs
                sig = stats_df[stats_df['fdr'] < 0.05]
                report['sig_positions'] = [int(p) for p in sig['POS']]
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
    return report

##########################
## Pipeline Functions   ##
##########################

def run_alignment(args):
    """Run the alignment portion of the pipeline."""
    # Temporary working folder, removed when alignment finishes
    wdir = tempfile.mkdtemp(prefix='methyl_aligner_')
    
    # Process reference FASTA (built-in NHIP_pro_forward_HY unless -r is given)
    if args.fasta:
        with open(args.fasta, 'r') as fh:
            next(fh)  # Skip header
            ref_seq = ''.join(line.strip() for line in fh)
    else:
        ref_seq = DEFAULT_REFERENCE_SEQ
    with open(f"{wdir}/ref.fa", 'w') as of:
        of.write(">r\n")
        of.write(ref_seq + "\n")

    # Warn if the reference is not the default: SNP calling does not apply
    if ref_seq.upper() != DEFAULT_REFERENCE_SEQ:
        sys.stderr.write(f"\nWARNING: using a custom reference ({args.fasta}) instead of the default\n"
                         f"NHIP_pro_forward_HY reference. You are proceeding at your own risk.\n"
                         f"SNP genotyping is skipped for this reference.\n\n")

    # Process FASTQ and check for duplicates
    read_counts = defaultdict(int)
    opener = gzip.open if args.fastq.endswith('.gz') else open
    with opener(args.fastq, 'rt') as fq:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline().strip().upper()
            fq.readline()  # +
            fq.readline()  # quality
            read_counts[seq] += 1
    
    duplicates = sum(1 for cnt in read_counts.values() if cnt > 1)
    if duplicates:
        sys.stderr.write(f"{duplicates} duplicates found and skipped\n")
    
    # Create bisulfite scoring matrix: reference C to read T scores as a match
    methyl_matrix = """    A   C   G   T   N
A   1  -1  -1  -1   0
C  -1   1  -1   1   0
G  -1  -1   1  -1   0
T  -1  -1  -1   1   0
N   0   0   0   0   0
"""
    
    with open(f"{wdir}/mat", 'w') as mat:
        mat.write(methyl_matrix)
    
    # Set up threading
    query_queue = Queue()
    result_queue = Queue()
    
    for seq in read_counts.keys():
        query_queue.put(seq)
    
    workers = []
    for _ in range(args.threads):
        t = threading.Thread(target=align_worker, args=(query_queue, result_queue, 
                                               f"{wdir}/ref.fa", wdir, 
                                               f"{wdir}/mat"))
        t.start()
        workers.append(t)
    
    for t in workers:
        t.join()
    
    # Process results
    all_alignments = []
    skipped = 0
    while not result_queue.empty():
        alignment = result_queue.get()
        if 100 * alignment['pct'] < args.percent:
            skipped += 1
            continue
        all_alignments.append(alignment)
    
    sys.stderr.write(f"{skipped} alignments below {args.percent}% skipped ({len(all_alignments)} total)\n")
    
    if not all_alignments:
        sys.stderr.write("no alignments to process\n")
        shutil.rmtree(wdir)
        return None, ref_seq
    
    # Process aligned sequences
    aseqs = []
    for alignment in all_alignments:
        aseq = []
        for r, s in zip(alignment['ref'], alignment['seq']):
            if r != '-':
                aseq.append(s)
        aseqs.append(''.join(aseq))
    
    # Create position information dataframe
    nts = ['A', 'C', 'G', 'T']
    data = []
    for i in range(len(ref_seq)):
        counts = {'POS': i+1, 'REF': ref_seq[i]}
        for nt in nts:
            counts[nt] = 0
        
        for s in aseqs:
            snt = s[i] if i < len(s) else ' '
            if snt in counts:
                counts[snt] += 1
        data.append(counts)
    
    df = pd.DataFrame(data)
    
    # Clean up the temporary folder
    shutil.rmtree(wdir)
    
    return df, ref_seq

def run_visualization(df, sample_name, outdir="results"):
    """Run the visualization portion of the pipeline."""
    # Create output directory based on sample name
    output_dir = os.path.join(outdir, 'samples', sample_name)
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter CpG sites
    df_ref_c = filter_cpg_sites(df)
    
    # Save data
    df.to_csv(f'{output_dir}/{sample_name}_matrix.csv', index=False)
    df_ref_c.to_csv(f'{output_dir}/{sample_name}_matrix_ref_c.csv', index=False)
    print(f"Data saved to {output_dir}/{sample_name}_matrix.csv and {output_dir}/{sample_name}_matrix_ref_c.csv")
    
    # Create plots
    create_plot(df, f'Nucleotide Counts by Position - {sample_name} (All Data)', 
                f'{output_dir}/{sample_name}_line_plot_all_data.pdf', calculate_percentages(df))
    
    create_plot(df_ref_c, f'Nucleotide Counts by Position - {sample_name} (CpG Only)', 
                f'{output_dir}/{sample_name}_line_plot_cpg_only.pdf', calculate_percentages(df_ref_c))
    
    return output_dir

##########################
## Commands             ##
##########################

def check_needle():
    """Stop with a helpful message if EMBOSS needle is not installed."""
    if shutil.which('needle') is None:
        sys.stderr.write("You must install 'emboss' and have 'needle' in your path.\n"
                         "If it is in a conda environment, activate it first, for example: conda activate emboss\n")
        sys.exit(1)

def align_samples(args, fastq_files):
    """Align each FASTQ and return {fastq: matrix CSV path} for samples that aligned."""
    check_needle()
    matrices = {}
    for n, fastq in enumerate(fastq_files, 1):
        sample_name = get_sample_name(fastq)
        print(f"\n[{n}/{len(fastq_files)}] Aligning {sample_name}")
        sample_args = argparse.Namespace(**vars(args))
        sample_args.fastq = fastq
        df, _ = run_alignment(sample_args)
        if df is None:
            print(f"Skipping {sample_name}: no reads passed the {args.percent}% identity filter")
            continue
        output_dir = run_visualization(df, sample_name, args.outdir)
        matrices[fastq] = os.path.join(output_dir, f'{sample_name}_matrix.csv')
    return matrices

def command_align(args):
    matrices = align_samples(args, args.fastq)
    if matrices:
        print(f"\nDone. Results for {len(matrices)} sample(s) are in {os.path.join(args.outdir, 'samples')}/")
        print("Next, make the figures with plot, for example:")
        print("  python3 methyl_aligner.py plot --sample " + " ".join(matrices.values()))

def analysis_pdf_path(name, outdir):
    """Figures and tables for an analysis go to <outdir>/analyses/<name>/<name>.pdf and friends."""
    pdf = os.path.join(outdir, 'analyses', name, f'{name}.pdf')
    os.makedirs(os.path.dirname(pdf), exist_ok=True)
    return pdf

def uses_default_reference(matrix_csv):
    refs = pd.read_csv(matrix_csv, usecols=['POS', 'REF']).sort_values('POS')['REF'].astype(str)
    return ''.join(refs).upper() == DEFAULT_REFERENCE_SEQ

def write_summary(path, args, report, case_files, control_files, sample_mode):
    """Write a plain text record of the command, inputs and results."""
    lines = [f"methyl_aligner analysis: {display_name(report['files'][0])}",
             f"Date: {datetime.datetime.now():%Y-%m-%d %H:%M}",
             "Command: " + shlex.join([os.path.basename(sys.argv[0])] + sys.argv[1:]), ""]
    if getattr(args, 'percent', None) is not None:
        reference = args.fasta if args.fasta else 'built-in NHIP_pro_forward_HY'
        lines += ["Alignment", f"  reference: {reference}",
                  f"  percent identity minimum: {args.percent}%", ""]
    groups = [('Sample', case_files)] if sample_mode else [('Case', case_files), ('Control', control_files)]
    lines.append("Inputs")
    for group, files in groups:
        lines.append(f"  {group} (n = {len(files)})")
        for f in files:
            ref = 'default reference' if uses_default_reference(f) else 'CUSTOM reference'
            lines.append(f"    {display_name(f)}: {f} ({ref})")
    lines.append("")
    if report['snps']:
        lines.append(f"SNP genotypes ({report['snp_sample']}), binomial test vs 50:50, P < 0.05 means homozygous")
        for rsid, pos, ref, alt, r, a, genotype, p_val in report['snps']:
            p_txt = f"P = {p_val:.2g}" if p_val is not None else "not called"
            lines.append(f"  {rsid} {ref}>{alt} (position {pos}): {genotype}, {ref} = {r}, {alt} = {a}, {p_txt}")
        lines.append("")
    region = report['region']
    if region:
        lines += [f"Differential methylation ({region['test']}, n = {region['n_case']} vs {region['n_control']})",
                  f"  overall: case is {region['delta']:+.1f}% vs control, P = {region['p_value']:.2g}",
                  f"  significant CpGs (FDR < 0.05): {region['n_sig']} of {region['n_tested']}"
                  + (f", at positions {', '.join(map(str, report['sig_positions']))}" if report['sig_positions'] else ""),
                  ""]
    lines.append("Files")
    lines += [f"  {f}" for f in report['files'] + [path]]
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

def make_analysis(args, name, case_files, control_files, sample_mode):
    """Make one figure with its tables and summary in results/analyses/<name>/."""
    print(f"\nAnalysis: {name}")
    output_pdf = analysis_pdf_path(name, args.outdir)
    report = create_methylation_plot(
        case_files=case_files,
        control_files=control_files,
        output_pdf=output_pdf,
        add_snp_lines=sample_mode,
        use_colors=args.color,
        sample_mode=sample_mode
    )
    if report:
        summary = output_pdf.rsplit('.', 1)[0] + '_summary.txt'
        write_summary(summary, args, report, case_files, control_files, sample_mode)
        print(f"Summary saved to {summary}")
        return os.path.dirname(output_pdf)
    return None

def command_plot(args):
    if args.sample:
        folders = [make_analysis(args, args.name or display_name(f), [f], [], sample_mode=True)
                   for f in args.sample]
    else:
        folders = [make_analysis(args, args.name or 'case_vs_control', args.case, args.control,
                                 sample_mode=False)]
    folders = [f for f in folders if f]
    if folders:
        print("\nDone. Results are in:")
        for folder in folders:
            print(f"  {folder}/")

def command_run(args):
    fastqs = args.sample or (args.case + args.control)
    matrices = align_samples(args, fastqs)
    missing = [f for f in fastqs if f not in matrices]
    if missing:
        sys.stderr.write("Cannot make figures: these samples had no aligned reads: " + ", ".join(missing) + "\n")
        sys.exit(1)
    if args.sample:
        args.sample = [matrices[f] for f in args.sample]
    else:
        args.case = [matrices[f] for f in args.case]
        args.control = [matrices[f] for f in args.control]
    command_plot(args)

###################
## Main Function ##
###################

def main(argv=None):
    args = parse_arguments(argv)
    {'align': command_align, 'plot': command_plot, 'run': command_run}[args.command](args)

if __name__ == '__main__':
    main()
