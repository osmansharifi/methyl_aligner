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
import shlex
import shutil
import subprocess
import sys
import tempfile
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from itertools import zip_longest

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

# Bisulfite scoring matrix for needle: reference C read as T scores as a match
BISULFITE_MATRIX = """    A   C   G   T   N
A   1  -1  -1  -1   0
C  -1   1  -1   1   0
G  -1  -1   1  -1   0
T  -1  -1  -1   1   0
N   0   0   0   0   0
"""

COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')
MIN_READS = 10  # reads needed to call a SNP or test a CpG

#######################
## Argument Parsing  ##
#######################

EXAMPLES = """
Examples:
  # One or more samples: methylation figure and SNP genotypes for each
  python3 methyl_aligner.py run --sample sample1.fastq sample2.fastq

  # Case vs control: differential methylation, plus SNP genotypes for every sample
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
        """Print usage, the error and how to get help, then exit."""
        self.print_usage(sys.stderr)
        sys.stderr.write(f"\nError: {message}\nFor help: python3 {self.prog} -h\n")
        sys.exit(2)

def parse_arguments(argv=None):
    """Parse and check all command line arguments."""
    argv = sys.argv[1:] if argv is None else list(argv)
    parser = FriendlyParser(
        prog='methyl_aligner.py',
        description='Bisulfite amplicon long read alignment, CpG methylation figures, '
                    'SNP genotyping and case vs control differential methylation.',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EXAMPLES)
    commands = parser.add_subparsers(dest='command', metavar='command', title='commands')

    def add_command(name, help_text, description, samples=None, align=True, figures=True):
        """Add one command (run, align or plot) with the options it uses."""
        p = commands.add_parser(name, help=help_text, description=description, epilog=EXAMPLES,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
        if samples:
            group = p.add_argument_group('samples (use --sample, or --case and --control)')
            group.add_argument('--sample', nargs='+', metavar=samples,
                               help='one or more samples: a methylation figure and SNP genotypes '
                                    '(rs4823981, rs4823982) for each')
            group.add_argument('--case', nargs='+', metavar=samples,
                               help='case group, compared with --control for differential methylation; '
                                    'each sample also gets its own SNP genotypes')
            group.add_argument('--control', nargs='+', metavar=samples,
                               help='control group, compared with --case')
        else:
            p.add_argument('fastq', nargs='+', metavar='FASTQ', help='one or more FASTQ files')
        if align:
            group = p.add_argument_group('alignment options')
            group.add_argument('-r', '--reference', dest='fasta', metavar='FASTA',
                               help='custom reference (default: built-in NHIP_pro_forward_HY). '
                                    'At your own risk: SNP genotyping is skipped for a custom reference.')
            group.add_argument('-p', '--percent', type=int, default=70, metavar='N',
                               help='minimum percent identity to keep a read (default: 70)')
            group.add_argument('-t', '--threads', type=int, default=4, metavar='N',
                               help='number of threads (default: 4)')
        group = p.add_argument_group('output options')
        if figures:
            group.add_argument('-o', '--name',
                               help='name of the analysis folder (default: the sample name for --sample, '
                                    'case_vs_control for --case/--control; per sample SNP folders '
                                    'always use the sample name)')
            group.add_argument('--color', action='store_true',
                               help='give each case and control sample its own color')
        group.add_argument('--outdir', default='results', metavar='DIR',
                           help='folder for all results (default: results)')
        return p

    subparsers = {
        'run': add_command('run', 'start here: align FASTQ files and make the figures in one step',
                           'Align FASTQ files, then make the figures. Use --sample for SNP genotypes, '
                           'or --case and --control for differential methylation.', samples='FASTQ'),
        'align': add_command('align', 'align FASTQ files only (step 1 of 2)',
                             'Align FASTQ files (.fastq, .fq, or gzipped) to the reference. '
                             'Writes results/samples/<sample>/ for each file.', figures=False),
        'plot': add_command('plot', 'make the figures from aligned samples (step 2 of 2)',
                            'Make figures from the matrix CSVs that "align" wrote to results/samples/. '
                            'Use --sample for SNP genotypes, or --case and --control for differential '
                            'methylation.', samples='CSV', align=False),
    }

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
    csvs = [f for f in inputs if f.lower().endswith('.csv')]
    if args.command == 'plot' and len(csvs) < len(inputs):
        p.error('plot takes the matrix CSVs made by align (results/samples/<sample>/<sample>_matrix.csv). '
                'To start from FASTQ files, use run instead.')
    if args.command != 'plot' and csvs:
        p.error(f'{args.command} takes FASTQ files. To use matrix CSVs from an earlier run, use plot instead.')
    if getattr(args, 'fasta', None) and not os.path.isfile(args.fasta):
        p.error(f'reference not found: {args.fasta}')
    return args

def sample_name(path):
    """Sample name from a FASTQ or matrix file: no folder, extension, or _matrix suffix."""
    name = os.path.basename(path)
    for ext in ('.gz', '.fastq', '.fq', '.csv'):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name[:-len('_matrix')] if name.endswith('_matrix') else name

#########################
## Alignment Functions ##
#########################

def align_read(seq, ref_file, mat_file):
    """Align a read and its reverse complement with needle; keep the strand with higher identity.
    Returns (reference alignment, read alignment, A/G percent identity)."""
    results = []
    for strand in (seq, seq.translate(COMPLEMENT)[::-1]):
        out = subprocess.run(['needle', '-asequence', ref_file, '-bsequence', f'asis::{strand}',
                              '-gapopen', '3', '-gapextend', '1', '-datafile', mat_file,
                              '-endweight', '-outfile', 'stdout', '-auto'],
                             capture_output=True, text=True, check=True).stdout
        ref_aln, read_aln = '', ''
        for line in out.splitlines():
            parts = line.split()
            if not line.startswith('#') and len(parts) == 4 and parts[1].isdigit():
                if parts[0] == 'r':
                    ref_aln += parts[2]
                else:
                    read_aln += parts[2]
        # Percent identity from A and G only: bisulfite conversion changes C to T
        pairs = [(r, s) for r, s in zip(ref_aln, read_aln) if r in 'AG' and s != '-']
        identity = sum(r == s for r, s in pairs) / len(pairs) if pairs else 0
        results.append((ref_aln, read_aln, identity))
    forward, reverse = results
    return forward if forward[2] > reverse[2] else reverse

def align_fastq(fastq, ref_seq, percent, threads):
    """Align every read in a FASTQ and return base counts at each reference position."""
    opener = gzip.open if fastq.endswith('.gz') else open
    with opener(fastq, 'rt') as fh:
        reads = [line.strip().upper() for i, line in enumerate(fh) if i % 4 == 1]
    unique_reads = list(dict.fromkeys(reads))  # identical reads are aligned and counted once
    if len(unique_reads) < len(reads):
        sys.stderr.write(f"{sum(c > 1 for c in Counter(reads).values())} duplicates found and skipped\n")

    with tempfile.TemporaryDirectory(prefix='methyl_aligner_') as wdir:
        ref_file, mat_file = os.path.join(wdir, 'ref.fa'), os.path.join(wdir, 'mat')
        with open(ref_file, 'w') as fh:
            fh.write(f">r\n{ref_seq}\n")
        with open(mat_file, 'w') as fh:
            fh.write(BISULFITE_MATRIX)
        with ThreadPoolExecutor(max_workers=threads) as pool:
            alignments = list(pool.map(lambda s: align_read(s, ref_file, mat_file), unique_reads))

    kept = [a for a in alignments if 100 * a[2] >= percent]
    sys.stderr.write(f"{len(alignments) - len(kept)} alignments below {percent}% skipped ({len(kept)} total)\n")
    if not kept:
        return None
    # Read base at each reference position (insertions in the read are dropped)
    bases = [''.join(s for r, s in zip(ref_aln, read_aln) if r != '-') for ref_aln, read_aln, _ in kept]
    columns = [Counter(col) for col in zip_longest(*bases, fillvalue=' ')][:len(ref_seq)]
    columns += [Counter()] * (len(ref_seq) - len(columns))
    df = pd.DataFrame({'POS': range(1, len(ref_seq) + 1), 'REF': list(ref_seq)})
    for nt in 'ACGT':
        df[nt] = [c[nt] for c in columns]
    return df

#############################
## Visualization Functions ##
#############################

def cpg_sites(df):
    """Rows of a count matrix where the reference has C followed by G."""
    return df[(df['REF'] == 'C') & (df['REF'].shift(-1) == 'G')]

################################
## Methylation Plot Functions ##
################################

def is_default_reference(df):
    """True if a count matrix was aligned to the built-in reference."""
    return ''.join(df.sort_values('POS')['REF'].astype(str)).upper() == DEFAULT_REFERENCE_SEQ

def call_snps(df):
    """Genotype each SNP from REF and ALT read counts with a binomial test against 50:50."""
    calls = []
    for rsid, pos, ref, alt in SNP_POSITIONS:
        row = df[df['POS'] == pos]
        if row.empty:
            calls.append((rsid, pos, ref, alt, 0, 0, "Not found", None))
            continue
        r, a = row[[ref, alt]].values[0]
        if r + a < MIN_READS:
            calls.append((rsid, pos, ref, alt, r, a, "Low coverage", None))
            continue
        p_val = binomtest(a, n=r + a, p=0.5).pvalue
        # Not significantly different from 50:50 means heterozygous
        genotype = f"{ref}/{alt}" if p_val >= 0.05 else (f"{alt}/{alt}" if a > r else f"{ref}/{ref}")
        calls.append((rsid, pos, ref, alt, r, a, genotype, p_val))
    return calls

def differential_methylation(case_cpgs, control_cpgs):
    """Test case vs control methylation at each CpG and across the whole region.

    With 2 or more samples per group, a Welch t-test on per-sample percent methylation
    (samples are the unit, so PCR duplicates do not inflate p values). With 1 sample in
    either group, Fisher's exact test on pooled C and T read counts.
    Per-CpG p values are Benjamini-Hochberg corrected.
    """
    def counts(cpgs):
        """C and T counts per CpG for each sample; CpGs below MIN_READS set to missing."""
        C = pd.concat([d.set_index('POS')['C'] for d in cpgs], axis=1, keys=range(len(cpgs)))
        T = pd.concat([d.set_index('POS')['T'] for d in cpgs], axis=1, keys=range(len(cpgs)))
        low = (C + T) < MIN_READS
        return C.mask(low), T.mask(low)

    case_C, case_T = counts(case_cpgs)
    ctrl_C, ctrl_T = counts(control_cpgs)
    case_pct = case_C / (case_C + case_T) * 100
    ctrl_pct = ctrl_C / (ctrl_C + ctrl_T) * 100
    replicates = len(case_cpgs) >= 2 and len(control_cpgs) >= 2

    def test(case_values, ctrl_values, case_counts, ctrl_counts):
        """P value for case vs control: Welch t-test with replicates, otherwise Fisher's exact on reads."""
        if replicates:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)  # identical values at a CpG
                p_val = ttest_ind(case_values, ctrl_values, equal_var=False).pvalue
        else:
            p_val = fisher_exact([case_counts, ctrl_counts]).pvalue
        return 1.0 if np.isnan(p_val) else p_val  # no variance in either group

    rows = []
    for pos in case_pct.index.intersection(ctrl_pct.index):
        a, b = case_pct.loc[pos].dropna(), ctrl_pct.loc[pos].dropna()
        if min(len(a), len(b)) < (2 if replicates else 1):
            continue
        p_val = test(a, b, [case_C.loc[pos].sum(), case_T.loc[pos].sum()],
                     [ctrl_C.loc[pos].sum(), ctrl_T.loc[pos].sum()])
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
    case_region, ctrl_region = case_pct.loc[tested].mean(), ctrl_pct.loc[tested].mean()
    region_p = test(case_region, ctrl_region,
                    [case_C.loc[tested].sum().sum(), case_T.loc[tested].sum().sum()],
                    [ctrl_C.loc[tested].sum().sum(), ctrl_T.loc[tested].sum().sum()])
    region = {'test': 'Welch t-test' if replicates else "Fisher's exact (reads, no replicates)",
              'n_case': len(case_cpgs), 'n_control': len(control_cpgs),
              'delta': round(case_region.mean() - ctrl_region.mean(), 1) + 0.0, 'p_value': region_p,
              'n_sig': int((stats_df['fdr'] < 0.05).sum()), 'n_tested': len(stats_df)}
    return stats_df, region

##########################
## Pipeline Functions   ##
##########################

def align_samples(args, fastq_files):
    """Align each FASTQ, save its matrices and plots, and return {fastq: matrix CSV path}."""
    if shutil.which('needle') is None:
        sys.stderr.write("You must install 'emboss' and have 'needle' in your path.\n"
                         "If it is in a conda environment, activate it first, for example: conda activate emboss\n")
        sys.exit(1)
    ref_seq = DEFAULT_REFERENCE_SEQ
    if args.fasta:
        with open(args.fasta) as fh:
            ref_seq = ''.join(line.strip() for line in fh if not line.startswith('>'))
    if ref_seq.upper() != DEFAULT_REFERENCE_SEQ:
        sys.stderr.write(f"\nWARNING: using a custom reference ({args.fasta}) instead of the default\n"
                         f"NHIP_pro_forward_HY reference. You are proceeding at your own risk.\n"
                         f"SNP genotyping is skipped for this reference.\n\n")

    colors = {'A': '#009E73', 'C': '#0072B2', 'G': '#000000', 'T': '#D55E00'}  # Okabe-Ito
    matrices = {}
    for n, fastq in enumerate(fastq_files, 1):
        name = sample_name(fastq)
        print(f"\n[{n}/{len(fastq_files)}] Aligning {name}")
        df = align_fastq(fastq, ref_seq, args.percent, args.threads)
        if df is None:
            print(f"Skipping {name}: no reads passed the {args.percent}% identity filter")
            continue

        folder = os.path.join(args.outdir, 'samples', name)
        os.makedirs(folder, exist_ok=True)
        df.to_csv(f'{folder}/{name}_matrix.csv', index=False)
        cpg_sites(df).to_csv(f'{folder}/{name}_matrix_ref_c.csv', index=False)
        print(f"Data saved to {folder}/{name}_matrix.csv and {folder}/{name}_matrix_ref_c.csv")

        # Nucleotide count plots: all positions and CpG sites only
        for data, label, suffix in ((df, 'All Data', 'all_data'), (cpg_sites(df), 'CpG Only', 'cpg_only')):
            total = data[list('ACGT')].sum().sum()
            fig, ax = plt.subplots(figsize=(7.2, 3.0))
            for nt in 'ACGT':
                ax.plot(data['POS'], data[nt], label=nt, color=colors[nt], linewidth=0.7)
            ax.set_xlabel('Position in amplicon (bp)')
            ax.set_ylabel('Read count')
            ax.set_title(f'Nucleotide Counts by Position - {name} ({label})', loc='left')
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
            ax.legend(title='Base', bbox_to_anchor=(1.01, 1), loc='upper left', alignment='left')
            ax.text(1.02, 0.42, f"C: {data['C'].sum() / total * 100:.1f}%\nT: {data['T'].sum() / total * 100:.1f}%",
                    transform=ax.transAxes, va='top', ha='left')
            fig.tight_layout()
            path = f'{folder}/{name}_line_plot_{suffix}.pdf'
            fig.savefig(path, bbox_inches='tight')
            plt.close(fig)
            print(f"Plot saved as {path}")
        matrices[fastq] = f'{folder}/{name}_matrix.csv'
    return matrices

def make_analysis(args, name, groups):
    """Make one analysis in <outdir>/analyses/<name>/: figure, SNP or statistics table, summary.
    groups is [('Sample', [csv])] for SNP genotyping or [('Case', csvs), ('Control', csvs)]."""
    print(f"\nAnalysis: {name}")
    folder = os.path.join(args.outdir, 'analyses', name)
    os.makedirs(folder, exist_ok=True)
    base = os.path.join(folder, name)
    sample_mode = groups[0][0] == 'Sample'
    files, notes, snps, region, sig_positions = [f'{base}.pdf'], [], [], None, []

    # Load each matrix; the methylation line is percent C at each CpG
    loaded = []
    for group, paths in groups:
        print(f"Reading {group.lower()} files...")
        for path in paths:
            df = pd.read_csv(path)
            cpgs = cpg_sites(df.sort_values('POS')).copy()
            cpgs['percent_C'] = cpgs['C'] / cpgs[list('ACGT')].sum(axis=1) * 100
            loaded.append((group, path, df, cpgs))

    fig, ax = plt.subplots(figsize=(7.2, 3.2))  # double-column journal width
    sigma = 1.5  # smoothing
    # Okabe-Ito colorblind-safe palette, used in this fixed order for --color mode
    palette = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#000000']
    group_colors = {'Case': '#D55E00', 'Control': '#0072B2', 'Sample': '#000000'}
    group_sizes = Counter(group for group, *_ in loaded)
    seen = Counter()
    for i, (group, path, df, cpgs) in enumerate(loaded):
        smoothed = gaussian_filter1d(cpgs['percent_C'], sigma)
        if args.color:
            ax.plot(cpgs['POS'], smoothed, color=palette[i % len(palette)], linewidth=1.2,
                    label=f'{group}: {sample_name(path)}')
            continue
        # One line per sample, one legend entry per group
        n = group_sizes[group]
        label = f'{group}: {sample_name(path)}' if n == 1 else (f'{group} (n = {n})' if not seen[group] else None)
        seen[group] += 1
        ax.plot(cpgs['POS'], smoothed, color=group_colors[group], linewidth=1.0, alpha=0.8, label=label)

    # CpG sites as tick marks along the x axis
    cpg_positions = pd.concat([cpgs['POS'] for *_, cpgs in loaded]).unique()
    ax.plot(cpg_positions, np.full(len(cpg_positions), 1.5), linestyle='none', marker='|',
            markersize=4, markeredgewidth=0.8, color='0.55', label='CpG site')

    if sample_mode:
        # SNP genotypes; positions only apply to the built-in reference
        _, path, df, _ = loaded[0]
        if is_default_reference(df):
            snps = call_snps(df)
        else:
            print(f"WARNING: {path} was not aligned to the default NHIP_pro_forward_HY reference.\n"
                  f"SNP genotyping is skipped for a custom reference.")
        snp_notes = []
        for rsid, pos, ref, alt, r, a, genotype, p_val in snps:
            print(f"{rsid} ({ref}>{alt}) position {pos}: {genotype}  {ref}={r} {alt}={a}")
            if p_val is None:
                snp_notes.append(f"{rsid} {genotype.lower()}")
                continue
            ax.axvline(x=pos, color='0.45', linestyle='--', linewidth=0.8, zorder=0)
            ax.text(pos, 1.01, f"{rsid}\n{genotype}", transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=7)
            snp_notes.append(f"{rsid} {genotype} (P = {p_val:.2g})")
        if snps:
            notes.append(f"Genotypes ({sample_name(path)}): " + ", ".join(snp_notes) + "; P < 0.05 means homozygous")
            pd.DataFrame([(sample_name(path),) + tuple(s) for s in snps],
                         columns=['sample', 'rsid', 'position', 'ref', 'alt', 'ref_reads', 'alt_reads',
                                  'genotype', 'p_value']).to_csv(f'{base}_snp_genotypes.csv', index=False)
            print(f"SNP genotypes saved to {base}_snp_genotypes.csv")
            files.append(f'{base}_snp_genotypes.csv')
    else:
        # Differential methylation: case vs control
        stats_df, region = differential_methylation([c for g, *_, c in loaded if g == 'Case'],
                                                    [c for g, *_, c in loaded if g == 'Control'])
        if region:
            stats_df.to_csv(f'{base}_differential_methylation.csv', index=False)
            print(f"Differential methylation ({region['test']}): mean difference "
                  f"{region['delta']:+.1f}%, p = {region['p_value']:.3g}, "
                  f"{region['n_sig']} of {region['n_tested']} CpGs FDR < 0.05")
            print(f"Per-CpG statistics saved to {base}_differential_methylation.csv")
            files.append(f'{base}_differential_methylation.csv')
            sig_positions = [int(p) for p in stats_df.loc[stats_df['fdr'] < 0.05, 'POS']]
            for pos in sig_positions:
                ax.text(pos, 94, '*', ha='center', va='center', fontsize=12, fontweight='bold')
            test = 't-test' if region['test'] == 'Welch t-test' else "Fisher's exact on reads"
            notes.append(f"* Significant CpG (FDR < 0.05): {region['n_sig']} of {region['n_tested']}. "
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
    fig.savefig(f'{base}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to {base}.pdf")

    # Summary: command, inputs, settings and results
    lines = [f"methyl_aligner analysis: {name}",
             f"Date: {datetime.datetime.now():%Y-%m-%d %H:%M}",
             "Command: " + shlex.join([os.path.basename(sys.argv[0])] + sys.argv[1:]), ""]
    if getattr(args, 'percent', None) is not None:
        lines += ["Alignment", f"  reference: {args.fasta or 'built-in NHIP_pro_forward_HY'}",
                  f"  percent identity minimum: {args.percent}%", ""]
    lines.append("Inputs")
    for group, paths in groups:
        lines.append(f"  {group} (n = {len(paths)})")
        for path in paths:
            df = next(d for _, p, d, _ in loaded if p == path)
            ref = 'default reference' if is_default_reference(df) else 'CUSTOM reference'
            lines.append(f"    {sample_name(path)}: {path} ({ref})")
    lines.append("")
    if snps:
        lines.append(f"SNP genotypes ({sample_name(loaded[0][1])}), binomial test vs 50:50, P < 0.05 means homozygous")
        for rsid, pos, ref, alt, r, a, genotype, p_val in snps:
            p_txt = f"P = {p_val:.2g}" if p_val is not None else "not called"
            lines.append(f"  {rsid} {ref}>{alt} (position {pos}): {genotype}, {ref} = {r}, {alt} = {a}, {p_txt}")
        lines.append("")
    if region:
        lines += [f"Differential methylation ({region['test']}, n = {region['n_case']} vs {region['n_control']})",
                  f"  overall: case is {region['delta']:+.1f}% vs control, P = {region['p_value']:.2g}",
                  f"  significant CpGs (FDR < 0.05): {region['n_sig']} of {region['n_tested']}"
                  + (f", at positions {', '.join(map(str, sig_positions))}" if sig_positions else ""), ""]
    lines.append("Files")
    lines += [f"  {f}" for f in files + [f'{base}_summary.txt']]
    with open(f'{base}_summary.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"Summary saved to {base}_summary.txt")
    return folder

###################
## Main Function ##
###################

def main(argv=None):
    """Run the chosen command: align, plot, or run (align then plot)."""
    args = parse_arguments(argv)

    if args.command in ('align', 'run'):
        fastqs = args.fastq if args.command == 'align' else (args.sample or args.case + args.control)
        matrices = align_samples(args, fastqs)
        if args.command == 'align':
            if matrices:
                print(f"\nDone. Results for {len(matrices)} sample(s) are in {os.path.join(args.outdir, 'samples')}/")
                print("Next, make the figures with plot, for example:")
                print("  python3 methyl_aligner.py plot --sample " + " ".join(matrices.values()))
            return
        missing = [f for f in fastqs if f not in matrices]
        if missing:
            sys.stderr.write("Cannot make figures: these samples had no aligned reads: " + ", ".join(missing) + "\n")
            sys.exit(1)
        args.sample = [matrices[f] for f in args.sample or []]
        args.case = [matrices[f] for f in args.case or []]
        args.control = [matrices[f] for f in args.control or []]

    if args.sample:
        folders = [make_analysis(args, args.name or sample_name(f), [('Sample', [f])]) for f in args.sample]
    else:
        folders = [make_analysis(args, args.name or 'case_vs_control',
                                 [('Case', args.case), ('Control', args.control)])]
        # Every case and control sample also gets its own SNP genotypes
        folders += [make_analysis(args, sample_name(f), [('Sample', [f])]) for f in args.case + args.control]
    print("\nDone. Results are in:")
    for folder in folders:
        print(f"  {folder}/")

if __name__ == '__main__':
    main()
