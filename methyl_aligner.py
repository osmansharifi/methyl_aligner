#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
import threading
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from queue import Queue
from collections import defaultdict

#######################
## Argument Parsing  ##
#######################

def parse_arguments():
    """Parse and return all command line arguments."""
    parser = argparse.ArgumentParser(description='Align FASTQ reads to a FASTA reference and visualize results.')
    
    # Required arguments
    parser.add_argument('fasta', help='Reference FASTA file')
    parser.add_argument('fastq', help='Query FASTQ file')
    
    # Alignment options
    parser.add_argument('-p', '--percent', type=int, default=70,
                       help='Percent identity minimum (default: 70)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads (default: 4)')
    parser.add_argument('-w', '--wrap', type=int, default=50,
                       help='Wrap length for output (default: 50)')
    parser.add_argument('-d', '--tempdir', 
                       help='Create and keep a named working directory for temp files')
    parser.add_argument('-m', '--methyl', action='store_true',
                       help='Methyl mode: uses asymmetric scoring matrix with C:T matching')
    parser.add_argument('-f', '--force', action='store_true',
                       help='Force rewrite over named working directory')
    
    return parser.parse_args()

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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[nt] for nt in seq[::-1]])

def run_needle(ref_file, seq, strand, wdir, tid, methyl_mode, mat_file):
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
    
    match, total = 0, 0
    for r, s in zip(ra, sa):
        if r == '-' or s == '-':
            continue
        if methyl_mode:
            if r in ['A', 'G']:
                if r == s:
                    match += 1
                total += 1
        else:
            if r == s:
                match += 1
            total += 1
    
    return {'ref': ra, 'seq': sa, 'pct': match/total if total > 0 else 0}

def align_worker(q, results, ref_file, wdir, methyl_mode, mat_file):
    """Worker function for alignment threads."""
    tid = threading.get_ident()
    alignments = []
    while True:
        try:
            plus = q.get_nowait()
        except:
            break
        
        anti_seq = anti(plus)
        a1 = run_needle(ref_file, plus, f"{tid}+", wdir, tid, methyl_mode, mat_file)
        a2 = run_needle(ref_file, anti_seq, f"{tid}-", wdir, tid, methyl_mode, mat_file)
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
    plt.figure(figsize=(10, 6))
    for nucleotide in ['A', 'C', 'G', 'T']:
        plt.plot(df['POS'], df[nucleotide], label=nucleotide, alpha=0.7)

    plt.xlabel('Position')
    plt.ylabel('Count')
    plt.title(title)

    percent_C, percent_T = percentages
    legend = plt.legend(title='Nucleotide', bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.text(1.02, 0.75, f'C: {percent_C:.2f}%\nT: {percent_T:.2f}%', 
             transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', 
             horizontalalignment='left', bbox=dict(facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Plot saved as {output_path}")

##########################
## Pipeline Functions   ##
##########################

def run_alignment(args):
    """Run the alignment portion of the pipeline."""
    # Check for needle
    try:
        subprocess.run(['which', 'needle'], check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        sys.stderr.write("You must install 'emboss' and have 'needle' in your path\n")
        sys.exit(1)
    
    # Set up working directory
    if args.tempdir:
        if os.path.exists(args.tempdir):
            if not args.force:
                sys.stderr.write(f"{args.tempdir} already exists, will not overwrite\n")
                sys.exit(1)
        wdir = args.tempdir
    else:
        wdir = f"temp.{os.getpid()}.methalign"
    
    os.makedirs(wdir, exist_ok=True)
    
    # Process reference FASTA
    with open(args.fasta, 'r') as fh, open(f"{wdir}/ref.fa", 'w') as of:
        of.write(">r\n")
        next(fh)  # Skip header
        ref_seq = ''.join(line.strip() for line in fh)
        of.write(ref_seq + "\n")
    
    # Process FASTQ and check for duplicates
    read_counts = defaultdict(int)
    with open(args.fastq, 'r') as fq:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline().strip()
            fq.readline()  # +
            fq.readline()  # quality
            read_counts[seq] += 1
    
    duplicates = sum(1 for cnt in read_counts.values() if cnt > 1)
    if duplicates:
        sys.stderr.write(f"{duplicates} duplicates found and skipped\n")
    
    # Create scoring matrix
    normal_matrix = """    A   C   G   T
A   1  -1  -1  -1
C  -1   1  -1   1
G  -1  -1   1  -1
T  -1  -1  -1   1
"""
    
    methyl_matrix = """    A   C   G   T
A   1  -1  -1  -1
C  -1   1  -1   1
G  -1  -1   1  -1
T  -1  -1  -1   1
"""
    
    with open(f"{wdir}/mat", 'w') as mat:
        mat.write(methyl_matrix if args.methyl else normal_matrix)
    
    # Set up threading
    query_queue = Queue()
    result_queue = Queue()
    
    for seq in read_counts.keys():
        query_queue.put(seq)
    
    workers = []
    for _ in range(args.threads):
        t = threading.Thread(target=align_worker, args=(query_queue, result_queue, 
                                               f"{wdir}/ref.fa", wdir, 
                                               args.methyl, f"{wdir}/mat"))
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
        if not args.tempdir:
            shutil.rmtree(wdir)
        sys.exit(0)
    
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
    
    # Clean up if using temporary directory
    if not args.tempdir:
        shutil.rmtree(wdir)
    
    return df, ref_seq

def run_visualization(df, sample_name, ref_seq):
    """Run the visualization portion of the pipeline."""
    # Create output directory based on sample name
    output_dir = os.path.join("results", sample_name)
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
    
    # Save reference sequence
    with open(f'{output_dir}/{sample_name}_reference.txt', 'w') as f:
        f.write(ref_seq)
    
    return output_dir

###################
## Main Function ##
###################
if __name__ == '__main__':
        # Parse all arguments first
    args = parse_arguments()
    
    # Get sample name from FASTQ filename
    sample_name = get_sample_name(args.fastq)
    
    # Run alignment pipeline
    df, ref_seq = run_alignment(args)
    
    # Run visualization pipeline
    output_dir = run_visualization(df, sample_name, ref_seq)
    
    print(f"\nAll outputs saved to: {output_dir}")