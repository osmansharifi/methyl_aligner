# methyl_aligner

methyl_aligner aligns bisulfite-converted long reads from the NHIP promoter
amplicon, plots CpG methylation, genotypes two promoter SNPs (rs4823981 and
rs4823982), and tests case vs control differential methylation.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/pipeline_schematic_dark.png">
  <img alt="methyl_aligner pipeline" src="docs/images/pipeline_schematic.png">
</picture>

## Installation ##

You need EMBOSS (for `needle`) and a few Python packages. The easiest way is
conda. Install miniconda (or anaconda), clone this repository, and from the
repository folder type:

```
conda env create -f environment.yml
conda activate methyl_aligner_env
```

On an Apple silicon Mac, if conda cannot find `emboss`, create the environment
for Intel Macs instead (it runs through Rosetta):

```
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate methyl_aligner_env
conda config --env --set subdir osx-64
```

Remember to `conda activate methyl_aligner_env` before each session.

## How to use it ##

Everything starts from your FASTQ files (plain or gzipped) and one command,
`run`. Pick the line that matches your question:

**Methylation and SNP genotypes for each sample**

```
python3 methyl_aligner.py run --sample sample1.fastq sample2.fastq
```

**Is methylation different between two groups?**

```
python3 methyl_aligner.py run --case case1.fastq case2.fastq case3.fastq \
                              --control ctrl1.fastq ctrl2.fastq ctrl3.fastq
```

That's it. Results go into `results/` (see [Outputs](#outputs)), and the
program tells you where each file was saved.

**Aligning once, plotting many times.** Alignment is the slow part. `run` is
simply `align` followed by `plot`, so you can do the two steps yourself and
re-plot later without realigning, for example to try different groupings:

```
python3 methyl_aligner.py align sample1.fastq sample2.fastq sample3.fastq sample4.fastq
python3 methyl_aligner.py plot --case    results/samples/sample1/sample1_matrix.csv results/samples/sample2/sample2_matrix.csv \
                               --control results/samples/sample3/sample3_matrix.csv results/samples/sample4/sample4_matrix.csv
```

`align` takes FASTQ files; `plot` takes the matrix CSVs that `align` saved.

**Help is always one command away:** `python3 methyl_aligner.py -h`, or
`python3 methyl_aligner.py run -h` for one command.

## Test everything ##

After installing, run every feature on the included example data (about a
minute):

```
bash example/run_examples.sh
```

The SNP calls for the example sample should be T/T at rs4823981 and A/A at
rs4823982, and the case vs control comparison should find 2 significant CpGs,
at positions 118 and 160.

## Options ##

```
samples (run and plot; use --sample, or --case and --control)
  --sample FILE ...     one or more samples: a figure and SNP genotypes for each
  --case FILE ...       case group      } compared for differential methylation
  --control FILE ...    control group   }

alignment options (run and align)
  -r FASTA              custom reference (default: built-in NHIP_pro_forward_HY)
                          at your own risk: SNP genotyping is skipped
  -p N                  minimum percent identity to keep a read [70]
  -t N                  number of threads [4]; more is faster

output options
  -o NAME               name of the analysis folder (run and plot)
                          [the sample name for --sample, case_vs_control for --case/--control]
  --color               give each case and control sample its own color (run and plot)
  --outdir DIR          folder for all results [results]
```

## Outputs ##

Everything goes into one `results/` folder: one subfolder per sample, and one
per analysis.

```
results/
├── samples/
│   └── sample1/
│       ├── sample1_matrix.csv              base counts at every reference position
│       ├── sample1_matrix_ref_c.csv        the same, CpG sites only
│       ├── sample1_line_plot_all_data.pdf  nucleotide counts, all positions
│       └── sample1_line_plot_cpg_only.pdf  nucleotide counts, CpG sites only
└── analyses/
    ├── sample1/                            from --sample
    │   ├── sample1.pdf                     methylation figure with SNP genotypes
    │   ├── sample1_snp_genotypes.csv       SNP calls with read counts and P values
    │   └── sample1_summary.txt             command, inputs, settings, results
    └── case_vs_control/                    from --case and --control
        ├── case_vs_control.pdf                          methylation figure, * = FDR < 0.05
        ├── case_vs_control_differential_methylation.csv statistics for every CpG
        └── case_vs_control_summary.txt                  command, inputs, settings, results
```

Each `summary.txt` records exactly how the analysis was made (command, date,
reference, alignment settings, every input file) and its key results, so any
figure can be traced back to its data.

## Methods ##

**Alignment.** Each read and its reverse complement are aligned to the
reference with EMBOSS `needle` (global alignment, gap open 3, gap extend 1),
using a bisulfite scoring matrix in which a reference C read as T scores as a
match. The strand with the higher identity is kept. Percent identity ignores
gaps and uses only the reference A and G positions, which bisulfite conversion
does not change. Reads below the `-p` identity are dropped. Identical reads are
aligned and counted once; with long reads, exact duplicates are rare. Bases
are then counted at every reference position.

**Reference.** The NHIP_pro_forward_HY reference (801 bp) is built into the
program, so no reference file is needed. With `-r` you can use your own, at
your own risk: the SNP positions are coordinates on NHIP_pro_forward_HY, so
SNP genotyping is skipped (with a warning) for any other reference. Global
alignment is slow and memory hungry, so a custom reference should be just large
enough to cover the amplicon.

**CpG methylation.** At each CpG, percent methylation is the percent of reads
with C (a C protected from conversion is methylated).

**SNP genotyping** (`--sample`). For rs4823981 (G>T, position 238) and
rs4823982 (G>A, position 333), the REF and ALT read counts are tested against
the 50:50 split expected for a heterozygote (binomial test). P < 0.05 means the
reads are too uneven for a heterozygote, so the site is called homozygous for
the allele with more reads; P >= 0.05 is called heterozygous. Sites with fewer
than 10 REF + ALT reads are not called. Genotypes can differ between samples,
so each sample is genotyped on its own, and case vs control figures do not
show SNPs.

**Differential methylation** (`--case` and `--control`). Every CpG with at
least 10 C + T reads in each sample is tested, along with the whole amplicon
(mean methylation of all tested CpGs). P values across CpGs are corrected with
the Benjamini-Hochberg FDR.

* 2 or more samples per group: Welch t-test on per-sample percent methylation.
  Samples are the unit of replication, so PCR duplicates do not inflate P values.
* 1 sample in either group: Fisher's exact test on pooled C and T reads. This
  treats reads as independent and ignores sample-to-sample variation, so treat
  results as exploratory and use biological replicates for real comparisons.

## Example figures ##

### Nucleotide counts, all positions
![Nucleotide counts, all positions](docs/images/test_sample_line_plot_all_data.png)

### Nucleotide counts, CpG sites only
![Nucleotide counts, CpG sites only](docs/images/test_sample_line_plot_cpg_only.png)

### CpG methylation with SNP genotypes (`--sample`)
![CpG methylation with SNP genotypes](docs/images/test_sample_percent_mC.png)

### Differential methylation (`--case` and `--control`)
![Differential methylation, case vs control](docs/images/example_differential_methylation.png)

Example with 3 case and 3 control samples simulated from `example/test_sample.fq`
(the matrices are in `example/matrices/`), where cases were made 30 points more
methylated at 5 CpGs near positions 116 to 166. Each line is one sample, and an
asterisk (*) marks CpGs with FDR < 0.05.

## Repository layout ##

```
methyl_aligner.py           the program (run, align, plot)
environment.yml             conda environment with every dependency
example/
├── test_sample.fq          example bisulfite long reads
├── matrices/               3 simulated case and 3 simulated control matrices
└── run_examples.sh         runs every feature on the example data
docs/images/                pipeline schematic (light and dark) and example figures
```
