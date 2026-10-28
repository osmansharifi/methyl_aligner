#!/usr/bin/env bash
# Runs every methyl_aligner feature on the example data.
# From the repository folder, with the conda environment active:
#   bash example/run_examples.sh
# Takes about a minute. Everything is written to results/.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== 1/5  run one sample: align, methylation figure, SNP genotypes"
python3 methyl_aligner.py run --sample example/test_sample.fq

echo
echo "== 2/5  plot several samples at once from saved matrices (no realigning)"
python3 methyl_aligner.py plot --sample example/matrices/sim_case1_matrix.csv \
                                        example/matrices/sim_control1_matrix.csv

echo
echo "== 3/5  case vs control differential methylation (3 simulated cases, 3 simulated controls)"
python3 methyl_aligner.py plot \
    --case example/matrices/sim_case1_matrix.csv example/matrices/sim_case2_matrix.csv example/matrices/sim_case3_matrix.csv \
    --control example/matrices/sim_control1_matrix.csv example/matrices/sim_control2_matrix.csv example/matrices/sim_control3_matrix.csv

echo
echo "== 4/5  same comparison with one color per sample, saved under a custom name"
python3 methyl_aligner.py plot --color -o case_vs_control_colors \
    --case example/matrices/sim_case1_matrix.csv example/matrices/sim_case2_matrix.csv example/matrices/sim_case3_matrix.csv \
    --control example/matrices/sim_control1_matrix.csv example/matrices/sim_control2_matrix.csv example/matrices/sim_control3_matrix.csv

echo
echo "== 5/5  align only (step 1 of 2), with 8 threads"
python3 methyl_aligner.py align -t 8 example/test_sample.fq

echo
echo "All examples finished. Outputs:"
find results -type f | sort
