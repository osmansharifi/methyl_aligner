#!/usr/bin/env bash
# Runs every methyl_aligner feature on the example data.
# From the repository folder, with the conda environment active:
#   bash example/run_examples.sh
# Takes about a minute. Everything is written to results/.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== 1/6  align: example reads -> results/samples/test_sample/"
python3 methyl_aligner.py align example/test_sample.fq

echo
echo "== 2/6  plot one sample with SNP genotypes -> results/analyses/test_sample_snps/"
python3 methyl_aligner.py plot test_sample_snps \
    --sample results/samples/test_sample/test_sample_matrix.csv -snp

echo
echo "== 3/6  case vs control differential methylation -> results/analyses/case_vs_control/"
echo "        (3 simulated cases and 3 simulated controls from example/matrices/)"
python3 methyl_aligner.py plot case_vs_control \
    --case example/matrices/sim_case1_matrix.csv example/matrices/sim_case2_matrix.csv example/matrices/sim_case3_matrix.csv \
    --control example/matrices/sim_control1_matrix.csv example/matrices/sim_control2_matrix.csv example/matrices/sim_control3_matrix.csv

echo
echo "== 4/6  same comparison, one color per sample -> results/analyses/case_vs_control_colors/"
python3 methyl_aligner.py plot case_vs_control_colors --color \
    --case example/matrices/sim_case1_matrix.csv example/matrices/sim_case2_matrix.csv example/matrices/sim_case3_matrix.csv \
    --control example/matrices/sim_control1_matrix.csv example/matrices/sim_control2_matrix.csv example/matrices/sim_control3_matrix.csv

echo
echo "== 5/6  run: align and plot in one step -> results/analyses/run_test_sample/"
python3 methyl_aligner.py run --sample example/test_sample.fq -snp -o run_test_sample

echo
echo "== 6/6  older command still works (percent_mC.py) -> results/analyses/legacy_check/"
python3 percent_mC.py legacy_check --sample results/samples/test_sample/test_sample_matrix.csv -snp

echo
echo "All examples finished. Outputs:"
find results -type f | sort
