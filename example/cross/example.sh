#!/bin/bash
set -euo pipefail

# Test workflow with biparental crosses

# Create library from all parents (high-density individuals in masked.ped)
python ../../pythonhmm-script.py -createlib -out lib -ped masked.ped -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 -ncorrect

# Impute using only the parents (in the haplotype library lib.ped) for each low-density individual
python ../../pythonhmm-script.py -impute -out imputed -library lib.ped -ped masked.ped -founders founders.txt -n_haplotypes 20 -seed 42

# Check accuracy - should be 0.966 (with seed 42)
python ../../pythonhmm-script.py -decode true.ped imputed.ped
python accuracy.py

