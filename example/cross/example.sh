#!/bin/bash
set -euo pipefail

# Test workflow with biparental crosses

# Create library from all parents (high-density individuals in masked.ped)
AlphaPlantImpute2 -createlib -out lib -ped masked.ped -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 

# Impute using only the parents (in the haplotype library lib.ped) for each low-density individual
AlphaPlantImpute2 -impute -out imputed -library lib.ped -ped masked.ped -founders founders.txt -n_haplotypes 20 -seed 42

# Check accuracy - should be 0.966 (with seed 42)
AlphaPlantImpute2 -decode true.ped imputed.ped
python accuracy.py

