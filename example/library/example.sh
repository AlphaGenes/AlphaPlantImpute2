#!/bin/bash
set -euo pipefail

# Test workflow with library update 

# Create library from half the high-density individuals
AlphaPlantImpute2 -createlib -out lib1 -ped hd_1.ped -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 

# Update the library with second half of the high-density individuals
AlphaPlantImpute2 -createlib -library lib1.ped -out lib2 -ped hd_2.ped -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 

# Impute using this library
AlphaPlantImpute2 -impute -out imputed -library lib2.ped -ped masked.ped -n_haplotypes 20 -seed 42

# Check accuracy - should be 1.0 for HD and 0.9638 for LD (with seed 42)
AlphaPlantImpute2 -decode true.ped imputed.ped
python ../accuracy.py

