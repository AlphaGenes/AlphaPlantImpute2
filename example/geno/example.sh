#!/bin/bash
set -euo pipefail

# Test workflow with genotype files

# Create library from .genotypes (only using HD genotypes in .ped)
AlphaPlantImpute2 -createlib -out lib -genotypes masked.genotypes -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 

# Impute using this library
AlphaPlantImpute2 -impute -out imputed -libphase lib.phase -genotypes masked.genotypes -n_haplotypes 20 -seed 42

# Check accuracy - should be 1.0 for HD and 0.9679 for LD (with seed 42)
python ../accuracy.py

