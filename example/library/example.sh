#!/bin/bash
set -euo pipefail

# Test workflow with library update 

# Create library from half the high-density individuals
python ../../pythonhmm-script.py -createlib -out lib1 -ped hd_1.ped -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 -ncorrect

# Update the library with second half of the high-density individuals
python ../../pythonhmm-script.py -createlib -library lib1.ped -out lib2 -ped hd_2.ped -n_haplotypes 20 -n_sample_rounds 5 -hd_threshold 0.9 -seed 42 -ncorrect

# Impute using this library
python ../../pythonhmm-script.py -impute -out imputed -library lib2.ped -ped masked.ped -n_haplotypes 20 -seed 42

# Check accuracy - should be 0.9993 for HD and 0.9648 for LD (with seed 42)
python ../../pythonhmm-script.py -decode true.ped imputed.ped
python ../accuracy.py

