Example data
============

Example data consists of genotypes from 100 simulated individuals with 1000 markers on one chromosome.
These individuals are randomly crossed to give 100 descendants. The descendant genotypes are masked
to simulate a low-density panel with 20 % of the original markers.

Files:
  masked.genotypes - genotypes for the 100 parents and 100 descendants, descendants have 20 % of the markers
  true.genotypes   - as above but with all 1000 markers retained - for imputation accuracy calculation

Imputation with AlphaPlantImpute2
=================================

To impute the genotypes with AlphaPlantImpute2 on the example data, run:

AlphaPlantImpute2 -out imputed -genotypes masked.genotypes -n_sample_rounds 5 -n_haplotypes 20

Imputation accuracy
===================

To report the imputation accuracy for the high-density individuals (parents) and low-density 
individuals (descendants), run:

python accuracy.py
