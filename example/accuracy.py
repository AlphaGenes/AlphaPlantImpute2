import pandas as pd
import numpy as np

# Simple script to check imputation accuracy

# Load files as pandas dataframes
true_genotypes = pd.read_csv('true.genotypes', sep=' ', header=None, index_col=0)
imputed_genotypes = pd.read_csv('imputed.genotypes', sep=' ', header=None, index_col=0)

# High density individuals are the first 100 individuals, low-density the second 100
high_density = np.arange(len(true_genotypes)) < 100
low_density = np.arange(len(true_genotypes)) >= 100

# Simple accuracy measure is correlation between true genotypes and imputed genotypes
accuracy = true_genotypes.corrwith(imputed_genotypes, axis=1)

# Report average accuracy for low- and high-density individuals
print('Mean imputation accuracy')
print(f'  high-density individuals {accuracy[high_density].mean():.4}')
print(f'  low-density individuals  {accuracy[low_density].mean():.4}')
