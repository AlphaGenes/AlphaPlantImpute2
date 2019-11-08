import pandas as pd
import numpy as np

# Load files as pandas dataframes
masked_genotypes = pd.read_csv('masked.genotypes', sep=' ', header=None, index_col=0)
true_genotypes = pd.read_csv('true.genotypes', sep=' ', header=None, index_col=0)
imputed_dosages = pd.read_csv('imputed.dosages', sep=' ', header=None, index_col=0)

# Determine which individuals are genotyped at low- and high-density
low_density = np.mean(masked_genotypes != 9, axis=1) < 0.9
high_density = ~low_density

# Accuracy is correlation between true genotypes and imputed dosages
accuracy = true_genotypes.corrwith(imputed_dosages, axis=1)

# Report average accuracy for low- and high-density individuals
print('Mean imputation accuracy')
print(f'  high-density individuals {accuracy[high_density].mean():.4}')
print(f'  low-density individuals  {accuracy[low_density].mean():.4}')
