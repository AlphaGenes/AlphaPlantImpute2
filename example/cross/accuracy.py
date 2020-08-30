import pandas as pd
import numpy as np

# Simple script to check imputation accuracy

# Load files as pandas dataframes
true_genotypes = pd.read_csv('true.genotypes', sep=' ', header=None, index_col=0)
imputed_genotypes = pd.read_csv('imputed.genotypes', sep=' ', header=None, index_col=0)

# Simple accuracy measure is correlation between true genotypes and imputed genotypes
accuracy = true_genotypes.corrwith(imputed_genotypes, axis=1)

# Report average accuracy for low- and high-density individuals
print(f'Mean imputation accuracy: {accuracy.mean():.4}')
