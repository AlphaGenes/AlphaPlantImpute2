"""Simple haploid HMM implementation"""

import argparse
import numpy as np
from .tinyhouse import Pedigree
from .tinyhouse import InputOutput
from .tinyhouse import BasicHMM


# Create dummy profile decorator if not defined
try:
    profile
except NameError as error:
    def profile(dummy):
        """Dummy decorator"""
        return dummy


def get_args():
    """Get and parse command-line arguments"""
    parser = argparse.ArgumentParser(description='')
    core_parser = parser.add_argument_group("Core arguments")
    core_parser.add_argument('-out', required=True, type=str, help='The output file prefix.')
    InputOutput.addInputFileParser(parser)

    pythonhmm_parser = parser.add_argument_group('pythonhmm Arguments')
    pythonhmm_parser.add_argument('-method', default='dosages', required=False, type=str, choices=['dosages', 'sampler'],
                                  help="Imputation method.")
    pythonhmm_parser.add_argument('-hdthreshold', default=0.9, required=False, type=float,
                                  help='Threshold for determining which haplotypes to include in the haplotype library.')
    pythonhmm_parser.add_argument('-nsamples', default=10, required=False, type=int,
                                  help="Number of samples to to use with the 'sampler' method.")

    return InputOutput.parseArgs("pythonhmm", parser)


def get_haplotype_library(pedigree, threshold):
    """Get haplotype library from high density individuals in a pedigree"""

    haplotype_library = []
    for individual in pedigree:
        for haplotype in individual.haplotypes:
            is_high_density = np.mean(haplotype != 9) >= threshold
            if is_high_density:
                haplotype_library += [haplotype]

    return haplotype_library


def do_imputation(args, haplotype_library, pedigree):
    """Impute individuals using haploidHMM and a haplotype library, setting the dosages for each individual
    The method can be 'dosages' which uses the forward-backward probabilities or 'sampler' which samples haplotypes 
    from the forward and backward probabilities"""
    
    n_loci = pedigree.nLoci
    for individual in pedigree:
        dosages = np.full(n_loci, 0., dtype=np.float32)
        for haplotype in individual.haplotypes:
            dosages += BasicHMM.haploidHMM(haplotype, haplotype_library,
                                           error=0.01, recombinationRate=1/n_loci,
                                           n_samples=args.nsamples,
                                           callingMethod=args.method)
        individual.dosages = dosages


@profile
def main():
    """Main execution code"""

    args = get_args()

    pedigree = Pedigree.Pedigree()
    InputOutput.readInPedigreeFromInputs(pedigree, args, genotypes=True, haps=True, reads=False)

    haplotype_library = get_haplotype_library(pedigree, args.hdthreshold)
    do_imputation(args, haplotype_library, pedigree)

    pedigree.writeDosages(args.out + '.dosages')

if __name__ == "__main__":
    main()
    