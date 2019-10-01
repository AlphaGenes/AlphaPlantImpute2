"""Simple haploid HMM implementation"""

import argparse
import concurrent.futures
from itertools import repeat
from numba import jit
import numpy as np
import random
from .tinyhouse import BasicHMM, InputOutput, Pedigree, HaplotypeLibrary

# Create dummy profile decorator if not defined
try:
    profile
except NameError as error:
    def profile(dummy):
        """Dummy decorator"""
        return dummy


def getargs():
    """Get and parse command-line arguments"""
    parser = argparse.ArgumentParser(description='')
    core_parser = parser.add_argument_group("Core arguments")
    core_parser.add_argument('-out', required=True, type=str, help='The output file prefix.')
    InputOutput.addInputFileParser(parser)

    pythonhmm_parser = parser.add_argument_group('PythonHMM Arguments')
    pythonhmm_parser.add_argument('-hdthreshold', default=0.9, required=False, type=float,
                                  help='Threshold for determining which haplotypes to include in the haplotype library.')
    pythonhmm_parser.add_argument('-nsamples', default=10, required=False, type=int,
                                  help="Number of samples to use with the 'sampler' method.")
    pythonhmm_parser.add_argument('-nhaplotypes', default=200, required=False, type=int,
                                  help="Number of haplotypes to use in the haplotype library.")
    pythonhmm_parser.add_argument('-nrounds', default=20, required=False, type=int,
                                  help="Number of rounds of HMM sampling.")
    pythonhmm_parser.add_argument('-maxthreads', default=1, required=False, type=int,
                                  help='Number of threads to use. Default: 1.')

    return InputOutput.parseArgs("pythonhmm", parser)


def high_density_individuals(pedigree):
    """Return a list of high density individuals in a pedigree
    Note: should probably go in Pedigree()"""
    return [individual
            for individual in pedigree
            if individual.is_high_density]


@jit(nopython=True)
def generate_haplotypes(genotype, maf):
    """Randomly generate a pair of (paternal, maternal) haplotypes from a genotype
    Homozygous loci:    de-facto phased: 0 -> (0,0), 2 -> (1,1)
    Heterozygous loci:  randomly assigned: 50% (0,1) and 50% (1,0)
    Missing loci:       randomly assigned according to minor allele frequency (maf):
                        allele 0 with probability 1-maf, allele 1 with probability maf
    Inputs:
      genotype:   array of genotypes {0,1,2,9} at each locus
      maf:        array of minor allele frequencies at each locus
    Returns (p_hap, m_hap):
      p_hap:      array of paternal haplotypes {0,1} at each locus (dtype np.int8)
      m_hap:      array of maternal haplotypes {0,1} at each locus (dtype np.int8)"""

    n_loci = len(genotype)
    p_hap = np.empty(n_loci, dtype=np.int8)
    m_hap = np.empty(n_loci, dtype=np.int8)
    for i in range(n_loci):
        if genotype[i] == 0:         # homozygous
            p_hap[i] = m_hap[i] = 0
        elif genotype[i] == 2:       # homozygous
            p_hap[i] = m_hap[i] = 1
        elif genotype[i] == 1:       # heterozygous
            p_hap[i] = np.random.randint(0, 2)
            m_hap[i] = 1 - p_hap[i]
        elif genotype[i] == 9:       # missing
            p_hap[i] = np.random.binomial(1, p=maf[i])
            m_hap[i] = np.random.binomial(1, p=maf[i])
        else:
            raise ValueError('Genotype not in {0,1,2,9}')

    return p_hap, m_hap


def create_haplotype_library(individuals, maf):  # can/should this be a member of HaplotypeLibrary() or pedigree()?
    """Create a haplotype library from list of individuals
    The population's minor allele frequency (maf) is used to randomly create alleles at missing loci"""

    n_loci = len(maf)
    haplotype_library = HaplotypeLibrary.HaplotypeLibrary(n_loci=n_loci)
    for individual in individuals:
        paternal_haplotype, maternal_haplotype = generate_haplotypes(individual.genotypes, maf)
        haplotype_library.append(paternal_haplotype, identifier=individual.idx)
        haplotype_library.append(maternal_haplotype, identifier=individual.idx)

    haplotype_library.freeze()
    return haplotype_library


@jit(nopython=True)
def correct_haplotypes(paternal_hap, maternal_hap, true_genotype, maf):
    """Correct paternal and maternal haplotype pair so that they match the true genotype"""

    # Mask to select just the mismatched loci (ignoring missing genotypes)
    mask = ((paternal_hap + maternal_hap) - true_genotype) != 0
    mask &= true_genotype != 9                   # ignore missing genotypes (& is bitwise and)
    assert(np.sum(true_genotype[mask]==9) == 0)  # double check there are no missing genotypes - remove in production code

    # Generate (matching) haplotypes at just these loci
    p, m = generate_haplotypes(true_genotype[mask], maf[mask])

    # Update
    paternal_hap[mask] = p
    maternal_hap[mask] = m


@jit(nopython=True, nogil=True)
def sample_haplotype_pair(genotype, haplotype_library, recombination_rate, error, maf):
    """Sample a pair of haplotypes for an individual with genotype 'genotypes'
    against a haplotype library 'haplotype_library'
    Note: the supplied haplotype_library should not contain a copy of the individual's haplotypes"""

    # Pass missing haplotypes (all 9) to getDiploidPointEstimates(), so that the genotypes are used directly
    n_loci = len(genotype)
    haplotypes = np.full((2, n_loci), 9, dtype=np.int8)
    
    n_pat = haplotype_library.shape[0]
    n_mat = haplotype_library.shape[0]
    point_estimate = np.empty((n_loci, n_pat, n_mat), dtype=np.float32)

    BasicHMM.getDiploidPointEstimates_geno(genotype, haplotype_library, haplotype_library, error, point_estimate)
    forward_probs = BasicHMM.diploidForward(point_estimate, recombination_rate)
    haplotypes = BasicHMM.diploidSampleHaplotypes(forward_probs, recombination_rate,
                                                  haplotype_library, haplotype_library)

    correct_haplotypes(haplotypes[0], haplotypes[1], genotype, maf)

    return haplotypes


@profile
def refine_library(args, individuals, haplotype_library, maf, recombination_rate, error):
    """Refine haplotype library"""

    print(f'Refining haplotype library, {args.nrounds} iterations')

    # List of genotypes and identifiers to iterate over
    genotypes = [individual.genotypes for individual in individuals]
    identifiers = [individual.idx for individual in individuals]

    # Loop
    for iteration in range(args.nrounds):
        print('  Iteration', iteration)

        # Choose random sample of haplotypes for each iteration
        sample = haplotype_library.sample(args.nhaplotypes) # sample is a HaplotypeLibrary()

        # Generator of haplotype libraries for ThreadPoolExecutor.map()
        # each subsequent library has the corresponding individual's haplotypes masked out
        haplotype_libraries = (sample.masked(individual.idx) for individual in individuals)

        # Sample haplotypes for all individuals in the library
        if args.maxthreads == 1:
            # Single threaded
            results = map(sample_haplotype_pair, genotypes, haplotype_libraries,
                          repeat(recombination_rate), repeat(error), repeat(maf))
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(sample_haplotype_pair, genotypes, haplotype_libraries,
                                       repeat(recombination_rate), repeat(error), repeat(maf))

        # Update library
        for haplotypes, identifier in zip(results, identifiers):
            haplotype_library.update_pair(haplotypes[0], haplotypes[1], identifier)


@jit(nopython=True, nogil=True)
def get_dosages(genotype, haplotype_library, recombination_rate, error):
    """"""
    # Pass missing haplotypes (all 9) to getDiploidPointEstimates(), so that the genotypes are used directly
    
    n_loci = len(genotype)
    nPat = haplotype_library.shape[0]
    nMat = haplotype_library.shape[0]

    # Get dosages - should this be a helper function - its like diploidHMM(..., callingMethod='dosages')
    pointEst = np.full((n_loci, nPat, nMat), 1, dtype = np.float32)

    point_estimate = BasicHMM.getDiploidPointEstimates_geno(genotype, haplotype_library, haplotype_library, error, pointEst)
    total_probs = BasicHMM.diploidForwardBackward(point_estimate, recombination_rate)
    dosages = BasicHMM.getDiploidDosages(total_probs, haplotype_library, haplotype_library)

    return dosages


def impute_individuals(args, pedigree, haplotype_library, recombination_rate, error): # haplotype_library_full -> haplotype_library
    """Impute all individuals in the pedigree"""

    n_iterations = 5
    n_loci = pedigree.nLoci
    print(f'Imputing individuals, {n_iterations} iterations')

    # List of genotypes to iterate over
    genotypes = [individual.genotypes for individual in pedigree]

    # Set all dosages to zero, so they can be incrementally added to
    for individual in pedigree:
        individual.dosages = np.full(n_loci, 0., dtype=np.float32)

    # Loop
    for iteration in range(n_iterations):
        print('  Iteration', iteration)

        # Sample the haplotype library for each iteration
        haplotype_library_sample = haplotype_library.sample(args.nhaplotypes)._haplotypes  # UGLY, have a return_library=False option?

        # Get dosages for all individuals
        #   Does it make sense to impute high density individuals? Yes to fill in missingness, but not if it introduces additional errors
        #   Handle LD and HD separately?
        if args.maxthreads == 1:
            # Single threaded
            results = map(get_dosages, genotypes, repeat(haplotype_library_sample),
                          repeat(recombination_rate), repeat(error))
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(get_dosages, genotypes, repeat(haplotype_library_sample),
                                       repeat(recombination_rate), repeat(error))

        # Update dosages from results
        for dosages, individual in zip(results, pedigree):
            individual.dosages += dosages

    # Normalize dosages and generate genotypes
    for individual in pedigree:
        individual.dosages /= n_iterations
        individual.genotypes = np.int8(np.round(individual.dosages))


def set_seed(args):
    """Set random seed for np.random and random and numba equivalents"""
    if args.seed and args.maxthreads > 1:
        print('The random seed option requires maxthreads=1: setting maxthreads to 1')
        args.maxthreads = 1
    if args.seed:
        print(f'Setting random seed to {args.seed}')
        InputOutput.setNumbaSeeds(args.seed)
        np.random.seed(args.seed)
        random.seed(args.seed)


@profile
def main():
    """Main execution code"""

    args = getargs()

    # Random seed
    set_seed(args)

    # Read data
    pedigree = Pedigree.Pedigree()
    InputOutput.readInPedigreeFromInputs(pedigree, args, genotypes=True, haps=False, reads=False)
    n_loci = pedigree.nLoci

    # Calculate MAF and determine high density individuals
    pedigree.setMaf()
    pedigree.high_density_threshold = args.hdthreshold
    pedigree.set_high_density()

    # Various parameters
    error = np.full(n_loci, 0.01, dtype=np.float32)
    recombination_rate = np.full(n_loci, 1/n_loci, dtype=np.float32)

    # Library
    individuals = high_density_individuals(pedigree)
    print('# HD individuals', len(individuals))
    haplotype_library = create_haplotype_library(individuals, pedigree.maf)
    refine_library(args, individuals, haplotype_library, pedigree.maf, recombination_rate, error)

   # Imputation
    impute_individuals(args, pedigree, haplotype_library, recombination_rate, error)

    # Output
    pedigree.writeDosages(args.out + '.dosages')
    pedigree.writeGenotypes(args.out + '.genotypes')


if __name__ == "__main__":
    main()
