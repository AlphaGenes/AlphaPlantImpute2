"""PythonHMM: a program for imputation"""

import argparse
import concurrent.futures
import random
from itertools import repeat
import numpy as np
from numba import jit
from pkg_resources import get_distribution, DistributionNotFound
from .tinyhouse import BasicHMM, InputOutput, Pedigree, HaplotypeLibrary


try:
    __version__ = get_distribution('pythonhmm').version
except DistributionNotFound:
    # Package not installed
    __version__ = None


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


def create_haplotype_library(individuals, maf):
    """Create a haplotype library from list of individuals
    The population's minor allele frequency (maf) is used to
    randomly create alleles at missing loci"""

    n_loci = len(maf)
    haplotype_library = HaplotypeLibrary.HaplotypeLibrary(n_loci=n_loci)
    for individual in individuals:
        paternal_haplotype, maternal_haplotype = generate_haplotypes(individual.genotypes, maf)
        if individual.inbred:
            # Only append one haplotype for an inbred/double haploid individual
            haplotype_library.append(paternal_haplotype, identifier=individual.idx)
        else:
            haplotype_library.append(paternal_haplotype, identifier=individual.idx)
            haplotype_library.append(maternal_haplotype, identifier=individual.idx)

    haplotype_library.freeze()
    return haplotype_library


@jit(nopython=True)
def correct_haplotypes(paternal_haplotype, maternal_haplotype, true_genotype, maf):
    """Correct paternal and maternal haplotype pair so that they match the true genotype"""

    # Mask to select just the mismatched loci (ignoring missing genotypes)
    mask = ((paternal_haplotype + maternal_haplotype) - true_genotype) != 0
    # Ignore missing genotypes (& is bitwise and)
    mask &= (true_genotype != 9)
    # Double check there are no missing genotypes
    assert np.sum(true_genotype[mask] == 9) == 0

    # Generate (matching) haplotypes at just these loci
    paternal, maternal = generate_haplotypes(true_genotype[mask], maf[mask])

    # Update
    paternal_haplotype[mask] = paternal
    maternal_haplotype[mask] = maternal


def correct_haplotypes_inbred(haplotype, true_haplotype, maf):
    
    true_genotype = np.full_like(true_haplotype, 9, dtype=np.int8)
    non_missing_loci = (true_haplotype != 9)
    true_genotype[non_missing_loci] = true_haplotype[non_missing_loci] * 2

    correct_haplotypes(haplotype, haplotype, true_genotype)
    


@jit(nopython=True, nogil=True)
def sample_haplotype_inbred(true_haplotype, true_genotype, haplotype_library, recombination_rate, error, maf):
    """Sample a haplotype for an inbred/double haploid individual using haplotype library 'haplotype_library'
    Returns:
      A single haplotype
    Note: the supplied haplotype_library should NOT contain a copy of the individual's haplotypes"""

    # This block should be (mostly) in BasicHMM through haploidHMM() interface
    point_estimates = BasicHMM.getHaploidPointEstimates(true_haplotype, haplotype_library, error) # cf. getDiploidPointEstimates_GENO
    forward_probs = BasicHMM.haploidForward(point_estimates, recombination_rate)
    haplotype = BasicHMM.haploidSampleHaplotype(forward_probs, haplotype_library, recombination_rate)

    correct_haplotypes(haplotype, haplotype, true_genotype, maf)
    return haplotype


@jit(nopython=True, nogil=True)
def sample_haplotypes_outbred(true_genotype, haplotype_library, recombination_rate, error, maf):
    # This block should be (mostly) in BasicHMM through haploidHMM() interface
    n_loci = len(true_genotype)
    haplotypes = np.empty((2, n_loci), dtype=np.int8)

    n_pat = n_mat = haplotype_library.shape[0]
    point_estimate = np.empty((n_loci, n_pat, n_mat), dtype=np.float32)
    BasicHMM.getDiploidPointEstimates_geno(true_genotype, haplotype_library, haplotype_library,
                                           error, point_estimate)
    forward_probs = BasicHMM.diploidForward(point_estimate, recombination_rate)
    haplotypes = BasicHMM.diploidSampleHaplotypes(forward_probs, recombination_rate,
                                                  haplotype_library, haplotype_library)
    correct_haplotypes(haplotypes[0], haplotypes[1], true_genotype, maf)
    return haplotypes


def sample_haplotypes(individual, haplotype_library, recombination_rate, error, maf):
    """Sample haplotypes for an individual using haplotype library 'haplotype_library'
    Outbreds return a pair of haplotypes as a 2d array of shape (2, n_loci)
    Inbreds return a single haplotype as a 1d array of shape (n_loci,)
    Note: the supplied haplotype_library should not contain a copy of the individual's haplotypes"""

    n_loci = haplotype_library.shape[1]

    if individual.inbred:
        haplotype = individual.haplotypes
        genotype = individual.genotypes
        haplotypes = sample_haplotype_inbred(haplotype, genotype, haplotype_library, recombination_rate, error, maf)
    else:
        genotype = individual.genotypes
        haplotypes = sample_haplotypes_outbred(genotype, haplotype_library, recombination_rate, error, maf)

    return haplotypes


def refine_library(args, individuals, haplotype_library, maf, recombination_rate, error):
    """Refine haplotype library"""

    print(f'Refining haplotype library, {args.nrounds} iterations')

    # Loop over iterations
    for iteration in range(args.nrounds):
        print('  Iteration', iteration)

        # Generator of subsampled haplotype libraries for ThreadPoolExecutor.map()
        # each library in the generator has the corresponding individual's haplotypes masked out
        haplotype_libraries = (haplotype_library.exclude_identifiers_and_sample(individual.idx, args.nhaplotypes)
                               for individual in individuals)

        # Sample haplotypes for all individuals in the library
        if args.maxthreads == 1:
            # Single threaded
            results = map(sample_haplotypes, individuals, haplotype_libraries,
                          repeat(recombination_rate), repeat(error), repeat(maf))
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(sample_haplotypes, individuals, haplotype_libraries,
                                       repeat(recombination_rate), repeat(error), repeat(maf))

        # Update library
        identifiers = [individual.idx for individual in individuals]
        for haplotypes, identifier in zip(results, identifiers):
            haplotype_library.update(haplotypes, identifier)


@jit(nopython=True, nogil=True)
def get_dosages(genotype, haplotype_library, recombination_rate, error):
    """Get dosages for an individual's genotype"""

    n_loci = len(genotype)
    n_pat = n_mat = haplotype_library.shape[0]
    point_estimate = np.ones((n_loci, n_pat, n_mat), dtype=np.float32)

    point_estimate = BasicHMM.getDiploidPointEstimates_geno(genotype, haplotype_library, haplotype_library,
                                                            error, point_estimate)
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
        haplotype_library_sample = haplotype_library.sample(args.nhaplotypes)

        # Get dosages for all individuals
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


def print_boilerplate():
    """Print software name, version etc."""
    name = 'AlphaPlantImpute2'  # get from __name__?
    description = 'Software for phasing and imputing genotypes in plant populations'
    author = 'AlphaGenes Group (http://alphagenes.roslin.ed.ac.uk)'
    version = f'Version: {__version__}'
    box_text = ['', name, '', version, '']
    width = 80
    print('*' * width)
    for line in box_text:
        print(f'*{line:^{width-2}}*')
    print('*' * width)
    extra_text = ['', description, '', author]
    for line in extra_text:
        print(f'{line:<{width}}')
    print('-' * width)
    print(' ' * width)


def handle_inbreds(pedigree):
    """Handle any inbred/double haploid individuals: set any heterozygous loci to missing
    and warn"""
    inbred_individuals = [individual for individual in pedigree if individual.inbred]
    threshold = 0.0  # fraction of heterozygous loci
    # Set heterozygous loci to missing
    for individual in inbred_individuals:
        heterozygous_loci = (individual.genotypes == 1)
        heterozygosity = np.mean(heterozygous_loci)
        if heterozygosity > threshold:
            print(f"Inbred individual '{individual.idx}' has heterozygosity of {heterozygosity}; setting heterozygous loci to missing")
            individual.genotypes[heterozygous_loci] = 9
    # Create haplotype from genotype
    for individual in inbred_individuals:
        n_loci = len(individual.genotypes)
        haplotype = np.full(n_loci, 9, dtype=np.int8)
        non_missing_loci = (individual.genotypes != 9)
        haplotype[non_missing_loci] = individual.genotypes[non_missing_loci] // 2
        individual.haplotypes = haplotype


@profile
def main():
    """Main execution code"""

    print_boilerplate()

    # Handle command-line arguments
    args = getargs()

    # Set random seed
    set_seed(args)

    # Read data
    pedigree = Pedigree.Pedigree()
    InputOutput.readInPedigreeFromInputs(pedigree, args, genotypes=True, haps=False, reads=False)
    n_loci = pedigree.nLoci

    # Handle any inbred/double haploid individuals
    handle_inbreds(pedigree)

    # Calculate minor allele frequency and determine high density individuals
    pedigree.setMaf()
    pedigree.high_density_threshold = args.hdthreshold
    pedigree.set_high_density()
    individuals = high_density_individuals(pedigree)
    print('# HD individuals', len(individuals))

    # Various parameters
    error = np.full(n_loci, 0.01, dtype=np.float32)
    recombination_rate = np.full(n_loci, 1/n_loci, dtype=np.float32)

    # Library
    haplotype_library = create_haplotype_library(individuals, pedigree.maf)
    refine_library(args, individuals, haplotype_library, pedigree.maf, recombination_rate, error)

   # Imputation
    impute_individuals(args, pedigree, haplotype_library, recombination_rate, error)

    # Output
    pedigree.writeDosages(args.out + '.dosages')
    pedigree.writeGenotypes(args.out + '.genotypes')


if __name__ == "__main__":
    main()
