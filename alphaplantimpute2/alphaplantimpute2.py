"""AlphaPlantImpute2: software for phasing and imputing genotypes in plant populations """

import argparse
import concurrent.futures
import random
import sys
from itertools import repeat
import numpy as np
from numba import jit
from pkg_resources import get_distribution, DistributionNotFound
from . tinyhouse import CombinedHMM, InputOutput, Pedigree, HaplotypeLibrary, ProbMath


# Get version
try:
    __version__ = get_distribution('alphaplantimpute2').version
except DistributionNotFound:
    # Package has not been installed
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

    # Core options
    core_parser = parser.add_argument_group('Core Options')
    core_parser.add_argument('-out', required=True, type=str, help='The output file prefix.')

    # Input options
    input_parser = parser.add_argument_group('Input Options')
    InputOutput.add_arguments_from_dictionary(input_parser, InputOutput.get_input_options(), options=['genotypes', 'pedigree', 'startsnp', 'stopsnp', 'seed'])

    # Algorithm options
    algorithm_parser = parser.add_argument_group('Algorithm Options')
    algorithm_parser.add_argument('-hd_threshold', default=0.9, required=False, type=float,
                                  help='Fraction of non-missing markers required to classify an individual as high-density. '
                                  'Only high-density individuals are used to build the haplotype library. Default: 0.9.')
    algorithm_parser.add_argument('-n_haplotypes', default=100, required=False, type=int,
                                  help='Number of haplotypes to sample from the haplotype library. Default: 100.')
    algorithm_parser.add_argument('-n_sample_rounds', default=10, required=False, type=int,
                                  help='Number of rounds of library refinement. Default: 10.')
    algorithm_parser.add_argument('-n_impute_rounds', default=1, required=False, type=int,
                                  help='Number of rounds of imputation. Default: 1.')
    algorithm_parser.add_argument('-n_bins', default=5, required=False, type=int,
                                  help='Number of bins for targeted haplotype sampling. Default: 5.')
    InputOutput.add_arguments_from_dictionary(algorithm_parser, InputOutput.get_probability_options(), options=['error', 'recombination'])
    
    # Library options
    library_parser = parser.add_argument_group('Haplotype Library Options')
    library_parser.add_argument('-library', required=False, type=str, help='A haplotype library file in [TBD] format. '
                                'Read the haplotype library from file rather than building it from high-density individuals.')

    # Multithreading options
    multithread_parser = parser.add_argument_group('Multithreading Options')
    InputOutput.add_arguments_from_dictionary(multithread_parser, InputOutput.get_multithread_options(), options=['maxthreads', 'iothreads'])

    return InputOutput.parseArgs('alphaplantimpute2', parser)


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


def sample_haplotypes(model, individual, haplotype_library, recombination_rate, error_rate):
    """Sample haplotypes for an individual using haplotype library 'haplotype_library'
    Note: the supplied haplotype_library should *not* contain a copy of the individual's haplotypes"""

    if individual.inbred:
        HaploidHMM.haploidHMM(individual, haplotype_library,
                              error_rate, recombination_rate, calling_method='sample')
    else:
        genotype_probabilities = model.run_HMM(individual=individual, haplotype_library=haplotype_library, algorithm='sample')
        values = ProbMath.call_genotype_probs(genotype_probabilities, calling_threshold=0.1)
        individual.imputed_haplotypes = np.array(values.haplotypes)
    return individual


def refine_library(model, args, individuals, haplotype_library, maf, recombination_rate, error_rate):
    """Refine haplotype library"""

    rounds_str = 'rounds' if args.n_sample_rounds > 1 else 'round'
    print(f'Refining haplotype library: {args.n_sample_rounds} {rounds_str}, {args.n_haplotypes} haplotype samples per round')

    # Loop over iterations
    for iteration in range(args.n_sample_rounds):
        print(f'  Round {iteration}')

        # Generator of subsampled haplotype libraries for ThreadPoolExecutor.map()
        # each library in the generator has the corresponding individual's haplotypes masked out
        haplotype_libraries = (haplotype_library.exclude_identifiers_and_sample(individual.idx, args.n_haplotypes)
                               for individual in individuals)
        # Arguments to pass to sample_haplotypes() via map() and ThreadPoolExecutor.map()
        sample_haplotypes_args = (repeat(model), individuals, haplotype_libraries, repeat(recombination_rate), repeat(error_rate))

        # Sample haplotypes for all individuals in the library
        if args.maxthreads == 1:
            # Single threaded
            results = map(sample_haplotypes, *sample_haplotypes_args)
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(sample_haplotypes, *sample_haplotypes_args)

        # Update library
        for individual in results:
            # Use the input genotype (as read in from data) as the 'true' genotype to correct the haplotypes
            haplotypes = individual.imputed_haplotypes
            if individual.inbred:
                correct_haplotypes(haplotypes, haplotypes, individual.genotypes, maf)
            else:
                correct_haplotypes(haplotypes[0], haplotypes[1], individual.genotypes, maf)
            haplotype_library.update(haplotypes, individual.idx)


def get_dosages(model, individual, haplotype_library, recombination_rate, error_rate):
    """Get dosages for an individual
    HaploidHMM.haploidHMM() and DiploidHMM.diploidHMM() set the individual's dosage member variable"""
    if individual.inbred:
        HaploidHMM.haploidHMM(individual, haplotype_library, error_rate, recombination_rate, calling_method='dosages')
    else:
        genotype_probabilities = model.run_HMM(individual=individual, haplotype_library=haplotype_library, algorithm='marginalize')
        ProbMath.set_from_genotype_probs(individual, geno_probs=genotype_probabilities,
                                         calling_threshold=0.1, set_genotypes=True, set_dosages=True, set_haplotypes=True)
    return individual


def get_phase(individual, haplotype_library, recombination_rate, error_rate):
    """Get phase for an individual
    HaploidHMM.haploidHMM() and DiploidHMM.diploidHMM() set the individual's imputed_haplotypes member variable"""
    if individual.inbred:
        HaploidHMM.haploidHMM(individual, haplotype_library, error_rate, recombination_rate, calling_method='Viterbi')
    else:
        DiploidHMM.diploidHMM(individual, haplotype_library, haplotype_library, error_rate, recombination_rate,
                              calling_method='Viterbi', use_called_haps=False)
    return individual


def impute_individuals(model, args, pedigree, haplotype_library, recombination_rate, error_rate):
    """Impute all individuals in the pedigree"""

    n_loci = pedigree.nLoci
    rounds_str = 'rounds' if args.n_impute_rounds > 1 else 'round'
    print(f'Imputing individuals: {args.n_impute_rounds} {rounds_str}, {args.n_haplotypes} haplotype samples per round')

    # Iterate over all individuals in the Pedigree() object
    individuals = pedigree

    # Set all dosages to zero, so they can be incrementally added to
    dosages = np.zeros((len(individuals), n_loci), dtype=np.float32)

    # Loop over rounds
    for iteration in range(args.n_impute_rounds):
        print(f'  Round {iteration}')

        # Sample the haplotype library for each iteration
        haplotype_library_sample = (haplotype_library.sample_targeted(args.n_haplotypes, individual.genotypes, args.n_bins)
                                    for individual in individuals)
        # Arguments to pass to get_dosages() via map() and ThreadPoolExecutor.map()
        get_dosages_args = (repeat(model), individuals, haplotype_library_sample, repeat(recombination_rate), repeat(error_rate))

        # Get dosages for all individuals
        if args.maxthreads == 1:
            # Single threaded
            results = map(get_dosages, *get_dosages_args)
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(get_dosages, *get_dosages_args)

        # Construct average dosages from results
        for i, individual in enumerate(results):
            dosages[i] += individual.dosages

    # Normalize dosages and generate genotypes
    for i, individual in enumerate(individuals):
        dosages[i] /= args.n_impute_rounds
        if individual.inbred:
            # The calculated dosage for an inbred/DH individual is for its single haplotype
            #   round dosage and then multiply by 2 to get genotype
            #   this prevents inbred/double haploids having loci imputed as heterozygous
            individual.genotypes = 2*np.int8(np.round(dosages[i]))
            individual.dosages = 2*dosages[i]
        else:
            # individual.genotypes = np.int8(np.round(dosages[i]))
            individual.dosages = dosages[i]


def phase_individuals(args, pedigree, haplotype_library, maf, recombination_rate, error_rate):
    """Phase all individuals in the pedigree"""

    print(f'Phasing individuals: {args.n_haplotypes} haplotype samples')

    # Iterate over all individuals in the Pedigree() object
    individuals = pedigree

    # Sample the haplotype library
    haplotype_library_sample = (haplotype_library.sample_targeted(args.n_haplotypes, individual.genotypes, args.n_bins)
                                for individual in individuals)
    # Arguments to pass to get_phase() via map() and ThreadPoolExecutor.map()
    get_phase_args = (individuals, haplotype_library_sample, repeat(recombination_rate), repeat(error_rate))

    # Get dosages for all individuals
    if args.maxthreads == 1:
        # Single threaded
        results = map(get_phase, *get_phase_args)
    else:
        # Multithreaded
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
            results = executor.map(get_phase, *get_phase_args)

    # Loop over results
    for individual in results:
        # (After imputation, `individual.genotype` is the imputed genotype.)
        # Use the imputed genotype to correct the haplotypes, so that imputed phase and genotype are consistent
        haplotypes = individual.imputed_haplotypes
        if individual.inbred:
            correct_haplotypes(haplotypes, haplotypes, individual.genotypes, maf)
            haplotypes = np.vstack([haplotypes, haplotypes])
        else:
            correct_haplotypes(haplotypes[0], haplotypes[1], individual.genotypes, maf)
        # Copy phased haplotypes into `haplotypes` member variable ready for
        # writing out with Pedigree.writePhase()
        individual.haplotypes = haplotypes


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

    # Print boilerplate text
    name = 'AlphaPlantImpute2'
    InputOutput.print_boilerplate(name, __version__)

    # Handle command-line arguments
    args = getargs()

    # Set random seed
    set_seed(args)

    # Read data
    pedigree = Pedigree.Pedigree()
    InputOutput.readInPedigreeFromInputs(pedigree, args, genotypes=True, haps=False, reads=False)
    n_loci = pedigree.nLoci

    # Construct empty haplotypes/genotypes if not read in from file
    for individual in pedigree:
        individual.constructInfo(n_loci, genotypes=True, haps=True)

    # Handle any inbred/double haploid individuals
    handle_inbreds(pedigree)

    # Calculate minor allele frequency and determine high density individuals
    pedigree.setMaf()
    pedigree.high_density_threshold = args.hd_threshold
    pedigree.set_high_density()
    individuals = high_density_individuals(pedigree)
    print(f'Read in {len(pedigree)} individuals, {len(individuals)} are at high-density (threshold {args.hd_threshold})')

    # Error if no high-density individuals
    if len(individuals) == 0:
        print(f'ERROR: zero individuals genotyped at high-density\n'
        'Try reducing the high-density threshold (-hd_threshold), so that more individuals are classed as high-density\n'
        'Exiting...')
        sys.exit(2)

    # Probabilistic rates
    error_rate = np.full(n_loci, args.error, dtype=np.float32)
    recombination_rate = np.full(n_loci, args.recomb/n_loci, dtype=np.float32)
 
    # Model
    model = CombinedHMM.DiploidMarkovModel(n_loci, error_rate, recombination_rate)
    #model = CombinedHMM.HaploidMarkovModel(n_loci, error_rate, recombination_rate)


    # Library
    if args.library:
        # Load haplotype library from file
        print(f'Loading haplotype library from: {args.library}')
        haplotype_library = HaplotypeLibrary.load(args.library)
    else:
        # Create haplotype library from high-density genotypes
        haplotype_library = create_haplotype_library(individuals, pedigree.maf)
        refine_library(model, args, individuals, haplotype_library, pedigree.maf, recombination_rate, error_rate)
        filepath = 'haplotype_library.pkl'
        print(f'Writing haplotype library to: {filepath}')
        HaplotypeLibrary.save(filepath, haplotype_library)

    # Imputation
    impute_individuals(model, args, pedigree, haplotype_library, recombination_rate, error_rate)

    # Phasing
    # phase_individuals(args, pedigree, haplotype_library, pedigree.maf, recombination_rate, error_rate)

    # Output
    pedigree.writeDosages(args.out + '.dosages')
    pedigree.writeGenotypes(args.out + '.genotypes')
    # pedigree.writePhase(args.out + '.phase')

if __name__ == "__main__":
    main()
