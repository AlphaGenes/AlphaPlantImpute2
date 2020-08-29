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
    core_parser.add_argument('-createlib', action='store_true', required=False, help='Create a haplotype library.')
    core_parser.add_argument('-impute', action='store_true', required=False, help='Impute individuals.')
    core_parser.add_argument('-decodeped', action='store_true', required=False, help='Decode .ped files.') # remove?
    core_parser.add_argument('-out', required=True, type=str, help='The output file prefix.')

    # Input options
    input_parser = parser.add_argument_group('Input Options')
    InputOutput.add_arguments_from_dictionary(input_parser, InputOutput.get_input_options(), options=['genotypes', 'pedigree', 'startsnp', 'stopsnp', 'seed'])
    input_parser.add_argument('-ped', default=None, required=False, type=str, nargs='*', help='A file in PLINK plain text format (.ped)')
    # input_parser.add_argument('-bim', default=None, required=False, type=str, nargs='*', help='A allele coding file in PLINK plain text format (.bim)')
    input_parser.add_argument('-libped', default=None, required=False, type=str, help='A haplotype library file in PLINK plain text format (.ped & .bim)')
    input_parser.add_argument('-libphase', default=None, required=False, type=str, help='A haplotype library file in AlphaImpute phase format (.phase)')
    input_parser.add_argument('-founders', default=None, required=False, type=str, help='A file that gives the founder individuals for each individual.')
    input_parser.add_argument('-ncorrect', action='store_true', required=False, 
                              help='When building a haplotype library, print the average number of loci that had to be corrected.')



    # Algorithm options
    algorithm_parser = parser.add_argument_group('Algorithm Options')
    algorithm_parser.add_argument('-hd_threshold', default=0.0, required=False, type=float,
                                  help='Fraction of non-missing markers required to classify an individual as high-density. '
                                  'Only high-density individuals are used to build the haplotype library. Default: 0.0.')
    algorithm_parser.add_argument('-n_haplotypes', default=100, required=False, type=int,
                                  help='Number of haplotypes to sample from the haplotype library. Default: 100.')
    algorithm_parser.add_argument('-haploid', action='store_true', required=False, help='Run using a haploid HMM instead of the default diploid HMM.')
    algorithm_parser.add_argument('-joint', action='store_true', required=False, help='Run using a joint HMM instead of the default diploid HMM.')
    algorithm_parser.add_argument('-calling_threshold', default=0.1, required=False, type=float, help='Genotype calling threshold. '
                                  'Use a value less than 0.25 for best-guess genotypes. Default: 0.1.')
    algorithm_parser.add_argument('-n_sample_rounds', default=10, required=False, type=int,
                                  help='Number of rounds of library refinement. Default: 10.')
    algorithm_parser.add_argument('-n_impute_rounds', default=1, required=False, type=int,
                                  help='Number of rounds of imputation. Default: 1.')
    algorithm_parser.add_argument('-n_bins', default=5, required=False, type=int,
                                  help='Number of bins for targeted haplotype sampling. Default: 5.')
    InputOutput.add_arguments_from_dictionary(algorithm_parser, InputOutput.get_probability_options(), options=['error', 'recombination'])
    
    # Library options
    # library_parser = parser.add_argument_group('Haplotype Library Options')
    # library_parser.add_argument('-library', required=False, type=str, help='A haplotype library file in [TBD] format. '
    #                             'Read the haplotype library from file rather than building it from high-density individuals.')

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


def append_random_haplotypes(haplotype_library, args, individuals, maf):
    """Append randomly phased haplotypes to haplotype_library"""
    for individual in individuals:
        paternal_haplotype, maternal_haplotype = generate_haplotypes(individual.genotypes, maf)
        if args.haploid:
            # Only append one haplotype for haploid model - doesn't matter which
            haplotype_library.append(maternal_haplotype, identifier=individual.idx)
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

    # Return number of loci that were corrected
    return np.sum(mask)


def sample_haplotypes(model, individual, haplotype_library, calling_threshold):
    """Sample haplotypes for an individual using haplotype library 'haplotype_library'
    Note: the supplied haplotype_library should *not* contain a copy of the individual's haplotypes"""

    genotype_probabilities = model.run_HMM(individual=individual, haplotype_library=haplotype_library, algorithm='sample')
    ProbMath.set_from_genotype_probs(individual, geno_probs=genotype_probabilities,
                                     calling_threshold=calling_threshold,
                                     set_genotypes=False, set_dosages=False, set_haplotypes=True)
    return individual


def refine_library(model, args, individuals, haplotype_library, maf):
    """Refine haplotype library"""

    rounds_str = 'rounds' if args.n_sample_rounds > 1 else 'round'
    print(f'Building haplotype library: {args.n_sample_rounds} {rounds_str}, {args.n_haplotypes} haplotype samples per round')

    # Loop over iterations
    for iteration in range(args.n_sample_rounds):
        print(f'  Round {iteration}')

        # Generator of subsampled haplotype libraries for ThreadPoolExecutor.map()
        # each library in the generator has the corresponding individual's haplotypes masked out
        haplotype_libraries = (haplotype_library.exclude_identifiers_and_sample(individual.idx, args.n_haplotypes)
                               for individual in individuals)
        # Arguments to pass to sample_haplotypes() via map() and ThreadPoolExecutor.map()
        function_args = (repeat(model), individuals, haplotype_libraries, repeat(args.calling_threshold))

        # Sample haplotypes for all individuals in the library
        if args.maxthreads == 1:
            # Single threaded
            results = map(sample_haplotypes, *function_args)
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(sample_haplotypes, *function_args)

        # Update library
        n_corrected = 0.0
        for individual in results:
            # Use the input genotype (as read in from data) as the 'true' genotype to correct the haplotypes
            haplotypes = individual.haplotypes
            if args.haploid:
                n_corrected += correct_haplotypes(haplotypes[0], haplotypes[0], individual.genotypes, maf)
                haplotype_library.update(haplotypes[0], individual.idx)
            else:
                n_corrected += correct_haplotypes(haplotypes[0], haplotypes[1], individual.genotypes, maf)
                haplotype_library.update(np.array(haplotypes), individual.idx)
        if args.ncorrect:
            print(f'  corrected loci: {n_corrected/len(individuals):.3g}')


def get_dosages(model, individual, haplotype_library, calling_threshold):
    """Get dosages for an individual"""
    genotype_probabilities = model.run_HMM(individual=individual, haplotype_library=haplotype_library, algorithm='marginalize')
    ProbMath.set_from_genotype_probs(individual, geno_probs=genotype_probabilities,
                                     calling_threshold=calling_threshold,
                                     set_genotypes=True, set_dosages=True, set_haplotypes=True)
    return individual


def impute_individuals(model, args, pedigree, haplotype_library):
    """Impute all individuals in the pedigree"""

    n_loci = pedigree.nLoci
    rounds_str = 'rounds' if args.n_impute_rounds > 1 else 'round'
    print(f'Imputing individuals: {args.n_impute_rounds} {rounds_str}, {args.n_haplotypes} haplotype samples per round')

    # Iterate over all individuals in the Pedigree()
    individuals = pedigree

    # Set all dosages to zero, so they can be incrementally added to
    dosages = np.zeros((len(individuals), n_loci), dtype=np.float32)

    # Loop over rounds
    for iteration in range(args.n_impute_rounds):
        print(f'  Round {iteration}')

        # Sample the haplotype library for each iteration
        if args.founders:
            haplotype_library_sample = (haplotype_library.sample_only_identifiers(individual.founders)
                                        for individual in individuals)
        else:
            haplotype_library_sample = (haplotype_library.sample_targeted(args.n_haplotypes,
                                                                          individual.genotypes,
                                                                          args.n_bins)
                                        for individual in individuals)


        # Arguments to pass to get_dosages() via map() and ThreadPoolExecutor.map()
        function_args = (repeat(model), individuals, haplotype_library_sample, repeat(args.calling_threshold))

        # Get dosages for all individuals
        if args.maxthreads == 1:
            # Single threaded
            results = map(get_dosages, *function_args)
        else:
            # Multithreaded
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
                results = executor.map(get_dosages, *function_args)

        # Construct average dosages from results
        for i, individual in enumerate(results):
            dosages[i] += individual.dosages

    # Normalize dosages and generate genotypes
    # NOTE: the following block is temporary
    # it would be better to average over genotype_probabilities and then
    # use ProbMath.set_from_genotype_probs() to call genotypes
    for i, individual in enumerate(individuals):
        dosages[i] /= args.n_impute_rounds
        if not np.allclose(individual.dosages, dosages[i]):
            print(individual.idx, np.sum(np.abs(individual.dosages - dosages[i])))
        individual.dosages = dosages[i]


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
    threshold = 0.0  # fraction of heterozygous loci
    # Set heterozygous loci to missing
    for individual in pedigree: #inbred_individuals:
        heterozygous_loci = (individual.genotypes == 1)
        heterozygosity = np.mean(heterozygous_loci)
        if heterozygosity > threshold:
            print(f"Individual '{individual.idx}' has heterozygosity of {heterozygosity}: setting heterozygous loci to missing")
            individual.genotypes[heterozygous_loci] = 9


def read_library(args):
    """Read in a haplotype library. Returns a HaplotypeLibrary() and allele coding array"""
    assert args.libped or args.libphase
    filename = args.libped if args.libped else args.libphase
    print(f'Reading haplotype library from: {filename}')
    library = Pedigree.Pedigree()
    if args.libped:
        library.readInPed(args.libped, args.startsnp, args.stopsnp, haps=True, get_coding=True)
    elif args.libphase:
        library.readInPhase(args.libphase, args.startsnp, args.stopsnp)
    else:
        # This shouldn't happen
        raise ValueError('No library specified')

    print(f'Haplotype library contains {len(library)} individuals with {library.nLoci} markers')
    haplotype_library = HaplotypeLibrary.HaplotypeLibrary(library.nLoci)
    for individual in library:
        for haplotype in individual.haplotypes:
            haplotype_library.append(haplotype, individual.idx)
    haplotype_library.freeze()

    return haplotype_library, library.allele_coding


def write_library(args, library, allele_coding):
    """Write out haplotype library, first converting to a Pedigree object"""
    # Construct pedigree
    pedigree = Pedigree.Pedigree()
    pedigree.nLoci = library._n_loci
    pedigree.allele_coding = allele_coding
    for idx, haplotype in library:
        individual = pedigree.getIndividual(idx)
        if individual.haplotypes is None:
            individual.haplotypes = haplotype
        else:
            individual.haplotypes = np.vstack([individual.haplotypes, haplotype])
    # Write out
    if args.file_format == 'PLINK':
        print(f'Writing haplotype library to {args.out}.ped')
        pedigree.writePhasePed(f'{args.out}.ped')
    else:
        print(f'Writing haplotype library to {args.out}.phase')
        pedigree.writePhase(f'{args.out}.phase')


def imputation(args, model, genotypes, library, coding): # pass coding via genotypes.allele_coding
    """Handle imputation tasks"""
    assert np.alltrue(genotypes.allele_coding == coding)  # not needed if coding passed via genotypes
    # Impute
    impute_individuals(model, args, genotypes, library)

    # Write out
    print(f'Writing out {args.out}.dosages')
    genotypes.writeDosages(f'{args.out}.dosages')
    if args.file_format == 'PLINK':
        genotypes.allele_coding = coding  # use the library coding
        print(f'Writing out {args.out}.ped')
        genotypes.writeGenotypesPed(f'{args.out}.ped')
    else:
        print(f'Writing out {args.out}.genotypes')
        genotypes.writeGenotypes(f'{args.out}.genotypes')


def create_library(args, model, genotypes, haplotype_library): # pedigree -> genotypes
    """Handle library creation tasks"""
    n_loci = genotypes.nLoci

    # Calculate minor allele frequency and determine high density individuals
    genotypes.setMaf()
    genotypes.high_density_threshold = args.hd_threshold
    genotypes.set_high_density()
    hd_individuals = [individual for individual in genotypes if individual.is_high_density]

    # Error if there are no high-density individuals
    if len(hd_individuals) == 0:
        print(f'ERROR: zero individuals genotyped at high-density\n'
                'Try reducing the high-density threshold (-hd_threshold), so that more individuals are classed as high-density\n'
                'Exiting...')
        sys.exit(2)
    print(f'{len(hd_individuals)} individuals are genotyped at high-density (threshold {args.hd_threshold}) '
          f'and will be used to build the haplotype library')

    # Update an existing library with new genotypes
    if args.libped or args.libphase:
        print(f"Updating haplotype library with {len(hd_individuals)} "
              f'high-density genotypes')

        # Handle any genotypes already present in the haplotype library
        duplicate_identifiers = (set(haplotype_library.get_identifiers()) &
                                 set(individual.idx for individual in hd_individuals))
        if len(duplicate_identifiers) != 0:
            print(f'WARNING: {len(duplicate_identifiers)} genotypes are already present '
                  f'in haplotype library: these will be ignored')
            print(f"  Ignored genotype identifiers are:\n"
                  f"  {' '.join(sorted(duplicate_identifiers))}")
            # Ignore the genotypes
            hd_individuals = [individual for individual in hd_individuals
                             if individual.idx not in duplicate_identifiers]

        haplotype_library.unfreeze()
    else:
        # Create empty haplotype library
        haplotype_library = HaplotypeLibrary.HaplotypeLibrary(n_loci=n_loci)

    append_random_haplotypes(haplotype_library, args, hd_individuals, genotypes.maf)
    refine_library(model, args, hd_individuals, haplotype_library, genotypes.maf)
    write_library(args, haplotype_library, genotypes.allele_coding)


def read_genotypes(args, allele_coding=None):
    """Read genotypes from file(s)"""
    print('Reading genotypes...')
    genotypes = Pedigree.Pedigree(constructor=Pedigree.PlantImputeIndividual)

    # Use provided coding if given
    get_coding = True
    if allele_coding is not None:
        genotypes.allele_coding = allele_coding
        get_coding = False

    InputOutput.readInPedigreeFromInputs(genotypes, args, genotypes=True, haps=False, reads=False, get_coding=get_coding)
    if len(genotypes) == 0:
        print('ERROR: no genotypes supplied. Please supply them with the -genotypes or -ped/-bim options\nExiting...')
        sys.exit(2)
    print(f'Read in {len(genotypes)} individuals with {genotypes.nLoci} markers')

    # Construct empty haplotypes member variables if not read in from file
    for individual in genotypes:
        individual.constructInfo(genotypes.nLoci, genotypes=False, haps=True)

    # Handle inbred/double haploid individuals
    if args.haploid:
        handle_inbreds(genotypes)

    return genotypes


def check_arguments_consistent(args):
    """Checks that a consistent set of arguments has been supplied"""

    # Check one and only one of the major arguments are supplied
    n = np.sum(args.createlib or args.impute or args.decodeped)
    if n == 0:
        print('ERROR: one of -createlib, -impute or -decodeped needs to be specified\nExiting...')
        sys.exit(2)
    elif n > 1:
        print('ERROR: only one of -createlib, -impute or -decodeped should be specified\nExiting...')
        sys.exit(2)

    # Temp
    if args.decodeped:
        return

    if args.libped and args.libphase:
        print('ERROR: only one of -libped or -libphase should be specified\nExiting...')
        sys.exit(2)

    if args.genotypes and (args.ped or args.bim):
        print('ERROR: -genotypes and -ped/-bim cannot be specified together\n'
              'Please only use one of -genotypes or -ped\nExiting...')
        sys.exit(2)

    if args.libped and args.genotypes:
        print('ERROR: -libped and -genotypes cannot be specified together\n'
              'Please use a consistent format, either -libped and -ped or -libphase and -genotypes\nExiting...')
        sys.exit(2)
    elif args.libphase and args.ped:
        print('ERROR: -libphase and -ped cannot be specified together\n'
              'Please use a consistent format, either -libped and -ped or -libphase and -genotypes\nExiting...')
        sys.exit(2)

    if args.impute:
        print('Imputing genotypes with haplotype library...\n')
        # Check requirements
        if args.libped is None and args.libphase is None:
            print('ERROR: no haplotype library supplied. Please supply a library with the -lib option\nExiting...')
            sys.exit(2)
    if args.createlib:
        if args.libped or args.libphase:  # This would be best moved to the create library function
            print(f'Updating haplotype library FILENAME .phase or .ped/.bim from genotypes...\n')
        else:
            print('Creating haplotype library from genotypes...\n')
        if args.founders:
            print('ERROR: -founders can only be specified with -impute, not -createlib\nExiting...')
            sys.exit(2)


def read_in_founders(filename, pedigree, haplotype_library):
    """Read in founders from file modifying the supplied pedigree with that information"""
    print(f'Reading founders from {filename}...')
    focal_ids = []
    founder_ids = set()
    with open(filename) as f:
        lines = f.readlines()

    for split_line in [line.split() for line in lines]:
        identifier = split_line[0]
        try:
            individual = pedigree[identifier]
            focal_ids += [individual.idx]
        except KeyError:
            print(f'ERROR: individual {identifier} not found in genotype '
                  f'(files read in with -ped or -genotypes)\nExiting...')
            sys.exit(2)
        founder_ids.update(split_line[1:])
        individual.founders = split_line[1:]
        if len(individual.founders) == 0:
            print(f'ERROR: individual {individual.idx} has no founders\nExiting...')
            sys.exit(2)

    # Remove non focal ids from pedigree - a bit ugly but simplifies a lot of subsequent code
    ids_to_remove = set(i.idx for i in pedigree) - set(focal_ids)
    for idx in ids_to_remove:
        del pedigree.individuals[idx]
    pedigree.setUpGenerations()
    print(f'Read in {len(focal_ids)} focal individuals, which will be imputed')

    # Check all founders are in the haplotype library
    diff = founder_ids - set(haplotype_library._identifiers)
    if len(diff) > 0:
        print(f'ERROR: focal individuals: {list(diff)}\n'
              f'not found in haplotype library\nExiting...')
        sys.exit(2)

@profile
def main():
    """Main execution code"""

    # Print boilerplate text
    name = 'AlphaPlantImpute2'
    InputOutput.print_boilerplate(name, __version__)

    # Parse command-line arguments
    args = getargs()

    # Argument checks
    check_arguments_consistent(args)

    # Set random seed
    set_seed(args)

    if args.decodeped:
        true = Pedigree.Pedigree(constructor=Pedigree.PlantImputeIndividual)
        true.readInPed('true.ped', args.startsnp, args.stopsnp, haps=False, get_coding=True)
        true.writeGenotypes('true.recoded')
        masked = Pedigree.Pedigree(constructor=Pedigree.PlantImputeIndividual)
        masked.allele_coding = true.allele_coding
        masked.readInPed('masked.ped', args.startsnp, args.stopsnp, haps=False, get_coding=False)
        masked.writeGenotypes('masked.recoded')
        imputed = Pedigree.Pedigree(constructor=Pedigree.PlantImputeIndividual)
        imputed.allele_coding = true.allele_coding
        imputed.readInPed('imputed.ped', args.startsnp, args.stopsnp, haps=False, get_coding=False)
        imputed.writeGenotypes('imputed.recoded')
        sys.exit(2)

    # Are we using PLINK plain text or AlphaImpute format
    args.file_format = None  # Monkey patch args - makes it easy to pass to functions
    if args.genotypes or args.libphase:
        args.file_format = 'AlphaImpute'
    elif args.ped or args.libped:
        args.file_format = 'PLINK'

    # Read library if one has been provided
    haplotype_library, library_coding = None, None
    if args.libped or args.libphase:
        haplotype_library, library_coding = read_library(args)

    # Read genotypes
    genotypes = read_genotypes(args, library_coding)
    n_loci = genotypes.nLoci

    # Check library and genotypes alleles are consistent

    # Check library and genotypes allele codings are the same
    if genotypes.allele_coding is not None and library_coding is not None:
        if not np.alltrue(genotypes.allele_coding == library_coding):
            print(genotypes.allele_coding[:, :10])
            print(library_coding)
            print((genotypes.allele_coding != library_coding).sum())
            print('WARNING: Library and genotype allele codings are different\n'
                  'WARNING: The software will not work as expected if the input files are inconsistently coded\n'
                  '         It is recommended to only provide a coding for the haplotype library')

    # Read founders if supplied
    # sample_target not required for this - just use all founders; undo modifications to sample_target
    focal_ids = None
    if args.founders:
         read_in_founders(args.founders, genotypes, haplotype_library)

    # Probabilistic rates
    error_rate = np.full(n_loci, args.error, dtype=np.float32)
    recombination_rate = np.full(n_loci, args.recomb/n_loci, dtype=np.float32)

    # Choose hidden Markov model
    if args.haploid:
        model = CombinedHMM.HaploidMarkovModel(n_loci, error_rate, recombination_rate)
    elif args.joint:
        model = CombinedHMM.JointMarkovModel(n_loci, error_rate, recombination_rate)
    else:
        model = CombinedHMM.DiploidMarkovModel(n_loci, error_rate, recombination_rate)

    # Create or update haplotype library from high-density genotypes
    if args.createlib:
        create_library(args, model, genotypes, haplotype_library)

    # Imputation
    if args.impute:
        imputation(args, model, genotypes, haplotype_library, library_coding) # just set genotypes allele coding?

if __name__ == "__main__":
    main()
