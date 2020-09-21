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
    core_parser.add_argument('-out', required=False, type=str, help='The output file prefix.')

    # Input options
    input_parser = parser.add_argument_group('Input Options')
    input_parser.add_argument('-ped', default=None, required=False, type=str, nargs='*',
                              help='A genotypes file in PLINK plain text format (.ped)')
    input_parser.add_argument('-library', default=None, required=False, type=str,
                              help='A haplotype library file in PLINK plain text format (.ped)')
    input_parser.add_argument('-founders', default=None, required=False, type=str,
                              help='A file that gives the founder individuals for each individual.')
    InputOutput.add_arguments_from_dictionary(input_parser, InputOutput.get_input_options(), options=['startsnp', 'stopsnp', 'seed', 'genotypes'])
    input_parser.add_argument('-libphase', default=None, required=False, type=str,
                              help='A haplotype library file in AlphaGenes phase format (.phase)')
    input_parser.add_argument('-decode', default=None, required=False, type=str, nargs='*',
                              help='Decode .ped files to AlphaGenes format.') # remove?

    # Algorithm options
    algorithm_parser = parser.add_argument_group('Algorithm Options')
    algorithm_parser.add_argument('-hd_threshold', default=0.0, required=False, type=float,
                                  help='Fraction of non-missing markers required to classify an individual as high-density. '
                                  'Only high-density individuals are used to build the haplotype library. Default: 0.0.')
    algorithm_parser.add_argument('-n_haplotypes', default=100, required=False, type=int,
                                  help='Number of haplotypes to sample from the haplotype library. Default: 100.')
    algorithm_parser.add_argument('-haploid', action='store_true', required=False,
                                  help='Run using a haploid HMM instead of the default diploid HMM.')
    algorithm_parser.add_argument('-joint', action='store_true', required=False,
                                  help='Run using a joint HMM instead of the default diploid HMM.')
    algorithm_parser.add_argument('-calling_threshold', default=0.1, required=False, type=float,
                                  help='Genotype calling threshold. '
                                  'Controls whether uncertain (imputed) loci are marked as missing. Default: 0.1.')
    algorithm_parser.add_argument('-n_sample_rounds', default=10, required=False, type=int,
                                  help='Number of rounds of library refinement. Default: 10.')
    algorithm_parser.add_argument('-n_windows', default=5, required=False, type=int,
                                  help='Number of windows for targeted haplotype sampling. Default: 5.')
    InputOutput.add_arguments_from_dictionary(algorithm_parser, InputOutput.get_probability_options(),
                                              options=['error', 'recombination'])
    algorithm_parser.add_argument('-incorrect_loci', action='store_true', required=False,
                                  help='When building a haplotype library, print the average number of loci '
                                  'that were incorrectly sampled from the hidden Markov model.')
    algorithm_parser.add_argument('-overwrite', action='store_true', required=False,
                                  help='Overwrite imputed genotypes with original values at non-missing loci. Default: False')

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
            print(i)
            print(genotype[i])
            print(genotype.dtype)
            print(genotype)
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


def correct_genotypes(true, imputed):
    """Set imputed genotypes to true genotypes at non-missing loci"""
    non_missing = true != 9
    n_corrected = (imputed[non_missing] != true[non_missing]).sum()
    imputed[non_missing] = true[non_missing]
    return n_corrected


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
        if args.incorrect_loci:
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
    print(f'Imputing {len(pedigree)} individuals...')

    # Iterate over all individuals in the Pedigree()
    individuals = pedigree

    # Store input genotypes for later comparison with imputed ones
    if args.overwrite:
        for individual in individuals:
            individual.input_genotypes = individual.genotypes.copy()

    # Sample the haplotype library for each individual
    if args.founders:
        haplotype_library_sample = (haplotype_library.sample_only_identifiers(individual.founders)
                                    for individual in individuals)
    else:
        haplotype_library_sample = (haplotype_library.sample_targeted(args.n_haplotypes,
                                                                      individual.genotypes,
                                                                      args.n_windows)
                                    for individual in individuals)

    # Arguments to pass to get_dosages() via map() and ThreadPoolExecutor.map()
    function_args = (repeat(model), individuals, haplotype_library_sample, repeat(args.calling_threshold))

    # Handle single and multithreaded
    if args.maxthreads == 1:
        # Single threaded
        results = map(get_dosages, *function_args)
    else:
        # Multithreaded
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.maxthreads) as executor:
            results = executor.map(get_dosages, *function_args)

    # Run the imputation - get_dosages is not called until `results` is enumerated
    for _ in results:
        pass

    # If requested, correct imputed loci at non-missing markers
    n_corrected = 0
    if args.overwrite:
        for individual in individuals:
            n_corrected += correct_genotypes(individual.input_genotypes, individual.genotypes)
        print(f'  average number of corrected loci: {n_corrected/len(individuals):.3g}')


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
    assert args.library or args.libphase
    filename = args.library if args.library else args.libphase
    print(f'Reading haplotype library from: {filename}')
    library = Pedigree.Pedigree()
    if args.library:
        library.readInPed(args.library, args.startsnp, args.stopsnp, haps=True, update_coding=True)
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
            individual.haplotypes = haplotype  # first haplotype
        else:
            individual.haplotypes = np.vstack([individual.haplotypes, haplotype])  # second haplotype

    # Write out
    if args.file_format == 'PLINK':
        print(f'Writing haplotype library to {args.out}.ped')
        pedigree.writePhasePed(f'{args.out}.ped')
    else:
        print(f'Writing haplotype library to {args.out}.phase')
        pedigree.writePhase(f'{args.out}.phase')


def imputation(args, model, genotypes, library):
    """Handle imputation tasks"""
    # Impute
    impute_individuals(model, args, genotypes, library)

    # Write out
    print(f'Writing genotype dosages to {args.out}.dosages')
    genotypes.writeDosages(f'{args.out}.dosages')
    if args.file_format == 'PLINK':
        print(f'Writing genotypes to {args.out}.ped')
        genotypes.writeGenotypesPed(f'{args.out}.ped')
    else:
        print(f'Writing genotypes to {args.out}.genotypes')
        genotypes.writeGenotypes(f'{args.out}.genotypes')


def create_library(args, model, genotypes, haplotype_library):
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
    if args.library or args.libphase:
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
            # Get the HD individuals to ignore duplicates
            hd_individuals = [individual for individual in hd_individuals
                             if individual.idx not in duplicate_identifiers]

        haplotype_library.unfreeze()
    else:
        # Create empty haplotype library
        haplotype_library = HaplotypeLibrary.HaplotypeLibrary(n_loci=n_loci)

    # Append (randomly phased) haplotypes to either an empty library or the library being updated
    append_random_haplotypes(haplotype_library, args, hd_individuals, genotypes.maf)
    # Phase the newly added haplotypes
    refine_library(model, args, hd_individuals, haplotype_library, genotypes.maf)
    # Output
    write_library(args, haplotype_library, genotypes.allele_coding)


def read_genotypes(args, allele_coding=None):
    """Read genotypes from file(s)"""
    print('Reading genotypes...')
    genotypes = Pedigree.Pedigree(constructor=Pedigree.PlantImputeIndividual)
    genotypes.allele_coding = allele_coding

    InputOutput.readInPedigreeFromInputs(genotypes, args, genotypes=True, haps=False, reads=False, update_coding=True)
    if len(genotypes) == 0:
        print('ERROR: no genotypes supplied. Please supply them with the -ped or -genotypes options\nExiting...')
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
    decode = args.decode is not None
    n = args.createlib + args.impute + decode
    if n == 0:
        print('ERROR: one of -createlib, -impute or -decode needs to be specified\nExiting...')
        sys.exit(2)
    elif n > 1:
        print('ERROR: only one of -createlib, -impute or -decode should be specified\nExiting...')
        sys.exit(2)

    # If -decode is chosen, return as no more checks are required
    if args.decode:
        print('Decoding PLINK .ped files to AlphaGenes format...\n')
        return

    if not args.out:
        print('ERROR: missing output prefix. Please specify one with -out\nExiting...')
        sys.exit(2)

    if args.library and args.libphase:
        print('ERROR: only one of -library or -libphase should be specified\nExiting...')
        sys.exit(2)

    if args.genotypes and args.ped:
        print('ERROR: -genotypes and -ped cannot be specified together\n'
              'Please only use one of -genotypes or -ped\nExiting...')
        sys.exit(2)

    if args.library and args.genotypes:
        print('ERROR: -library and -genotypes cannot be specified together\n'
              'Please use a consistent format, either -library and -ped or -libphase and -genotypes\nExiting...')
        sys.exit(2)
    elif args.libphase and args.ped:
        print('ERROR: -libphase and -ped cannot be specified together\n'
              'Please use a consistent format, either -library and -ped or -libphase and -genotypes\nExiting...')
        sys.exit(2)

    if args.impute:
        print('Imputing genotypes with haplotype library...\n')
        # Check requirements
        if args.library is None and args.libphase is None:
            print('ERROR: no haplotype library supplied. Please supply a library with the -lib option\nExiting...')
            sys.exit(2)
    if args.createlib:

        if args.library or args.libphase:
            libfile = args.library if args.library else args.libphase
            print(f'Updating haplotype library {libfile} from genotypes...\n')
        else:
            print('Creating haplotype library from genotypes...\n')
        if args.founders:
            print('ERROR: -founders can only be specified with -impute, not -createlib\nExiting...')
            sys.exit(2)


def read_founders(filename, pedigree, haplotype_library):
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
        except KeyError:
            print(f'WARNING: individual {identifier} not found in genotype '
                  f'files (as read in with -ped or -genotypes) and will be ignored')
            continue
        founders = split_line[1:]
        if len(founders) == 0 or set(founders) == {'0'}:
            print(f'WARNING: individual {individual.idx} has no founders and will be ignored')
        else:
            focal_ids += [individual.idx]
            founder_ids.update(split_line[1:])
            individual.founders = split_line[1:]

    # Remove non focal ids from pedigree - a bit ugly but simplifies a lot of subsequent code
    ids_to_remove = set(i.idx for i in pedigree) - set(focal_ids)
    for idx in ids_to_remove:
        del pedigree.individuals[idx]
    pedigree.setUpGenerations()
    print(f'Read in {len(focal_ids)} focal individuals, which will be imputed')

    # Check all founders are in the haplotype library
    #founder_ids = founder_ids - {'0'}  # remove any unknown founders (coded as '0')
    diff = founder_ids - set(haplotype_library._identifiers)
    if len(diff) > 0:
        print(f'ERROR: focal individuals: {list(diff)}\n'
              f'not found in haplotype library\nExiting...')
        sys.exit(2)


def decode(args):
    """Read in a number of .ped files and decode to AlphaGenes format (.genotypes)
    The allele coding is read (interpreted) from the first .ped file"""
    print(f"Decoding: {', '.join(args.decode)}")
    allele_coding = None
    for pedfile in args.decode:
        basename = pedfile.split('.')[0]
        pedigree = Pedigree.Pedigree(constructor=Pedigree.PlantImputeIndividual)
        pedigree.allele_coding = allele_coding
        print(f'Reading and decoding {pedfile}...')
        pedigree.readInPed(pedfile, args.startsnp, args.stopsnp, haps=True, update_coding=True)
        print(f'Writing decoded file to {basename}.genotypes')
        pedigree.writeGenotypes(f'{basename}.genotypes')

        # Use previously read allele coding for each subsequent file
        allele_coding = pedigree.allele_coding


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

    # PLINK .ped decode to AlphaGenes format
    if args.decode:
        decode(args)
        sys.exit(0)

    # Are we using PLINK plain text or AlphaGenes format
    args.file_format = None  # Monkey patch args - makes it easy to pass to functions
    if args.genotypes or args.libphase:
        args.file_format = 'AlphaGenes'
    elif args.ped or args.library:
        args.file_format = 'PLINK'

    # Read library if one has been provided
    haplotype_library, library_coding = None, None
    if args.library or args.libphase:
        haplotype_library, library_coding = read_library(args)

    # Read genotypes (using library coding if provided)
    genotypes = read_genotypes(args, library_coding)
    n_loci = genotypes.nLoci

    # Read founders if supplied
    focal_ids = None
    if args.founders:
         read_founders(args.founders, genotypes, haplotype_library)

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
        imputation(args, model, genotypes, haplotype_library)

if __name__ == "__main__":
    main()
