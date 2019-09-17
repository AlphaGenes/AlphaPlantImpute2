"""Simple haploid HMM implementation"""

import argparse
import concurrent.futures
from itertools import repeat
from numba import jit
import numpy as np
from .tinyhouse import BasicHMM, InputOutput, Pedigree, HaplotypeLibrary

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

    pythonhmm_parser = parser.add_argument_group('PythonHMM Arguments')
#    pythonhmm_parser.add_argument('-method', default='dosages', required=False, type=str, choices=['dosages', 'sampler'],
#                                  help="Imputation method.")
    pythonhmm_parser.add_argument('-hdthreshold', default=0.9, required=False, type=float,
                                  help='Threshold for determining which haplotypes to include in the haplotype library.')
    pythonhmm_parser.add_argument('-nsamples', default=10, required=False, type=int,
                                  help="Number of samples to use with the 'sampler' method.")
    pythonhmm_parser.add_argument('-nhaplotypes', default=200, required=False, type=int,
                                  help="Number of haplotypes to use in the haplotype library.")
    pythonhmm_parser.add_argument('-nrounds', default=20, required=False, type=int,
                                  help="Number of rounds of HMM sampling.")

    return InputOutput.parseArgs("pythonhmm", parser)


def sample_individuals(pedigree, n_individuals):
    """Return a randomly sampled number of individuals from a pedigree 
    or equivalently a list of individuals
    Note: should probably go in Pedigree()"""
    if n_individuals > len(pedigree):
        n_individuals = len(pedigree)
    indices = np.random.choice(len(pedigree), size=n_individuals, replace=False)
    return [individual
            for i, individual in enumerate(pedigree)
            if i in indices]


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


def create_haplotype_library(individuals, n_loci, maf):  # can/should this be a member of HaplotypeLibrary() or pedigree()?
    """Create a haplotype library from list of individuals
    The population's minor allele frequency (maf) is used to randomly create alleles at missing loci"""
    
    haplotype_library = HaplotypeLibrary.HaplotypeLibrary2(n_loci=n_loci)
    for individual in individuals:
        paternal_haplotype, maternal_haplotype = generate_haplotypes(individual.genotypes, maf)
        haplotype_library.append(paternal_haplotype, identifier=individual.idx)
        haplotype_library.append(maternal_haplotype, identifier=individual.idx)

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
    Note: haplotype_library should not contain a copy of the individual's haplotypes"""
    # n_loci needs to be a global (or passed in) so we don't have to check with each call to this function
    n_loci = len(genotype)

    # Pass missing haplotypes (all 9) to getDiploidPointEstimates(), so that the genotypes are used directly
    haplotypes = np.full((2,n_loci), 9, dtype=np.int8)

    point_estimate = BasicHMM.getDiploidPointEstimates(genotype, haplotypes[0], haplotypes[1],
                                                       haplotype_library, haplotype_library, error)
    forward_probs = BasicHMM.diploidForward(point_estimate, recombination_rate)
    haplotypes = BasicHMM.diploidSampleHaplotypes(forward_probs, recombination_rate,
                                                  haplotype_library, haplotype_library)

    correct_haplotypes(haplotypes[0], haplotypes[1], genotype, maf)
    
    return haplotypes


def precompute_masks(library, individuals):  # better name and possibly move to HaplotypeLibrary! - needs individuals -> identifiers
    """Calculate haplotype library masks for all individuals
    The masks are a list of row indices into the haplotype library array that select all haplotypes except 
    for the given individual (identifier)
    Returns a generator of 1D NumPy arrays"""
    # Note: ugly use of the 'private' member: HaplotypeLibrary._identifiers()
    return ((library._identifiers != individual.idx) for individual in individuals)


def refine_library(args, individuals, haplotype_library, maf, recombination_rate, error):
    """Refine haplotype library"""

    genotypes = (individual.genotypes for individual in individuals) # generator

    for iteration in range(args.nrounds):
        print('Iteration', iteration)

        # Choose random sample of haplotypes for each iteration
        sample = haplotype_library.sample(args.nhaplotypes) # haplotype_library_sample is a HaplotypeLibrary()

        # Generator of haplotype libraries for ThreadPoolExecutor.map()
        # each subsequent library has the corresponding individual's haplotypes masked out
        masks = precompute_masks(sample, individuals)
        # Both following work, which is faster?
#        haplotype_libraries = (sample._haplotypes[mask] for mask in masks)
        haplotype_libraries = (sample.masked(individual.idx) for individual in individuals)

        # Parallel loop over individuals to sample haplotypes
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:    # NEED to get workers from command line
            results = executor.map(sample_haplotype_pair, genotypes, haplotype_libraries,
                                   repeat(recombination_rate), repeat(error), repeat(maf))

        # Update library
        identifiers = [individual.idx for individual in individuals]
        for i, haplotypes in enumerate(results):
            haplotype_library.update_pair(haplotypes[0], haplotypes[1], identifiers[i])
            
    return haplotype_library # this probably? creates a copy of the library. Need to update it in place


def print_mismatches(individual, paternal_hap, maternal_hap):
    # Use np.min([a,b]) to handle swapped haplotypes
    hap0 = individual.haplotypes[0]
    hap1 = individual.haplotypes[1]
    paternal_mismatches = np.min([np.sum((hap0 - paternal_hap) != 0), np.sum((hap1 - paternal_hap) != 0)])
    maternal_mismatches = np.min([np.sum((hap0 - maternal_hap) != 0), np.sum((hap1 - maternal_hap) != 0)])
    print(f'p, m, g, corr', paternal_mismatches, maternal_mismatches, np.sum(((hap0+hap1) - (paternal_hap+maternal_hap)) != 0),
          np.corrcoef(hap0+hap1, paternal_hap+maternal_hap)[0,1])


def impute_individuals(args, pedigree, haplotype_library_full, recombination_rate, error):
    """Impute all individuals in the pedigree"""
    
    print('Imputation')
    print(haplotype_library_full._haplotypes.shape)
    n_loci = pedigree.nLoci
    paternal_hap = np.full(n_loci, 9, dtype=np.int8)
    maternal_hap = np.full(n_loci, 9, dtype=np.int8)

    for individual in pedigree:
        individual.dosages = np.full(n_loci, 0., dtype=np.float32)
    
    n_iterations = 5
    for iteration in range(n_iterations):
        print('Iteration', iteration)
        
        haplotype_library = haplotype_library_full.sample(args.nhaplotypes)._haplotypes  # UGLY 

        # Does it make sense to impute high density individuals? Yes to fill in missingness, but not if it introduces additional errors 
        for individual in pedigree:

            print(f'  Individual idn,idx: {individual.idn},{individual.idx}')
            # Want to impute actual genotypes, not dosages
            
            paternal_hap[:] = 9
            maternal_hap[:] = 9
            point_estimate = BasicHMM.getDiploidPointEstimates(individual.genotypes, paternal_hap, maternal_hap,  
                                                               haplotype_library, haplotype_library, error)
            #forward_probs = BasicHMM.diploidForward(point_estimate, recombination_rate)   
            #paternal_hap, maternal_hap = BasicHMM.diploidSampleHaplotypes(forward_probs, recombination_rate, 
            #                                                              haplotype_library, haplotype_library)

            total_probs = BasicHMM.diploidForwardBackward(point_estimate, recombination_rate)
            individual.dosages += BasicHMM.getDiploidDosages(total_probs, haplotype_library, haplotype_library)
            
    for individual in pedigree:
        individual.dosages /= n_iterations
        individual.genotypes = np.int8(np.round(individual.dosages))
    

@profile
def main():
    """Main execution code"""

    args = get_args()

    pedigree = Pedigree.Pedigree()
    InputOutput.readInPedigreeFromInputs(pedigree, args, genotypes=True, haps=False, reads=False)

    # Calculate MAF and determine high density individuals
    pedigree.setMaf()
    pedigree.high_density_threshold = args.hdthreshold
    pedigree.set_high_density()
    
    # Various parameters
    n_loci = pedigree.nLoci
    error = np.full(n_loci, 0.01, dtype=np.float32)
    recombination_rate = np.full(n_loci, 1/n_loci, dtype=np.float32)
   
    # Subsample
    
    # Library
    individuals = high_density_individuals(pedigree)
    print('# HD individuals', len(individuals))
    haplotype_library = create_haplotype_library(individuals, n_loci, pedigree.maf)
    haplotype_library = refine_library(args, individuals, haplotype_library, pedigree.maf, recombination_rate, error)

    # Imputation
    impute_individuals(args, pedigree, haplotype_library, recombination_rate, error)

    # Output
    pedigree.writeDosages(args.out + '.dosages')
    pedigree.writeGenotypes(args.out + '.genotypes')
    

if __name__ == "__main__":
    main()