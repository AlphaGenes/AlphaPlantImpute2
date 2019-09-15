"""Simple haploid HMM implementation"""

import argparse
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
    

def refine_library(args, individuals, haplotype_library_full, maf, recombination_rate, error):
    """Refine haplotype library"""
    
    #haplotype_library_updated = np.empty_like(haplotype_library, dtype=np.int8)

    n_loci = haplotype_library_full._n_loci   # too hacky
    paternal_hap = np.full(n_loci, 9, dtype=np.int8)
    maternal_hap = np.full(n_loci, 9, dtype=np.int8)

    for iteration in range(args.nrounds):
        print('Iteration', iteration)

        # Choose random sample for each iteration
        haplotype_library = haplotype_library_full.sample(args.nhaplotypes)
        # haplotype_library is a HaplotypeLibrary()
        
        for individual in individuals:
            print(f'  Individual idn,idx: {individual.idn},{individual.idx}')#, end=', ')

            haplotype_library_masked = haplotype_library.masked(individual.idx) 
            # haplotype_library_masked is a numpy array
            # This is confusing

            # Pass missing haplotypes to getDiploidPointEstimates(), so that the genotypes are used directly
            paternal_hap[:] = 9
            maternal_hap[:] = 9

            # Sample
            point_estimate = BasicHMM.getDiploidPointEstimates(individual.genotypes, paternal_hap, maternal_hap,  
                                                               haplotype_library_masked, haplotype_library_masked, error)
            forward_probs = BasicHMM.diploidForward(point_estimate, recombination_rate)   
            paternal_hap, maternal_hap = BasicHMM.diploidSampleHaplotypes(forward_probs, recombination_rate, 
                                                                          haplotype_library_masked, haplotype_library_masked)

            # Correct paternal and maternal haplotype pair so that they match the true genotype
            # Could do this for the first n iterations (if iteraction < 5)        
            if True: 
                correct_haplotypes(paternal_hap, maternal_hap, individual.genotypes, maf)
            #print_mismatches(individual, paternal_hap, maternal_hap)

            # Save updated haplotypes in the full library
            haplotype_library_full.update_pair(paternal_hap, maternal_hap, individual.idx)
            
    return haplotype_library_full # this probably creates a copy of the library. Need to update it in place


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