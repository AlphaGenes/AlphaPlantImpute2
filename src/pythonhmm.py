"""Simple haploid HMM implementation"""

import argparse
from numba import jit
import numpy as np
from .tinyhouse import BasicHMM, HaplotypeLibrary, InputOutput, Pedigree

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


def create_haplotype_library(individuals, maf):  # can/should this be a member of HaplotypeLibrary() or pedigree()?
    """Create a haplotype library from list of individuals
    The population's minor allele frequency (maf) is used to randomly create alleles at missing loci"""
    
    haplotype_library = HaplotypeLibrary.HaplotypeLibrary()
    for individual in individuals:
        paternal_haplotype, maternal_haplotype = generate_haplotypes(individual.genotypes, maf)
        haplotype_library.append(paternal_haplotype)
        haplotype_library.append(maternal_haplotype)

    # Return a NumPy array - might want to keep as a HaplotypeLibrary() instance once this class has more functionality
    return haplotype_library.asMatrix()


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
    

def refine_library(individuals, haplotype_library, maf, recombination_rate, error, n_iterations=20):
    """Refine haplotype library
    Note: individuals - for the loop and individual.genotypes - just need the genotypes (could live in HaplotypeLibrary()?)"""
    
    assert(2*len(individuals) == haplotype_library.shape[0])
    n_loci = haplotype_library.shape[1]

    haplotype_library_updated = np.empty_like(haplotype_library, dtype=np.int8)
    mask = np.ones(len(haplotype_library), dtype=np.bool)
    paternal_hap = np.full(n_loci, 9, dtype=np.int8)
    maternal_hap = np.full(n_loci, 9, dtype=np.int8)
    
    for iteration in range(n_iterations):
        print('Iteration', iteration)
        for i, individual in enumerate(individuals):

            print(f'  Individual idn,idx: {individual.idn},{individual.idx}')#, end=', ')

            # Mask haplotype library to exclude this individual
            # Masked arrays might make this easier (and/or a 'mask' function in HaplotypeLibrary()) and faster 
            # see: https://stackoverflow.com/questions/7429118/how-do-i-get-all-the-values-from-a-numpy-array-excluding-a-certain-index
            mask[:] = True
            paternal_hap_idx = i*2 + 0  # index into haplotype_library
            maternal_hap_idx = i*2 + 1  # index into haplotype_library 
            mask[paternal_hap_idx] = mask[maternal_hap_idx] = False
            haplotype_library_masked = haplotype_library[mask]

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

            # Save updated haplotypes 
            haplotype_library_updated[paternal_hap_idx] = paternal_hap
            haplotype_library_updated[maternal_hap_idx] = maternal_hap

        # Swap haplotype libraries, so that we use the updated haplotypes in the next iteration
        haplotype_library, haplotype_library_updated = haplotype_library_updated, haplotype_library

    return haplotype_library_updated


def print_mismatches(individual, paternal_hap, maternal_hap):
    # Use np.min([a,b]) to handle swapped haplotypes
    hap0 = individual.haplotypes[0]
    hap1 = individual.haplotypes[1]
    paternal_mismatches = np.min([np.sum((hap0 - paternal_hap) != 0), np.sum((hap1 - paternal_hap) != 0)])
    maternal_mismatches = np.min([np.sum((hap0 - maternal_hap) != 0), np.sum((hap1 - maternal_hap) != 0)])
    print(f'p, m, g, corr', paternal_mismatches, maternal_mismatches, np.sum(((hap0+hap1) - (paternal_hap+maternal_hap)) != 0),
          np.corrcoef(hap0+hap1, paternal_hap+maternal_hap)[0,1])


def impute_individuals(pedigree, haplotype_library, recombination_rate, error):
    """Impute all individuals in the pedigree"""
    
    print('Imputation')
    n_loci = pedigree.nLoci
    paternal_hap = np.full(n_loci, 9, dtype=np.int8)
    maternal_hap = np.full(n_loci, 9, dtype=np.int8)

    for individual in pedigree:
        
        print(f'  Individual idn,idx: {individual.idn},{individual.idx}')
        # Want to impute actual genotypes, not dosages

        paternal_hap[:] = 9
        maternal_hap[:] = 9
        point_estimate = BasicHMM.getDiploidPointEstimates(individual.genotypes, paternal_hap, maternal_hap,  
                                                           haplotype_library, haplotype_library, error)
        forward_probs = BasicHMM.diploidForward(point_estimate, recombination_rate)   
        paternal_hap, maternal_hap = BasicHMM.diploidSampleHaplotypes(forward_probs, recombination_rate, 
                                                                      haplotype_library, haplotype_library)

        total_probs = BasicHMM.diploidForwardBackward(point_estimate, recombination_rate)
        individual.dosages = BasicHMM.getDiploidDosages(total_probs, haplotype_library, haplotype_library)
        individual.genotypes = np.int8(np.round(individual.dosages))

        #imputed_genotype = paternal_hap + maternal_hap
        #true_genotype = individual.haplotypes[0] + individual.haplotypes[1]
        #print(f'Individual {individual.idn}, sample {np.corrcoef(imputed_genotype, true_genotype)[0,1]:.3f}, dosage {np.corrcoef(np.round(individual.dosages), true_genotype)[0,1]:.3f}')


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
    n_haplotypes = args.nhaplotypes
    assert(np.mod(n_haplotypes, 2) == 0) # number of haplotype should be divisible by 2
    individuals = high_density_individuals(pedigree)
    individuals = sample_individuals(individuals, n_haplotypes//2)
    
    # Library
    n_iterations = args.nrounds
    haplotype_library = create_haplotype_library(individuals, pedigree.maf)
    haplotype_library = refine_library(individuals, haplotype_library, pedigree.maf, recombination_rate, error, n_iterations)
    
    # Imputation
    impute_individuals(pedigree, haplotype_library, recombination_rate, error)

    # Output
    pedigree.writeDosages(args.out + '.dosages')
    pedigree.writeGenotypes(args.out + '.genotypes')
    

if __name__ == "__main__":
    main()