.. AlphaPlantImpute2 documentation master file, created by
   sphinx-quickstart on Thu Nov  7 14:07:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|program|
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. highlight:: none


Introduction
~~~~~~~~~~~~


|program| is a genotype imputation algorithm for haploid and diploid plants. |program| supports individuals genotyped with SNP array data. 

Please report any issues to `alphagenes@roslin.ed.ac.uk <alphagenes@roslin.ed.ac.uk>`_ or `steve.thorn@roslin.ed.ac.uk <steve.thorn@roslin.ed.ac.uk>`_.

Conditions of use
-----------------
|program| is part of a suite of software that the AlphaGenes group has developed. It is currently only availble for use by KWS SAAT SE & Co. KGaA. All other interested parties are requested to contact John Hickey (`John.Hickey@roslin.ed.ac.uk <John.Hickey@roslin.ed.ac.uk>`_) to discuss the terms of use.

Authorship
----------

|program| is part of a body of imputation software developed by the AlphaGenes group under John Hickey. It is based on the hidden Markov model (HMM) developed by Roberto Antolín in AlphaImpute (Antolín et al. 2017), and is similar to the algorithm used in MaCH (Li et al. 2011). |program| was written by Steve Thorn and Andrew Whalen, and is currently being supported by Steve Thorn.

Disclaimer
----------

While every effort has been made to ensure that |program| does what it claims to do, there is absolutely no guarantee that the results provided are correct. Use of |program| is at your own risk.


Program options
~~~~~~~~~~~~~~~~~~~~~~~~~~

|program| takes a number of command line arguments to control its behaviour. To display a list of arguments, use either ``AlphaPlantImpute2`` or ``AlphaPlantImpute2 -h``. 


Core options 
--------------
::
  
      -out OUT              The output file prefix.

The ``-out`` argument specifies the output file prefix. For example, with ``-out prefix``, |program| outputs the imputed genotypes to ``prefix.genotypes`` and genotype dosages to ``prefix.dosages``.


Input options 
----------------
::

      -genotypes [GENOTYPES [GENOTYPES ...]]
                            A file in AlphaGenes format.
      -pedigree [PEDIGREE [PEDIGREE ...]]
                            A pedigree file in AlphaGenes format.
      -startsnp STARTSNP    The first marker to consider. The first marker in the
                            file is marker "1".
      -stopsnp STOPSNP      The last marker to consider.
      -seed SEED            A random seed to use for debugging.
      
|program| requires one or more genotype files and, optionally, a pedigree file that specifies whether individuals are inbred/double haploid or outbred.

|program| supports genotype files in the AlphaGenes format as specified by the ``-genotypes`` option. A pedigree file can be supplied using the ``-pedigree`` option. 

 
Algorithm options 
------------------------
::

      -hd_threshold HD_THRESHOLD
                            Fraction of non-missing markers required to classify
                            an individual as high-density. Only high-density
                            individuals are used to build the haplotype library.
                            Default: 0.9.
      -n_haplotypes N_HAPLOTYPES
                            Number of haplotypes to sample from the haplotype
                            library. Default: 100.
      -n_sample_rounds N_SAMPLE_ROUNDS
                            Number of rounds of library refinement. Default: 10.
      -n_impute_rounds N_IMPUTE_ROUNDS
                            Number of rounds of imputation. Default: 1.
      -n_bins N_BINS        Number of bins for targeted haplotype sampling.
                            Default: 5.
      -error ERROR          Genotyping error rate. Default: 0.01.
      -recomb RECOMB        Recombination rate per chromosome. Default: 1.


Haplotype library options
-------------------------
Note: haplotype library functionality is under development. These options are likely to change considerably in future versions.
::

      -library LIBRARY      A haplotype library file in [TBD] format. Read the
                            haplotype library from file rather than building it
                            from high-density individuals.


Multithreading options
------------------------
::

      -maxthreads MAXTHREADS
                            Maximum number of threads to use for analysis.
                            Default: 1.
      -iothreads IOTHREADS  Number of threads to use for input and output.
                            Default: 1.

|program| supports multithreading for running the analysis, and for reading in and writing out (large amounts of) data. The parameter ``-maxthreads`` controls the number of threads used by the analysis. The parameter ``-iothreads`` controls the number of threads used for input and output — setting this option to a value greater than 1 is only recommended for very large files (>10,000 individuals).


Algorithm summary
~~~~~~~~~~~~~~~~~

|program| performs the following steps:  

Create a haplotype library
    Pairs of haplotypes are generated from the genotypes of high-density individuals. High-density individuals are those with a fraction of non-missing markers greater than a given threshold (``-hd_threshold``). Homozygous loci are de-facto phased; heterozygous loci are randomly assigned with equal probability and missing loci are randomly assigned according to minor allele frequency.

Refine the haplotype library
    Each haplotype in the library is refined using a hidden Markov model. The hidden states (at each locus) are either haplotypes in the library (for inbred/double haploid individuals), or pairs of haplotypes (for outbred, diploid individuals). The model randomly generates a haplotype (inbred/double haploid) or pair of haplotypes (outbred, diploid) according to the HMM probabilities. The number of haplotypes considered by the HMM is reduced by randomly sampling a number of haplotypes (``-n_haplotypes``). This number is a trade-off between higher imputation accuracy (higher numbers of haplotypes) and faster computation time (lower numbers). Once each haplotype has been generated, the process iterates for a number of rounds (``-n_sample_rounds``) to refine the haplotypes. A small number of rounds (e.g. 10) is usually being sufficient to accurately estimate the haplotypes. This step, in effect, phases the high-density individuals. 

Impute individuals
    Each individual's genotype is imputed using the same hidden Markov Model with haplotypes from the refined library as hidden states. The number of haplotypes considered (``-n_haplotypes``) is reduced by targeted sampling of haplotypes from the library. Targeted sampling chooses haplotypes that have the fewest opposite homozygous markers compared to the focal individual within a bin of markers. The number of bins the chromosome is divided into is specified by ``-n_bins``. Haplotypes are randomly sampled in the case where many haplotypes have the same number of opposite homozygous markers. The imputation process repeats for a number of rounds (``-n_impute_rounds``), updating the average genotype dosages each time. Finally, the imputed genotypes are calculated by taking the integer-rounded values of the dosages.

Phase individuals
    Each individual is phased using the same hidden Markov Model and the Viterbi algorithm. The number of haplotypes considered (``-n_haplotypes``) is reduced by targeted sampling, but only one round is necessary as the Viterbi algorithm finds the *most likely* phased haplotype. Any inconsistencies between imputed genotype and imputed phase are corrected by modifying the haplotypes.

Input file formats
~~~~~~~~~~~~~~~~~~

Genotype file 
-------------

Genotype files contain the input genotypes for each individual. The first value in each line is the individual's identifier. The remaining values are the genotypes of the individual at each locus, either 0, 1, or 2 (or 9 if missing). The following example gives the genotypes for four individuals genotyped on four markers each.

Example: ::

  id1 0 2 9 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0

Pedigree file
-------------

The pedigree file is only used to specify which individuals are inbred/double haploid or outbred. Each line of the pedigree file has three or four values: the individual's identifier, their father's identifier and their mother's identifier (both '0' in this case, where '0' represents an unknown identifier). The, optional, fourth column specifies whether the individual is inbred or a double haploid ('inbred' or 'dh', case insensitive) or not ('outbred', case insensitive; or an empty field). The algorithm treats inbred and double haploid individuals identically — any residual heterozygous markers in the genotype file are set to missing.

The following example shows four individuals, two double haploids and two outbred individuals.

Example: ::

  id1 0 0 DH
  id2 0 0
  id3 0 0 DH
  id4 0 0
 

Output file formats
~~~~~~~~~~~~~~~~~~~

Genotype file
---------------

The genotype file gives the imputed genotypes (either 0, 1, or 2) for each individual at each locus.

Example: ::

  id1 0 2 2 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0
  
Dosage file
-----------

The dosage file gives the expected allele dosage for the alternative (or minor) allele for each individual. The first value in each line is the individual's identifier. The remaining values are the allele dosages at each locus. If there is a lot of uncertainty in the underlying genotype calls, this file can often prove more accurate than the called genotype values.

Example: ::

  id1    0.0503    2.0000    1.7890    0.0001
  id2    1.0000    1.0000    1.0000    1.0000
  id3    2.0000    0.0000    2.0000    0.0001
  id4    0.0000    2.0000    1.0000    0.0000

Phase file
---------------

The phase file gives the imputed phases (either 0, or 1) for each individual's haplotypes at each locus.

Example: ::

  id1 0 1 1 0 
  id1 0 1 1 0
  id2 0 0 0 1 
  id2 1 1 1 0


.. |program| replace:: AlphaPlantImpute2