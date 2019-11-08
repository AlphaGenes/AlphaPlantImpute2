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


|program| is a genotype imputation and phasing algorithm for diploid plants. |program| supports individuals genotyped with SNP array data. 

Please report any issues to `John.Hickey@roslin.ed.ac.uk <John.Hickey@roslin.ed.ac.uk>`_ or `steve.thorn@roslin.ed.ac.uk <steve.thorn@roslin.ed.ac.uk>`_.

Conditions of use
-----------------
|program| is part of a suite of software that our group has developed. It is currently only availble for use by KWS SAAT SE & Co. KGaA. All other interested commercial organizations are requested to contact John Hickey (`John.Hickey@roslin.ed.ac.uk <John.Hickey@roslin.ed.ac.uk>`_) to discuss the terms of use.

Citation and Authorship
-----------------------

|program| is part of a body of imputation software developed by the AlphaGenes group under Professor John Hickey. It is based on the hidden Markov model (HMM) used in AlphaImpute. |program| was written by Steve Thorn and Andrew Whalen, and is currently being supported by Steve Thorn.

Citation:
A citation could go here.


Disclaimer
----------

While every effort has been made to ensure that |program| does what it claims to do, there is absolutely no guarantee that the results provided are correct. Use of |program| is at your own risk.


Program Options
~~~~~~~~~~~~~~~~~~~~~~~~~~

|program| takes in a number of command line arguments to control the program's behaviour. To view a list of arguments, run |program| without any command line arguments, i.e. ``|program|`` or ``|program| -h``. 


Core Options 
--------------
::
  
    Core Options:
      -out OUT              The output file prefix.

The ``-out`` argument gives the output file prefix for the outputs of |program|. By default, |program| outputs a file with imputed genotypes, ``prefix.genotypes``, phased haplotypes ``prefix.phase``, and genotype dosages ``prefix.dosages``. For more information on which files are created, see "Output Arguments", below.


Input Options 
----------------
::

    Input Options:
      -genotypes [GENOTYPES [GENOTYPES ...]]
                            A file in AlphaGenes format.
      -pedigree [PEDIGREE [PEDIGREE ...]]
                            A pedigree file in AlphaGenes format.
      -startsnp STARTSNP    The first marker to consider. The first marker in the
                            file is marker "1".
      -stopsnp STOPSNP      The last marker to consider.
      -seed SEED            A random seed to use for debugging.
      
|program| requires one or more genotype files and, optionally, a pedigree file to run the analysis.

[How pedigree specifies inbred/dh status]

|program| supports binary plink files, ``-bfile`` and genotype files in the AlphaGenes format, ``-genotypes``. A pedigree file can also be supplied using the ``-pedigree`` option. 


Multithreading Options 
------------------------
::

    Multithreading Options:
      -maxthreads MAXTHREADS
                            Maximum number of threads to use for analysis.
                            Default: 1.
      -iothreads IOTHREADS  Number of threads to use for input and output.
                            Default: 1.

|program| supports multithreading for running the analysis, and for reading in and writing out (large amounts of) data. The parameter ``-maxthreads`` controls the number of threads used by the algorithm. The parameter ``-iothreads`` controls the number of threads used for input and output. Setting this option to a value greater than 1 is only recommended for very large files (i.e. >10,000 individuals).


Algorithm Options 
------------------------
::

    Algorithm Options:
      -hd_threshold HD_THRESHOLD
                            Percentage of non-missing markers to classify an
                            individual as high-density. Only high-density
                            individuals make up haplotype library. Default: 0.9.
      -n_haplotypes N_HAPLOTYPES
                            Number of haplotypes to sample from the haplotype
                            library in each HMM round. Default: 100.
      -n_sample_rounds N_SAMPLE_ROUNDS
                            Number of rounds of library refinement. Default: 10.
      -n_impute_rounds N_IMPUTE_ROUNDS
                            Number of rounds of imputation. Default: 5.

|program| performs the following steps:  

Create haplotype library
    Putative haplotypes are generated from the genotypes of high-density individuals.
    Homozygous loci are de-facto phased; heterozygous loci are randomly assigned with equal probability and missing loci are randomly assigned according to minor allele frequency. [Lots of passive voice]

Refine haplotype library
    Each haplotype in the library is refined using a hidden Markov model. The hidden states (at each locus) are the other haplotypes in the library (excluding the haplotype being considered). The model randomly generates a new haplotype according to the HMM's probabilities considering other haplotypes in the library, but *excluding* the haplotype being considered. The number of haplotypes considered is also reduced by randomly sampling a (user-defined) number of haplotypes. This number is a tradeoff between high imputation accuracy (high numbers) and faster computation time (low numbers). Once one set of new haplotypes has been generated, this process repeats using these newly generated haplotypes. A (user-specified) number of rounds is used to refine the haplotypes. This step, in effect, phases the high-density individuals. A small number of rounds (e.g. 10) is usually sufficient to converge to an accurate solution.

Impute individuals
    Each individual's genotype is imputed using a HMM, which uses all pairs of haplotypes in the refined library as hidden states. Again, the number of haplotypes considered is reduced by randomly sampling and the procsess repeats a (user-defined) number of times updating the average genotype dosages at each step. The imputed genotypes are simply the integer rounded values of the dosages.

This is a tradeoff between high imputation accuracy (high numbers) and faster computation time (low numbers)

Input file formats
~~~~~~~~~~~~~~~~~~

Genotype file 
-------------

Genotype files contain the input genotypes for each individual. The first value in each line is the individual's identifier. The remaining values are the genotypes of the individual at each locus, either 0, 1, or 2 (or 9 if missing). The following examples gives the genotypes for four individuals genotyped on four markers each.

Example: ::

  id1 0 2 9 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0

Pedigree file
-------------

Each line of a pedigree file has three values, the individual's identifier, their father's identifier, and their mother's identifier. '0' represents an unknown identifier.

Example: ::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2


Output file formats
~~~~~~~~~~~~~~~~~~~

Genotype file
-----------

The genotype file gives the imputed genotypes (either 0, 1, or 2) for each individual at each locus.

Example: ::

  id1 0 2 2 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0
  
Phase file
-----------

The phase file gives the phased haplotypes (either 0 or 1) for each [high-density?] individual in two lines.

Example: ::

  id1 0 1 1 0
  id1 0 1 1 0
  id2 1 1 1 0
  id2 0 0 0 1
  id3 1 0 1 0
  id3 1 0 1 0 
  id4 0 1 0 0
  id4 0 1 1 0

Dosage file
-----------

The dosage file gives the expected allele dosage for the alternative (or minor) allele for each individual. The first value in each line is the individual's identifier. The remaining values are the allele dosages at each loci. If there is a lot of uncertainty in the underlying genotype calls, this file can often prove more accurate than the called genotype values. [add some values like 1.5, 0.7 1.9 etc.]

Example: ::

  id1    0.0003    2.0000    2.0000    0.0001
  id2    1.0000    1.0000    1.0000    1.0000
  id3    2.0003    0.0000    2.0000    0.0001
  id4    0.0000    2.0000    1.0000    0.0000


.. |program| replace:: AlphaPlantImpute2