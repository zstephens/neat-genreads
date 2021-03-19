# NEAT moved
- The active development of NEAT will proceed on the NCSA project space. You can find the active repository here: https://github.com/ncsa/NEAT. We appreciate all the help Zach has provided in understanding the code and maintaining it.

# NEAT v3.0
- NEAT gen_reads now runs in Python 3 exclusively. The previous, Python 2 version is stored in the repo as v2.0, but will not be undergoing active development.
- Converted sequence objects to Biopython Sequence objects to take advantage of the Biopython library
- Converted cigar strings to lists. Now there are simply a functions that convert a generic cigar string to a list and vice versa.
- Tried to take better advantage of some Biopython libraries, such as for parsing fastas.
- For now, we've eliminated the "Jobs" option and the merge jobs function. We plan to implement multi-threading instead as a way to speed up NEAT's simulation process.
- Added a basic bacterial wrapper that will simulate multiple generations of bacteria based on an input fasta, mutate them and then produce the fastqs, bams, and vcfs for the resultant bacterial population.


## TODOs for v3.1
NEAT is still undergoing active development with many exciting upgrades planned. We also plan to bring the code up to full production scale and will continue to improve the following features (if you would like to [Contribute](CONTRIBUTING.md))
- Using Python's multithreading libraries, speed up NEAT's gen_reads tool significantly.
- Take advantage of pandas library for reading in bed files and other files.
- Code optimization for all gen_reads files (in source folder)
- Further cleanup to PEP8 standards
- Refactor the code to integrate NEAT's utilities into the package
- Improvements and standardization for the utilities, across the board
- VCF compare has some nice features and output, but is very slow to use. Can we improve this utility?

For improvements, we have a lot of great ideas for general improvements aimed at better simulating bacteria, but we believe this same improvements will have applications in other species as well. 
- Multiploidy - all right this has nothing to do with bacteria specifically, but it is a feature we would like to implement into gen_reads.
- Structural Variants - model large scale structural variants with an eye toward intergenic SVs.
- Transposable Elements - model transposons within the sequence
- Repeat regions - will bring a variety of interesting applications
