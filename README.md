# neat-genreads
NEAT-genReads is a fine-grained read simulator. GenReads simulates real-looking data using models learned from specific datasets. There are several supporting utilities for generating models used for simulation.

This is an in-progress v2.0 of the software. For a previous stable release please see: [genReads1](https://github.com/zstephens/genReads1)


## Requirements

* Python 2.7
* Numpy 1.9.1+

## Usage
Here's the simplest invocation of genReads using default parameters. This command produces a single ended fastq file with reads of length 100, ploidy 2, coverage 10X, using the default sequencing substitution, GC% bias, and mutation rate models.

```
python genReads.py -r ref.fa -R 100 -o simulated_data
``` 

The most commonly added options are --pe, --bam, --vcf, and -c. 


Option           |  Description
------           |:----------
-h, --help       |  Displays usage information
-r <str>         |  Reference sequence file in fasta format. Required.
-R <int>         |  Read length. Required. 
-o <str>         |  Output prefix. Use this option to specify where and what to call output files. Required
-c <float>       |  Average coverage across the entire dataset. Default: 10
-e <str>         |  Sequencing error model.
-E <float>       |  Average sequencing error rate. The sequencing error rate model is rescaled to make this the average value. 
-p <int>         |  ploidy [2]
-t <str>         |  bed file containing targeted regions; default coverage for targeted regions is 98% of -c option; default coverage outside targeted regions is 2% of -c option
-m <str>         |  mutation model directory
-M <float>       |  Average mutation rate. The mutation rate model is rescaled to make this the average value.
-s <str>         |  input sample model
-v <str>         |  Input VCF file. Variants from this VCF will be inserted into the simulated sequence.
--pe <int> <int> |  Paired-end fragment length mean and standard deviation. To produce paired end data, one of --pe or --pe-model must be specified.
--pe-model <str> |  Empirical fragment length distribution. Can be generated using [computeFraglen.py](#computefraglenpy). To produce paired end data, one of --pe or --pe-model must be specified.
--cancer         |  Produce tumor/normal datasets. See [Simulating matched tumor/normal pairs](#simulating-matched-tumornormal-pairs)
-cm <str>        |  Cancer mutation model directory. See [Simulating matched tumor/normal pairs](#simulating-matched-tumornormal-pairs)
-cp <float>      |  Tumor sample purity. Default 0.8. See [Simulating matched tumor/normal pairs](#simulating-matched-tumornormal-pairs)
--gc-model <str> |  Empirical GC coverage bias distribution.  Can be generated using [computeGC.py](#computegcpy)
--job <int> <int>|  Jobs IDs for generating reads in parallel
--bam            |  Output golden BAM file
--vcf            |  Output golden VCF file
--rng <int>      |  rng seed value 
--gz             |  Gzip output FQ and VCF


## Functionality

![Diagram describing the way that genReads simulates datasets](docs/flow3.png "Diagram describing the way that genReads simulates datasets")

NEAT genReads produces simulated sequencing datasets. It creates FASTQ files with reads sampled from a provided reference genome, using sequencing error rates and mutation rates learned from real sequencing data. The strength of genReads lies in the ability for the user to customize many sequencing parameters, produce 'golden', true positive datasets, and produce types of data that other simulators cannot (e.g. tumour/normal data).

Features:

- Simulate single-end and paired-end reads 
- Custom read length
- Can introduce random mutations and/or mutations from a VCF file
  - Supported mutation types include SNPs, indels (of any length), inversions, translocations, duplications
  - Can emulate multi-ploid heterozygosity for SNPs and small indels
- Can simulate targeted sequencing via BED input specifying regions to sample from
- Specify simple fragment length model with mean and standard deviation or an empirically learned fragment distribution using utilities/computeFraglen.py
- Simulates quality scores using either the default model or empirically learned quality scores using utilities/fastq_to_qscoreModel.py
- Introduces sequencing substitution errors using either the default model or empirically learned from utilities/
- Accounts for GC% coverage bias using model learned from utilities/computeGC.py
- Output a VCF file with the 'golden' set of true positive variants. These can be compared to bioinformatics workflow output (includes coverage and allele balance information).
- Output a BAM file with the 'golden' set of aligned reads. These indicate where each read originated and how it should be aligned with the reference
- Create paired tumour/normal datasets using characteristics learned from real tumour data
- Parallelized. Can run multiple "partial" simulations in parallel and merge results
- Low memory footprint. Constant (proportional to the size of the reference sequence)

## Examples

The following commands are examples for common types of data to be generated. The simulation uses a reference genome in fasta format to generate reads of 126 bases with default 10X coverage. Outputs paired fastq files, a BAM file and a VCF file. The random variants inserted into the sequence will be present in the VCF and all of the reads will show their proper alignment in the BAM. Unless specified, the simulator will also insert some "sequencing error" -- random variants in some reads that represents false positive results from sequencing.

### Whole genome simulation
Simulate whole genome dataset with random variants inserted according to the default model. 

```
python genReads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300,30
```

### Targeted region simulation
Simulate a targeted region of a genome, e.g. exome, with -t.

```
python genReads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300,30                 \
        -t hg19_exome.bed
```

### Insert specific variants
Simulate a whole genome dataset with only the variants in the provided VCF file using -v and -M.

```
python genReads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300,30                 \
        -v NA12878.vcf              \
        -M 0
```

### Single end reads
Simulate single-end reads by omitting the --pe option.

```
python genReads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       
```



### Simulating matched tumor/normal pairs
Simulate cancer sequencing data using the --cancer option.

```
python genReads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --cancer
```


# Utilities	
Several scripts are distributed with genReads that are used to generate the models used for simulation.

## computeGC.py

Takes .genomecov files produced by BEDtools genomeCov (with -d option).


## computeFraglen.py

Takes SAM file via stdin:

    ./samtools view toy.bam | python computeFraglen.py

and creates fraglen.p model in working directory.


## FindNucleotideContextOnReference.healthy.pl
This script takes in VCF files and generates variant frequency models for NEAT. Coordinates for each variant are located within the HG19 human reference. The corresponding trinucleotide context around that location on the reference is returned into a new column. 

### Running the Script
The script requires 5 arguments to be entered after the full path to FindNucleotideContextOnReference.healthy.pl

```
1. Full path to Fastahack
2. Full path to Reference Genome
3. Full path to input VCF
4. Full path to output file
5. Full path to human GFF

```

### Generating Trinucleotide Frequency Matrices
16 variant frequency matrices are generated based on the trinucleotide context of each variant. The script creates a hash of hashes as it runs through the input VCF. A hash is created for each trinucleotide context, when the middle allele is removed (N_N). Each trinucleotide hash contains hashes of the original mutated from allele, which each contain keys of the mutated to alleles. As the script runs through the input VCF, each key is tallied up.

After reading through the input VCF and tallying up the total of each key, in reference to the mutated from allele and trinucleotide context, each value is divided by the total number of it corresponding mutated from allele in the trinucleotide context. This produces a new matrix for each trinucleotide context, where the probabilities of each row add up to 1. Each matrix file is named "ContextN-N.trinuc" corresponding to the trinucleotide context.

### Generating Insertion and Deletion Length Probabilities
As the script is reading through the input VCF, it locates insertions and deletions and returns their lengths using a hash. For example, the insertion length hash has lengths of the insertions as its key and the total number of each length as the values.

Insertion length totals are then divided by the total number of insertions to calculate the probabity of each insertion length occurring in the data set. These probabilities are written to the files "SSM_insLength.prob" and "SSM_delLength.prob".

### Generating Overall Mutation Type Probabilities
The script counts the number of insertions and deletions found in the input VCF, as well as the total number of mutations. Insertion and deletion totals are divided by the total number of mutations to calculate their probability in the data set. These probabilities are written to the file "SSM_overall.prob".

### ANNOVAR VCFs
ANNOVAR output VCF files can be read by the script. Frequencies of each annotation are calculated and written out to the file "annofreq.prob". Frequencies of exon variant consequences are written out to the file "exonic_consequences.prob".

### Creating Variant Frequency per Region BED File from Human GFF/GTF
While reading through a human GTF or GFF file, coordinates of each variant from the input VCF are compared to each region found in the GTF. Regions are defined as the START and END columns of the input GTF. When a variant's location falls within a GTF region, the region key is incremented. Variant frequency is then calculated out of the length of each region. Variants per Region frequencies are then written to an output BED file with the following format:

```
START	END	Variant_Frequency

```

### Computing the germline-tumor Allele Mismatch Information
Each tumor_genotype is split and compared to the reference allele from the references genome. If the original germline allele does not match the reference genome allele, the middle allele of the trinucleotide context is replaced with the original germline allele.

### Controlled Data and Germline-Reference Allele Mismatch Information
ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.




