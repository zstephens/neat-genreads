# neat-genreads
NEAT-genReads is a fine-grained read simulator. GenReads simulates real-looking data using models learned from specific datasets. There are several supporting utilities for generating models used for simulation.

This is an in-progress v2.0 of the software. For a previous stable release please see: [genReads1](https://github.com/zstephens/genReads1)


Table of Contents
=================

  * [neat-genreads](#neat-genreads)
  * [Table of Contents](#table-of-contents)
    * [Requirements](#requirements)
    * [Usage](#usage)
    * [Functionality](#functionality)
    * [Examples](#examples)
      * [Whole genome simulation](#whole-genome-simulation)
      * [Targeted region simulation](#targeted-region-simulation)
      * [Insert specific variants](#insert-specific-variants)
      * [Single end reads](#single-end-reads)
      * [Parallelizing simulation](#parallelizing-simulation)
  * [Utilities](#utilities)
    * [computeGC.py](#computegcpy)
    * [computeFraglen.py](#computefraglenpy)
    * [genMutModel.py](#genmutmodelpy)
    * [genSeqErrorModel.py](#genseqerrormodelpy)
    * [plotMutModel.py](#plotmutmodelpy)
    * [vcf_compare_OLD.py](#vcf_compare_oldpy)
      * [Note on Sensitive Patient Data](#note-on-sensitive-patient-data)




## Requirements

* Python 2.7
* Numpy 1.9.1+

## Usage
Here's the simplest invocation of genReads using default parameters. This command produces a single ended fastq file with reads of length 101, ploidy 2, coverage 10X, using the default sequencing substitution, GC% bias, and mutation rate models.

```
python genReads.py -r ref.fa -R 101 -o simulated_data
``` 

The most commonly added options are --pe, --bam, --vcf, and -c. 


Option           |  Description
------           |:----------
-h, --help       |  Displays usage information
-r <str>         |  Reference sequence file in fasta format. A reference index (.fai) will be created if one is not found in the directory of the reference as [reference filename].fai. Required. The index can be created using samtools faidx.
-R <int>         |  Read length. Required. 
-o <str>         |  Output prefix. Use this option to specify where and what to call output files. Required
-c <float>       |  Average coverage across the entire dataset. Default: 10
-e <str>         |  Sequencing error model pickle file
-E <float>       |  Average sequencing error rate. The sequencing error rate model is rescaled to make this the average value. 
-p <int>         |  ploidy [2]
-t <str>         |  bed file containing targeted regions; default coverage for targeted regions is 98% of -c option; default coverage outside targeted regions is 2% of -c option
-to <float>      |  off-target coverage scalar [0.02]
-m <str>         |  mutation model pickle file
-M <float>       |  Average mutation rate. The mutation rate model is rescaled to make this the average value. Must be between 0 and 0.3. These random mutations are inserted in addition to the once specified in the -v option.
-s <str>         |  input sample model
-v <str>         |  Input VCF file. Variants from this VCF will be inserted into the simulated sequence with 100% certainty.
--pe <int> <int> |  Paired-end fragment length mean and standard deviation. To produce paired end data, one of --pe or --pe-model must be specified.
--pe-model <str> |  Empirical fragment length distribution. Can be generated using [computeFraglen.py](#computefraglenpy). To produce paired end data, one of --pe or --pe-model must be specified.
--gc-model <str> |  Empirical GC coverage bias distribution.  Can be generated using [computeGC.py](#computegcpy)
--job <int> <int>|  Jobs IDs for generating reads in parallel
--nnr            |  save non-N ref regions (for parallel jobs)
--bam            |  Output golden BAM file
--vcf            |  Output golden VCF file
--rng <int>      |  rng seed value; identical RNG value should produce identical runs of the program, so things like read locations, variant positions, error positions, etc, should all be the same.
--gz             |  Gzip output FQ and VCF
--no-fastq       |  Bypass generation of FASTQ read files


## Functionality

![Diagram describing the way that genReads simulates datasets](docs/flow_new.png "Diagram describing the way that genReads simulates datasets")

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
        --pe 300 30
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
        --pe 300 30                 \
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
        --pe 300 30                 \
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



### Parallelizing simulation
When possible, simulation can be done in parallel via multiple executions with different --job options. The resultant files will then need to be merged using utilities/mergeJobs.py. The following example shows splitting a simulation into 4 separate jobs (which can be run independently):

```
python genReads.py  -r hg19.fa -R 126 -o /home/me/simulated_reads --bam --vcf --job 1 4
python genReads.py  -r hg19.fa -R 126 -o /home/me/simulated_reads --bam --vcf --job 2 4
python genReads.py  -r hg19.fa -R 126 -o /home/me/simulated_reads --bam --vcf --job 3 4
python genReads.py  -r hg19.fa -R 126 -o /home/me/simulated_reads --bam --vcf --job 4 4

python mergeJobs.py -i /home/me/simulated_reads -o /home/me/simulated_reads_merged -s /path/to/samtools
```

In future revisions the dependence on SAMtools will be removed.

To simulate human WGS 50X, try 50 chunks or less.

# Utilities	
Several scripts are distributed with genReads that are used to generate the models used for simulation.

## computeGC.py

Computes GC% coverage bias distribution from sample (bedrolls genomecov) data.
Takes .genomecov files produced by BEDtools genomeCov (with -d option).

```
bedtools genomecov
        -d                          \
        -ibam normal.bam            \
        -g reference.fa
```

```
python computeGC.py                 \
        -r reference.fa             \
        -i genomecovfile            \
        -w [sliding window length]  \
        -o /path/to/model.p
```

## computeFraglen.py

Computes empirical fragment length distribution from sample data.
Takes SAM file via stdin:

    ./samtools view toy.bam | python computeFraglen.py

and creates fraglen.p model in working directory.

## genMutModel.py

Takes references genome and TSV file to generate mutation models:

```
python genMutModel.py               \
        -r hg19.fa                  \
        -m inputVariants.tsv        \
        -o /home/me/models.p
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.

## genSeqErrorModel.py

Generates sequence error model for genReads.py -e option.
This script needs revision, to improve the quality-score model eventually, and to include code to learn sequencing errors from pileup data.

```
python genSeqErrorModel.py                            \
        -i input_read1.fq (.gz) / input_read1.sam     \
        -o output.p                                   \
        -i2 input_read2.fq (.gz) / input_read2.sam    \
        -p input_alignment.pileup                     \
        -q quality score offset [33]                  \
        -Q maximum quality score [41]                 \
        -n maximum number of reads to process [all]   \
        -s number of simulation iterations [1000000]  \
        --plot perform some optional plotting
```
## plotMutModel.py

Performs plotting and comparison of mutation models generated from genMutModel.py.

```
python plotMutModel.py                                        \
        -i model1.p [model2.p] [model3.p]...                  \
        -l legend_label1 [legend_label2] [legend_label3]...   \
        -o path/to/pdf_plot_prefix
```

## vcf_compare_OLD.py

Tool for comparing VCF files.

```
python vcf_compare_OLD.py
        --version          show program's version number and exit      \
        -h, --help         show this help message and exit             \
        -r <ref.fa>        * Reference Fasta                           \
        -g <golden.vcf>    * Golden VCF                                \
        -w <workflow.vcf>  * Workflow VCF                              \
        -o <prefix>        * Output Prefix                             \
        -m <track.bed>     Mappability Track                           \
        -M <int>           Maptrack Min Len                            \
        -t <regions.bed>   Targetted Regions                           \
        -T <int>           Min Region Len                              \
        -c <int>           Coverage Filter Threshold [15]              \
        -a <float>         Allele Freq Filter Threshold [0.3]          \
        --vcf-out          Output Match/FN/FP variants [False]         \
        --no-plot          No plotting [False]                         \
        --incl-homs        Include homozygous ref calls [False]        \
        --incl-fail        Include calls that failed filters [False]   \
        --fast             No equivalent variant detection [False]
```
Mappability track examples: https://github.com/zstephens/neat-repeat/tree/master/example_mappabilityTracks

### Note on Sensitive Patient Data
ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.




