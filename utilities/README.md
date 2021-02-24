# compute_gc.py

Takes .genomecov files produced by BEDtools genomeCov (with -d option).

```
bedtools genomecov 
	-d				\
	-ibam normal.bam		\
        -g reference.fa             
```

```
python compute_gc.py			\
        -r reference.fasta         	\
        -i genomecov            	\
        -w [sliding window length]  	\
        -o /path/to/output
```

The main function in this file processes the inputs (reference.fasta, genome.cov, window length), and outputs a GC count for the sequence in the form of a pickle file at the location and with the name from the path the user provides with the -o command.



# compute_fraglen.py

Takes SAM or BAM files and uses console commands for processing:

```
python compute_fralgen.py			\
	-i path to sam file		\
	-o path/to/output
```

The main function in this file will save a pickle (.p) in the location and with the name from the path the user provides with the -o command.

**Please be aware that pysam is not usable on Windows, so any BAM file will need to be turned into a SAM file using samtools beforehand.**

To convert a BAM file to a SAM file using samtools, use the following command:

```
samtools view nameof.bam > nameof.sam
```




# genMutModel.py

Takes references genome and TSV file to generate mutation models:

```
python gen_mut_model.py               \
        -r hg19.fa                  \
        -m inputVariants.tsv        \
        -o /home/me/models.p
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.

# genSeqErrorModel.py

Generates sequence error model for genReads.py -e option.

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

# plotMutModel.py

Performs plotting and comparison of mutation models generated from genMutModel.py.

```
python plotMutModel.py                                        \
        -i model1.p [model2.p] [model3.p]...                  \
        -l legend_label1 [legend_label2] [legend_label3]...   \
        -o path/to/pdf_plot_prefix
```

# vcf_compare_OLD.py

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

## Controlled Data and Germline-Reference Allele Mismatch Information
ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.

