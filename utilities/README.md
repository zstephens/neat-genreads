# computeGC.py

Takes .genomecov files produced by BEDtools genomeCov (with -d option).


# computeFraglen.py

Takes SAM file via stdin:

./samtools view toy.bam | python computeFraglen.py

and creates fraglen.p model in working directory.


# computeTrinucStats.py

Takes references genome and TSV file to generate mutation models:

```
python computeTrinucStats.py        \
        -r hg19.fa                  \
        -m inputVariants.tsv        \
        -o /home/me/models.p
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.


# FindNucleotideContextOnReference.healthy.pl
This script takes in VCF files and generates variant frequency models for NEAT. Coordinates for each variant are located within the HG19 human reference. The corresponding trinucleotide context around that location on the reference is returned into a new column.

## Running the Script
The script requires 5 arguments to be entered after the full path to FindNucleotideContextOnReference.healthy.pl

```
1. Full path to Fastahack
2. Full path to Reference Genome
3. Full path to input VCF
4. Full path to output file
5. Full path to human GFF
```

## Generating Trinucleotide Frequency Matrices
16 variant frequency matrices are generated based on the trinucleotide context of each variant. The script creates a hash of hashes as it runs through the input VCF. A hash is created for each trinucleotide context, when the middle allele is removed (N_N). Each trinucleotide hash contains hashes of the original mutated from allele, which each contain keys of the mutated to alleles. As the script runs through the input VCF, each key is tallied up.

After reading through the input VCF and tallying up the total of each key, in reference to the mutated from allele and trinucleotide context, each value is divided by the total number of it corresponding mutated from allele in the trinucleotide context. This produces a new matrix for each trinucleotide context, where the probabilities of each row add up to 1. Each matrix file is named "ContextN-N.trinuc" corresponding to the trinucleotide context.

## Generating Insertion and Deletion Length Probabilities
As the script is reading through the input VCF, it locates insertions and deletions and returns their lengths using a hash. For example, the insertion length hash has lengths of the insertions as its key and the total number of each length as the values.

Insertion length totals are then divided by the total number of insertions to calculate the probabity of each insertion length occurring in the data set. These probabilities are written to the files "SSM_insLength.prob" and "SSM_delLength.prob".

## Generating Overall Mutation Type Probabilities
The script counts the number of insertions and deletions found in the input VCF, as well as the total number of mutations. Insertion and deletion totals are divided by the total number of mutations to calculate their probability in the data set. These probabilities are written to the file "SSM_overall.prob".

## ANNOVAR VCFs
ANNOVAR output VCF files can be read by the script. Frequencies of each annotation are calculated and written out to the file "annofreq.prob". Frequencies of exon variant consequences are written out to the file "exonic_consequences.prob".

## Creating Variant Frequency per Region BED File from Human GFF/GTF
While reading through a human GTF or GFF file, coordinates of each variant from the input VCF are compared to each region found in the GTF. Regions are defined as the START and END columns of the input GTF. When a variant's location falls within a GTF region, the region key is incremented. Variant frequency is then calculated out of the length of each region. Variants per Region frequencies are then written to an output BED file with the following format:

```
START   END     Variant_Frequency

```

## Computing the germline-tumor Allele Mismatch Information
Each tumor_genotype is split and compared to the reference allele from the references genome. If the original germline allele does not match the reference genome allele, the middle allele of the trinucleotide context is replaced with the original germline allele.

## Controlled Data and Germline-Reference Allele Mismatch Information
ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.

