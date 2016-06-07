# computeGC.py

Takes .genomecov files produced by BEDtools genomeCov (with -d option).


# computeFraglen.py

Takes SAM file via stdin:

./samtools view toy.bam | python computeFraglen.py

and creates fraglen.p model in working directory.


# genMutModel.py

Takes references genome and TSV file to generate mutation models:

```
python genMutModel.py               \
        -r hg19.fa                  \
        -m inputVariants.tsv        \
        -o /home/me/models.p
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.

## Controlled Data and Germline-Reference Allele Mismatch Information
ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.

