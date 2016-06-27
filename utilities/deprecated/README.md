#Deprecated Perl Scripts
These scripts were updated and rewritten in python to improve ease of use and speed. Usage and a quick description of the deprecated scripts can be found below. Please use genMutModel.py to generate mutation models.

##FindNucleotideContextOnReference.pl
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
