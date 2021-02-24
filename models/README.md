# models
Used by genReads.py for simulating various characteristics of NGS datasets.


## mutation models

* Used via '-m' parameter

**MutModel_NA12878.p** - Mimic mutations via statistics derived from NA12878 germline variant calls.
(Note that due to filesize restrictions, this one will need to be unzipped first)
**MutModel_BRCA_US_ICGC.p** - Aggregate breast cancer mutation statistics from deidentified ICGC data.
**MutModel_CLLE-ES_ICGC.p** - Aggregate leukemia mutation statistics from deidentified ICGC data.
**MutModel_SKCM-US_ICGC.p** - Aggregate melanoma mutation statistics from deidentified ICGC data.


## sequencing error models

* Used via '-e' parameter

**errorModel_toy.p** - Sequencing error statistics derived from in-house NGS data.


## paired-end fragment length distribution model

* Used via '--pe-model' parameter

**fraglenModel_toy.p** - Fragment length statistics derived from in-house NGS data.


## GC% coverage bias model

* Used via '--gc-model' parameter

**gcBias_toy.p** - GC% coverage bias statistics derived from in-house NGS data.
