#!/usr/bin/env python

# Python 3 ready

import os
import re
import pickle
import argparse
import numpy as np
from Bio import SeqIO
import pandas as pd


#########################################################
#				VARIOUS HELPER FUNCTIONS				#
#########################################################


def cluster_list(list_to_cluster: list, delta: float) -> list:
    """
    Clusters a sorted list
    :param list_to_cluster: a sorted list
    :param delta: the value to compare list items to
    :return: a clustered list of values
    """
    out_list = [[list_to_cluster[0]]]
    previous_value = list_to_cluster[0]
    current_index = 0
    for item in list_to_cluster[1:]:
        if item - previous_value <= delta:
            out_list[current_index].append(item)
        else:
            current_index += 1
            out_list.append([])
            out_list[current_index].append(item)
        previous_value = item
    return out_list


#####################################
#				main()				#
#####################################


def main():
    # Some constants we'll need later
    REF_WHITELIST = [str(n) for n in range(1, 30)] + ['x', 'y', 'X', 'Y', 'mt', 'Mt', 'MT']
    REF_WHITELIST += ['chr' + n for n in REF_WHITELIST]
    VALID_NUCL = ['A', 'C', 'G', 'T']
    VALID_TRINUC = [VALID_NUCL[i] + VALID_NUCL[j] + VALID_NUCL[k] for i in range(len(VALID_NUCL)) for j in
                    range(len(VALID_NUCL)) for k in range(len(VALID_NUCL))]
    # if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
    VCF_DEFAULT_POP_FREQ = 0.00001

    parser = argparse.ArgumentParser(description='gen_mut_model.py')
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-m', type=str, required=True, metavar='/path/to/mutations.vcf',
                        help="Mutation file for organism in VCF format")
    parser.add_argument('-o', type=str, required=True, metavar='/path/to/output/and/prefix',
                        help="Name of output file (final model will append \'.p\')")
    parser.add_argument('-bi', type=str, required=False, metavar='Bed file of regions to include '
                                                                 '(use bedtools complement if you have a '
                                                                 'bed of exclusion areas)', default=None,
                        help="only_use_these_regions.bed")
    parser.add_argument('--save-trinuc', required=False, action='store_true', default=False,
                        help='save trinucleotide counts for reference')
    parser.add_argument('--human-sample', required=False, action='store_true', default=False,
                        help='To skip unnumbered scaffolds in human references')
    parser.add_argument('--skip-common', required=False, action='store_true', default=False,
                        help='Do not save common snps + high mut regions')
    args = parser.parse_args()

    (ref, vcf, out_pickle, save_trinuc, skip_common) = (
        args.r, args.m, args.o, args.save_trinuc, args.skip_common)

    is_human = args.human_sample

    # how many times do we observe each trinucleotide in the reference (and input bed region, if present)?
    TRINUC_REF_COUNT = {}
    TRINUC_BED_COUNT = {}
    # [(trinuc_a, trinuc_b)] = # of times we observed a mutation from trinuc_a into trinuc_b
    TRINUC_TRANSITION_COUNT = {}
    # total count of SNPs
    SNP_COUNT = 0
    # overall SNP transition probabilities
    SNP_TRANSITION_COUNT = {}
    # total count of indels, indexed by length
    INDEL_COUNT = {}
    # tabulate how much non-N reference sequence we've eaten through
    TOTAL_REFLEN = 0
    # detect variants that occur in a significant percentage of the input samples (pos,ref,alt,pop_fraction)
    COMMON_VARIANTS = []
    # tabulate how many unique donors we've encountered (this is useful for identifying common variants)
    TOTAL_DONORS = {}
    # identify regions that have significantly higher local mutation rates than the average
    HIGH_MUT_REGIONS = []

    # Process bed file,
    is_bed = False
    mybed = None
    if args.bi is not None:
        print('Processing bed file...')
        try:
            mybed = pd.read_csv(args.bi, sep='\t', header=None, index_col=None)
            is_bed = True
        except ValueError:
            print('Problem parsing bed file. Ensure bed file is tab separated, standard bed format')
            exit(1)
        mybed = mybed.rename(columns={0: 'chrom', 1: 'start', 2: 'end'})
        # Adding a couple of columns we'll need for later calculations
        mybed['coords'] = list(zip(mybed.start, mybed.end))
        mybed['track_len'] = mybed.end - mybed.start + 1

    # Process reference file
    print('Processing reference...')
    try:
        reference = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    except ValueError:
        print("Problems parsing reference file. Ensure reference is in proper fasta format")
        exit(1)

    # simplify naming and filter out actual human genomes from scaffolding
    ref_dict = {}
    for key in reference.keys():
        key_split = key.split("|")
        if is_human:
            if key_split[0] in REF_WHITELIST:
                ref_dict[key_split[0]] = reference[key]
            else:
                continue
        else:
            ref_dict[key_split[0]] = reference[key]

    ref_list = list(ref_dict.keys())

    # Process VCF file. First check if it's been entered as a TSV
    if vcf[-3:] == 'tsv':
        print("Warning! TSV file must follow VCF specifications.")

    # Pre-parsing to find all the matching chromosomes between ref and vcf
    print('Processing VCF file...')
    try:
        variants = pd.read_csv(vcf, sep='\t', comment='#', index_col=None, header=None)
    except ValueError:
        print("VCF must be in standard VCF format with tab-separated columns")
        exit(1)

    # Narrow chromosomes to those matching the reference
    # This is in part to make sure the names match
    variant_chroms = variants[0].to_list()
    variant_chroms = list(set(variant_chroms))
    matching_chromosomes = []
    for ref_name in ref_list:
        if ref_name not in variant_chroms:
            continue
        else:
            matching_chromosomes.append(ref_name)

    # Check to make sure there are some matches
    if not matching_chromosomes:
        print("Found no chromosomes in common between VCF and Fasta. Please fix the chromosome names and try again")
        exit(1)

    # Double check that there are matches
    try:
        matching_variants = variants[variants[0].isin(matching_chromosomes)]
    except ValueError:
        print("Problem matching variants with reference.")
        exit(1)
    if matching_variants.empty:
        print("There is no overlap between reference and variant file. This could be a chromosome naming problem")
        exit(1)


    # Rename header in dataframe for processing
    matching_variants = matching_variants.rename(columns={0: "CHROM", 1: 'chr_start', 2: 'ID', 3: 'REF', 4: 'ALT',
                                                          5: 'QUAL', 6: 'FILTER', 7: 'INFO'})

    # Change the indexing by -1 to match python format indexing (0-based)
    matching_variants['chr_start'] = matching_variants['chr_start'] - 1
    matching_variants['chr_end'] = matching_variants['chr_start']

    # Process the variant table
    indices_to_indels = \
        matching_variants.loc[matching_variants.ALT.apply(len) != matching_variants.REF.apply(len)].index

    # indels in vcf don't include the preserved first nucleotide, so lets trim the vcf alleles
    ref_values_to_change = matching_variants.loc[indices_to_indels, 'REF'].copy().str[1:]
    alt_values_to_change = matching_variants.loc[indices_to_indels, 'ALT'].copy().str[1:]
    matching_variants.loc[indices_to_indels, "REF"] = ref_values_to_change
    matching_variants.loc[indices_to_indels, "ALT"] = alt_values_to_change
    matching_variants.replace('', '-', inplace=True)

    # If multi-alternate alleles are present, lets just ignore this variant. I may come back and improve this later
    indices_to_ignore = matching_variants[matching_variants['ALT'].str.contains(',')].index
    matching_variants = matching_variants.drop(indices_to_ignore)

    # if we encounter a multi-np (i.e. 3 nucl --> 3 different nucl), let's skip it for now...

    # Alt and Ref contain no dashes
    no_dashes = matching_variants[
        ~matching_variants['REF'].str.contains('-') & ~matching_variants['ALT'].str.contains('-')].index
    # Alt and Ref lengths are greater than 1
    long_variants = matching_variants[
        (matching_variants['REF'].apply(len) > 1) & (matching_variants['ALT'].apply(len) > 1)].index
    complex_variants = list(set(no_dashes) & set(long_variants))
    matching_variants = matching_variants.drop(complex_variants)

    # This is solely to make regex easier later, since we can't predict where in the line a string will be
    new_info = ';' + matching_variants['INFO'].copy() + ';'
    matching_variants['INFO'] = new_info

    # Now we check that the bed and vcf have matching regions
    # This also checks that the vcf and bed have the same naming conventions and cuts out scaffolding.
    if is_bed:
        bed_chroms = list(set(mybed['chrom']))
        matching_bed_keys = list(set(bed_chroms) & set(variant_chroms))
        try:
            matching_bed = mybed[mybed['chrom'].isin(matching_bed_keys)]
        except ValueError:
            print('Problem matching bed chromosomes to variant file.')
            exit(1)
        if matching_bed.empty:
            print("There is no overlap between bed and variant file. "
                  "This could be a chromosome naming problem")
            exit(1)

    # Count Trinucleotides in reference, based on bed or not
    print('Counting trinucleotides in reference...')
    if is_bed:
        print("since you're using a bed input, we have to count trinucs in bed region even if\n"
              "you already have a trinuc count file for the reference...")
        sub_regions = matching_bed['coords'].to_list()
        for sr in sub_regions:
            for ref_name in matching_chromosomes:
                for i in range(sr[0], sr[1] - 2):
                    trinuc = str(ref_dict[ref_name][i:i + 3].seq)
                    # skip if trinuc contains invalid characters
                    if not trinuc in VALID_TRINUC:
                        continue
                    if trinuc not in TRINUC_REF_COUNT:
                        TRINUC_REF_COUNT[trinuc] = 0
                    TRINUC_REF_COUNT[trinuc] += 1

    elif not os.path.isfile(ref + '.trinucCounts'):
        for ref_name in matching_chromosomes:
            for i in range(len(ref_dict[ref_name]) - 2):
                trinuc = str(ref_dict[ref_name][i:i + 3].seq)
                # skip if trinuc contains invalid characters
                if not trinuc in VALID_TRINUC:
                    continue
                if trinuc not in TRINUC_REF_COUNT:
                    TRINUC_REF_COUNT[trinuc] = 0
                TRINUC_REF_COUNT[trinuc] += 1
    else:
        print('Found trinucCounts file, using that.')

    # Load and process variants in each reference sequence individually, for memory reasons...
    print('Creating mutational model...')
    for ref_name in matching_chromosomes:
        # Count the number of non-N nucleotides for the reference
        TOTAL_REFLEN += len(ref_dict[ref_name].seq) - ref_dict[ref_name].seq.count('N')

        # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
        VDAT_COMMON = []

        # Create a view that narrows variants list to current ref
        variants_to_process = matching_variants[matching_variants["CHROM"] == ref_name].copy()
        bed_to_process = matching_bed[matching_bed['chrom'] == ref_name].copy()
        ref_sequence = str(ref_dict[ref_name].seq)

        # we want only snps
        # so, no '-' characters allowed, and chrStart must be same as chrEnd
        snp_df = variants_to_process[~variants_to_process.index.isin(indices_to_indels)]
        snp_df = snp_df.loc[snp_df['chr_start'] == snp_df['chr_end']]
        if not snp_df.empty:
            # only consider positions where ref allele in vcf matches the nucleotide in our reference
            for index, row in snp_df.iterrows():
                trinuc_to_analyze = str(ref_sequence[row.chr_start - 1: row.chr_start + 2])
                if trinuc_to_analyze not in VALID_TRINUC:
                    continue
                if row.REF == trinuc_to_analyze[1]:
                    trinuc_ref = trinuc_to_analyze
                    trinuc_alt = trinuc_to_analyze[0] + snp_df.loc[index, 'ALT'] + trinuc_to_analyze[2]
                    if trinuc_alt not in VALID_TRINUC:
                        continue
                    key = (trinuc_ref, trinuc_alt)
                    if key not in TRINUC_TRANSITION_COUNT:
                        TRINUC_TRANSITION_COUNT[key] = 0
                    TRINUC_TRANSITION_COUNT[key] += 1
                    SNP_COUNT += 1
                    key2 = (str(row.REF), str(row.ALT))
                    if key2 not in SNP_TRANSITION_COUNT:
                        SNP_TRANSITION_COUNT[key2] = 0
                    SNP_TRANSITION_COUNT[key2] += 1

                    my_pop_freq = VCF_DEFAULT_POP_FREQ
                    if ';CAF=' in snp_df.loc[index, 'INFO']:
                        caf_str = re.findall(r";CAF=.*?(?=;)", row.INFO)[0]
                        if ',' in caf_str:
                            my_pop_freq = float(caf_str[5:].split(',')[1])
                    VDAT_COMMON.append(
                        (row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))
                else:
                    print('\nError: ref allele in variant call does not match reference.\n')
                    exit(1)

        # now let's look for indels...
        indel_df = variants_to_process[variants_to_process.index.isin(indices_to_indels)]
        if not indel_df.empty:
            for index, row in indel_df.iterrows():
                if "-" in row.REF:
                    len_ref = 0
                else:
                    len_ref = len(row.REF)
                if "-" in row.ALT:
                    len_alt = 0
                else:
                    len_alt = len(row.ALT)
                if len_ref != len_alt:
                    indel_len = len_alt - len_ref
                    if indel_len not in INDEL_COUNT:
                        INDEL_COUNT[indel_len] = 0
                    INDEL_COUNT[indel_len] += 1

                    my_pop_freq = VCF_DEFAULT_POP_FREQ
                    if ';CAF=' in row.INFO:
                        caf_str = re.findall(r";CAF=.*?(?=;)", row.INFO)[0]
                        if ',' in caf_str:
                            my_pop_freq = float(caf_str[5:].split(',')[1])
                    VDAT_COMMON.append((row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))

        # if we didn't find anything, skip ahead along to the next reference sequence
        if not len(VDAT_COMMON):
            print('Found no variants for this reference.')
            continue

        # identify common mutations
        percentile_var = 95
        min_value = np.percentile([n[4] for n in VDAT_COMMON], percentile_var)
        for k in sorted(VDAT_COMMON):
            if k[4] >= min_value:
                COMMON_VARIANTS.append((ref_name, k[0], k[1], k[3], k[4]))
        VDAT_COMMON = {(n[0], n[1], n[2], n[3]): n[4] for n in VDAT_COMMON}

        # identify areas that have contained significantly higher random mutation rates
        dist_thresh = 2000
        percentile_clust = 97
        scaler = 1000
        # identify regions with disproportionately more variants in them
        VARIANT_POS = sorted([n[0] for n in VDAT_COMMON.keys()])
        clustered_pos = cluster_list(VARIANT_POS, dist_thresh)
        by_len = [(len(clustered_pos[i]), min(clustered_pos[i]), max(clustered_pos[i]), i) for i in
                 range(len(clustered_pos))]
        # Not sure what this was intended to do or why it is commented out. Leaving it here for now.
        # by_len  = sorted(by_len,reverse=True)
        # minLen = int(np.percentile([n[0] for n in by_len],percentile_clust))
        # by_len  = [n for n in by_len if n[0] >= minLen]
        candidate_regions = []
        for n in by_len:
            bi = int((n[1] - dist_thresh) / float(scaler)) * scaler
            bf = int((n[2] + dist_thresh) / float(scaler)) * scaler
            candidate_regions.append((n[0] / float(bf - bi), max([0, bi]), min([len(ref_dict[ref_name]), bf])))
        minimum_value = np.percentile([n[0] for n in candidate_regions], percentile_clust)
        for n in candidate_regions:
            if n[0] >= minimum_value:
                HIGH_MUT_REGIONS.append((ref_name, n[1], n[2], n[0]))
        # collapse overlapping regions
        for i in range(len(HIGH_MUT_REGIONS) - 1, 0, -1):
            if HIGH_MUT_REGIONS[i - 1][2] >= HIGH_MUT_REGIONS[i][1] and HIGH_MUT_REGIONS[i - 1][0] == \
                    HIGH_MUT_REGIONS[i][0]:
                # Might need to research a more accurate way to get the mutation rate for this region
                avgMutRate = 0.5 * HIGH_MUT_REGIONS[i - 1][3] + 0.5 * HIGH_MUT_REGIONS[i][
                    3]
                HIGH_MUT_REGIONS[i - 1] = (
                    HIGH_MUT_REGIONS[i - 1][0], HIGH_MUT_REGIONS[i - 1][1], HIGH_MUT_REGIONS[i][2], avgMutRate)
                del HIGH_MUT_REGIONS[i]

    # if we didn't count ref trinucs because we found file, read in ref counts from file now
    if os.path.isfile(ref + '.trinucCounts'):
        print('reading pre-computed trinuc counts...')
        f = open(ref + '.trinucCounts', 'r')
        for line in f:
            splt = line.strip().split('\t')
            TRINUC_REF_COUNT[splt[0]] = int(splt[1])
        f.close()
    # otherwise, save trinuc counts to file, if desired
    elif save_trinuc:
        if is_bed:
            print('unable to save trinuc counts to file because using input bed region...')
        else:
            print('saving trinuc counts to file...')
            f = open(ref + '.trinucCounts', 'w')
            for trinuc in sorted(TRINUC_REF_COUNT.keys()):
                f.write(trinuc + '\t' + str(TRINUC_REF_COUNT[trinuc]) + '\n')
            f.close()

    # if for some reason we didn't find any valid input variants, exit gracefully...
    total_var = SNP_COUNT + sum(INDEL_COUNT.values())
    if total_var == 0:
        print(
            '\nError: No valid variants were found, model could not be created. (Are you using the correct reference?)\n')
        exit(1)


    ###	COMPUTE PROBABILITIES

    # frequency that each trinuc mutated into anything else
    TRINUC_MUT_PROB = {}
    # frequency that a trinuc mutates into another trinuc, given that it mutated
    TRINUC_TRANS_PROBS = {}
    # frequency of snp transitions, given a snp occurs.
    SNP_TRANS_FREQ = {}

    for trinuc in sorted(TRINUC_REF_COUNT.keys()):
        myCount = 0
        for k in sorted(TRINUC_TRANSITION_COUNT.keys()):
            if k[0] == trinuc:
                myCount += TRINUC_TRANSITION_COUNT[k]
        TRINUC_MUT_PROB[trinuc] = myCount / float(TRINUC_REF_COUNT[trinuc])
        for k in sorted(TRINUC_TRANSITION_COUNT.keys()):
            if k[0] == trinuc:
                TRINUC_TRANS_PROBS[k] = TRINUC_TRANSITION_COUNT[k] / float(myCount)

    for n1 in VALID_NUCL:
        rollingTot = sum([SNP_TRANSITION_COUNT[(n1, n2)] for n2 in VALID_NUCL if (n1, n2) in SNP_TRANSITION_COUNT])
        for n2 in VALID_NUCL:
            key2 = (n1, n2)
            if key2 in SNP_TRANSITION_COUNT:
                SNP_TRANS_FREQ[key2] = SNP_TRANSITION_COUNT[key2] / float(rollingTot)

    # compute average snp and indel frequencies
    SNP_FREQ = SNP_COUNT / float(total_var)
    AVG_INDEL_FREQ = 1. - SNP_FREQ
    INDEL_FREQ = {k: (INDEL_COUNT[k] / float(total_var)) / AVG_INDEL_FREQ for k in INDEL_COUNT.keys()}

    if is_bed:
        track_sum = float(mybed['track_len'].sum())
        AVG_MUT_RATE = total_var / track_sum
    else:
        AVG_MUT_RATE = total_var / float(TOTAL_REFLEN)

    #	if values weren't found in data, appropriately append null entries
    print_trinuc_warning = False
    for trinuc in VALID_TRINUC:
        trinuc_mut = [trinuc[0] + n + trinuc[2] for n in VALID_NUCL if n != trinuc[1]]
        if trinuc not in TRINUC_MUT_PROB:
            TRINUC_MUT_PROB[trinuc] = 0.
            print_trinuc_warning = True
        for trinuc2 in trinuc_mut:
            if (trinuc, trinuc2) not in TRINUC_TRANS_PROBS:
                TRINUC_TRANS_PROBS[(trinuc, trinuc2)] = 0.
                print_trinuc_warning = True
    if print_trinuc_warning:
        print(
            'Warning: Some trinucleotides transitions were not encountered in the input dataset, '
            'probabilities of 0.0 have been assigned to these events.')

    #
    #	print some stuff
    #
    for k in sorted(TRINUC_MUT_PROB.keys()):
        print('p(' + k + ' mutates) =', TRINUC_MUT_PROB[k])

    for k in sorted(TRINUC_TRANS_PROBS.keys()):
        print('p(' + k[0] + ' --> ' + k[1] + ' | ' + k[0] + ' mutates) =', TRINUC_TRANS_PROBS[k])

    for k in sorted(INDEL_FREQ.keys()):
        if k > 0:
            print('p(ins length = ' + str(abs(k)) + ' | indel occurs) =', INDEL_FREQ[k])
        else:
            print('p(del length = ' + str(abs(k)) + ' | indel occurs) =', INDEL_FREQ[k])

    for k in sorted(SNP_TRANS_FREQ.keys()):
        print('p(' + k[0] + ' --> ' + k[1] + ' | SNP occurs) =', SNP_TRANS_FREQ[k])


    print('p(snp)   =', SNP_FREQ)
    print('p(indel) =', AVG_INDEL_FREQ)
    print('overall average mut rate:', AVG_MUT_RATE)
    print('total variants processed:', total_var)

    #
    # save variables to file
    #
    if skip_common:
        OUT_DICT = {'AVG_MUT_RATE': AVG_MUT_RATE,
                    'SNP_FREQ': SNP_FREQ,
                    'SNP_TRANS_FREQ': SNP_TRANS_FREQ,
                    'INDEL_FREQ': INDEL_FREQ,
                    'TRINUC_MUT_PROB': TRINUC_MUT_PROB,
                    'TRINUC_TRANS_PROBS': TRINUC_TRANS_PROBS}
    else:
        OUT_DICT = {'AVG_MUT_RATE': AVG_MUT_RATE,
                    'SNP_FREQ': SNP_FREQ,
                    'SNP_TRANS_FREQ': SNP_TRANS_FREQ,
                    'INDEL_FREQ': INDEL_FREQ,
                    'TRINUC_MUT_PROB': TRINUC_MUT_PROB,
                    'TRINUC_TRANS_PROBS': TRINUC_TRANS_PROBS,
                    'COMMON_VARIANTS': COMMON_VARIANTS,
                    'HIGH_MUT_REGIONS': HIGH_MUT_REGIONS}
    pickle.dump(OUT_DICT, open(out_pickle, "wb"))


if __name__ == "__main__":
    main()
