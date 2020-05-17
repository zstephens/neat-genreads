#!/usr/bin/env python

# Python 3 ready

import sys
import os
import re
import bisect
import pickle
import argparse
import numpy as np
from Bio import SeqIO
import pandas as pd
#
# SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-2])
# sys.path.append(SIM_PATH)
# print(sys.path)


#########################################################
#				VARIOUS HELPER FUNCTIONS				#
#########################################################


# cluster a sorted list
def cluster_list(l: list, delta: float) -> list:
    outList = [[l[0]]]
    prevVal = l[0]
    currentInd = 0
    for n in l[1:]:
        if n - prevVal <= delta:
            outList[currentInd].append(n)
        else:
            currentInd += 1
            outList.append([])
            outList[currentInd].append(n)
        prevVal = n
    return outList


def list_2_countDict(l: list) -> dict:
    cDict = {}
    for n in l:
        if n not in cDict:
            cDict[n] = 0
        cDict[n] += 1
    return cDict


def getBedTracks(fn: str) -> dict:
    f = open(fn, 'r')
    trackDict = {}
    for line in f:
        splt = line.strip().split('\t')
        if splt[0] not in trackDict:
            trackDict[splt[0]] = []
        trackDict[splt[0]].extend([int(splt[1]), int(splt[2])])
    f.close()
    return trackDict


def getTrackLen(trackDict: dict) -> float:
    totSum = 0
    for k in trackDict.keys():
        for i in range(0, len(trackDict[k]), 2):
            totSum += trackDict[k][i + 1] - trackDict[k][i] + 1
    return totSum


def isInBed(track, ind):
    myInd = bisect.bisect(track, ind)
    if myInd & 1:
        return True
    if myInd < len(track):
        if track[myInd - 1] == ind:
            return True
    return False


#####################################
#				main()				#c
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

    parser = argparse.ArgumentParser(description='genMutModel.py')
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-m', type=str, required=True, metavar='/path/to/mutations.vcf',
                        help="Mutation file for organism in VCF format")
    parser.add_argument('-o', type=str, required=True, metavar='/path/to/output/and/prefix',
                        help="Name of output file (final model will append \'.p\')")
    parser.add_argument('-bi', type=str, required=False, metavar='Bed file of regions to include', default=None,
                        help="only_use_these_regions.bed")
    parser.add_argument('-be', type=str, required=False, metavar='Bed file of regions to exclued', default=None,
                        help="exclude_these_regions.bed")
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
    print_bed_warning = True
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

    mybed = None
    if args.bi is not None:
        print('only considering variants in specified bed regions...')
        mybed = (getBedTracks(args.bi), True)
    elif args.be is not None:
        print('only considering variants outside of specified bed regions...')
        mybed = (getBedTracks(args.be), False)

    if vcf[-4:] == '.tsv':
        print("Warning! TSV file must follow VCF specifications.")

    reference = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))

    # simplify naming
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

    # Pre-parsing
    print('reading input variants...')
    header = ["CHROM", 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    variants = pd.read_csv(vcf, sep='\t', comment='#', names=header)
    # hard-code index values based on expected columns in vcf
    (c1, c2, c3, m1, m2, m3) = (0, 1, 1, 3, 3, 4)
    variant_chroms = variants['CHROM'].to_list()
    matching_chromosomes = []
    for ref_name in ref_list:
        if ref_name not in variant_chroms:
            continue
        else:
            matching_chromosomes.append(ref_name)

    # Check to make sure there are some matches
    if not matching_chromosomes:
        print("Found no chromosomes in common between VCF and Fasta. Please fix the chromosome names and try again")
        exit(-1)

    matching_variants = variants[variants['CHROM'].isin(matching_chromosomes)]

    if matching_variants.empty:
        print("No variants matched with the Reference. This may be a chromosome naming problem.")
        exit(-1)

    # Rename variants dataframe for processing
    matching_variants = matching_variants.rename(columns={'POS': 'chr_start'})
    # Change the indexing by -1 to match VCF format indexing
    matching_variants['chr_start'] = matching_variants['chr_start'] - 1
    matching_variants['chr_end'] = matching_variants['chr_start']

    # load and process variants in each reference sequence individually, for memory reasons...
    for ref_name in matching_chromosomes:
        # Count the number of non-N nucleotides for the reference
        TOTAL_REFLEN += len(ref_dict[ref_name].seq) - ref_dict[ref_name].seq.count('N')

        # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
        VDAT_COMMON = []

        # narrow variants list to current ref
        df_to_process = matching_variants[matching_variants["CHROM"] == ref_name]

        """ ##########################################################################
        ###						COUNT TRINUCLEOTIDES IN ref						   ###
        ########################################################################## """

        if mybed is not None:
            if print_bed_warning:
                print("since you're using a bed input, "
                      "we have to count trinucs in bed region even if "
                      "you specified a trinuc count file for the reference...")
                print_bed_warning = False
            if ref_name in mybed[0]:
                ref_key = ref_name
            elif ('chr' in ref_name) and (ref_name not in mybed[0]) and (ref_name[3:] in mybed[0]):
                ref_key = ref_name[3:]
            elif ('chr' not in ref_name) and (ref_name not in mybed[0]) and ('chr' + ref_name in mybed[0]):
                ref_key = 'chr' + ref_name
            if ref_key in mybed[0]:
                subRegions = [(mybed[0][ref_key][n], mybed[0][ref_key][n + 1]) for n in
                              range(0, len(mybed[0][ref_key]), 2)]
                for sr in subRegions:
                    for i in range(sr[0], sr[1] + 1 - 2):
                        trinuc = ref_dict[ref_name][i:i + 3].seq
                        if not trinuc in VALID_TRINUC:
                            continue  # skip if trinuc contains invalid characters, or not in specified bed region
                        if trinuc not in TRINUC_BED_COUNT:
                            TRINUC_BED_COUNT[trinuc] = 0
                        TRINUC_BED_COUNT[trinuc] += 1

        if not os.path.isfile(ref + '.trinucCounts'):
            print('counting trinucleotides in reference...')
            for i in range(len(ref_dict[ref_name]) - 2):
                if i % 1000000 == 0 and i > 0:
                    print(i, '/', len(ref_dict[ref_name]))
                    # break
                trinuc = ref_dict[ref_name][i:i + 3].seq
                if not str(trinuc) in VALID_TRINUC:
                    continue  # skip if trinuc contains invalid characters
                if str(trinuc) not in TRINUC_REF_COUNT:
                    TRINUC_REF_COUNT[trinuc] = 0
                TRINUC_REF_COUNT[trinuc] += 1
        else:
            print('skipping trinuc counts (for whole reference) because we found a file...')

        """ ##########################################################################
        ###							READ INPUT VARIANTS							   ###
        ########################################################################## """

        print('reading input variants...')

        # TODO change all this to pandas

        if len(allele_ref) != len(allele_alternate):
            # indels in vcf don't include the preserved first nucleotide, so lets trim the vcf alleles
            [allele_ref, allele_normal, allele_alternate] = [allele_ref[1:], allele_normal[1:], allele_alternate[1:]]
        if not allele_ref:
            allele_ref = '-'
        if not allele_normal:
            allele_normal = '-'
        if not allele_alternate:
            allele_alternate = '-'
        # if alternate alleles are present, lets just ignore this variant. I may come back and improve this later
        if ',' in allele_alternate:
            continue
        vcf_info = ';' + splt[7] + ';'

        # if we encounter a multi-np (i.e. 3 nucl --> 3 different nucl), let's skip it for now...
        if ('-' not in allele_normal and '-' not in allele_alternate) and (
                len(allele_normal) > 1 or len(allele_alternate) > 1):
            print('skipping a complex variant...')
            continue

        # to deal with '1' vs 'chr1' references, manually change names. this is hacky and bad.
        if 'chr' not in chrName:
            chrName = 'chr' + chrName
        if 'chr' not in ref_name:
            ref_name = 'chr' + ref_name
        # skip irrelevant variants
        if chrName != ref_name:
            continue

        # if variant is outside the regions we're interested in (if specified), skip it...
        if mybed is not None:
            ref_key = ref_name
            if not ref_key in mybed[0] and ref_key[3:] in mybed[0]:  # account for 1 vs chr1, again...
                ref_key = ref_key[3:]
            if ref_key not in mybed[0]:
                inBed = False
            else:
                inBed = isInBed(mybed[0][ref_key], chrStart)
            if inBed != mybed[1]:
                continue

        # we want only snps
        # so, no '-' characters allowed, and chrStart must be same as chrEnd
        if '-' not in allele_normal and '-' not in allele_alternate and chrStart == chrEnd:
            trinuc_ref = ref_dict[ref_name][chrStart - 1:chrStart + 2]
            if not trinuc_ref in VALID_TRINUC:
                continue  # skip ref trinuc with invalid characters
            # only consider positions where ref allele in vcf matches the nucleotide in our reference
            if allele_ref == trinuc_ref[1]:
                trinuc_normal = ref_dict[ref_name][chrStart - 1] + allele_normal + ref_dict[ref_name][chrStart + 1]
                trinuc_tumor = ref_dict[ref_name][chrStart - 1] + allele_alternate + ref_dict[ref_name][chrStart + 1]
                if not trinuc_normal in VALID_TRINUC or not trinuc_tumor in VALID_TRINUC:
                    continue  # skip if mutation contains invalid char
                key = (trinuc_normal, trinuc_tumor)
                if key not in TRINUC_TRANSITION_COUNT:
                    TRINUC_TRANSITION_COUNT[key] = 0
                TRINUC_TRANSITION_COUNT[key] += 1
                SNP_COUNT += 1
                key2 = (allele_normal, allele_alternate)
                if key2 not in SNP_TRANSITION_COUNT:
                    SNP_TRANSITION_COUNT[key2] = 0
                SNP_TRANSITION_COUNT[key2] += 1

                if is_vcf:
                    myPopFreq = VCF_DEFAULT_POP_FREQ
                    if ';CAF=' in vcf_info:
                        cafStr = re.findall(r";CAF=.*?(?=;)", vcf_info)[0]
                        if ',' in cafStr:
                            myPopFreq = float(cafStr[5:].split(',')[1])
                    VDAT_COMMON.append((chrStart, allele_ref, allele_normal, allele_alternate, myPopFreq))
                else:
                    VDAT_COMMON.append((chrStart, allele_ref, allele_normal, allele_alternate))
                    TOTAL_DONORS[donor_id] = True
            else:
                print('\nError: ref allele in variant call does not match reference.\n')
                exit(1)

        # now let's look for indels...
        if '-' in allele_normal:
            len_normal = 0
        else:
            len_normal = len(allele_normal)
        if '-' in allele_alternate:
            len_tumor = 0
        else:
            len_tumor = len(allele_alternate)
        if len_normal != len_tumor:
            indel_len = len_tumor - len_normal
            if indel_len not in INDEL_COUNT:
                INDEL_COUNT[indel_len] = 0
            INDEL_COUNT[indel_len] += 1

            if is_vcf:
                myPopFreq = VCF_DEFAULT_POP_FREQ
                if ';CAF=' in vcf_info:
                    cafStr = re.findall(r";CAF=.*?(?=;)", vcf_info)[0]
                    if ',' in cafStr:
                        myPopFreq = float(cafStr[5:].split(',')[1])
                VDAT_COMMON.append((chrStart, allele_ref, allele_normal, allele_alternate, myPopFreq))
            else:
                VDAT_COMMON.append((chrStart, allele_ref, allele_normal, allele_alternate))
                TOTAL_DONORS[donor_id] = True

        # if we didn't find anything, skip ahead along to the next reference sequence
        if not len(VDAT_COMMON):
            print('Found no variants for this reference, moving along...')
            continue

        #
        # identify common mutations
        #
        percentile_var = 95
        if is_vcf:
            minVal = np.percentile([n[4] for n in VDAT_COMMON], percentile_var)
            for k in sorted(VDAT_COMMON):
                if k[4] >= minVal:
                    COMMON_VARIANTS.append((ref_name, k[0], k[1], k[3], k[4]))
            VDAT_COMMON = {(n[0], n[1], n[2], n[3]): n[4] for n in VDAT_COMMON}
        else:
            N_DONORS = len(TOTAL_DONORS)
            VDAT_COMMON = list_2_countDict(VDAT_COMMON)
            minVal = int(np.percentile(VDAT_COMMON.values(), percentile_var))
            for k in sorted(VDAT_COMMON.keys()):
                if VDAT_COMMON[k] >= minVal:
                    COMMON_VARIANTS.append((ref_name, k[0], k[1], k[3], VDAT_COMMON[k] / float(N_DONORS)))

        #
        # identify areas that have contained significantly higher random mutation rates
        #
        dist_thresh = 2000
        percentile_clust = 97
        qptn = 1000
        # identify regions with disproportionately more variants in them
        VARIANT_POS = sorted([n[0] for n in VDAT_COMMON.keys()])
        clustered_pos = cluster_list(VARIANT_POS, dist_thresh)
        byLen = [(len(clustered_pos[i]), min(clustered_pos[i]), max(clustered_pos[i]), i) for i in
                 range(len(clustered_pos))]
        # byLen  = sorted(byLen,reverse=True)
        # minLen = int(np.percentile([n[0] for n in byLen],percentile_clust))
        # byLen  = [n for n in byLen if n[0] >= minLen]
        candidate_regions = []
        for n in byLen:
            bi = int((n[1] - dist_thresh) / float(qptn)) * qptn
            bf = int((n[2] + dist_thresh) / float(qptn)) * qptn
            candidate_regions.append((n[0] / float(bf - bi), max([0, bi]), min([len(ref_dict[ref_name]), bf])))
        minVal = np.percentile([n[0] for n in candidate_regions], percentile_clust)
        for n in candidate_regions:
            if n[0] >= minVal:
                HIGH_MUT_REGIONS.append((ref_name, n[1], n[2], n[0]))
        # collapse overlapping regions
        for i in range(len(HIGH_MUT_REGIONS) - 1, 0, -1):
            if HIGH_MUT_REGIONS[i - 1][2] >= HIGH_MUT_REGIONS[i][1] and HIGH_MUT_REGIONS[i - 1][0] == \
                    HIGH_MUT_REGIONS[i][0]:
                avgMutRate = 0.5 * HIGH_MUT_REGIONS[i - 1][3] + 0.5 * HIGH_MUT_REGIONS[i][
                    3]  # not accurate, but I'm lazy
                HIGH_MUT_REGIONS[i - 1] = (
                    HIGH_MUT_REGIONS[i - 1][0], HIGH_MUT_REGIONS[i - 1][1], HIGH_MUT_REGIONS[i][2], avgMutRate)
                del HIGH_MUT_REGIONS[i]

    #
    # if we didn't count ref trinucs because we found file, read in ref counts from file now
    #
    if os.path.isfile(ref + '.trinucCounts'):
        print('reading pre-computed trinuc counts...')
        f = open(ref + '.trinucCounts', 'r')
        for line in f:
            splt = line.strip().split('\t')
            TRINUC_REF_COUNT[splt[0]] = int(splt[1])
        f.close()
    # otherwise, save trinuc counts to file, if desired
    elif save_trinuc:
        if mybed is not None:
            print('unable to save trinuc counts to file because using input bed region...')
        else:
            print('saving trinuc counts to file...')
            f = open(ref + '.trinucCounts', 'w')
            for trinuc in sorted(TRINUC_REF_COUNT.keys()):
                f.write(trinuc + '\t' + str(TRINUC_REF_COUNT[trinuc]) + '\n')
            f.close()

    #
    # if using an input bed region, make necessary adjustments to trinuc ref counts based on the bed region trinuc counts
    #
    if mybed is not None:
        if mybed[1] == True:  # we are restricting our attention to bed regions, so ONLY use bed region trinuc counts
            TRINUC_REF_COUNT = TRINUC_BED_COUNT
        else:  # we are only looking outside bed regions, so subtract bed region trinucs from entire reference trinucs
            for k in TRINUC_REF_COUNT.keys():
                if k in TRINUC_BED_COUNT:
                    TRINUC_REF_COUNT[k] -= TRINUC_BED_COUNT[k]

    # if for some reason we didn't find any valid input variants, exit gracefully...
    totalVar = SNP_COUNT + sum(INDEL_COUNT.values())
    if totalVar == 0:
        print(
            '\nError: No valid variants were found, model could not be created. (Are you using the correct reference?)\n')
        exit(1)

    """ ##########################################################################
    ###							COMPUTE PROBABILITIES						   ###
    ########################################################################## """

    # for k in sorted(TRINUC_REF_COUNT.keys()):
    #		print k, TRINUC_REF_COUNT[k]
    #
    # for k in sorted(TRINUC_TRANSITION_COUNT.keys()):
    #	print k, TRINUC_TRANSITION_COUNT[k]

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
    SNP_FREQ = SNP_COUNT / float(totalVar)
    AVG_INDEL_FREQ = 1. - SNP_FREQ
    INDEL_FREQ = {k: (INDEL_COUNT[k] / float(totalVar)) / AVG_INDEL_FREQ for k in INDEL_COUNT.keys()}
    if mybed is not None:
        if mybed[1] == True:
            AVG_MUT_RATE = totalVar / float(getTrackLen(mybed[0]))
        else:
            AVG_MUT_RATE = totalVar / float(TOTAL_REFLEN - getTrackLen(mybed[0]))
    else:
        AVG_MUT_RATE = totalVar / float(TOTAL_REFLEN)

    #
    #	if values weren't found in data, appropriately append null entries
    #
    printTrinucWarning = False
    for trinuc in VALID_TRINUC:
        trinuc_mut = [trinuc[0] + n + trinuc[2] for n in VALID_NUCL if n != trinuc[1]]
        if trinuc not in TRINUC_MUT_PROB:
            TRINUC_MUT_PROB[trinuc] = 0.
            printTrinucWarning = True
        for trinuc2 in trinuc_mut:
            if (trinuc, trinuc2) not in TRINUC_TRANS_PROBS:
                TRINUC_TRANS_PROBS[(trinuc, trinuc2)] = 0.
                printTrinucWarning = True
    if printTrinucWarning:
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

    # for n in COMMON_VARIANTS:
    #	print n

    # for n in HIGH_MUT_REGIONS:
    #	print n

    print('p(snp)   =', SNP_FREQ)
    print('p(indel) =', AVG_INDEL_FREQ)
    print('overall average mut rate:', AVG_MUT_RATE)
    print('total variants processed:', totalVar)

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
