import sys
import time
import os
import re
import random


def parse_line(vcf_line, col_dict, col_samp):
    #	check if we want to proceed..
    reference_allele = vcf_line[col_dict['REF']]
    alternate_allele = vcf_line[col_dict['ALT']]
    # enough columns?
    if len(vcf_line) != len(col_dict):
        return None
    # exclude homs / filtered?
    # There was an attempt I think to make these optional that was never implemented.
    if alternate_allele == '.' or alternate_allele == '' or alternate_allele == reference_allele:
        return None
    if vcf_line[col_dict['FILTER']] != 'PASS' and vcf_line[col_dict['FILTER']] != '.':
        return None

    #	default vals
    alt_alleles = [alternate_allele]
    alt_freqs = []

    gt_per_samp = []

    #	any alt alleles?
    alt_split = alternate_allele.split(',')
    if len(alt_split) > 1:
        alt_alleles = alt_split

    #	check INFO for AF
    af = None
    if 'INFO' in col_dict and ';AF=' in ';' + vcf_line[col_dict['INFO']]:
        info = vcf_line[col_dict['INFO']] + ';'
        af = re.findall(r"AF=.*?(?=;)", info)[0][3:]
    if af is not None:
        af_splt = af.split(',')
        while (len(af_splt) < len(alt_alleles)):  # are we lacking enough AF values for some reason?
            af_splt.append(af_splt[-1])  # phone it in.
        if len(af_splt) != 0 and af_splt[0] != '.' and af_splt[0] != '':  # missing data, yay
            alt_freqs = [float(n) for n in af_splt]
    else:
        alt_freqs = [None] * max([len(alt_alleles), 1])

    gt_per_samp = None
    #	if available (i.e. we simulated it) look for WP in info
    if len(col_samp) == 0 and 'INFO' in col_dict and 'WP=' in vcf_line[col_dict['INFO']]:
        info = vcf_line[col_dict['INFO']] + ';'
        gt_per_samp = [re.findall(r"WP=.*?(?=;)", info)[0][3:]]
    else:
        #	if no sample columns, check info for GT
        if len(col_samp) == 0 and 'INFO' in col_dict and 'GT=' in vcf_line[col_dict['INFO']]:
            info = vcf_line[col_dict['INFO']] + ';'
            gt_per_samp = [re.findall(r"GT=.*?(?=;)", info)[0][3:]]
        elif len(col_samp):
            fmt = ':' + vcf_line[col_dict['FORMAT']] + ':'
            if ':GT:' in fmt:
                gt_ind = fmt.split(':').index('GT')
                gt_per_samp = [vcf_line[col_samp[iii]].split(':')[gt_ind - 1] for iii in range(len(col_samp))]
                for i in range(len(gt_per_samp)):
                    gt_per_samp[i] = gt_per_samp[i].replace('.', '0')
        if gt_per_samp is None:
            gt_per_samp = [None] * max([len(col_samp), 1])

    return alt_alleles, alt_freqs, gt_per_samp


def parse_vcf(vcf_path, tumor_normal=False, ploidy=2):
    tt = time.time()
    print('--------------------------------')
    sys.stdout.write('reading input VCF...\n')
    sys.stdout.flush()

    col_dict = {}
    col_samp = []
    n_skipped = 0
    n_skipped_because_hash = 0
    all_vars = {}  # [ref][pos]
    samp_names = []
    printed_warning = False
    f = open(vcf_path, 'r')
    for line in f:

        if line[0] != '#':
            if len(col_dict) == 0:
                print('\n\nERROR: VCF has no header?\n' + vcf_path + '\n\n')
                f.close()
                exit(1)
            splt = line[:-1].split('\t')
            pl_out = parse_line(splt, col_dict, col_samp)
            if pl_out is None:
                n_skipped += 1
            else:
                (aa, af, gt) = pl_out

                # make sure at least one allele somewhere contains the variant
                if tumor_normal:
                    gt_eval = gt[:2]
                else:
                    gt_eval = gt[:1]
                # For some reason this had an additional "if True" inserted. I guess it was supposed to be an option
                # the user could set but was never implemented.
                if None in gt_eval:
                    if not printed_warning:
                        print('Warning: Found variants without a GT field, assuming heterozygous...')
                        printed_warning = True
                    for i in range(len(gt_eval)):
                        tmp = ['0'] * ploidy
                        tmp[random.randint(0, ploidy - 1)] = '1'
                        gt_eval[i] = '/'.join(tmp)
                non_reference = False
                for gtVal in gt_eval:
                    if gtVal is not None:
                        if '1' in gtVal:
                            non_reference = True
                if not non_reference:
                    # skip if no genotype actually contains this variant
                    n_skipped += 1
                    continue

                chrom = splt[0]
                pos = int(splt[1])
                ref = splt[3]
                # skip if position is <= 0
                if pos <= 0:
                    n_skipped += 1
                    continue

                # hash variants to avoid inserting duplicates (there are some messy VCFs out there...)
                if chrom not in all_vars:
                    all_vars[chrom] = {}
                if pos not in all_vars[chrom]:
                    all_vars[chrom][pos] = (pos, ref, aa, af, gt_eval)
                else:
                    n_skipped_because_hash += 1

        else:
            if line[1] != '#':
                cols = line[1:-1].split('\t')
                for i in range(len(cols)):
                    if 'FORMAT' in col_dict:
                        col_samp.append(i)
                    col_dict[cols[i]] = i
                if len(col_samp):
                    samp_names = cols[-len(col_samp):]
                    if len(col_samp) == 1:
                        pass
                    elif len(col_samp) == 2 and tumor_normal:
                        print('Detected 2 sample columns in input VCF, assuming tumor/normal.')
                    else:
                        print(
                            'Warning: Multiple sample columns present in input VCF. By default genReads uses '
                            'only the first column.')
                else:
                    samp_names = ['Unknown']
                if tumor_normal:
                    # tumorInd  = samp_names.index('TUMOR')
                    # normalInd = samp_names.index('NORMAL')
                    if 'NORMAL' not in samp_names or 'TUMOR' not in samp_names:
                        print('\n\nERROR: Input VCF must have a "NORMAL" and "TUMOR" column.\n')
    f.close()

    vars_out = {}
    for r in all_vars.keys():
        vars_out[r] = [list(all_vars[r][k]) for k in sorted(all_vars[r].keys())]
        # prune unnecessary sequence from ref/alt alleles
        for i in range(len(vars_out[r])):
            while len(vars_out[r][i][1]) > 1 and all([n[-1] == vars_out[r][i][1][-1] for n in vars_out[r][i][2]]):
                vars_out[r][i][1] = vars_out[r][i][1][:-1]
                vars_out[r][i][2] = [n[:-1] for n in vars_out[r][i][2]]
            vars_out[r][i] = tuple(vars_out[r][i])

    print('found', sum([len(n) for n in all_vars.values()]), 'valid variants in input vcf.')
    print(' *', n_skipped, 'variants skipped: (qual filtered / ref genotypes / invalid syntax)')
    print(' *', n_skipped_because_hash, 'variants skipped due to multiple variants found per position')
    print('--------------------------------')
    return samp_names, vars_out
