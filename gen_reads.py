#!/usr/bin/env python
# encoding: utf-8
""" ////////////////////////////////////////////////////////////////////////////////
   ///                                                                          ///
  ///       gen_reads.py                                                        ///
 ///        VERSION 2.0: HARDER, BETTER, FASTER, STRONGER!                    ///
///////                                                                      //////
   ///      Variant and read simulator for benchmarking NGS workflows          ///
  ///                                                                         ///
 ///        Written by:     Zach Stephens                                    ///
///////     For:            DEPEND Research Group, UIUC                     ///////
   ///      Date:           May 29, 2015                                       ///
  ///       Contact:        zstephe2@illinois.edu                             ///
 ///                                                                         ///
/////////////////////////////////////////////////////////////////////////////// """

import sys
import copy
import random
import re
import time
import bisect
import pickle
import numpy as np
import argparse
import pathlib

from py.input_checking import required_field, check_file_open, is_in_range
from py.ref_func import index_ref, read_ref, get_all_ref_regions, partition_ref_regions
from py.vcf_func import parse_vcf
from py.output_file_writer import OutputFileWriter, reverse_complement, sam_flag
from py.probability import DiscreteDistribution, mean_ind_of_weighted_list
from py.SequenceContainer import SequenceContainer, ReadContainer, parse_input_mutation_model

"""//////////////////////////////////////////////////
////////////    PARSE INPUT ARGUMENTS    ////////////
//////////////////////////////////////////////////"""


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='NEAT-genReads V2.0')
    parser.add_argument('-r', type=str, required=True, metavar='<str>', help="* ref.fa")
    parser.add_argument('-R', type=int, required=True, metavar='<int>', help="* read length")
    parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output prefix")
    parser.add_argument('-c', type=float, required=False, metavar='<float>', default=10., help="average coverage")
    parser.add_argument('-e', type=str, required=False, metavar='<str>', default=None, help="sequencing error model")
    parser.add_argument('-E', type=float, required=False, metavar='<float>', default=-1,
                        help="rescale avg sequencing error rate to this")
    parser.add_argument('-p', type=int, required=False, metavar='<int>', default=2, help="ploidy")
    parser.add_argument('-t', type=str, required=False, metavar='<str>', default=None,
                        help="bed file containing targeted regions")
    parser.add_argument('-d', type=str, required=False, metavar='<str>', default=None, help="discard_regions.bed")
    parser.add_argument('-to', type=float, required=False, metavar='<float>', default=0.00,
                        help="off-target coverage scalar")
    parser.add_argument('-m', type=str, required=False, metavar='<str>', default=None,
                        help="mutation model pickle file")
    parser.add_argument('-M', type=float, required=False, metavar='<float>', default=-1,
                        help="rescale avg mutation rate to this")
    parser.add_argument('-Mb', type=str, required=False, metavar='<str>', default=None,
                        help="bed file containing positional mut rates")
    parser.add_argument('-N', type=int, required=False, metavar='<int>', default=-1,
                        help="below this qual, replace base-calls with 'N's")
    parser.add_argument('-v', type=str, required=False, metavar='<str>', default=None, help="input VCF file")

    parser.add_argument('--pe', nargs=2, type=int, required=False, metavar=('<int>', '<int>'), default=(None, None),
                        help='paired-end fragment length mean and std')
    parser.add_argument('--pe-model', type=str, required=False, metavar='<str>', default=None,
                        help='empirical fragment length distribution')
    parser.add_argument('--gc-model', type=str, required=False, metavar='<str>', default=None,
                        help='empirical GC coverage bias distribution')
    parser.add_argument('--bam', required=False, action='store_true', default=False, help='output golden BAM file')
    parser.add_argument('--vcf', required=False, action='store_true', default=False, help='output golden VCF file')
    parser.add_argument('--fa', required=False, action='store_true', default=False,
                        help='output FASTA instead of FASTQ')
    parser.add_argument('--rng', type=int, required=False, metavar='<int>', default=-1,
                        help='rng seed value; identical RNG value should produce identical runs of the program, so '
                             'things like read locations, variant positions, error positions, etc, '
                             'should all be the same.')
    # TODO check if this argument does anything at all. Near as I can tell the results are ALWAYS gzipped.
    parser.add_argument('--gz', required=False, action='store_true', default=False, help='gzip output FQ and VCF')
    parser.add_argument('--no-fastq', required=False, action='store_true', default=False,
                        help='bypass fastq generation')
    parser.add_argument('--discard-offtarget', required=False, action='store_true', default=False,
                        help='discard reads outside of targeted regions')
    parser.add_argument('--force-coverage', required=False, action='store_true', default=False,
                        help='[debug] ignore fancy models, force coverage to be constant')
    parser.add_argument('--rescale-qual', required=False, action='store_true', default=False,
                        help='rescale quality scores to match -E input')
    args = parser.parse_args(raw_args)

    """
    Set variables for processing
    """
    # absolute path to this script
    sim_path = pathlib.Path(__file__).resolve().parent
    print(sim_path)

    # if coverage val for a given window/position is below this value, consider it effectively zero.
    low_cov_thresh = 50

    # required args
    (reference, read_len, out_prefix) = (args.r, args.R, args.o)
    # various dataset parameters
    (coverage, ploids, input_bed, discard_bed, se_model, se_rate, mut_model, mut_rate, mut_bed, input_vcf) = \
        (args.c, args.p, args.t, args.d, args.e, args.E, args.m, args.M, args.Mb, args.v)
    # cancer params (disabled currently)
    # (cancer, cancer_model, cancer_purity) = (args.cancer, args.cm, args.cp)
    (cancer, cancer_model, cancer_purity) = (False, None, 0.8)
    (off_target_scalar, off_target_discard, force_coverage, rescale_qual) = (args.to, args.discard_offtarget,
                                                                             args.force_coverage, args.rescale_qual)
    # important flags
    (save_bam, save_vcf, fasta_instead, gzipped_out, no_fastq) = \
        (args.bam, args.vcf, args.fa, args.gz, args.no_fastq)

    # sequencing model parameters
    (fragment_size, fragment_std) = args.pe
    (fraglen_model, gc_bias_model) = args.pe_model, args.gc_model
    n_max_qual = args.N

    rng_seed = args.rng

    """************************************************
    ****            INPUT ERROR CHECKING
    ************************************************"""
    if (fragment_size is None and fragment_std is not None) or (fragment_size is not None and fragment_std is None):
        print('\nError: --pe argument takes 2 space-separated arguments.\n')
        exit(1)

    only_vcf = no_fastq and save_vcf and not save_bam and not fasta_instead
    if only_vcf:
        print('Only producing VCF output, that should speed things up a bit...')

    check_file_open(reference, 'ERROR: could not open reference', required=True)
    check_file_open(input_vcf, 'ERROR: could not open input VCF', required=False)
    check_file_open(input_bed, 'ERROR: could not open input BED', required=False)
    required_field(out_prefix, 'ERROR: no output prefix provided')

    # if user specified mean/std, use artificial fragment length distribution, otherwise use
    # the empirical model specified. If neither, then we're doing single-end reads.
    if fragment_size is not None and fragment_std is not None:
        paired_end = True
        paired_end_artificial = True
    elif fraglen_model is not None:
        paired_end = True
        paired_end_artificial = False
    else:
        paired_end = False
        paired_end_artificial = False

    # if user specified no fastq, not fasta only, and no bam and no vcf, then print error and exit.
    if no_fastq and not fasta_instead and not save_bam and not save_vcf:
        print('\nError: No files would be written.\n')
        exit(1)

    if no_fastq:
        print('Bypassing FASTQ generation...')

    if rng_seed == -1:
        rng_seed = random.randint(1, 99999999)
    random.seed(rng_seed)

    """************************************************
    ****             LOAD INPUT MODELS
    ************************************************"""

    # mutation models
    mut_model = parse_input_mutation_model(mut_model, 1)
    if cancer:
        cancer_model = parse_input_mutation_model(cancer_model, 2)
    if mut_rate < 0.:
        mut_rate = None

    # sequencing error model
    if se_rate < 0.:
        se_rate = None
    if se_model is None:
        print('Using default sequencing error model.')
        se_model = sim_path / 'models/errorModel_toy.p'
        se_class = ReadContainer(read_len, se_model, se_rate)
    else:
        # probably need to do some sanity checking
        se_class = ReadContainer(read_len, se_model, se_rate)

    # GC-bias model
    if gc_bias_model is None:
        print('Using default gc-bias model.')
        gc_bias_model = sim_path / 'models/gcBias_toy.p'
        [gc_scale_count, gc_scale_val] = pickle.load(open(gc_bias_model, 'rb'))
        gc_window_size = gc_scale_count[-1]
    else:
        [gc_scale_count, gc_scale_val] = pickle.load(open(gc_bias_model, 'rb'))
        gc_window_size = gc_scale_count[-1]

    # fragment length distribution
    if paired_end and not paired_end_artificial:
        print('Using empirical fragment length distribution.')
        [potential_values, potential_prob] = pickle.load(open(fraglen_model, 'rb'))
        fraglen_vals = []
        fraglen_prob = []
        for i in range(len(potential_values)):
            if potential_values[i] > read_len:
                fraglen_vals.append(potential_values[i])
                fraglen_prob.append(potential_prob[i])
        # should probably add some validation and sanity-checking code here...
        fraglen_distribution = DiscreteDistribution(fraglen_prob, fraglen_vals)
        fragment_size = fraglen_vals[mean_ind_of_weighted_list(fraglen_prob)]

    # fragment length distribution: normal distribution that goes out to +- 6 standard deviations
    if paired_end and paired_end_artificial:
        print(
            'Using artificial fragment length distribution. mean=' + str(fragment_size) + ', std=' + str(fragment_std))
        if fragment_std == 0:
            fraglen_distribution = DiscreteDistribution([1], [fragment_size], degenerate_val=fragment_size)
        else:
            potential_values = range(fragment_size - 6 * fragment_std, fragment_size + 6 * fragment_std + 1)
            fraglen_vals = []
            for i in range(len(potential_values)):
                if potential_values[i] > read_len:
                    fraglen_vals.append(potential_values[i])
            fraglen_prob = [np.exp(-(((n - float(fragment_size)) ** 2) / (2 * (fragment_std ** 2)))) for n in
                            fraglen_vals]
            fraglen_distribution = DiscreteDistribution(fraglen_prob, fraglen_vals)

    """************************************************
    ****            HARD-CODED CONSTANTS
    ************************************************"""

    # target window size for read sampling. how many times bigger than read/frag length
    WINDOW_TARGET_SCALE = 100

    # allowed nucleotides
    ALLOWED_NUCL = ['A', 'C', 'G', 'T']

    """************************************************
    ****          MORE INPUT ERROR CHECKING
    ************************************************"""

    is_in_range(read_len, 10, 1000000, 'Error: -R must be between 10 and 1,000,000')
    is_in_range(coverage, 0, 1000000, 'Error: -c must be between 0 and 1,000,000')
    is_in_range(ploids, 1, 100, 'Error: -p must be between 1 and 100')
    is_in_range(off_target_scalar, 0, 1, 'Error: -to must be between 0 and 1')
    if mut_rate != -1 and mut_rate is not None:
        is_in_range(mut_rate, 0, 0.3, 'Error: -M must be between 0 and 0.3')
    if se_rate != -1 and se_rate is not None:
        is_in_range(se_rate, 0, 0.3, 'Error: -E must be between 0 and 0.3')
    if n_max_qual != -1:
        is_in_range(n_max_qual, 1, 40, 'Error: -N must be between 1 and 40')

    """************************************************
    ****   Process Inputs
    ************************************************"""

    # index reference
    ref_index = index_ref(reference)
    if paired_end:
        n_handling = ('random', fragment_size)
    else:
        n_handling = ('ignore', read_len)
    indices_by_ref_name = {ref_index[n][0]: n for n in range(len(ref_index))}

    # parse input variants, if present
    input_variants = []
    if input_vcf is not None:
        if cancer:
            (sample_names, input_variants) = parse_vcf(input_vcf, tumor_normal=True, ploidy=ploids)
            tumor_ind = sample_names.index('TUMOR')
            normal_ind = sample_names.index('NORMAL')
        else:
            (sample_names, input_variants) = parse_vcf(input_vcf, ploidy=ploids)
        for k in sorted(input_variants.keys()):
            input_variants[k].sort()

    # parse input targeted regions, if present
    input_regions = {}
    ref_list = [n[0] for n in ref_index]
    if input_bed is not None:
        with open(input_bed, 'r') as f:
            for line in f:
                [my_chr, pos1, pos2] = line.strip().split('\t')[:3]
                if my_chr not in input_regions:
                    input_regions[my_chr] = [-1]
                input_regions[my_chr].extend([int(pos1), int(pos2)])
        # some validation
        n_in_bed_only = 0
        n_in_ref_only = 0
        for k in ref_list:
            if k not in input_regions:
                n_in_ref_only += 1
        for k in input_regions.keys():
            if not k in ref_list:
                n_in_bed_only += 1
                del input_regions[k]
        if n_in_ref_only > 0:
            print('Warning: Reference contains sequences not found in targeted regions BED file.')
        if n_in_bed_only > 0:
            print(
                'Warning: Targeted regions BED file contains sequence names not found in reference (regions ignored).')
    # parse discard bed similarly
    discardRegions = {}
    if discard_bed != None:
        f = open(discard_bed, 'r')
        for line in f:
            [myChr, pos1, pos2] = line.strip().split('\t')[:3]
            if myChr not in discardRegions:
                discardRegions[myChr] = [-1]
            discardRegions[myChr].extend([int(pos1), int(pos2)])
        f.close()

    # parse input mutation rate rescaling regions, if present
    mut_rate_regions = {}
    mut_rate_values = {}
    if mut_bed is not None:
        with open(mut_bed, 'r') as f:
            for line in f:
                [my_chr, pos1, pos2, meta_data] = line.strip().split('\t')[:4]
                mut_str = re.findall(r"mut_rate=.*?(?=;)", meta_data + ';')
                (pos1, pos2) = (int(pos1), int(pos2))
                if len(mut_str) and (pos2 - pos1) > 1:
                    # mut_rate = #_mutations / length_of_region, let's bound it by a reasonable amount
                    mut_rate = max([0.0, min([float(mut_str[0][9:]), 0.3])])
                    if my_chr not in mut_rate_regions:
                        mut_rate_regions[my_chr] = [-1]
                        mut_rate_values[my_chr] = [0.0]
                    mut_rate_regions[my_chr].extend([pos1, pos2])
                    # TODO figure out what the next line is supposed to do and fix
                    mut_rate_values.extend([mut_rate * (pos2 - pos1)] * 2)

    # initialize output files (part I)
    bam_header = None
    if save_bam:
        bam_header = [copy.deepcopy(ref_index)]
    vcf_header = None
    if save_vcf:
        vcf_header = [reference]

    # initialize output files (part II)
    if cancer:
        OFW = OutputFileWriter(out_prefix + '_normal', paired=paired_end, bam_header=bam_header, vcf_header=vcf_header,
                               gzipped=gzipped_out, no_fastq=no_fastq, fasta_instead=fasta_instead)
        OFW_CANCER = OutputFileWriter(out_prefix + '_tumor', paired=paired_end, bam_header=bam_header,
                                      vcf_header=vcf_header, gzipped=gzipped_out,
                                      no_fastq=no_fastq, fasta_instead=fasta_instead)
    else:
        OFW = OutputFileWriter(out_prefix, paired=paired_end, bam_header=bam_header, vcf_header=vcf_header,
                               gzipped=gzipped_out, no_fastq=no_fastq,
                               fasta_instead=fasta_instead)
    OUT_PREFIX_NAME = out_prefix.split('/')[-1]

    """************************************************
    ****        LET'S GET THIS PARTY STARTED...
    ************************************************"""

    read_name_count = 1  # keep track of the number of reads we've sampled, for read-names
    unmapped_records = []

    for RI in range(len(ref_index)):

        # read in reference sequence and notate blocks of Ns
        (refSequence, n_regions) = read_ref(reference, ref_index[RI], n_handling)

        # count total bp we'll be spanning so we can get an idea of how far along we are
        # (for printing progress indicators)
        total_bp_span = sum([n[1] - n[0] for n in n_regions['non_N']])
        current_progress = 0
        current_percent = 0
        have_printed100 = False

        """prune invalid input variants, e.g variants that:
                - try to delete or alter any N characters
                - don't match the reference base at their specified position
                - any alt allele contains anything other than allowed characters"""
        valid_variants = []
        n_skipped = [0, 0, 0]
        if ref_index[RI][0] in input_variants:
            for n in input_variants[ref_index[RI][0]]:
                span = (n[0], n[0] + len(n[1]))
                r_seq = str(refSequence[span[0] - 1:span[1] - 1])  # -1 because going from VCF coords to array coords
                any_bad_chr = any((nn not in ALLOWED_NUCL) for nn in [item for sublist in n[2] for item in sublist])
                if r_seq != n[1]:
                    n_skipped[0] += 1
                    continue
                elif 'N' in r_seq:
                    n_skipped[1] += 1
                    continue
                elif any_bad_chr:
                    n_skipped[2] += 1
                    continue
                # if bisect.bisect(n_regions['big'],span[0])%2 or bisect.bisect(n_regions['big'],span[1])%2:
                # continue
                valid_variants.append(n)
            print('found', len(valid_variants), 'valid variants for ' + ref_index[RI][0] + ' in input VCF...')
            if any(n_skipped):
                print(sum(n_skipped), 'variants skipped...')
                print(' - [' + str(n_skipped[0]) + '] ref allele does not match reference')
                print(' - [' + str(n_skipped[1]) + '] attempting to insert into N-region')
                print(' - [' + str(n_skipped[2]) + '] alt allele contains non-ACGT characters')

        # add large random structural variants
        #
        #	TBD!!!

        # determine sampling windows based on read length, large N regions, and structural mutations.
        # in order to obtain uniform coverage, windows should overlap by:
        #		- read_len, if single-end reads
        #		- fragment_size (mean), if paired-end reads
        # ploidy is fixed per large sampling window,
        # coverage distributions due to GC% and targeted regions are specified within these windows
        sampling_windows = []
        ALL_VARIANTS_OUT = {}
        sequences = None
        if paired_end:
            target_size = WINDOW_TARGET_SCALE * fragment_size
            overlap = fragment_size
            overlap_min_window_size = max(fraglen_distribution.values) + 10
        else:
            target_size = WINDOW_TARGET_SCALE * read_len
            overlap = read_len
            overlap_min_window_size = read_len + 10

        print('--------------------------------')
        if only_vcf:
            print('generating vcf...')
        else:
            print('sampling reads...')
        tt = time.time()

        for i in range(len(n_regions['non_N'])):
            (pi, pf) = n_regions['non_N'][i]
            n_target_windows = max([1, (pf - pi) // target_size])
            bpd = int((pf - pi) / float(n_target_windows))
            # bpd += gc_window_size - bpd%gc_window_size

            # print len(refSequence), (pi,pf), n_target_windows
            # print structural_vars

            # if for some reason our region is too small to process, skip it! (sorry)
            if n_target_windows == 1 and (pf - pi) < overlap_min_window_size:
                # print 'Does this ever happen?'
                continue

            start = pi
            end = min([start + bpd, pf])
            # print '------------------RAWR:', (pi,pf), n_target_windows, bpd
            vars_from_prev_overlap = []
            vars_cancer_from_prev_overlap = []
            vind_from_prev = 0
            is_last_time = False
            have_printed100 = False

            while True:

                # which inserted variants are in this window?
                vars_in_window = []
                updated = False
                for j in range(vind_from_prev, len(valid_variants)):
                    v_pos = valid_variants[j][0]
                    if start < v_pos < end:  # update: changed >= to >, so variant cannot be inserted in first position
                        vars_in_window.append(tuple([v_pos - 1] + list(valid_variants[j][1:])))  # vcf --> array coords
                    if v_pos >= end - overlap - 1 and updated is False:
                        updated = True
                        vind_from_prev = j
                    if v_pos >= end:
                        break

                # determine which structural variants will affect our sampling window positions
                structural_vars = []
                for n in vars_in_window:
                    buffer_needed = max([max([abs(len(n[1]) - len(alt_allele)), 1]) for alt_allele in
                                         n[2]])  # change: added abs() so that insertions are also buffered.
                    structural_vars.append(
                        (n[0] - 1, buffer_needed))  # -1 because going from VCF coords to array coords

                # adjust end-position of window based on inserted structural mutations
                buffer_added = 0
                keep_going = True
                while keep_going:
                    keep_going = False
                    for n in structural_vars:
                        # adding "overlap" here to prevent SVs from being introduced in overlap regions
                        # (which can cause problems if random mutations from the previous window land on top of them)
                        delta = (end - 1) - (n[0] + n[1]) - 2 - overlap
                        if delta < 0:
                            # print 'DELTA:', delta, 'END:', end, '-->',
                            buffer_added = -delta
                            end += buffer_added
                            ####print end
                            keep_going = True
                            break
                next_start = end - overlap
                next_end = min([next_start + bpd, pf])
                if next_end - next_start < bpd:
                    end = next_end
                    is_last_time = True

                # print progress indicator
                ####print 'PROCESSING WINDOW:',(start,end), [buffer_added], 'next:', (next_start,next_end)
                current_progress += end - start
                new_percent = int((current_progress * 100) / float(total_bp_span))
                if new_percent > current_percent:
                    if new_percent <= 99 or (new_percent == 100 and not have_printed100):
                        sys.stdout.write(str(new_percent) + '% ')
                        sys.stdout.flush()
                    current_percent = new_percent
                    if current_percent == 100:
                        have_printed100 = True

                skip_this_window = False

                # compute coverage modifiers
                coverage_avg = None
                coverage_dat = [gc_window_size, gc_scale_val, []]
                target_hits = 0
                if input_bed is None:
                    coverage_dat[2] = [1.0] * (end - start)
                else:
                    if ref_index[RI][0] not in input_regions:
                        coverage_dat[2] = [off_target_scalar] * (end - start)
                    else:
                        for j in range(start, end):
                            if not (bisect.bisect(input_regions[ref_index[RI][0]], j) % 2):
                                coverage_dat[2].append(1.0)
                                target_hits += 1
                            else:
                                coverage_dat[2].append(off_target_scalar)

                # offtarget and we're not interested?
                if off_target_discard and target_hits <= read_len:
                    coverage_avg = 0.0
                    skip_this_window = True

                # print len(coverage_dat[2]), sum(coverage_dat[2])
                if sum(coverage_dat[2]) < low_cov_thresh:
                    coverage_avg = 0.0
                    skip_this_window = True

                # check for small window sizes
                if (end - start) < overlap_min_window_size:
                    skip_this_window = True

                if skip_this_window:
                    # skip window, save cpu time
                    start = next_start
                    end = next_end
                    if is_last_time:
                        break
                    if end >= pf:
                        is_last_time = True
                    vars_from_prev_overlap = []
                    continue

                # construct sequence data that we will sample reads from
                if sequences is None:
                    sequences = SequenceContainer(start, refSequence[start:end], ploids, overlap, read_len,
                                                  [mut_model] * ploids, mut_rate, only_vcf=only_vcf)
                else:
                    sequences.update(start, refSequence[start:end], ploids, overlap, read_len, [mut_model] * ploids,
                                     mut_rate)

                # insert variants
                sequences.insert_mutations(vars_from_prev_overlap + vars_in_window)
                all_inserted_variants = sequences.random_mutations()
                # print all_inserted_variants

                # init coverage
                if sum(coverage_dat[2]) >= low_cov_thresh:
                    if paired_end:
                        coverage_avg = sequences.init_coverage(tuple(coverage_dat), frag_dist=fraglen_distribution)
                    else:
                        coverage_avg = sequences.init_coverage(tuple(coverage_dat))

                # unused cancer stuff
                if cancer:
                    tumor_sequences = SequenceContainer(start, refSequence[start:end], ploids, overlap, read_len,
                                                        [cancer_model] * ploids, mut_rate, coverage_dat)
                    tumor_sequences.insert_mutations(vars_cancer_from_prev_overlap + all_inserted_variants)
                    all_cancer_variants = tumor_sequences.random_mutations()

                # which variants do we need to keep for next time (because of window overlap)?
                vars_from_prev_overlap = []
                vars_cancer_from_prev_overlap = []
                for n in all_inserted_variants:
                    if n[0] >= end - overlap - 1:
                        vars_from_prev_overlap.append(n)
                if cancer:
                    for n in all_cancer_variants:
                        if n[0] >= end - overlap - 1:
                            vars_cancer_from_prev_overlap.append(n)

                # if we're only producing VCF, no need to go through the hassle of generating reads
                if only_vcf:
                    pass
                else:
                    window_span = end - start
                    if paired_end:
                        if force_coverage:
                            reads_to_sample = int((window_span * float(coverage)) / (2 * read_len)) + 1
                        else:
                            reads_to_sample = int((window_span * float(coverage) * coverage_avg) / (2 * read_len)) + 1
                    else:
                        if force_coverage:
                            reads_to_sample = int((window_span * float(coverage)) / read_len) + 1
                        else:
                            reads_to_sample = int((window_span * float(coverage) * coverage_avg) / read_len) + 1

                    # if coverage is so low such that no reads are to be sampled, skip region
                    #      (i.e., remove buffer of +1 reads we add to every window)
                    if reads_to_sample == 1 and sum(coverage_dat[2]) < low_cov_thresh:
                        reads_to_sample = 0

                    # sample reads
                    ASDF2_TT = time.time()
                    for i in range(reads_to_sample):

                        is_unmapped = []
                        if paired_end:
                            my_fraglen = fraglen_distribution.sample()
                            my_read_data = sequences.sample_read(se_class, my_fraglen)
                            if my_read_data is None:  # skip if we failed to find a valid position to sample read
                                continue
                            if my_read_data[0][0] is None:
                                is_unmapped.append(True)
                            else:
                                is_unmapped.append(False)
                                my_read_data[0][0] += start  # adjust mapping position based on window start
                            if my_read_data[1][0] is None:
                                is_unmapped.append(True)
                            else:
                                is_unmapped.append(False)
                                my_read_data[1][0] += start
                        else:
                            my_read_data = sequences.sample_read(se_class)
                            if my_read_data is None:  # skip if we failed to find a valid position to sample read
                                continue
                            if my_read_data[0][0] is None:  # unmapped read (lives in large insertion)
                                is_unmapped = [True]
                            else:
                                is_unmapped = [False]
                                my_read_data[0][0] += start  # adjust mapping position based on window start

                        # are we discarding offtargets?
                        outside_boundaries = []
                        if off_target_discard and input_bed is not None:
                            outside_boundaries += [bisect.bisect(input_regions[ref_index[RI][0]], n[0]) % 2 for n
                                                   in my_read_data]
                            outside_boundaries += [
                                bisect.bisect(input_regions[ref_index[RI][0]], n[0] + len(n[2])) % 2 for n in
                                my_read_data]
                        if discard_bed is not None:
                            outside_boundaries += [bisect.bisect(discardRegions[ref_index[RI][0]], n[0]) % 2 for
                                                   n in my_read_data]
                            outside_boundaries += [
                                bisect.bisect(discardRegions[ref_index[RI][0]], n[0] + len(n[2])) % 2 for n in
                                my_read_data]
                        if len(outside_boundaries) and any(outside_boundaries):
                            continue

                        my_read_name = OUT_PREFIX_NAME + '-' + ref_index[RI][0] + '-' + str(read_name_count)
                        read_name_count += len(my_read_data)

                        # if desired, replace all low-quality bases with Ns
                        if n_max_qual > -1:
                            for j in range(len(my_read_data)):
                                my_read_string = [n for n in my_read_data[j][2]]
                                for k in range(len(my_read_data[j][3])):
                                    adjusted_qual = ord(my_read_data[j][3][k]) - se_class.offQ
                                    if adjusted_qual <= n_max_qual:
                                        my_read_string[k] = 'N'
                                my_read_data[j][2] = ''.join(my_read_string)

                        # flip a coin, are we forward or reverse strand?
                        is_forward = (random.random() < 0.5)

                        # if read (or read + mate for PE) are unmapped, put them at end of bam file
                        if all(is_unmapped):
                            if paired_end:
                                if is_forward:
                                    flag1 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'first', 'mate_reverse'])
                                    flag2 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'second', 'reverse'])
                                else:
                                    flag1 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'second', 'mate_reverse'])
                                    flag2 = sam_flag(['paired', 'unmapped', 'mate_unmapped', 'first', 'reverse'])
                                unmapped_records.append((my_read_name + '/1', my_read_data[0], flag1))
                                unmapped_records.append((my_read_name + '/2', my_read_data[1], flag2))
                            else:
                                flag1 = sam_flag(['unmapped'])
                                unmapped_records.append((my_read_name + '/1', my_read_data[0], flag1))

                        my_ref_index = indices_by_ref_name[ref_index[RI][0]]

                        # write SE output
                        if len(my_read_data) == 1:
                            if no_fastq is not True:
                                if is_forward:
                                    OFW.write_fastq_record(my_read_name, my_read_data[0][2], my_read_data[0][3])
                                else:
                                    OFW.write_fastq_record(my_read_name, reverse_complement(my_read_data[0][2]),
                                                           my_read_data[0][3][::-1])
                            if save_bam:
                                if is_unmapped[0] is False:
                                    if is_forward:
                                        flag1 = 0
                                        OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[0][0],
                                                             my_read_data[0][1], my_read_data[0][2], my_read_data[0][3],
                                                             sam_flag=flag1)
                                    else:
                                        flag1 = sam_flag(['reverse'])
                                        OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[0][0],
                                                             my_read_data[0][1], my_read_data[0][2], my_read_data[0][3],
                                                             sam_flag=flag1)
                        # write PE output
                        elif len(my_read_data) == 2:
                            if no_fastq is not True:
                                OFW.write_fastq_record(my_read_name, my_read_data[0][2], my_read_data[0][3],
                                                       read2=my_read_data[1][2], qual2=my_read_data[1][3],
                                                       orientation=is_forward)
                            if save_bam:
                                if is_unmapped[0] is False and is_unmapped[1] is False:
                                    if is_forward:
                                        flag1 = sam_flag(['paired', 'proper', 'first', 'mate_reverse'])
                                        flag2 = sam_flag(['paired', 'proper', 'second', 'reverse'])
                                    else:
                                        flag1 = sam_flag(['paired', 'proper', 'second', 'mate_reverse'])
                                        flag2 = sam_flag(['paired', 'proper', 'first', 'reverse'])
                                    OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[0][0],
                                                         my_read_data[0][1], my_read_data[0][2], my_read_data[0][3],
                                                         sam_flag=flag1,
                                                         mate_pos=my_read_data[1][0])
                                    OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[1][0],
                                                         my_read_data[1][1], my_read_data[1][2], my_read_data[1][3],
                                                         sam_flag=flag2, mate_pos=my_read_data[0][0])
                                elif is_unmapped[0] is False and is_unmapped[1] is True:
                                    if is_forward:
                                        flag1 = sam_flag(['paired', 'first', 'mate_unmapped', 'mate_reverse'])
                                        flag2 = sam_flag(['paired', 'second', 'unmapped', 'reverse'])
                                    else:
                                        flag1 = sam_flag(['paired', 'second', 'mate_unmapped', 'mate_reverse'])
                                        flag2 = sam_flag(['paired', 'first', 'unmapped', 'reverse'])
                                    OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[0][0],
                                                         my_read_data[0][1], my_read_data[0][2], my_read_data[0][3],
                                                         sam_flag=flag1, mate_pos=my_read_data[0][0])
                                    OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[0][0],
                                                         my_read_data[1][1], my_read_data[1][2], my_read_data[1][3],
                                                         sam_flag=flag2, mate_pos=my_read_data[0][0], aln_map_quality=0)
                                elif is_unmapped[0] is True and is_unmapped[1] is False:
                                    if is_forward:
                                        flag1 = sam_flag(['paired', 'first', 'unmapped', 'mate_reverse'])
                                        flag2 = sam_flag(['paired', 'second', 'mate_unmapped', 'reverse'])
                                    else:
                                        flag1 = sam_flag(['paired', 'second', 'unmapped', 'mate_reverse'])
                                        flag2 = sam_flag(['paired', 'first', 'mate_unmapped', 'reverse'])
                                    OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[1][0],
                                                         my_read_data[0][1], my_read_data[0][2], my_read_data[0][3],
                                                         sam_flag=flag1, mate_pos=my_read_data[1][0], aln_map_quality=0)
                                    OFW.write_bam_record(my_ref_index, my_read_name, my_read_data[1][0],
                                                         my_read_data[1][1], my_read_data[1][2], my_read_data[1][3],
                                                         sam_flag=flag2, mate_pos=my_read_data[1][0])
                        else:
                            print('\nError: Unexpected number of reads generated...\n')
                            exit(1)
                    # print 'READS:',time.time()-ASDF2_TT

                    if not is_last_time:
                        OFW.flush_buffers(bam_max=next_start)
                    else:
                        OFW.flush_buffers(bam_max=end + 1)

                # tally up all the variants that got successfully introduced
                for n in all_inserted_variants:
                    ALL_VARIANTS_OUT[n] = True

                # prepare indices of next window
                start = next_start
                end = next_end
                if is_last_time:
                    break
                if end >= pf:
                    is_last_time = True

        if current_percent != 100 and not have_printed100:
            print('100%')
        else:
            print('')
        if only_vcf:
            print('VCF generation completed in ', end='')
        else:
            print('Read sampling completed in ', end='')
        print(int(time.time() - tt), '(sec)')

        # write all output variants for this reference
        if save_vcf:
            print('Writing output VCF...')
            for k in sorted(ALL_VARIANTS_OUT.keys()):
                current_ref = ref_index[RI][0]
                my_id = '.'
                my_quality = '.'
                my_filter = 'PASS'
                # k[0] + 1 because we're going back to 1-based vcf coords
                OFW.write_vcf_record(current_ref, str(int(k[0]) + 1), my_id, k[1], k[2], my_quality, my_filter, k[4])

    # break

    # write unmapped reads to bam file
    if save_bam and len(unmapped_records):
        print('writing unmapped reads to bam file...')
        for umr in unmapped_records:
            if paired_end:
                OFW.write_bam_record(-1, umr[0], 0, umr[1][1], umr[1][2], umr[1][3], sam_flag=umr[2], mate_pos=0,
                                     aln_map_quality=0)
            else:
                OFW.write_bam_record(-1, umr[0], 0, umr[1][1], umr[1][2], umr[1][3], sam_flag=umr[2], aln_map_quality=0)

    # close output files
    OFW.close_files()
    if cancer:
        OFW_CANCER.close_files()


if __name__ == '__main__':
    main()
