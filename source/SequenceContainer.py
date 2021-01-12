import pdb
import random
import copy
import pathlib
import bisect
import pickle
import sys
from time import time

import numpy as np
from Bio.Seq import Seq

from source.neat_cigar_rework import CigarString
from source.probability import DiscreteDistribution, poisson_list
from source.neat_cigar import CigarString as CigarStringOld

"""
Constants needed for analysis
"""
MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3  # the maximum percentage of a window that can contain mutations

NUCL = ['A', 'C', 'G', 'T']
TRI_IND = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
           'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
NUC_IND = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALL_TRI = [NUCL[i] + NUCL[j] + NUCL[k] for i in range(len(NUCL)) for j in range(len(NUCL)) for k in range(len(NUCL))]
ALL_IND = {ALL_TRI[i]: i for i in range(len(ALL_TRI))}

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

"""
DEFAULT MUTATION MODELS
"""

DEFAULT_1_OVERALL_MUT_RATE = 0.001
DEFAULT_1_HOMOZYGOUS_FREQ = 0.010
DEFAULT_1_INDEL_FRACTION = 0.05
DEFAULT_1_INS_VS_DEL = 0.6
DEFAULT_1_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_1_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033]
DEFAULT_1_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_1_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_1 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_1_TRI_FREQS = [copy.deepcopy(example_matrix_1) for _ in range(16)]
DEFAULT_1_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
DEFAULT_MODEL_1 = [DEFAULT_1_OVERALL_MUT_RATE,
                   DEFAULT_1_HOMOZYGOUS_FREQ,
                   DEFAULT_1_INDEL_FRACTION,
                   DEFAULT_1_INS_VS_DEL,
                   DEFAULT_1_INS_LENGTH_VALUES,
                   DEFAULT_1_INS_LENGTH_WEIGHTS,
                   DEFAULT_1_DEL_LENGTH_VALUES,
                   DEFAULT_1_DEL_LENGTH_WEIGHTS,
                   DEFAULT_1_TRI_FREQS,
                   DEFAULT_1_TRINUC_BIAS]

DEFAULT_2_OVERALL_MUT_RATE = 0.002
DEFAULT_2_HOMOZYGOUS_FREQ = 0.200
DEFAULT_2_INDEL_FRACTION = 0.1
DEFAULT_2_INS_VS_DEL = 0.3
DEFAULT_2_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_2_INS_LENGTH_WEIGHTS = [0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
# noinspection DuplicatedCode
DEFAULT_2_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_2_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_2 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_2_TRI_FREQS = [copy.deepcopy(example_matrix_2) for _ in range(16)]
DEFAULT_2_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
DEFAULT_MODEL_2 = [DEFAULT_2_OVERALL_MUT_RATE,
                   DEFAULT_2_HOMOZYGOUS_FREQ,
                   DEFAULT_2_INDEL_FRACTION,
                   DEFAULT_2_INS_VS_DEL,
                   DEFAULT_2_INS_LENGTH_VALUES,
                   DEFAULT_2_INS_LENGTH_WEIGHTS,
                   DEFAULT_2_DEL_LENGTH_VALUES,
                   DEFAULT_2_DEL_LENGTH_WEIGHTS,
                   DEFAULT_2_TRI_FREQS,
                   DEFAULT_2_TRINUC_BIAS]


class SequenceContainer:
    """
    Container for reference sequences, applies mutations
    """

    def __init__(self, x_offset, sequence, ploidy, window_overlap, read_len, mut_models=None, mut_rate=None,
                 only_vcf=False):

        # initialize basic variables
        self.only_vcf = only_vcf
        self.x = x_offset
        self.ploidy = ploidy
        self.read_len = read_len
        self.sequences = [Seq(str(sequence)).tomutable() for _ in range(self.ploidy)]
        self.seq_len = len(sequence)
        self.indel_list = [[] for _ in range(self.ploidy)]
        self.snp_list = [[] for _ in range(self.ploidy)]
        self.all_cigar = [[] for _ in range(self.ploidy)]
        self.all_cigar2 = [[] for _ in range(self.ploidy)]
        self.fm_pos = [[] for _ in range(self.ploidy)]
        self.fm_span = [[] for _ in range(self.ploidy)]

        # Blacklist explanation:
        # black_list[ploid][pos] = 0		safe to insert variant here
        # black_list[ploid][pos] = 1		indel inserted here
        # black_list[ploid][pos] = 2		snp inserted here
        # black_list[ploid][pos] = 3		invalid position for various processing reasons
        self.black_list = [np.zeros(self.seq_len, dtype='<i4') for _ in range(self.ploidy)]

        # disallow mutations to occur on window overlap points
        self.win_buffer = window_overlap
        for p in range(self.ploidy):
            self.black_list[p][-self.win_buffer] = 3
            self.black_list[p][-self.win_buffer - 1] = 3

        # initialize mutation models
        if not mut_models:
            default_model = [copy.deepcopy(DEFAULT_MODEL_1) for _ in range(self.ploidy)]
            self.model_data = default_model[:self.ploidy]
        else:
            if len(mut_models) != self.ploidy:
                print('\nError: Number of mutation models received is not equal to specified ploidy\n')
                sys.exit(1)
            self.model_data = copy.deepcopy(mut_models)

        # do we need to rescale mutation frequencies?
        mut_rate_sum = sum([n[0] for n in self.model_data])
        self.mut_rescale = mut_rate
        if self.mut_rescale is None:
            self.mut_scalar = 1.0
        else:
            self.mut_scalar = float(self.mut_rescale) // (mut_rate_sum / float(len(self.model_data)))

        # how are mutations spread to each ploid, based on their specified mut rates?
        self.ploid_mut_frac = [float(n[0]) / mut_rate_sum for n in self.model_data]
        self.ploid_mut_prior = DiscreteDistribution(self.ploid_mut_frac, range(self.ploidy))

        # init mutation models
        #
        # self.models[ploid][0] = average mutation rate
        # self.models[ploid][1] = p(mut is homozygous | mutation occurs)
        # self.models[ploid][2] = p(mut is indel | mut occurs)
        # self.models[ploid][3] = p(insertion | indel occurs)
        # self.models[ploid][4] = distribution of insertion lengths
        # self.models[ploid][5] = distribution of deletion lengths
        # self.models[ploid][6] = distribution of trinucleotide SNP transitions
        # self.models[ploid][7] = p(trinuc mutates)
        self.models = []
        for n in self.model_data:
            self.models.append([self.mut_scalar * n[0], n[1], n[2], n[3], DiscreteDistribution(n[5], n[4]),
                                DiscreteDistribution(n[7], n[6]), []])
            for m in n[8]:
                # noinspection PyTypeChecker
                self.models[-1][6].append([DiscreteDistribution(m[0], NUCL), DiscreteDistribution(m[1], NUCL),
                                           DiscreteDistribution(m[2], NUCL), DiscreteDistribution(m[3], NUCL)])
            self.models[-1].append([m for m in n[9]])

        # initialize poisson attributes
        self.indel_poisson, self.snp_poisson = self.init_poisson()

        # sample the number of variants that will be inserted into each ploid
        self.indels_to_add = [n.sample() for n in self.indel_poisson]
        self.snps_to_add = [n.sample() for n in self.snp_poisson]

        # initialize trinuc snp bias
        # compute mutation positional bias given trinucleotide strings of the sequence (ONLY AFFECTS SNPs)
        #
        # note: since indels are added before snps, it's possible these positional biases aren't correctly utilized
        #       at positions affected by indels. At the moment I'm going to consider this negligible.
        trinuc_snp_bias = [[0. for _ in range(self.seq_len)] for _ in range(self.ploidy)]
        self.trinuc_bias = [None for _ in range(self.ploidy)]
        for p in range(self.ploidy):
            for i in range(self.win_buffer + 1, self.seq_len - 1):
                trinuc_snp_bias[p][i] = self.models[p][7][ALL_IND[str(self.sequences[p][i - 1:i + 2])]]
            self.trinuc_bias[p] = DiscreteDistribution(trinuc_snp_bias[p][self.win_buffer + 1:self.seq_len - 1],
                                                       range(self.win_buffer + 1, self.seq_len - 1))

        # initialize coverage attributes
        self.window_size = None
        self.coverage_distribution = None
        self.fraglen_ind_map = None

    def update_basic_vars(self, x_offset, sequence, ploidy, window_overlap, read_len):
        self.x = x_offset
        self.ploidy = ploidy
        self.read_len = read_len
        self.sequences = [Seq(str(sequence)).tomutable() for _ in range(self.ploidy)]
        self.seq_len = len(sequence)
        self.indel_list = [[] for _ in range(self.ploidy)]
        self.snp_list = [[] for _ in range(self.ploidy)]
        self.all_cigar = [[] for _ in range(self.ploidy)]
        self.all_cigar2 = [[] for _ in range(self.ploidy)]
        self.fm_pos = [[] for _ in range(self.ploidy)]
        self.fm_span = [[] for _ in range(self.ploidy)]
        self.black_list = [np.zeros(self.seq_len, dtype='<i4') for _ in range(self.ploidy)]

        # disallow mutations to occur on window overlap points
        self.win_buffer = window_overlap
        for p in range(self.ploidy):
            self.black_list[p][-self.win_buffer] = 3
            self.black_list[p][-self.win_buffer - 1] = 3

    def update_mut_models(self, mut_models, mut_rate):
        if not mut_models:
            default_model = [copy.deepcopy(DEFAULT_MODEL_1) for _ in range(self.ploidy)]
            self.model_data = default_model[:self.ploidy]
        else:
            if len(mut_models) != self.ploidy:
                print('\nError: Number of mutation models received is not equal to specified ploidy\n')
                sys.exit(1)
            self.model_data = copy.deepcopy(mut_models)

        # do we need to rescale mutation frequencies?
        mut_rate_sum = sum([n[0] for n in self.model_data])
        self.mut_rescale = mut_rate
        if self.mut_rescale is None:
            self.mut_scalar = 1.0
        else:
            self.mut_scalar = float(self.mut_rescale) // (mut_rate_sum / float(len(self.model_data)))

        # how are mutations spread to each ploid, based on their specified mut rates?
        self.ploid_mut_frac = [float(n[0]) / mut_rate_sum for n in self.model_data]
        self.ploid_mut_prior = DiscreteDistribution(self.ploid_mut_frac, range(self.ploidy))
        self.models = []
        for n in self.model_data:
            self.models.append([self.mut_scalar * n[0], n[1], n[2], n[3], DiscreteDistribution(n[5], n[4]),
                                DiscreteDistribution(n[7], n[6]), []])
            for m in n[8]:
                # noinspection PyTypeChecker
                self.models[-1][6].append([DiscreteDistribution(m[0], NUCL), DiscreteDistribution(m[1], NUCL),
                                           DiscreteDistribution(m[2], NUCL), DiscreteDistribution(m[3], NUCL)])
            self.models[-1].append([m for m in n[9]])

    def update_trinuc_bias(self):
        trinuc_snp_bias = [[0. for _ in range(self.seq_len)] for _ in range(self.ploidy)]
        self.trinuc_bias = [None for _ in range(self.ploidy)]
        for p in range(self.ploidy):
            for i in range(self.win_buffer + 1, self.seq_len - 1):
                trinuc_snp_bias[p][i] = self.models[p][7][ALL_IND[str(self.sequences[p][i - 1:i + 2])]]
            self.trinuc_bias[p] = DiscreteDistribution(trinuc_snp_bias[p][self.win_buffer + 1:self.seq_len - 1],
                                                       range(self.win_buffer + 1, self.seq_len - 1))

    def init_coverage(self, coverage_data, frag_dist=None):
        """
        Initializes coverage for the sequence container. Only makes changes if we are not in vcf-only mode.

        :param coverage_data: A tuple containing the window size, gc scalars and target coverage values.
        :param frag_dist: A probability distribution of the fragment size.
        :return: Mean coverage value
        """

        # If we're only creating a vcf, skip some expensive initialization related to coverage depth
        if not self.only_vcf:
            (self.window_size, gc_scalars, target_cov_vals) = coverage_data
            gc_cov_vals = [[] for _ in self.sequences]
            tr_cov_vals = [[] for _ in self.sequences]
            avg_out = []
            self.coverage_distribution = []
            for i in range(len(self.sequences)):
                max_coord = min([len(self.sequences[i]) - self.read_len, len(self.all_cigar[i]) - self.read_len])
                # Trying to fix a problem wherein the above line gives a negative answer
                if max_coord <= 0:
                    max_coord = min([len(self.sequences[i]), len(self.all_cigar[i])])
                # compute gc-bias
                j = 0
                while j + self.window_size < len(self.sequences[i]):
                    gc_c = self.sequences[i][j:j + self.window_size].count('G') + \
                           self.sequences[i][j:j + self.window_size].count('C')
                    gc_cov_vals[i].extend([gc_scalars[gc_c]] * self.window_size)
                    j += self.window_size
                gc_c = self.sequences[i][-self.window_size:].count('G') + \
                       self.sequences[i][-self.window_size:].count('C')
                gc_cov_vals[i].extend([gc_scalars[gc_c]] * (len(self.sequences[i]) - len(gc_cov_vals[i])))
                #
                tr_cov_vals[i].append(target_cov_vals[0])
                prev_val = self.fm_pos[i][0]
                for j in range(1, max_coord):
                    if self.fm_pos[i][j] is None:
                        tr_cov_vals[i].append(target_cov_vals[prev_val])
                    elif self.fm_span[i][j] - self.fm_pos[i][j] <= 1:
                        tr_cov_vals[i].append(target_cov_vals[prev_val])
                    else:
                        tr_cov_vals[i].append(sum(target_cov_vals[self.fm_pos[i][j]:self.fm_span[i][j]]) / float(
                            self.fm_span[i][j] - self.fm_pos[i][j]))
                        prev_val = self.fm_pos[i][j]

                # shift by half of read length
                if len(tr_cov_vals[i]) > int(self.read_len / 2.):
                    tr_cov_vals[i] = [0.0] * int(self.read_len // 2) + tr_cov_vals[i][:-int(self.read_len / 2.)]
                # fill in missing indices
                tr_cov_vals[i].extend([0.0] * (len(self.sequences[i]) - len(tr_cov_vals[i])))

                #
                coverage_vector = np.cumsum([tr_cov_vals[i][nnn] *
                                             gc_cov_vals[i][nnn] for nnn in range(len(tr_cov_vals[i]))])
                coverage_vals = []
                for j in range(0, max_coord):
                    coverage_vals.append(coverage_vector[j + self.read_len] - coverage_vector[j])
                avg_out.append(np.mean(coverage_vals) / float(self.read_len))

                if frag_dist is None:
                    self.coverage_distribution.append(DiscreteDistribution(coverage_vals, range(len(coverage_vals))))

                # fragment length nightmare
                else:
                    current_thresh = 0.
                    index_list = [0]
                    for j in range(len(frag_dist.cum_prob)):
                        if frag_dist.cum_prob[j] >= current_thresh + COV_FRAGLEN_PERCENTILE / 100.0:
                            current_thresh = frag_dist.cum_prob[j]
                            index_list.append(j)
                    flq = [frag_dist.values[nnn] for nnn in index_list]
                    if frag_dist.values[-1] not in flq:
                        flq.append(frag_dist.values[-1])
                    flq.append(LARGE_NUMBER)

                    self.fraglen_ind_map = {}
                    for j in frag_dist.values:
                        b_ind = bisect.bisect(flq, j)
                        if abs(flq[b_ind - 1] - j) <= abs(flq[b_ind] - j):
                            self.fraglen_ind_map[j] = flq[b_ind - 1]
                        else:
                            self.fraglen_ind_map[j] = flq[b_ind]

                    self.coverage_distribution.append({})
                    for flv in sorted(list(set(self.fraglen_ind_map.values()))):
                        buffer_val = self.read_len
                        for j in frag_dist.values:
                            if self.fraglen_ind_map[j] == flv and j > buffer_val:
                                buffer_val = j
                        max_coord = min([len(self.sequences[i]) - buffer_val - 1,
                                         len(self.all_cigar[i]) - buffer_val + self.read_len - 2])
                        # print 'BEFORE:', len(self.sequences[i])-buffer_val
                        # print 'AFTER: ', len(self.all_cigar[i])-buffer_val+self.read_len-2
                        # print 'AFTER2:', max_coord
                        coverage_vals = []
                        for j in range(0, max_coord):
                            coverage_vals.append(
                                coverage_vector[j + self.read_len] - coverage_vector[j] + coverage_vector[j + flv] -
                                coverage_vector[
                                    j + flv - self.read_len])

                        # EXPERIMENTAL
                        # quantized_covVals = quantize_list(coverage_vals)
                        # self.coverage_distribution[i][flv] = DiscreteDistribution([n[2] for n in quantized_covVals],[(n[0],n[1]) for n in quantized_covVals])

                        # TESTING
                        # import matplotlib.pyplot as mpl
                        # print len(coverage_vals),'-->',len(quantized_covVals)
                        # mpl.figure(0)
                        # mpl.plot(range(len(coverage_vals)), coverage_vals)
                        # for qcv in quantized_covVals:
                        #	mpl.plot([qcv[0], qcv[1]+1], [qcv[2],qcv[2]], 'r')
                        # mpl.show()
                        # sys.exit(1)

                        self.coverage_distribution[i][flv] = DiscreteDistribution(coverage_vals,
                                                                                  range(len(coverage_vals)))

            return np.mean(avg_out)

    def init_poisson(self):
        ind_l_list = [self.seq_len * self.models[i][0] * self.models[i][2] * self.ploid_mut_frac[i] for i in
                      range(len(self.models))]
        snp_l_list = [self.seq_len * self.models[i][0] * (1. - self.models[i][2]) * self.ploid_mut_frac[i] for i in
                      range(len(self.models))]
        k_range = range(int(self.seq_len * MAX_MUTFRAC))
        # return (indel_poisson, snp_poisson)
        return [poisson_list(k_range, ind_l_list[n]) for n in range(len(self.models))], \
               [poisson_list(k_range, snp_l_list[n]) for n in range(len(self.models))]

    def update(self, x_offset, sequence, ploidy, window_overlap, read_len, mut_models=None, mut_rate=None):
        # if mutation model is changed, we have to reinitialize it...
        if ploidy != self.ploidy or mut_rate != self.mut_rescale or mut_models is not None:
            self.ploidy = ploidy
            self.mut_rescale = mut_rate
            self.update_mut_models(mut_models, mut_rate)
        # if sequence length is different than previous window, we have to redo snp/indel poissons
        if len(sequence) != self.seq_len:
            self.seq_len = len(sequence)
            self.indel_poisson, self.snp_poisson = self.init_poisson()
        # basic vars
        self.update_basic_vars(x_offset, sequence, ploidy, window_overlap, read_len)
        self.indels_to_add = [n.sample() for n in self.indel_poisson]
        self.snps_to_add = [n.sample() for n in self.snp_poisson]
        # initialize trinuc snp bias
        if not IGNORE_TRINUC:
            self.update_trinuc_bias()

    def insert_mutations(self, input_list):
        for input_variable in input_list:
            which_ploid = []
            wps = input_variable[4][0]

            # if no genotype given, assume heterozygous and choose a single ploid based on their mut rates
            if wps is None:
                which_ploid.append(self.ploid_mut_prior.sample())
                which_alt = [0]
            else:
                if '/' in wps or '|' in wps:
                    if '/' in wps:
                        splt = wps.split('/')
                    else:
                        splt = wps.split('|')
                    which_ploid = []
                    for i in range(len(splt)):
                        if splt[i] == '1':
                            which_ploid.append(i)
                    # assume we're just using first alt for inserted variants?
                    which_alt = [0 for _ in which_ploid]
                # otherwise assume monoploidy
                else:
                    which_ploid = [0]
                    which_alt = [0]

            # ignore invalid ploids
            for i in range(len(which_ploid) - 1, -1, -1):
                if which_ploid[i] >= self.ploidy:
                    del which_ploid[i]

            for i in range(len(which_ploid)):
                p = which_ploid[i]
                my_alt = input_variable[2][which_alt[i]]
                my_var = (input_variable[0] - self.x, input_variable[1], my_alt)
                in_len = max([len(input_variable[1]), len(my_alt)])

                if my_var[0] < 0 or my_var[0] >= len(self.black_list[p]):
                    print('\nError: Attempting to insert variant out of window bounds:')
                    print(my_var, '--> blackList[0:' + str(len(self.black_list[p])) + ']\n')
                    sys.exit(1)
                if len(input_variable[1]) == 1 and len(my_alt) == 1:
                    if self.black_list[p][my_var[0]]:
                        continue
                    self.snp_list[p].append(my_var)
                    self.black_list[p][my_var[0]] = 2
                else:
                    indel_failed = False
                    for k in range(my_var[0], my_var[0] + in_len + 1):
                        if k >= len(self.black_list[p]):
                            indel_failed = True
                            continue
                        if self.black_list[p][k]:
                            indel_failed = True
                            continue
                    if indel_failed:
                        continue
                    for k in range(my_var[0], my_var[0] + in_len + 1):
                        self.black_list[p][k] = 1
                    self.indel_list[p].append(my_var)

    def random_mutations(self):

        # add random indels
        all_indels = [[] for _ in self.sequences]
        for i in range(self.ploidy):
            for j in range(self.indels_to_add[i]):
                # insert homozygous indel
                if random.random() <= self.models[i][1]:
                    which_ploid = range(self.ploidy)
                # insert heterozygous indel
                else:
                    which_ploid = [self.ploid_mut_prior.sample()]

                # try to find suitable places to insert indels
                event_pos = -1
                for attempt in range(MAX_ATTEMPTS):
                    event_pos = random.randint(self.win_buffer, self.seq_len - 1)
                    for p in which_ploid:
                        if self.black_list[p][event_pos]:
                            event_pos = -1
                    if event_pos != -1:
                        break
                if event_pos == -1:
                    continue

                # insertion
                if random.random() <= self.models[i][3]:
                    in_len = self.models[i][4].sample()
                    # sequence content of random insertions is uniformly random (change this later, maybe)
                    in_seq = ''.join([random.choice(NUCL) for _ in range(in_len)])
                    ref_nucl = self.sequences[i][event_pos]
                    my_indel = (event_pos, ref_nucl, ref_nucl + in_seq)
                # deletion
                else:
                    in_len = self.models[i][5].sample()
                    # skip if deletion too close to boundary
                    if event_pos + in_len + 1 >= len(self.sequences[i]):
                        continue
                    if in_len == 1:
                        in_seq = self.sequences[i][event_pos + 1]
                    else:
                        in_seq = str(self.sequences[i][event_pos + 1:event_pos + in_len + 1])
                    ref_nucl = self.sequences[i][event_pos]
                    my_indel = (event_pos, ref_nucl + in_seq, ref_nucl)

                # if event too close to boundary, skip. if event conflicts with other indel, skip.
                skip_event = False
                if event_pos + len(my_indel[1]) >= self.seq_len - self.win_buffer - 1:
                    skip_event = True
                if skip_event:
                    continue
                for p in which_ploid:
                    for k in range(event_pos, event_pos + in_len + 1):
                        if self.black_list[p][k]:
                            skip_event = True
                if skip_event:
                    continue

                for p in which_ploid:
                    for k in range(event_pos, event_pos + in_len + 1):
                        self.black_list[p][k] = 1
                    all_indels[p].append(my_indel)

        # add random snps
        all_snps = [[] for _ in self.sequences]
        for i in range(self.ploidy):
            for j in range(self.snps_to_add[i]):
                # insert homozygous SNP
                if random.random() <= self.models[i][1]:
                    which_ploid = range(self.ploidy)
                # insert heterozygous SNP
                else:
                    which_ploid = [self.ploid_mut_prior.sample()]

                # try to find suitable places to insert snps
                event_pos = -1
                for attempt in range(MAX_ATTEMPTS):
                    # based on the mutation model for the specified ploid, choose a SNP location based on trinuc bias
                    # (if there are multiple ploids, choose one at random)
                    if IGNORE_TRINUC:
                        event_pos = random.randint(self.win_buffer + 1, self.seq_len - 2)
                    else:
                        ploid_to_use = which_ploid[random.randint(0, len(which_ploid) - 1)]
                        event_pos = self.trinuc_bias[ploid_to_use].sample()
                    for p in which_ploid:
                        if self.black_list[p][event_pos]:
                            event_pos = -1
                    if event_pos != -1:
                        break
                if event_pos == -1:
                    continue

                ref_nucl = self.sequences[i][event_pos]
                context = str(self.sequences[i][event_pos - 1]) + str(self.sequences[i][event_pos + 1])
                # sample from tri-nucleotide substitution matrices to get SNP alt allele
                new_nucl = self.models[i][6][TRI_IND[context]][NUC_IND[ref_nucl]].sample()
                my_snp = (event_pos, ref_nucl, new_nucl)

                for p in which_ploid:
                    all_snps[p].append(my_snp)
                    self.black_list[p][my_snp[0]] = 2

        # combine random snps with inserted snps, remove any snps that overlap indels
        for p in range(len(all_snps)):
            all_snps[p].extend(self.snp_list[p])
            all_snps[p] = [n for n in all_snps[p] if self.black_list[p][n[0]] != 1]

        # MODIFY REFERENCE STRING: SNPS
        for i in range(len(all_snps)):
            for j in range(len(all_snps[i])):
                v_pos = all_snps[i][j][0]

                if all_snps[i][j][1] != self.sequences[i][v_pos]:
                    print('\nError: Something went wrong!\n', all_snps[i][j], self.sequences[i][v_pos], '\n')
                    print(all_snps[i][j])
                    print(self.sequences[i][v_pos])
                    sys.exit(1)
                else:
                    self.sequences[i][v_pos] = all_snps[i][j][2]

        # organize the indels we want to insert
        for i in range(len(all_indels)):
            all_indels[i].extend(self.indel_list[i])
        all_indels_ins = [sorted([list(m) for m in n]) for n in all_indels]

        # MODIFY REFERENCE STRING: INDELS
        for i in range(len(all_indels_ins)):
            rolling_adj = 0
            temp_symbol_string = CigarString(str(len(self.sequences[i])) + "M")
            # TODO Delete commented out lines once CigarString works 100%
            temp_symbol_string2 = ['M' for _ in self.sequences[i]]
            assert(temp_symbol_string.cigar == CigarStringOld(list_in=temp_symbol_string2).get_string())

            for j in range(len(all_indels_ins[i])):
                v_pos = all_indels_ins[i][j][0] + rolling_adj
                v_pos2 = v_pos + len(all_indels_ins[i][j][1])
                indel_length = len(all_indels_ins[i][j][2]) - len(all_indels_ins[i][j][1])
                rolling_adj += indel_length

                if all_indels_ins[i][j][1] != str(self.sequences[i][v_pos:v_pos2]):
                    print('\nError: Something went wrong!\n', all_indels_ins[i][j], [v_pos, v_pos2],
                          str(self.sequences[i][v_pos:v_pos2]), '\n')
                    sys.exit(1)
                else:
                    # alter reference sequence
                    self.sequences[i] = self.sequences[i][:v_pos] + Seq(all_indels_ins[i][j][2]).tomutable() + \
                                        self.sequences[i][v_pos2:]
                    # notate indel positions for cigar computation
                    if indel_length > 0:
                        cigar_to_insert = CigarString(str(indel_length) + 'I')
                        temp_symbol_string.insert_cigar_element(v_pos + 1, cigar_to_insert,
                                                                len(all_indels_ins[i][j][1]))
                        # TODO Delete commented out lines once CigarString works 100%
                        check = temp_symbol_string2[:v_pos + 1]
                        temp_symbol_string2 = temp_symbol_string2[:v_pos + 1] + \
                                              ['I'] * indel_length + temp_symbol_string2[v_pos2 + 1:]
                        assert(temp_symbol_string.cigar == CigarStringOld(list_in=temp_symbol_string2).get_string())
                    elif indel_length < 0:
                        cigar_to_insert = CigarString(str(abs(indel_length)) + 'D1M')
                        temp_symbol_string.insert_cigar_element(v_pos + 1, cigar_to_insert)
                        # TODO Delete commented out lines once CigarString works 100%
                        check = temp_symbol_string2[v_pos + 1]
                        temp_symbol_string2[v_pos + 1] = 'D' * abs(indel_length) + 'M'
                        assert (temp_symbol_string.cigar == CigarStringOld(list_in=temp_symbol_string2).get_string())
                        if temp_symbol_string.cigar == '28662M4D4821M':
                            print("check this one")
                            pdb.set_trace()


            # pre-compute cigar strings
            for j in range(len(temp_symbol_string) - self.read_len):
                self.all_cigar[i].append(temp_symbol_string.get_cigar_fragment(j, j + self.read_len))
                # TODO Delete commented out lines once CigarString works 100%
                self.all_cigar2[i].append(CigarStringOld(list_in=temp_symbol_string2[j:j + self.read_len]).get_string())
                if temp_symbol_string.get_cigar_fragment(j, j + self.read_len).cigar != CigarStringOld(list_in=temp_symbol_string2[j:j + self.read_len]).get_string():
                    print("new " + str(temp_symbol_string.get_cigar_fragment(j, j + self.read_len)))
                    print("old " + CigarStringOld(list_in=temp_symbol_string2[j:j + self.read_len]).get_string())
                    pdb.set_trace()

            # create some data structures we will need later:
            # --- self.FM_pos[ploid][pos]: position of the left-most matching base (IN REFERENCE COORDINATES, i.e.
            #       corresponding to the unmodified reference genome)
            # --- self.FM_span[ploid][pos]: number of reference positions spanned by a read originating from
            #       this coordinate
            md_so_far = 0
            temp_symbol_string_list = temp_symbol_string.string_to_list()
            for j in range(len(temp_symbol_string_list)):
                self.fm_pos[i].append(md_so_far)
                # fix an edge case with deletions
                if temp_symbol_string_list[j] == 'D':
                    self.fm_pos[i][-1] += temp_symbol_string_list[j].count('D')
                # compute number of ref matches for each read
                span_dif = len([n for n in temp_symbol_string_list[j:j + self.read_len] if 'M' in n])
                self.fm_span[i].append(self.fm_pos[i][-1] + span_dif)
                md_so_far += temp_symbol_string_list[j].count('M') + temp_symbol_string_list[j].count('D')

        # tally up all the variants we handled...
        count_dict = {}
        all_variants = [sorted(all_snps[i] + all_indels[i]) for i in range(self.ploidy)]
        for i in range(len(all_variants)):
            for j in range(len(all_variants[i])):
                all_variants[i][j] = tuple([all_variants[i][j][0] + self.x]) + all_variants[i][j][1:]
                t = tuple(all_variants[i][j])
                if t not in count_dict:
                    count_dict[t] = []
                count_dict[t].append(i)

        # TODO: combine multiple variants that happened to occur at same position into single vcf entry?

        output_variants = []
        for k in sorted(count_dict.keys()):
            output_variants.append(k + tuple([len(count_dict[k]) / float(self.ploidy)]))
            ploid_string = ['0' for _ in range(self.ploidy)]
            for k2 in [n for n in count_dict[k]]:
                ploid_string[k2] = '1'
            output_variants[-1] += tuple(['WP=' + '/'.join(ploid_string)])
        return output_variants

    def sample_read(self, sequencing_model, frag_len=None):

        # choose a ploid
        my_ploid = random.randint(0, self.ploidy - 1)

        # stop attempting to find a valid position if we fail enough times
        MAX_READPOS_ATTEMPTS = 100
        attempts_thus_far = 0

        # choose a random position within the ploid, and generate quality scores / sequencing errors
        reads_to_sample = []
        if frag_len is None:
            r_pos = self.coverage_distribution[my_ploid].sample()

            # sample read position and call function to compute quality scores / sequencing errors
            r_dat = self.sequences[my_ploid][r_pos:r_pos + self.read_len]
            (my_qual, my_errors) = sequencing_model.get_sequencing_errors(r_dat)
            reads_to_sample.append([r_pos, my_qual, my_errors, r_dat])

        else:
            r_pos1 = self.coverage_distribution[my_ploid][self.fraglen_ind_map[frag_len]].sample()

            # EXPERIMENTAL
            # coords_to_select_from = self.coverage_distribution[my_ploid][self.fraglens_ind_map[frag_len]].sample()
            # r_pos1 = random.randint(coords_to_select_from[0],coords_to_select_from[1])

            r_pos2 = r_pos1 + frag_len - self.read_len
            r_dat1 = self.sequences[my_ploid][r_pos1:r_pos1 + self.read_len]
            r_dat2 = self.sequences[my_ploid][r_pos2:r_pos2 + self.read_len]
            (my_qual1, my_errors1) = sequencing_model.get_sequencing_errors(r_dat1)
            (my_qual2, my_errors2) = sequencing_model.get_sequencing_errors(r_dat2, is_reverse_strand=True)
            reads_to_sample.append([r_pos1, my_qual1, my_errors1, r_dat1])
            reads_to_sample.append([r_pos2, my_qual2, my_errors2, r_dat2])

        # error format:
        # myError[i] = (type, len, pos, ref, alt)

        """
        examine sequencing errors to-be-inserted.
            - remove deletions that don't have enough bordering sequence content to "fill in"
            if error is valid, make the changes to the read data
        """
        read_out = []
        for read in reads_to_sample:
            try:
                my_cigar = self.all_cigar[my_ploid][read[0]]
            except IndexError:
                print('Index error when attempting to find cigar string.')
                print(my_ploid, len(self.all_cigar[my_ploid]), read[0])
                if frag_len is not None:
                    print((r_pos1, r_pos2))
                    print(frag_len, self.fraglen_ind_map[frag_len])
                sys.exit(1)
            total_d = sum([error[1] for error in read[2] if error[0] == 'D'])
            total_i = sum([error[1] for error in read[2] if error[0] == 'I'])
            avail_b = len(self.sequences[my_ploid]) - read[0] - self.read_len - 1

            if total_d:
                print("there is a deletion error")
                pdb.set_trace()
            # add buffer sequence to fill in positions that get deleted
            read[3] += self.sequences[my_ploid][read[0] + self.read_len:read[0] + self.read_len + total_d]
            expanded_cigar = []
            adj = 0
            sse_adj = [0 for _ in range(self.read_len + max(sequencing_model.err_p[3]))]
            any_indel_err = False

            # sort by letter (D > I > S) such that we introduce all indel errors before substitution errors
            # secondarily, sort by index
            arranged_errors = {'D': [], 'I': [], 'S': []}
            for error in read[2]:
                arranged_errors[error[0]].append((error[2], error))
            sorted_errors = []
            for k in sorted(arranged_errors.keys()):
                sorted_errors.extend([n[1] for n in sorted(arranged_errors[k])])

            skip_indels = False

            # FIXED TdB 05JUN2018
            # Moved this outside the for error loop, since it messes up the CIGAR string when
            # more than one deletion is in the same read
            extra_cigar_val = []
            # END FIXED TdB

            for error in sorted_errors:
                e_len = error[1]
                e_pos = error[2]
                if error[0] == 'D' or error[0] == 'I':
                    any_indel_err = True

                    # FIXED TdB 05JUN2018
                    # Moved this OUTSIDE the for error loop, since it messes up the CIGAR string
                    # when more than one deletion is in the same read
                    # extra_cigar_val = []
                    # END FIXED TdB

                    if total_d > avail_b:  # if not enough bases to fill-in deletions, skip all indel erors
                        continue
                    if not expanded_cigar:
                        expanded_cigar = my_cigar.string_to_list()
                        # TODO delete these lines using old CigarString once it is working 100%
                        expanded_cigar2 = CigarStringOld(string_in=my_cigar).get_list()
                        assert(expanded_cigar == expanded_cigar2)
                        fill_to_go = total_d - total_i + 1
                        if fill_to_go > 0:
                            try:
                                # TODO delete these lines using old CigarString once it is working 100%
                                extra_cigar_val2 = CigarStringOld(string_in=self.all_cigar2[my_ploid][read[0]
                                                               + fill_to_go]).get_list()[-fill_to_go:]
                                extra_cigar_val = self.all_cigar[my_ploid][read[0]
                                                                           + fill_to_go].string_to_list()[-fill_to_go:]

                            except IndexError:
                                # Applying the deletions we want requires going beyond region boundaries.
                                # Skip all indel errors
                                skip_indels = True

                            assert(extra_cigar_val == extra_cigar_val2)

                    if skip_indels:
                        continue

                    # insert deletion error into read and update cigar string accordingly
                    if error[0] == 'D':
                        my_adj = sse_adj[e_pos]
                        pi = e_pos + my_adj
                        pf = e_pos + my_adj + e_len + 1
                        if str(read[3][pi:pf]) == str(error[3]):
                            read[3] = read[3][:pi + 1] + read[3][pf:]
                            expanded_cigar = expanded_cigar[:pi + 1] + expanded_cigar[pf:]
                            # weird edge case with del at very end of region. Make a guess and add a "M"
                            if pi + 1 == len(expanded_cigar):
                                expanded_cigar.append('M')
                            try:
                                expanded_cigar[pi + 1] = 'D' * e_len + expanded_cigar[pi + 1]
                            except IndexError:
                                print("Bug!! Index error on expanded cigar")
                                pdb.set_trace()
                                sys.exit(1)

                        else:
                            print('\nError, ref does not match alt while attempting to insert deletion error!\n')
                            sys.exit(1)
                        adj -= e_len
                        for i in range(e_pos, len(sse_adj)):
                            sse_adj[i] -= e_len

                    # insert insertion error into read and update cigar string accordingly
                    else:
                        my_adj = sse_adj[e_pos]
                        if str(read[3][e_pos + my_adj]) == error[3]:
                            read[3] = read[3][:e_pos + my_adj] + error[4] + read[3][e_pos + my_adj + 1:]
                            expanded_cigar = expanded_cigar[:e_pos + my_adj] + ['I'] * e_len + expanded_cigar[
                                                                                               e_pos + my_adj:]
                        else:
                            print('\nError, ref does not match alt while attempting to insert insertion error!\n')
                            print('---', chr(read[3][e_pos + my_adj]), '!=', error[3])
                            sys.exit(1)
                        adj += e_len
                        for i in range(e_pos, len(sse_adj)):
                            sse_adj[i] += e_len

                else:  # substitution errors, much easier by comparison...
                    if str(read[3][e_pos + sse_adj[e_pos]]) == error[3]:
                        read[3][e_pos + sse_adj[e_pos]] = error[4]
                    else:
                        print('\nError, ref does not match alt while attempting to insert substitution error!\n')
                        sys.exit(1)

            if len(my_cigar) != 100:
                print(str(my_cigar) + " is not equal to 100.")
                pdb.set_trace()
            if any_indel_err:
                if len(expanded_cigar):
                    relevant_cigar = (expanded_cigar + extra_cigar_val)[:self.read_len]
                    my_cigar = CigarString(CigarString.list_to_string(relevant_cigar))

                    # debugging
                    if len(my_cigar) != 100:
                        print(my_cigar)
                        pdb.set_trace()

                    # TODO delete this line once new cigar is 100% working
                    my_cigar2 = CigarStringOld(list_in=relevant_cigar).get_string()

                assert(my_cigar == my_cigar2)
                read[3] = read[3][:self.read_len]

            read_out.append([self.fm_pos[my_ploid][read[0]], my_cigar, read[3], str(read[1])])

        # read_out[i] = (pos, cigar, read_string, qual_string)
        return read_out


class SequencingError:
    """
    Container to model sequencing errors: computes quality scores and positions to insert errors
    """

    def __init__(self, read_len, error_model, rescaled_error):

        self.read_len = read_len

        model_path = pathlib.Path(error_model)
        try:
            error_dat = pickle.load(open(model_path, 'rb'), encoding="bytes")
        except IOError:
            print("\nProblem opening the sequencing error model.\n")
            sys.exit(1)

        self.uniform = False

        # uniform-error SE reads (e.g., PacBio)
        if len(error_dat) == 4:
            self.uniform = True
            [q_scores, off_q, avg_error, error_params] = error_dat
            self.uniform_q_score = int(-10. * np.log10(avg_error) + 0.5)
            print('Using uniform sequencing error model. (q=' + str(self.uniform_q_score) + '+' + str(
                off_q) + ', p(err)={0:0.2f}%)'.format(100. * avg_error))

        # only 1 q-score model present, use same model for both strands
        elif len(error_dat) == 6:
            [init_q1, prob_q1, q_scores, off_q, avg_error, error_params] = error_dat
            self.pe_models = False

        # found a q-score model for both forward and reverse strands
        elif len(error_dat) == 8:
            [init_q1, prob_q1, init_q2, prob_q2, q_scores, off_q, avg_error, error_params] = error_dat
            self.pe_models = True
            if len(init_q1) != len(init_q2) or len(prob_q1) != len(prob_q2):
                print('\nError: R1 and R2 quality score models are of different length.\n')
                sys.exit(1)

        # This serves as a sanity check for the input model
        else:
            print('\nError: Something wrong with error model.\n')
            sys.exit(1)

        self.q_err_rate = [0.] * (max(q_scores) + 1)
        for q in q_scores:
            self.q_err_rate[q] = 10. ** (-q / 10.)
        self.off_q = off_q
        self.err_p = error_params
        # Selects a new nucleotide based on the error model
        self.err_sse = [DiscreteDistribution(n, NUCL) for n in self.err_p[0]]
        # allows for selection of indel length based on the parameters of the model
        self.err_sie = DiscreteDistribution(self.err_p[2], self.err_p[3])
        # allows for indel insertion based on the length above and the probability from the model
        self.err_sin = DiscreteDistribution(self.err_p[5], NUCL)

        # adjust sequencing error frequency to match desired rate
        if rescaled_error is None:
            self.error_scale = 1.0
        else:
            self.error_scale = rescaled_error / avg_error
            print('Warning: Quality scores no longer exactly representative of error probability. '
                  'Error model scaled by {0:.3f} to match desired rate...'.format(self.error_scale))

        if not self.uniform:
            # adjust length to match desired read length
            if self.read_len == len(init_q1):
                self.q_ind_remap = range(self.read_len)
            else:
                print('Warning: Read length of error model (' + str(len(init_q1)) + ') does not match -R value (' + str(
                    self.read_len) + '), rescaling model...')
                self.q_ind_remap = [max([1, len(init_q1) * n // read_len]) for n in range(read_len)]

            # initialize probability distributions
            self.init_dist_by_pos_1 = [DiscreteDistribution(init_q1[i], q_scores) for i in range(len(init_q1))]
            self.prob_dist_by_pos_by_prev_q1 = [None]
            for i in range(1, len(init_q1)):
                self.prob_dist_by_pos_by_prev_q1.append([])
                for j in range(len(init_q1[0])):
                    # if we don't have sufficient data for a transition, use the previous quality score
                    if np.sum(prob_q1[i][j]) <= 0.:
                        self.prob_dist_by_pos_by_prev_q1[-1].append(
                            DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
                    else:
                        self.prob_dist_by_pos_by_prev_q1[-1].append(DiscreteDistribution(prob_q1[i][j], q_scores))

            # If paired-end, initialize probability distributions for the other strand
            if self.pe_models:
                self.init_dist_by_pos_2 = [DiscreteDistribution(init_q2[i], q_scores) for i in range(len(init_q2))]
                self.prob_dist_by_pos_by_prev_q2 = [None]
                for i in range(1, len(init_q2)):
                    self.prob_dist_by_pos_by_prev_q2.append([])
                    for j in range(len(init_q2[0])):
                        if np.sum(prob_q2[i][
                                      j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                            self.prob_dist_by_pos_by_prev_q2[-1].append(
                                DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
                        else:
                            self.prob_dist_by_pos_by_prev_q2[-1].append(DiscreteDistribution(prob_q2[i][j], q_scores))

    def get_sequencing_errors(self, read_data, is_reverse_strand=False):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.

        :param read_data: sequence to insert errors into
        :param is_reverse_strand: whether to treat this as the reverse strand or not
        :return: modified sequence and associate quality scores
        """

        q_out = [0] * self.read_len
        s_err = []

        if self.uniform:
            my_q = [self.uniform_q_score + self.off_q] * self.read_len
            q_out = ''.join([chr(n) for n in my_q])
            for i in range(self.read_len):
                if random.random() < self.error_scale * self.q_err_rate[self.uniform_q_score]:
                    s_err.append(i)
        else:
            if self.pe_models and is_reverse_strand:
                my_q = self.init_dist_by_pos_2[0].sample()
            else:
                my_q = self.init_dist_by_pos_1[0].sample()
            q_out[0] = my_q

            for i in range(1, self.read_len):
                if self.pe_models and is_reverse_strand:
                    my_q = self.prob_dist_by_pos_by_prev_q2[self.q_ind_remap[i]][my_q].sample()
                else:
                    my_q = self.prob_dist_by_pos_by_prev_q1[self.q_ind_remap[i]][my_q].sample()
                q_out[i] = my_q

            if is_reverse_strand:
                q_out = q_out[::-1]

            for i in range(self.read_len):
                if random.random() < self.error_scale * self.q_err_rate[q_out[i]]:
                    s_err.append(i)

            q_out = ''.join([chr(n + self.off_q) for n in q_out])

        if self.error_scale == 0.0:
            return q_out, []

        s_out = []
        n_del_so_far = 0
        # don't allow indel errors to occur on subsequent positions
        prev_indel = -2
        # don't allow other sequencing errors to occur on bases removed by deletion errors
        del_blacklist = []

        for ind in s_err[::-1]:  # for each error that we're going to insert...

            # determine error type
            is_sub = True
            if ind != 0 and ind != self.read_len - 1 - max(self.err_p[3]) and abs(ind - prev_indel) > 1:
                if random.random() < self.err_p[1]:
                    is_sub = False

            # insert substitution error
            if is_sub:
                my_nucl = str(read_data[ind])
                new_nucl = self.err_sse[NUC_IND[my_nucl]].sample()
                s_out.append(('S', 1, ind, my_nucl, new_nucl))

            # insert indel error
            else:
                indel_len = self.err_sie.sample()

                # insertion error
                if random.random() < self.err_p[4]:
                    my_nucl = str(read_data[ind])
                    new_nucl = my_nucl + ''.join([self.err_sin.sample() for n in range(indel_len)])
                    s_out.append(('I', len(new_nucl) - 1, ind, my_nucl, new_nucl))

                # deletion error (prevent too many of them from stacking up)
                elif ind < self.read_len - 2 - n_del_so_far:
                    my_nucl = str(read_data[ind:ind + indel_len + 1])
                    new_nucl = str(read_data[ind])
                    n_del_so_far += len(my_nucl) - 1
                    s_out.append(('D', len(my_nucl) - 1, ind, my_nucl, new_nucl))
                    for i in range(ind + 1, ind + indel_len + 1):
                        del_blacklist.append(i)
                prev_indel = ind

        # remove blacklisted errors
        for i in range(len(s_out) - 1, -1, -1):
            if s_out[i][2] in del_blacklist:
                del s_out[i]

        return q_out, s_out


# parse mutation model pickle file
def parse_input_mutation_model(model=None, which_default=1):
    if which_default == 1:
        out_model = [copy.deepcopy(n) for n in DEFAULT_MODEL_1]
    elif which_default == 2:
        out_model = [copy.deepcopy(n) for n in DEFAULT_MODEL_2]
    else:
        print('\nError: Unknown default mutation model specified\n')
        sys.exit(1)

    if model is not None:
        pickle_dict = pickle.load(open(model, "rb"))
        out_model[0] = pickle_dict['AVG_MUT_RATE']
        out_model[2] = 1. - pickle_dict['SNP_FREQ']

        ins_list = pickle_dict['INDEL_FREQ']
        if len(ins_list):
            ins_count = sum([ins_list[k] for k in ins_list.keys() if k >= 1])
            del_count = sum([ins_list[k] for k in ins_list.keys() if k <= -1])
            ins_vals = [k for k in sorted(ins_list.keys()) if k >= 1]
            ins_wght = [ins_list[k] / float(ins_count) for k in ins_vals]
            del_vals = [k for k in sorted([abs(k) for k in ins_list.keys() if k <= -1])]
            del_wght = [ins_list[-k] / float(del_count) for k in del_vals]
        else:  # degenerate case where no indel stats are provided
            ins_count = 1
            del_count = 1
            ins_vals = [1]
            ins_wght = [1.0]
            del_vals = [1]
            del_wght = [1.0]
        out_model[3] = ins_count / float(ins_count + del_count)
        out_model[4] = ins_vals
        out_model[5] = ins_wght
        out_model[6] = del_vals
        out_model[7] = del_wght

        trinuc_trans_prob = pickle_dict['TRINUC_TRANS_PROBS']
        for k in sorted(trinuc_trans_prob.keys()):
            my_ind = TRI_IND[k[0][0] + k[0][2]]
            (k1, k2) = (NUC_IND[k[0][1]], NUC_IND[k[1][1]])
            out_model[8][my_ind][k1][k2] = trinuc_trans_prob[k]
        for i in range(len(out_model[8])):
            for j in range(len(out_model[8][i])):
                for l in range(len(out_model[8][i][j])):
                    # if trinuc not present in input mutation model, assign it uniform probability
                    if float(sum(out_model[8][i][j])) < 1e-12:
                        out_model[8][i][j] = [0.25, 0.25, 0.25, 0.25]
                    else:
                        out_model[8][i][j][l] /= float(sum(out_model[8][i][j]))

        trinuc_mut_prob = pickle_dict['TRINUC_MUT_PROB']
        which_have_we_seen = {n: False for n in ALL_TRI}
        trinuc_mean = np.mean(list(trinuc_mut_prob.values()))
        for trinuc in trinuc_mut_prob.keys():
            out_model[9][ALL_IND[trinuc]] = trinuc_mut_prob[trinuc]
            which_have_we_seen[trinuc] = True
        for trinuc in which_have_we_seen.keys():
            if not which_have_we_seen[trinuc]:
                out_model[9][ALL_IND[trinuc]] = trinuc_mean

    return out_model
