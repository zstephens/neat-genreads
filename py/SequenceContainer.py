import random
import copy
import os
import bisect
import pickle
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from py.probability import DiscreteDistribution, poisson_list, quantize_list
from py.neat_cigar import CigarString

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
DEFAULT_1_TRI_FREQS = [copy.deepcopy(example_matrix_1) for n in range(16)]
DEFAULT_1_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for n in ALL_TRI]
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
DEFAULT_2_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_2_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_2 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_2_TRI_FREQS = [copy.deepcopy(example_matrix_2) for n in range(16)]
DEFAULT_2_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for n in ALL_TRI]
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


#
#	Container for reference sequences, applies mutations
#
class SequenceContainer:
    def __init__(self, x_offset, sequence, ploidy, window_overlap, read_len, mut_models=[], mut_rate=None,
                 only_vcf=False):
        # initialize basic variables
        self.only_vcf = only_vcf
        self.init_basic_vars(x_offset, sequence, ploidy, window_overlap, read_len)
        # initialize mutation models
        self.init_mut_models(mut_models, mut_rate)
        # sample the number of variants that will be inserted into each ploid
        self.init_poisson()
        self.indels_to_add = [n.sample() for n in self.ind_pois]
        self.snps_to_add = [n.sample() for n in self.snp_pois]
        # initialize trinuc snp bias
        self.init_trinuc_bias()

    def init_basic_vars(self, x_offset, sequence, ploidy, window_overlap, read_len):
        self.x = x_offset
        self.ploidy = ploidy
        self.read_len = read_len
        self.sequences = [Seq(str(sequence), IUPAC.unambiguous_dna).tomutable() for n in range(self.ploidy)]
        self.seq_len = len(sequence)
        self.indel_list = [[] for n in range(self.ploidy)]
        self.snp_list = [[] for n in range(self.ploidy)]
        self.all_cigar = [[] for n in range(self.ploidy)]
        self.fm_pos = [[] for n in range(self.ploidy)]
        self.fm_span = [[] for n in range(self.ploidy)]
        # blackList[ploid][pos] = 0		safe to insert variant here
        # blackList[ploid][pos] = 1		indel inserted here
        # blackList[ploid][pos] = 2		snp inserted here
        # blackList[ploid][pos] = 3		invalid position for various processing reasons
        self.black_list = [np.zeros(self.seq_len, dtype='<i4') for n in range(self.ploidy)]

        # disallow mutations to occur on window overlap points
        self.win_buffer = window_overlap
        for p in range(self.ploidy):
            self.black_list[p][-self.win_buffer] = 3
            self.black_list[p][-self.win_buffer - 1] = 3

    def init_coverage(self, coverage_data, frag_dist=None):
        # if we're only creating a vcf, skip some expensive initialization related to coverage depth
        if not self.only_vcf:
            (self.window_size, gc_scalars, target_cov_vals) = coverage_data
            gc_cov_vals = [[] for n in self.sequences]
            tr_cov_vals = [[] for n in self.sequences]
            self.coverage_distribution = []
            avg_out = []
            for i in range(len(self.sequences)):
                max_coord = min([len(self.sequences[i]) - self.read_len, len(self.all_cigar[i]) - self.read_len])
                # Trying to fix a problem wherein the above line gives a negative answer
                if max_coord <= 0:
                    max_coord = min([len(self.sequences[i]), len(self.all_cigar[i])])
                # compute gc-bias
                j = 0
                while j + self.window_size < len(self.sequences[i]):
                    gc_c = self.sequences[i][j:j + self.window_size].count('G') + self.sequences[i][
                                                                                  j:j + self.window_size].count('C')
                    gc_cov_vals[i].extend([gc_scalars[gc_c]] * self.window_size)
                    j += self.window_size
                gc_c = self.sequences[i][-self.window_size:].count('G') + self.sequences[i][-self.window_size:].count(
                    'C')
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
                    # print (i,j), self.adj[i][j], self.all_cigar[i][j], self.FM_pos[i][j], self.FM_span[i][j]
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

                    self.fraglens_ind_map = {}
                    for j in frag_dist.values:
                        bInd = bisect.bisect(flq, j)
                        if abs(flq[bInd - 1] - j) <= abs(flq[bInd] - j):
                            self.fraglens_ind_map[j] = flq[bInd - 1]
                        else:
                            self.fraglens_ind_map[j] = flq[bInd]

                    self.coverage_distribution.append({})
                    for flv in sorted(list(set(self.fraglens_ind_map.values()))):
                        buffer_val = self.read_len
                        for j in frag_dist.values:
                            if self.fraglens_ind_map[j] == flv and j > buffer_val:
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
                        # exit(1)

                        self.coverage_distribution[i][flv] = DiscreteDistribution(coverage_vals,
                                                                                  range(len(coverage_vals)))

            return np.mean(avg_out)

    def init_mut_models(self, mut_models, mut_rate):
        if not mut_models:
            default_model = [copy.deepcopy(DEFAULT_MODEL_1) for n in range(self.ploidy)]
            self.model_data = default_model[:self.ploidy]
        else:
            if len(mut_models) != self.ploidy:
                print('\nError: Number of mutation models recieved is not equal to specified ploidy\n')
                exit(1)
            self.model_data = copy.deepcopy(mut_models)

        # do we need to rescale mutation frequencies?
        mut_rateSum = sum([n[0] for n in self.model_data])
        self.mut_rescale = mut_rate
        if self.mut_rescale == None:
            self.mut_scalar = 1.0
        else:
            self.mut_scalar = float(self.mut_rescale) // (mut_rateSum / float(len(self.model_data)))

        # how are mutations spread to each ploid, based on their specified mut rates?
        self.ploid_mut_frac = [float(n[0]) / mut_rateSum for n in self.model_data]
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
                self.models[-1][6].append([DiscreteDistribution(m[0], NUCL), DiscreteDistribution(m[1], NUCL),
                                           DiscreteDistribution(m[2], NUCL), DiscreteDistribution(m[3], NUCL)])
            self.models[-1].append([m for m in n[9]])

    def init_poisson(self):
        ind_l_list = [self.seq_len * self.models[i][0] * self.models[i][2] * self.ploid_mut_frac[i] for i in
                      range(len(self.models))]
        snp_l_list = [self.seq_len * self.models[i][0] * (1. - self.models[i][2]) * self.ploid_mut_frac[i] for i in
                      range(len(self.models))]
        k_range = range(int(self.seq_len * MAX_MUTFRAC))
        self.ind_pois = [poisson_list(k_range, ind_l_list[n]) for n in range(len(self.models))]
        self.snp_pois = [poisson_list(k_range, snp_l_list[n]) for n in range(len(self.models))]

    def init_trinuc_bias(self):
        # compute mutation positional bias given trinucleotide strings of the sequence (ONLY AFFECTS SNPs)
        #
        # note: since indels are added before snps, it's possible these positional biases aren't correctly utilized
        #       at positions affected by indels. At the moment I'm going to consider this negligible.
        trinuc_snp_bias = [[0. for n in range(self.seq_len)] for m in range(self.ploidy)]
        self.trinuc_bias = [None for n in range(self.ploidy)]
        for p in range(self.ploidy):
            for i in range(self.win_buffer + 1, self.seq_len - 1):
                trinuc_snp_bias[p][i] = self.models[p][7][ALL_IND[str(self.sequences[p][i - 1:i + 2])]]
            self.trinuc_bias[p] = DiscreteDistribution(trinuc_snp_bias[p][self.win_buffer + 1:self.seq_len - 1],
                                                       range(self.win_buffer + 1, self.seq_len - 1))

    def update(self, x_offset, sequence, ploidy, window_overlap, read_len, mut_models=[], mut_rate=None):
        # if mutation model is changed, we have to reinitialize it...
        if ploidy != self.ploidy or mut_rate != self.mut_rescale or mut_models != []:
            self.ploidy = ploidy
            self.mut_rescale = mut_rate
            self.init_mut_models(mut_models, mut_rate)
        # if sequence length is different than previous window, we have to redo snp/indel poissons
        if len(sequence) != self.seq_len:
            self.seq_len = len(sequence)
            self.init_poisson()
        # basic vars
        self.init_basic_vars(x_offset, sequence, ploidy, window_overlap, read_len)
        self.indels_to_add = [n.sample() for n in self.ind_pois]
        self.snps_to_add = [n.sample() for n in self.snp_pois]
        # initialize trinuc snp bias
        if not IGNORE_TRINUC:
            self.init_trinuc_bias()

    def insert_mutations(self, input_list):
        for inpV in input_list:
            which_ploid = []
            wps = inpV[4][0]
            if wps is None:  # if no genotype given, assume heterozygous and choose a single ploid based on their mut rates
                which_ploid.append(self.ploid_mut_prior.sample())
                which_alt = [0]
            else:
                if '/' in wps or '|' in wps:
                    if '/' in wps:
                        splt = wps.split('/')
                    else:
                        splt = wps.split('|')
                    which_ploid = []
                    which_alt = []
                    for i in range(len(splt)):
                        if splt[i] == '1':
                            which_ploid.append(i)
                    # assume we're just using first alt for inserted variants?
                    which_alt = [0 for n in which_ploid]
                else:  # otherwise assume monoploidy
                    which_ploid = [0]
                    which_alt = [0]

            # ignore invalid ploids
            for i in range(len(which_ploid) - 1, -1, -1):
                if which_ploid[i] >= self.ploidy:
                    del which_ploid[i]

            for i in range(len(which_ploid)):
                p = which_ploid[i]
                myAlt = inpV[2][which_alt[i]]
                myVar = (inpV[0] - self.x, inpV[1], myAlt)
                inLen = max([len(inpV[1]), len(myAlt)])

                if myVar[0] < 0 or myVar[0] >= len(self.black_list[p]):
                    print('\nError: Attempting to insert variant out of window bounds:')
                    print(myVar, '--> blackList[0:' + str(len(self.black_list[p])) + ']\n')
                    exit(1)
                if len(inpV[1]) == 1 and len(myAlt) == 1:
                    if self.black_list[p][myVar[0]]:
                        continue
                    self.snp_list[p].append(myVar)
                    self.black_list[p][myVar[0]] = 2
                else:
                    indel_failed = False
                    for k in range(myVar[0], myVar[0] + inLen + 1):
                        if k >= len(self.black_list[p]):
                            indel_failed = True
                            continue
                        if self.black_list[p][k]:
                            indel_failed = True
                            continue
                    if indel_failed:
                        continue
                    for k in range(myVar[0], myVar[0] + inLen + 1):
                        self.black_list[p][k] = 1
                    self.indel_list[p].append(myVar)

    def random_mutations(self):

        #	add random indels
        all_indels = [[] for n in self.sequences]
        for i in range(self.ploidy):
            for j in range(self.indels_to_add[i]):
                if random.random() <= self.models[i][1]:  # insert homozygous indel
                    whichPloid = range(self.ploidy)
                else:  # insert heterozygous indel
                    whichPloid = [self.ploid_mut_prior.sample()]

                # try to find suitable places to insert indels
                eventPos = -1
                for attempt in range(MAX_ATTEMPTS):
                    eventPos = random.randint(self.win_buffer, self.seq_len - 1)
                    for p in whichPloid:
                        if self.black_list[p][eventPos]:
                            eventPos = -1
                    if eventPos != -1:
                        break
                if eventPos == -1:
                    continue

                if random.random() <= self.models[i][3]:  # insertion
                    inLen = self.models[i][4].sample()
                    # sequence content of random insertions is uniformly random (change this later, maybe)
                    inSeq = ''.join([random.choice(NUCL) for n in range(inLen)])
                    ref_nucl = self.sequences[i][eventPos]
                    myIndel = (eventPos, ref_nucl, ref_nucl + inSeq)
                else:  # deletion
                    inLen = self.models[i][5].sample()
                    if eventPos + inLen + 1 >= len(self.sequences[i]):  # skip if deletion too close to boundary
                        continue
                    if inLen == 1:
                        inSeq = self.sequences[i][eventPos + 1]
                    else:
                        inSeq = str(self.sequences[i][eventPos + 1:eventPos + inLen + 1])
                    ref_nucl = self.sequences[i][eventPos]
                    myIndel = (eventPos, ref_nucl + inSeq, ref_nucl)

                # if event too close to boundary, skip. if event conflicts with other indel, skip.
                skipEvent = False
                if eventPos + len(myIndel[1]) >= self.seq_len - self.win_buffer - 1:
                    skipEvent = True
                if skipEvent:
                    continue
                for p in whichPloid:
                    for k in range(eventPos, eventPos + inLen + 1):
                        if self.black_list[p][k]:
                            skipEvent = True
                if skipEvent:
                    continue

                for p in whichPloid:
                    for k in range(eventPos, eventPos + inLen + 1):
                        self.black_list[p][k] = 1
                    all_indels[p].append(myIndel)

        #	add random snps
        all_snps = [[] for n in self.sequences]
        for i in range(self.ploidy):
            for j in range(self.snps_to_add[i]):
                if random.random() <= self.models[i][1]:  # insert homozygous SNP
                    whichPloid = range(self.ploidy)
                else:  # insert heterozygous SNP
                    whichPloid = [self.ploid_mut_prior.sample()]

                # try to find suitable places to insert snps
                eventPos = -1
                for attempt in range(MAX_ATTEMPTS):
                    # based on the mutation model for the specified ploid, choose a SNP location based on trinuc bias
                    # (if there are multiple ploids, choose one at random)
                    if IGNORE_TRINUC:
                        eventPos = random.randint(self.win_buffer + 1, self.seq_len - 2)
                    else:
                        ploid_to_use = whichPloid[random.randint(0, len(whichPloid) - 1)]
                        eventPos = self.trinuc_bias[ploid_to_use].sample()
                    for p in whichPloid:
                        if self.black_list[p][eventPos]:
                            eventPos = -1
                    if eventPos != -1:
                        break
                if eventPos == -1:
                    continue

                ref_nucl = self.sequences[i][eventPos]
                context = str(self.sequences[i][eventPos - 1]) + str(self.sequences[i][eventPos + 1])
                # sample from tri-nucleotide substitution matrices to get SNP alt allele
                new_nucl = self.models[i][6][TRI_IND[context]][NUC_IND[ref_nucl]].sample()
                mySNP = (eventPos, ref_nucl, new_nucl)

                for p in whichPloid:
                    all_snps[p].append(mySNP)
                    self.black_list[p][mySNP[0]] = 2

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
                    exit(1)
                else:
                    self.sequences[i][v_pos] = all_snps[i][j][2]

        # organize the indels we want to insert
        for i in range(len(all_indels)):
            all_indels[i].extend(self.indel_list[i])
        all_indels_ins = [sorted([list(m) for m in n]) for n in all_indels]

        # MODIFY REFERENCE STRING: INDELS
        adj_to_add = [[] for n in range(self.ploidy)]
        for i in range(len(all_indels_ins)):
            rolling_adj = 0
            temp_symbol_string = ['M' for n in self.sequences[i]]
            # there's an off-by-one error somewhere in the position sampling routines.. this might fix it
            # temp_symbol_string.append('M')
            for j in range(len(all_indels_ins[i])):
                v_pos = all_indels_ins[i][j][0] + rolling_adj
                v_pos2 = v_pos + len(all_indels_ins[i][j][1])
                rolling_adj += len(all_indels_ins[i][j][2]) - len(all_indels_ins[i][j][1])

                if all_indels_ins[i][j][1] != str(self.sequences[i][v_pos:v_pos2]):
                    print('\nError: Something went wrong!\n', all_indels_ins[i][j], [v_pos, v_pos2],
                          str(self.sequences[i][v_pos:v_pos2]), '\n')
                    exit(1)
                else:
                    # alter reference sequence
                    self.sequences[i] = self.sequences[i][:v_pos] + Seq(all_indels_ins[i][j][2],
                                                                        IUPAC.unambiguous_dna).tomutable() + \
                                        self.sequences[i][v_pos2:]
                    # notate indel positions for cigar computation
                    d = len(all_indels_ins[i][j][2]) - len(all_indels_ins[i][j][1])
                    if d > 0:
                        temp_symbol_string = temp_symbol_string[:v_pos + 1] + ['I'] * d + temp_symbol_string[
                                                                                          v_pos2 + 1:]
                    elif d < 0:
                        temp_symbol_string[v_pos + 1] = 'D' * abs(d) + 'M'

            # precompute cigar strings
            for j in range(len(temp_symbol_string) - self.read_len):
                self.all_cigar[i].append(CigarString(listIn=temp_symbol_string[j:j + self.read_len]).get_string())

            # create some data structures we will need later:
            # --- self.FM_pos[ploid][pos]: position of the left-most matching base (IN REFERENCE COORDINATES, i.e. corresponding to the unmodified reference genome)
            # --- self.FM_span[ploid][pos]: number of reference positions spanned by a read originating from this coordinate
            MD_so_far = 0
            for j in range(len(temp_symbol_string)):
                self.fm_pos[i].append(MD_so_far)
                # fix an edge case with deletions
                if 'D' in temp_symbol_string[j]:
                    self.fm_pos[i][-1] += temp_symbol_string[j].count('D')
                # compute number of ref matches for each read
                span_dif = len([n for n in temp_symbol_string[j:j + self.read_len] if 'M' in n])
                self.fm_span[i].append(self.fm_pos[i][-1] + span_dif)
                MD_so_far += temp_symbol_string[j].count('M') + temp_symbol_string[j].count('D')

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

        #
        #	TODO: combine multiple variants that happened to occur at same position into single vcf entry?
        #

        output_variants = []
        for k in sorted(count_dict.keys()):
            output_variants.append(k + tuple([len(count_dict[k]) / float(self.ploidy)]))
            ploid_string = ['0' for n in range(self.ploidy)]
            for k2 in [n for n in count_dict[k]]:
                ploid_string[k2] = '1'
            output_variants[-1] += tuple(['WP=' + '/'.join(ploid_string)])
        return output_variants

    def sample_read(self, sequencingModel, frag_len=None):

        # choose a ploid
        my_ploid = random.randint(0, self.ploidy - 1)

        # stop attempting to find a valid position if we fail enough times
        MAX_READPOS_ATTEMPTS = 100
        attempts_thus_far = 0

        # choose a random position within the ploid, and generate quality scores / sequencing errors
        readsToSample = []
        if frag_len == None:
            r_pos = self.coverage_distribution[my_ploid].sample()

            # sample read position and call function to compute quality scores / sequencing errors
            r_dat = self.sequences[my_ploid][r_pos:r_pos + self.read_len]
            (myQual, myErrors) = sequencingModel.get_sequencing_errors(r_dat)
            readsToSample.append([r_pos, myQual, myErrors, r_dat])

        else:
            r_pos1 = self.coverage_distribution[my_ploid][self.fraglens_ind_map[frag_len]].sample()

            # EXPERIMENTAL
            # coords_to_select_from = self.coverage_distribution[my_ploid][self.fraglens_ind_map[frag_len]].sample()
            # r_pos1 = random.randint(coords_to_select_from[0],coords_to_select_from[1])

            r_pos2 = r_pos1 + frag_len - self.read_len
            r_dat1 = self.sequences[my_ploid][r_pos1:r_pos1 + self.read_len]
            r_dat2 = self.sequences[my_ploid][r_pos2:r_pos2 + self.read_len]
            (myQual1, myErrors1) = sequencingModel.get_sequencing_errors(r_dat1)
            (myQual2, myErrors2) = sequencingModel.get_sequencing_errors(r_dat2, is_reverse_strand=True)
            readsToSample.append([r_pos1, myQual1, myErrors1, r_dat1])
            readsToSample.append([r_pos2, myQual2, myErrors2, r_dat2])

        # error format:
        # myError[i] = (type, len, pos, ref, alt)

        # examine sequencing errors to-be-inserted.
        #	- remove deletions that don't have enough bordering sequence content to "fill in"
        # if error is valid, make the changes to the read data
        r_out = []
        for r in readsToSample:
            try:
                myCigar = self.all_cigar[my_ploid][r[0]]
            except IndexError:
                print('Index error when attempting to find cigar string.')
                print(my_ploid, len(self.all_cigar[my_ploid]), r[0])
                if frag_len is not None:
                    print((r_pos1, r_pos2))
                    print(frag_len, self.fraglens_ind_map[frag_len])
                exit(1)
            total_d = sum([error[1] for error in r[2] if error[0] == 'D'])
            total_i = sum([error[1] for error in r[2] if error[0] == 'I'])
            availB = len(self.sequences[my_ploid]) - r[0] - self.read_len - 1
            # add buffer sequence to fill in positions that get deleted
            r[3] += self.sequences[my_ploid][r[0] + self.read_len:r[0] + self.read_len + total_d]
            expanded_cigar = []
            extraCigar = []
            adj = 0
            sse_adj = [0 for n in range(self.read_len + max(sequencingModel.errP[3]))]
            any_indel_err = False

            # sort by letter (D > I > S) such that we introduce all indel errors before substitution errors
            # secondarily, sort by index
            arrangedErrors = {'D': [], 'I': [], 'S': []}
            for error in r[2]:
                arrangedErrors[error[0]].append((error[2], error))
            sortedErrors = []
            for k in sorted(arrangedErrors.keys()):
                sortedErrors.extend([n[1] for n in sorted(arrangedErrors[k])])

            skip_indels = False

            # FIXED TdB 05JUN2018
            # Moved this outside the for error loop, since it messes up the CIGAR string when more than one deletion is in the same read
            extra_cigar_val = []
            # END FIXED TdB

            for error in sortedErrors:
                e_len = error[1]
                e_pos = error[2]
                if error[0] == 'D' or error[0] == 'I':
                    any_indel_err = True

                    # FIXED TdB 05JUN2018
                    # Moved this OUTSIDE the for error loop, since it messes up the CIGAR string when more than one deletion is in the same read
                    # extra_cigar_val = []
                    # END FIXED TdB

                    if total_d > availB:  # if not enough bases to fill-in deletions, skip all indel erors
                        continue
                    if expanded_cigar == []:
                        expanded_cigar = CigarString(stringIn=myCigar).get_list()
                        fill_to_go = total_d - total_i + 1
                        if fill_to_go > 0:
                            try:
                                extra_cigar_val = CigarString(
                                    stringIn=self.all_cigar[my_ploid][r[0] + fill_to_go]).get_list()[-fill_to_go:]
                            except IndexError:  # applying the deletions we want requires going beyond region boundaries. skip all indel errors
                                skip_indels = True

                    if skip_indels:
                        continue

                    # insert deletion error into read and update cigar string accordingly
                    if error[0] == 'D':
                        my_adj = sse_adj[e_pos]
                        pi = e_pos + my_adj
                        pf = e_pos + my_adj + e_len + 1
                        if str(r[3][pi:pf]) == str(error[3]):
                            r[3] = r[3][:pi + 1] + r[3][pf:]
                            expanded_cigar = expanded_cigar[:pi + 1] + expanded_cigar[pf:]
                            if pi + 1 == len(
                                    expanded_cigar):  # weird edge case with del at very end of region. Make a guess and add a "M"
                                expanded_cigar.append('M')
                            expanded_cigar[pi + 1] = 'D' * e_len + expanded_cigar[pi + 1]
                        else:
                            print('\nError, ref does not match alt while attempting to insert deletion error!\n')
                            exit(1)
                        adj -= e_len
                        for i in range(e_pos, len(sse_adj)):
                            sse_adj[i] -= e_len

                    # insert insertion error into read and update cigar string accordingly
                    else:
                        my_adj = sse_adj[e_pos]
                        if str(r[3][e_pos + my_adj]) == error[3]:
                            r[3] = r[3][:e_pos + my_adj] + error[4] + r[3][e_pos + my_adj + 1:]
                            expanded_cigar = expanded_cigar[:e_pos + my_adj] + ['I'] * e_len + expanded_cigar[
                                                                                               e_pos + my_adj:]
                        else:
                            print('\nError, ref does not match alt while attempting to insert insertion error!\n')
                            print('---', chr(r[3][e_pos + my_adj]), '!=', error[3])
                            exit(1)
                        adj += e_len
                        for i in range(e_pos, len(sse_adj)):
                            sse_adj[i] += e_len

                else:  # substitution errors, much easier by comparison...
                    if str(r[3][e_pos + sse_adj[e_pos]]) == error[3]:
                        r[3][e_pos + sse_adj[e_pos]] = error[4]
                    else:
                        print('\nError, ref does not match alt while attempting to insert substitution error!\n')
                        exit(1)

            if any_indel_err:
                if len(expanded_cigar):
                    relevantCigar = (expanded_cigar + extra_cigar_val)[:self.read_len]
                    myCigar = CigarString(listIn=relevantCigar).get_string()

                r[3] = r[3][:self.read_len]

            r_out.append([self.fm_pos[my_ploid][r[0]], myCigar, str(r[3]), str(r[1])])

        # r_out[i] = (pos, cigar, read_string, qual_string)
        return r_out


class ReadContainer:
    """
    Container for read data, computes quality scores and positions to insert errors
    """
    def __init__(self, read_len, error_model, rescaled_error):

        self.read_len = read_len

        error_dat = pickle.load(open(error_model, 'rb'), encoding="bytes")
        self.UNIFORM = False
        if len(error_dat) == 4:  # uniform-error SE reads (e.g. PacBio)
            self.UNIFORM = True
            [q_scores, off_q, avg_error, error_params] = error_dat
            self.uniform_q_score = int(-10. * np.log10(avg_error) + 0.5)
            print('Using uniform sequencing error model. (q=' + str(self.uniform_q_score) + '+' + str(
                off_q) + ', p(err)={0:0.2f}%)'.format(100. * avg_error))
        elif len(error_dat) == 6:  # only 1 q-score model present, use same model for both strands
            [init_q1, prob_q1, q_scores, off_q, avg_error, error_params] = error_dat
            self.PE_MODELS = False
        elif len(error_dat) == 8:  # found a q-score model for both forward and reverse strands
            [init_q1, prob_q1, init_q2, prob_q2, q_scores, off_q, avg_error, error_params] = error_dat
            self.PE_MODELS = True
            if len(init_q1) != len(init_q2) or len(prob_q1) != len(prob_q2):
                print('\nError: R1 and R2 quality score models are of different length.\n')
                exit(1)
        else:
            print('\nError: Something wrong with error model.\n')
            exit(1)
        print(q_scores)
        self.q_err_rate = [0.] * (max(q_scores) + 1)
        for q in q_scores:
            self.q_err_rate[q] = 10. ** (-q / 10.)
        self.off_q = off_q

        # error_params = [SSE_PROB, SIE_RATE, SIE_PROB, SIE_VAL, SIE_INS_FREQ, SIE_INS_NUCL]
        self.errP = error_params
        self.errSSE = [DiscreteDistribution(n, NUCL) for n in self.errP[0]]
        self.errSIE = DiscreteDistribution(self.errP[2], self.errP[3])
        self.errSIN = DiscreteDistribution(self.errP[5], NUCL)

        # adjust sequencing error frequency to match desired rate
        if rescaled_error is None:
            self.errorScale = 1.0
        else:
            self.errorScale = rescaled_error / avg_error
            print('Warning: Quality scores no longer exactly representative of error probability. '
                  'Error model scaled by {0:.3f} to match desired rate...'.format(self.errorScale))

        if not self.UNIFORM:
            # adjust length to match desired read length
            if self.read_len == len(init_q1):
                self.qIndRemap = range(self.read_len)
            else:
                print('Warning: Read length of error model (' + str(len(init_q1)) + ') does not match -R value (' + str(
                    self.read_len) + '), rescaling model...')
                self.qIndRemap = [max([1, len(init_q1) * n // read_len]) for n in range(read_len)]

            # initialize probability distributions
            self.initDistByPos1 = [DiscreteDistribution(init_q1[i], q_scores) for i in range(len(init_q1))]
            self.probDistByPosByPrevQ1 = [None]
            for i in range(1, len(init_q1)):
                self.probDistByPosByPrevQ1.append([])
                for j in range(len(init_q1[0])):
                    if np.sum(prob_q1[i][
                                  j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                        self.probDistByPosByPrevQ1[-1].append(
                            DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
                    else:
                        self.probDistByPosByPrevQ1[-1].append(DiscreteDistribution(prob_q1[i][j], q_scores))

            if self.PE_MODELS:
                self.initDistByPos2 = [DiscreteDistribution(init_q2[i], q_scores) for i in range(len(init_q2))]
                self.probDistByPosByPrevQ2 = [None]
                for i in range(1, len(init_q2)):
                    self.probDistByPosByPrevQ2.append([])
                    for j in range(len(init_q2[0])):
                        if np.sum(prob_q2[i][
                                      j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                            self.probDistByPosByPrevQ2[-1].append(
                                DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
                        else:
                            self.probDistByPosByPrevQ2[-1].append(DiscreteDistribution(prob_q2[i][j], q_scores))

    def get_sequencing_errors(self, read_data, is_reverse_strand=False):

        q_out = [0] * self.read_len
        s_err = []

        if self.UNIFORM:
            my_q = [self.uniform_q_score + self.off_q for n in range(self.read_len)]
            q_out = ''.join([chr(n) for n in my_q])
            for i in range(self.read_len):
                if random.random() < self.errorScale * self.q_err_rate[self.uniform_q_score]:
                    s_err.append(i)
        else:

            if self.PE_MODELS and is_reverse_strand:
                my_q = self.initDistByPos2[0].sample()
            else:
                my_q = self.initDistByPos1[0].sample()
            q_out[0] = my_q

            for i in range(1, self.read_len):
                if self.PE_MODELS and is_reverse_strand:
                    my_q = self.probDistByPosByPrevQ2[self.qIndRemap[i]][my_q].sample()
                else:
                    my_q = self.probDistByPosByPrevQ1[self.qIndRemap[i]][my_q].sample()
                q_out[i] = my_q

            if is_reverse_strand:
                q_out = q_out[::-1]

            for i in range(self.read_len):
                if random.random() < self.errorScale * self.q_err_rate[q_out[i]]:
                    s_err.append(i)

            q_out = ''.join([chr(n + self.off_q) for n in q_out])

        if self.errorScale == 0.0:
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
            if ind != 0 and ind != self.read_len - 1 - max(self.errP[3]) and abs(ind - prev_indel) > 1:
                if random.random() < self.errP[1]:
                    is_sub = False

            # error_out = (type, len, pos, ref, alt)

            if is_sub:  # insert substitution error
                my_nucl = str(read_data[ind])
                new_nucl = self.errSSE[NUC_IND[my_nucl]].sample()
                s_out.append(('S', 1, ind, my_nucl, new_nucl))
            else:  # insert indel error
                indel_len = self.errSIE.sample()
                if random.random() < self.errP[4]:  # insertion error
                    my_nucl = str(read_data[ind])
                    new_nucl = my_nucl + ''.join([self.errSIN.sample() for n in range(indel_len)])
                    s_out.append(('I', len(new_nucl) - 1, ind, my_nucl, new_nucl))
                elif ind < self.read_len - 2 - n_del_so_far:  # deletion error (prevent too many of them from stacking up)
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
        exit(1)

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

