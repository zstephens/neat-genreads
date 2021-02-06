#!/usr/bin/env source

#
#	a quick script for comparing mutation models
#
#	source plotMutModel.source -i model1.p [model2.p] [model3.p]... -l legend_label1 [legend_label2] [legend_label3]... -o path/to/pdf_plot_prefix
#
# Python 3 ready

import sys
import pickle
import bisect
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import argparse

# mpl.rc('text',usetex=True)
# mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

parser = argparse.ArgumentParser(description='Plot and compare mutation models from gen_mut_model.source Usage: '
                                             'source plotMutModel.source -i model1.p [model2.p] [model3.p]... '
                                             '-l legend_label1 [legend_label2] [legend_label3]... '
                                             '-o path/to/pdf_plot_prefix',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
parser.add_argument('-i', type=str, required=True, metavar='<str>', nargs='+',
                    help="* mutation_model_1.p [mutation_model_2.p] [mutation_model_3] ...")
parser.add_argument('-l', type=str, required=True, metavar='<str>', nargs='+',
                    help="* legend labels: model1_name [model2_name] [model3_name]...")
parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output pdf prefix")
args = parser.parse_args()


def get_color(i, N, colormap='jet'):
    cm = mpl.get_cmap(colormap)
    c_norm = colors.Normalize(vmin=0, vmax=N + 1)
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cm)
    color_val = scalar_map.to_rgba(i)
    return color_val


def is_in_bed(track, ind):
    if ind in track:
        return True
    elif bisect.bisect(track, ind) % 1 == 1:
        return True
    else:
        return False


def get_bed_overlap(track, ind_s, ind_e):
    if ind_s in track:
        my_ind = track.index(ind_s)
        return min([track[my_ind + 1] - ind_s + 1, ind_e - ind_s + 1])
    else:
        my_ind = bisect.bisect(track, ind_s)
        if my_ind % 1 and my_ind < len(track) - 1:
            return min([track[my_ind + 1] - ind_s + 1, ind_e - ind_s + 1])
    return 0


# a waaaaaaay slower version of the above function ^^
# def getTrackOverlap(track1,track2):
#	otrack = [0 for n in range(max(track1+track2)+1)]
#	for i in range(0,len(track1),2):
#		for j in range(track1[i],track1[i+1]+1):
#			otrack[j] = 1
#	ocount = 0
#	for i in range(0,len(track2),2):
#		for j in range(track2[i],track2[i+1]+1):
#			if otrack[j]:
#				ocount += 1
#	return ocount

OUP = args.o
LAB = args.l
# print LAB
INP = args.i

N_FILES = len(INP)

mpl.rcParams.update({'font.size': 13, 'font.weight': 'bold', 'lines.linewidth': 3})

#################################################
#
#	BASIC STATS
#
#################################################
mpl.figure(0, figsize=(12, 10))

mpl.subplot(2, 2, 1)
color_ind = 0
for fn in INP:
    my_col = get_color(color_ind, N_FILES)
    color_ind += 1
    DATA_DICT = pickle.load(open(fn, "rb"), encoding="utf-8")
    [AVG_MUT_RATE, SNP_FREQ, INDEL_FREQ] = [DATA_DICT['AVG_MUT_RATE'], DATA_DICT['SNP_FREQ'], DATA_DICT['INDEL_FREQ']]
    mpl.bar([color_ind - 1], [AVG_MUT_RATE], 1., color=my_col)
mpl.xlim([-1, N_FILES + 1])
mpl.grid()
mpl.xticks([], [])
mpl.ylabel('Frequency')
mpl.title('Overall mutation rate (1/bp)')

mpl.subplot(2, 2, 2)
color_ind = 0
for fn in INP:
    my_col = get_color(color_ind, N_FILES)
    color_ind += 1
    DATA_DICT = pickle.load(open(fn, "rb"), encoding='utf-8')
    [AVG_MUT_RATE, SNP_FREQ, INDEL_FREQ] = [DATA_DICT['AVG_MUT_RATE'], DATA_DICT['SNP_FREQ'], DATA_DICT['INDEL_FREQ']]
    mpl.bar([color_ind - 1], [SNP_FREQ], 1., color=my_col)
    mpl.bar([color_ind - 1], [1. - SNP_FREQ], 1., color=my_col, bottom=[SNP_FREQ], hatch='/')
mpl.axis([-1, N_FILES + 1, 0, 1.2])
mpl.grid()
mpl.xticks([], [])
mpl.yticks([0, .2, .4, .6, .8, 1.], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
mpl.ylabel('Frequency')
mpl.title('SNP freq [  ] & indel freq [//]')

mpl.subplot(2, 1, 2)
color_ind = 0
leg_text = LAB
for fn in INP:
    my_col = get_color(color_ind, N_FILES)
    color_ind += 1
    DATA_DICT = pickle.load(open(fn, "rb"))
    [AVG_MUT_RATE, SNP_FREQ, INDEL_FREQ] = [DATA_DICT['AVG_MUT_RATE'], DATA_DICT['SNP_FREQ'], DATA_DICT['INDEL_FREQ']]
    x = sorted(INDEL_FREQ.keys())
    y = [INDEL_FREQ[n] for n in x]
    mpl.plot(x, y, color=my_col)
# leg_text.append(fn)
mpl.grid()
mpl.xlabel('Indel size (bp)', fontweight='bold')
mpl.ylabel('Frequency')
mpl.title('Indel frequency by size (- deletion, + insertion)')
mpl.legend(leg_text)
# mpl.show()
mpl.savefig(OUP + '_plot1_mutRates.pdf')

#################################################
#
#	TRINUC PRIOR PROB
#
#################################################
mpl.figure(1, figsize=(14, 6))
color_ind = 0
leg_text = LAB
for fn in INP:
    my_col = get_color(color_ind, N_FILES)
    color_ind += 1
    DATA_DICT = pickle.load(open(fn, "rb"))
    TRINUC_MUT_PROB = DATA_DICT['TRINUC_MUT_PROB']

    x = range(color_ind - 1, len(TRINUC_MUT_PROB) * N_FILES, N_FILES)
    xt = sorted(TRINUC_MUT_PROB.keys())
    y = [TRINUC_MUT_PROB[k] for k in xt]
    markerline, stemlines, baseline = mpl.stem(x, y, '-.')
    mpl.setp(markerline, 'markerfacecolor', my_col)
    mpl.setp(markerline, 'markeredgecolor', my_col)
    mpl.setp(baseline, 'color', my_col, 'linewidth', 0)
    mpl.setp(stemlines, 'color', my_col, 'linewidth', 3)
    if color_ind == 1:
        mpl.xticks(x, xt, rotation=90)
# leg_text.append(fn)
mpl.grid()
mpl.ylabel('p(trinucleotide mutates)')
mpl.legend(leg_text)
# mpl.show()
mpl.savefig(OUP + '_plot2_trinucPriors.pdf')

#################################################
#
#	TRINUC TRANS PROB
#
#################################################
plot_num = 3
for fn in INP:
    fig = mpl.figure(plot_num, figsize=(12, 10))
    DATA_DICT = pickle.load(open(fn, "rb"))
    TRINUC_TRANS_PROBS = DATA_DICT['TRINUC_TRANS_PROBS']

    xt2 = [m[3] for m in sorted([(n[0], n[2], n[1], n) for n in xt])]
    reverse_dict = {xt2[i]: i for i in range(len(xt2))}
    Z = np.zeros((64, 64))
    L = [['' for n in range(64)] for m in range(64)]
    for k in TRINUC_TRANS_PROBS:
        i = reverse_dict[k[0]]
        j = reverse_dict[k[1]]
        Z[i][j] = TRINUC_TRANS_PROBS[k]

    HARDCODED_LABEL = ['A_A', 'A_C', 'A_G', 'A_T',
                       'C_A', 'C_C', 'C_G', 'C_T',
                       'G_A', 'G_C', 'G_G', 'G_T',
                       'T_A', 'T_C', 'T_G', 'T_T']

    for pi in range(16):
        mpl.subplot(4, 4, pi + 1)
        Z2 = Z[pi * 4:(pi + 1) * 4, pi * 4:(pi + 1) * 4]
        X, Y = np.meshgrid(range(0, len(Z2[0]) + 1), range(0, len(Z2) + 1))
        im = mpl.pcolormesh(X, Y, Z2[::-1, :], vmin=0.0, vmax=0.5)
        mpl.axis([0, 4, 0, 4])
        mpl.xticks([0.5, 1.5, 2.5, 3.5], ['A', 'C', 'G', 'T'])
        mpl.yticks([0.5, 1.5, 2.5, 3.5], ['T', 'G', 'C', 'A'])
        mpl.text(1.6, 1.8, HARDCODED_LABEL[pi], color='white')

    # colorbar haxx
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label(r"p(X$Y_1$Z->X$Y_2$Z | X_Z mutates)")

    # mpl.tight_layout()
    # mpl.figtext(0.24,0.94,'Trinucleotide Mutation Frequency',size=20)
    # mpl.show()
    mpl.savefig(OUP + '_plot' + str(plot_num) + '_trinucTrans.pdf')
    plot_num += 1

#################################################
#
#	HIGH MUT REGIONS
#
#################################################
track_byFile_byChr = [{} for n in INP]
bp_total_byFile = [0 for n in INP]
color_ind = 0
for fn in INP:
    DATA_DICT = pickle.load(open(fn, "rb"))
    HIGH_MUT_REGIONS = DATA_DICT['HIGH_MUT_REGIONS']
    for region in HIGH_MUT_REGIONS:
        if region[0] not in track_byFile_byChr[color_ind]:
            track_byFile_byChr[color_ind][region[0]] = []
        track_byFile_byChr[color_ind][region[0]].extend([region[1], region[2]])
        bp_total_byFile[color_ind] += region[2] - region[1] + 1
    color_ind += 1

bp_overlap_count = [[0 for m in INP] for n in INP]
for i in range(N_FILES):
    bp_overlap_count[i][i] = bp_total_byFile[i]
    for j in range(i + 1, N_FILES):
        for k in track_byFile_byChr[i].keys():
            if k in track_byFile_byChr[j]:
                for ii in range(len(track_byFile_byChr[i][k][::2])):
                    bp_overlap_count[i][j] += get_bed_overlap(track_byFile_byChr[j][k], track_byFile_byChr[i][k][ii * 2],
                                                            track_byFile_byChr[i][k][ii * 2 + 1])

print('')
print('HIGH_MUT_REGION OVERLAP BETWEEN ' + str(N_FILES) + ' MODELS...')
for i in range(N_FILES):
    for j in range(i, N_FILES):
        n_dissimilar = (bp_overlap_count[i][i] - bp_overlap_count[i][j]) + (
                    bp_overlap_count[j][j] - bp_overlap_count[i][j])
        if bp_overlap_count[i][j] == 0:
            percentage_v = 0.0
        else:
            percentage_v = bp_overlap_count[i][j] / float(bp_overlap_count[i][j] + n_dissimilar)
    print('overlap[' + str(i) + ',' + str(j) + '] = ' + str(bp_overlap_count[i][j]) + ' bp ({0:.3f}%)'.format(
        percentage_v * 100.))
print('')

#################################################
#
#	COMMON VARIANTS
#
#################################################
set_of_vars = [set([]) for n in INP]
color_ind = 0
for fn in INP:
    DATA_DICT = pickle.load(open(fn, "rb"))
    COMMON_VARIANTS = DATA_DICT['COMMON_VARIANTS']
    for n in COMMON_VARIANTS:
        set_of_vars[color_ind].add(n)
    color_ind += 1

print('')
print('COMMON_VARIANTS OVERLAP BETWEEN ' + str(N_FILES) + ' MODELS...')
for i in range(N_FILES):
    for j in range(i, N_FILES):
        overlap_count = len(set_of_vars[i].intersection(set_of_vars[j]))
        n_dissimilar = (len(set_of_vars[i]) - overlap_count) + (len(set_of_vars[j]) - overlap_count)
        if overlap_count == 0:
            percentage_v = 0.0
        else:
            percentage_v = overlap_count / float(overlap_count + n_dissimilar)
    print('overlap[' + str(i) + ',' + str(j) + '] = ' + str(overlap_count) + ' variants ({0:.3f}%)'.format(
        percentage_v * 100.))
print('')
