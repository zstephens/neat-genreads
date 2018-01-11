#!/usr/bin/env python

#
#	a quick script for comparing mutation models
#
#	python plotMutModel.py -i model1.p [model2.p] [model3.p]... -l legend_label1 [legend_label2] [legend_label3]... -o path/to/pdf_plot_prefix 
#

import sys
import pickle
import bisect
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import argparse

#mpl.rc('text',usetex=True)
#mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

parser = argparse.ArgumentParser(description='Plot and compare mutation models from genMutModel.py Usage: python plotMutModel.py -i model1.p [model2.p] [model3.p]... -l legend_label1 [legend_label2] [legend_label3]... -o path/to/pdf_plot_prefix')
parser.add_argument('-i',  type=str,   required=True,   metavar='<str>',   nargs='+',                help="* mutation_model_1.p [mutation_model_2.p] [mutation_model_3] ...")
parser.add_argument('-l',  type=str,   required=True,   metavar='<str>',   nargs='+',                help="* legend labels: model1_name [model2_name] [model3_name]...")
parser.add_argument('-o',  type=str,   required=True,   metavar='<str>',                             help="* output pdf prefix")
args = parser.parse_args()



def getColor(i,N,colormap='jet'):
	cm = mpl.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=0, vmax=N+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	colorVal = scalarMap.to_rgba(i)
	return colorVal

def isInBed(track,ind):
	if ind in track:
		return True
	elif bisect.bisect(track,ind)%1 == 1:
		return True
	else:
		return False

def getBedOverlap(track,ind_s,ind_e):
	if ind_s in track:
		myInd = track.index(ind_s)
		return min([track[myInd+1]-ind_s+1,ind_e-ind_s+1])
	else:
		myInd = bisect.bisect(track,ind_s)
		if myInd%1 and myInd < len(track)-1:
			return min([track[myInd+1]-ind_s+1,ind_e-ind_s+1])
	return 0

# a waaaaaaay slower version of the above function ^^
#def getTrackOverlap(track1,track2):
#	otrack = [0 for n in xrange(max(track1+track2)+1)]
#	for i in xrange(0,len(track1),2):
#		for j in xrange(track1[i],track1[i+1]+1):
#			otrack[j] = 1
#	ocount = 0
#	for i in xrange(0,len(track2),2):
#		for j in xrange(track2[i],track2[i+1]+1):
#			if otrack[j]:
#				ocount += 1
#	return ocount

OUP  = args.o
LAB = args.l
#print LAB
INP  = args.i

N_FILES = len(INP)

mpl.rcParams.update({'font.size': 13, 'font.weight':'bold', 'lines.linewidth': 3})

#################################################
#
#	BASIC STATS
#
#################################################
mpl.figure(0,figsize=(12,10))

mpl.subplot(2,2,1)
colorInd = 0
for fn in INP:
	myCol = getColor(colorInd,N_FILES)
	colorInd += 1
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	[AVG_MUT_RATE, SNP_FREQ, INDEL_FREQ] = [DATA_DICT['AVG_MUT_RATE'], DATA_DICT['SNP_FREQ'], DATA_DICT['INDEL_FREQ']]
	mpl.bar([colorInd-1],[AVG_MUT_RATE],1.,color=myCol)
mpl.xlim([-1,N_FILES+1])
mpl.grid()
mpl.xticks([],[])
mpl.ylabel('Frequency')
mpl.title('Overall mutation rate (1/bp)')

mpl.subplot(2,2,2)
colorInd = 0
for fn in INP:
	myCol = getColor(colorInd,N_FILES)
	colorInd += 1
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	[AVG_MUT_RATE, SNP_FREQ, INDEL_FREQ] = [DATA_DICT['AVG_MUT_RATE'], DATA_DICT['SNP_FREQ'], DATA_DICT['INDEL_FREQ']]
	mpl.bar([colorInd-1],[SNP_FREQ],1.,color=myCol)
	mpl.bar([colorInd-1],[1.-SNP_FREQ],1.,color=myCol,bottom=[SNP_FREQ],hatch='/')
mpl.axis([-1,N_FILES+1,0,1.2])
mpl.grid()
mpl.xticks([],[])
mpl.yticks([0,.2,.4,.6,.8,1.],[0,0.2,0.4,0.6,0.8,1.0])
mpl.ylabel('Frequency')
mpl.title('SNP freq [  ] & indel freq [//]')

mpl.subplot(2,1,2)
colorInd = 0
legText  = LAB
for fn in INP:
	myCol = getColor(colorInd,N_FILES)
	colorInd += 1
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	[AVG_MUT_RATE, SNP_FREQ, INDEL_FREQ] = [DATA_DICT['AVG_MUT_RATE'], DATA_DICT['SNP_FREQ'], DATA_DICT['INDEL_FREQ']]
	x = sorted(INDEL_FREQ.keys())
	y = [INDEL_FREQ[n] for n in x]
	mpl.plot(x,y,color=myCol)
	#legText.append(fn)
mpl.grid()
mpl.xlabel('Indel size (bp)', fontweight='bold')
mpl.ylabel('Frequency')
mpl.title('Indel frequency by size (- deletion, + insertion)')
mpl.legend(legText)
#mpl.show()
mpl.savefig(OUP+'_plot1_mutRates.pdf')

#################################################
#
#	TRINUC PRIOR PROB
#
#################################################
mpl.figure(1,figsize=(14,6))
colorInd = 0
legText  = LAB
for fn in INP:
	myCol = getColor(colorInd,N_FILES)
	colorInd += 1
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	TRINUC_MUT_PROB = DATA_DICT['TRINUC_MUT_PROB']

	x  = range(colorInd-1,len(TRINUC_MUT_PROB)*N_FILES,N_FILES)
	xt = sorted(TRINUC_MUT_PROB.keys())
	y  = [TRINUC_MUT_PROB[k] for k in xt]
	markerline, stemlines, baseline = mpl.stem(x,y,'-.')
	mpl.setp(markerline, 'markerfacecolor', myCol)
	mpl.setp(markerline, 'markeredgecolor', myCol)
	mpl.setp(baseline, 'color', myCol, 'linewidth', 0)
	mpl.setp(stemlines, 'color', myCol, 'linewidth', 3)
	if colorInd == 1:
		mpl.xticks(x,xt,rotation=90)
	#legText.append(fn)
mpl.grid()
mpl.ylabel('p(trinucleotide mutates)')
mpl.legend(legText)
#mpl.show()
mpl.savefig(OUP+'_plot2_trinucPriors.pdf')

#################################################
#
#	TRINUC TRANS PROB
#
#################################################
plotNum = 3
for fn in INP:
	fig = mpl.figure(plotNum,figsize=(12,10))
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	TRINUC_TRANS_PROBS = DATA_DICT['TRINUC_TRANS_PROBS']

	xt2 = [m[3] for m in sorted([(n[0],n[2],n[1],n) for n in xt])]
	reverse_dict = {xt2[i]:i for i in xrange(len(xt2))}
	Z = np.zeros((64,64))
	L = [['' for n in xrange(64)] for m in xrange(64)]
	for k in TRINUC_TRANS_PROBS:
		i = reverse_dict[k[0]]
		j = reverse_dict[k[1]]
		Z[i][j] = TRINUC_TRANS_PROBS[k]

	HARDCODED_LABEL = ['A_A','A_C','A_G','A_T',
	                   'C_A','C_C','C_G','C_T',
	                   'G_A','G_C','G_G','G_T',
	                   'T_A','T_C','T_G','T_T']

	for pi in xrange(16):
		mpl.subplot(4,4,pi+1)
		Z2 = Z[pi*4:(pi+1)*4,pi*4:(pi+1)*4]
		X, Y = np.meshgrid( range(0,len(Z2[0])+1), range(0,len(Z2)+1) )
		im = mpl.pcolormesh(X,Y,Z2[::-1,:],vmin=0.0,vmax=0.5)
		mpl.axis([0,4,0,4])
		mpl.xticks([0.5,1.5,2.5,3.5],['A','C','G','T'])
		mpl.yticks([0.5,1.5,2.5,3.5],['T','G','C','A'])
		mpl.text(1.6, 1.8, HARDCODED_LABEL[pi], color='white')

	# colorbar haxx
	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	cb = fig.colorbar(im,cax=cbar_ax)
	cb.set_label(r"p(X$Y_1$Z->X$Y_2$Z | X_Z mutates)")

	#mpl.tight_layout()
	#mpl.figtext(0.24,0.94,'Trinucleotide Mutation Frequency',size=20)
	#mpl.show()
	mpl.savefig(OUP+'_plot'+str(plotNum)+'_trinucTrans.pdf')
	plotNum += 1

#################################################
#
#	HIGH MUT REGIONS
#
#################################################
track_byFile_byChr = [{} for n in INP]
bp_total_byFile    = [0 for n in INP]
colorInd = 0
for fn in INP:
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	HIGH_MUT_REGIONS = DATA_DICT['HIGH_MUT_REGIONS']
	for region in HIGH_MUT_REGIONS:
		if region[0] not in track_byFile_byChr[colorInd]:
			track_byFile_byChr[colorInd][region[0]] = []
		track_byFile_byChr[colorInd][region[0]].extend([region[1],region[2]])
		bp_total_byFile[colorInd] += region[2]-region[1]+1
	colorInd += 1

bp_overlap_count = [[0 for m in INP] for n in INP]
for i in xrange(N_FILES):
	bp_overlap_count[i][i] = bp_total_byFile[i]
	for j in xrange(i+1,N_FILES):
		for k in track_byFile_byChr[i].keys():
			if k in track_byFile_byChr[j]:
				for ii in xrange(len(track_byFile_byChr[i][k][::2])):
					bp_overlap_count[i][j] += getBedOverlap(track_byFile_byChr[j][k],track_byFile_byChr[i][k][ii*2],track_byFile_byChr[i][k][ii*2+1])

print ''				
print 'HIGH_MUT_REGION OVERLAP BETWEEN '+str(N_FILES)+' MODELS...'
for i in xrange(N_FILES):
	for j in xrange(i,N_FILES):
		nDissimilar = (bp_overlap_count[i][i]-bp_overlap_count[i][j]) + (bp_overlap_count[j][j]-bp_overlap_count[i][j])
		if bp_overlap_count[i][j] == 0:
			percentageV = 0.0
		else:
			percentageV = bp_overlap_count[i][j]/float(bp_overlap_count[i][j]+nDissimilar)
		print 'overlap['+str(i)+','+str(j)+'] = '+str(bp_overlap_count[i][j])+' bp ({0:.3f}%)'.format(percentageV*100.)
print ''

#################################################
#
#	COMMON VARIANTS
#
#################################################
setofVars = [set([]) for n in INP]
colorInd = 0
for fn in INP:
	DATA_DICT = pickle.load( open( fn, "rb" ) )
	COMMON_VARIANTS = DATA_DICT['COMMON_VARIANTS']
	for n in COMMON_VARIANTS:
		setofVars[colorInd].add(n)
	colorInd += 1

print ''
print 'COMMON_VARIANTS OVERLAP BETWEEN '+str(N_FILES)+' MODELS...'
for i in xrange(N_FILES):
	for j in xrange(i,N_FILES):
		overlapCount = len(setofVars[i].intersection(setofVars[j]))
		nDissimilar  = (len(setofVars[i])-overlapCount) + (len(setofVars[j])-overlapCount)
		if overlapCount == 0:
			percentageV = 0.0
		else:
			percentageV = overlapCount/float(overlapCount+nDissimilar)
		print 'overlap['+str(i)+','+str(j)+'] = '+str(overlapCount)+' variants ({0:.3f}%)'.format(percentageV*100.)
print ''



