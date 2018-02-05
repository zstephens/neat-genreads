#!/usr/bin/env python
# encoding: utf-8
""" ////////////////////////////////////////////////////////////////////////////////
   ///                                                                          ///
  ///       genReads.py                                                        ///
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

import os
import sys
import copy
import random
import re
import time
import bisect
import cPickle as pickle
import numpy as np
#import matplotlib.pyplot as mpl
import argparse

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1])
sys.path.append(SIM_PATH+'/py/')

from inputChecking		import requiredField, checkFileOpen, checkDir, isInRange
from refFunc			import indexRef, readRef, getAllRefRegions, partitionRefRegions, ALLOWED_NUCL
from vcfFunc			import parseVCF
from OutputFileWriter	import OutputFileWriter
from probability		import DiscreteDistribution, mean_ind_of_weighted_list
from SequenceContainer	import SequenceContainer, ReadContainer, parseInputMutationModel

# if coverage val for a given window/position is below this value, consider it effectively zero.
LOW_COV_THRESH = 50

"""//////////////////////////////////////////////////
////////////    PARSE INPUT ARGUMENTS    ////////////
//////////////////////////////////////////////////"""


parser = argparse.ArgumentParser(description='NEAT-genReads V2.0')
parser.add_argument('-r', type=str,   required=True,  metavar='<str>',                  help="* ref.fa")
parser.add_argument('-R', type=int,   required=True,  metavar='<int>',                  help="* read length")
parser.add_argument('-o', type=str,   required=True,  metavar='<str>',                  help="* output prefix")
parser.add_argument('-c', type=float, required=False, metavar='<float>', default=10.,   help="average coverage")
parser.add_argument('-e', type=str,   required=False, metavar='<str>',   default=None,  help="sequencing error model")
parser.add_argument('-E', type=float, required=False, metavar='<float>', default=-1,    help="rescale avg sequencing error rate to this")
parser.add_argument('-p', type=int,   required=False, metavar='<int>',   default=2,     help="ploidy")
parser.add_argument('-t', type=str,   required=False, metavar='<str>',   default=None,  help="bed file containing targeted regions")
parser.add_argument('-to',type=float, required=False, metavar='<float>', default=0.00,  help="off-target coverage scalar")
parser.add_argument('-m', type=str,   required=False, metavar='<str>',   default=None,  help="mutation model pickle file")
parser.add_argument('-M', type=float, required=False, metavar='<float>', default=-1,    help="rescale avg mutation rate to this")
parser.add_argument('-Mb',type=str,   required=False, metavar='<str>',   default=None,  help="bed file containing positional mut rates")
parser.add_argument('-N', type=int,   required=False, metavar='<int>',   default=-1,    help="below this qual, replace base-calls with 'N's")
#parser.add_argument('-s', type=str,   required=False, metavar='<str>',   default=None,  help="input sample model")
parser.add_argument('-v', type=str,   required=False, metavar='<str>',   default=None,  help="input VCF file")

parser.add_argument('--pe',       nargs=2, type=int,   required=False, metavar=('<int>','<int>'), default=(None,None), help='paired-end fragment length mean and std')
parser.add_argument('--pe-model',          type=str,   required=False, metavar='<str>',           default=None,  help='empirical fragment length distribution')
#parser.add_argument('--cancer',                        required=False, action='store_true',       default=False, help='produce tumor/normal datasets')
#parser.add_argument('-cm',                 type=str,   required=False, metavar='<str>',           default=None,  help="cancer mutation model directory")
#parser.add_argument('-cp',                 type=float, required=False, metavar='<float>',         default=0.8,   help="tumor sample purity")
parser.add_argument('--gc-model',          type=str,   required=False, metavar='<str>',           default=None,  help='empirical GC coverage bias distribution')
parser.add_argument('--job',      nargs=2, type=int,   required=False, metavar=('<int>','<int>'), default=(0,0), help='jobs IDs for generating reads in parallel')
parser.add_argument('--nnr',                           required=False, action='store_true',       default=False, help='save non-N ref regions (for parallel jobs)')
parser.add_argument('--bam',                           required=False, action='store_true',       default=False, help='output golden BAM file')
parser.add_argument('--vcf',                           required=False, action='store_true',       default=False, help='output golden VCF file')
parser.add_argument('--rng',               type=int,   required=False, metavar='<int>',           default=-1,    help='rng seed value; identical RNG value should produce identical runs of the program, so things like read locations, variant positions, error positions, etc, should all be the same.')
parser.add_argument('--gz',                            required=False, action='store_true',       default=False, help='gzip output FQ and VCF')
parser.add_argument('--no-fastq',                      required=False, action='store_true',       default=False, help='bypass fastq generation')
args = parser.parse_args()

# required args
(REFERENCE, READLEN, OUT_PREFIX) = (args.r, args.R, args.o)
# various dataset parameters
(COVERAGE, PLOIDS, INPUT_BED, SE_MODEL, SE_RATE, MUT_MODEL, MUT_RATE, MUT_BED, INPUT_VCF) = (args.c, args.p, args.t, args.e, args.E, args.m, args.M, args.Mb, args.v)
# cancer params (disabled currently)
#(CANCER, CANCER_MODEL, CANCER_PURITY) = (args.cancer, args.cm, args.cp)
(CANCER, CANCER_MODEL, CANCER_PURITY) = (False, None, 0.8)
(OFFTARGET_SCALAR) = (args.to)
# important flags
(SAVE_BAM, SAVE_VCF, GZIPPED_OUT, NO_FASTQ) = (args.bam, args.vcf, args.gz, args.no_fastq)

ONLY_VCF = (NO_FASTQ and SAVE_VCF and not(SAVE_BAM))
if ONLY_VCF:
	print 'Only producing VCF output, that should speed things up a bit...'

# sequencing model parameters
(FRAGMENT_SIZE, FRAGMENT_STD) = args.pe
FRAGLEN_MODEL = args.pe_model
GC_BIAS_MODEL = args.gc_model
N_MAX_QUAL    = args.N

# if user specified no fastq, no bam, no vcf, then inform them of their wasteful ways and exit
if NO_FASTQ == True and SAVE_BAM == False and SAVE_VCF == False:
        print '\nError: No files will be written when --no-fastq is specified without --vcf or --bam.'
        exit(1)

# if user specified mean/std, use artificial fragment length distribution, otherwise use
# the empirical model specified. If neither, then we're doing single-end reads.
PAIRED_END = False
PAIRED_END_ARTIFICIAL = False
if FRAGMENT_SIZE != None and FRAGMENT_STD != None:
	PAIRED_END = True
	PAIRED_END_ARTIFICIAL = True
elif FRAGLEN_MODEL != None:
	PAIRED_END = True
	PAIRED_END_ARTIFICIAL = False

(MYJOB, NJOBS) = args.job
if MYJOB == 0:
	MYJOB = 1
	NJOBS = 1
SAVE_NON_N = args.nnr

RNG_SEED = args.rng
if RNG_SEED == -1:
	RNG_SEED = random.randint(1,99999999)
random.seed(RNG_SEED)


"""************************************************
****            INPUT ERROR CHECKING
************************************************"""


checkFileOpen(REFERENCE,'ERROR: could not open reference',required=True)
checkFileOpen(INPUT_VCF,'ERROR: could not open input VCF',required=False)
checkFileOpen(INPUT_BED,'ERROR: could not open input BED',required=False)
requiredField(OUT_PREFIX,'ERROR: no output prefix provided')
if (FRAGMENT_SIZE == None and FRAGMENT_STD != None) or (FRAGMENT_SIZE != None and FRAGMENT_STD == None):
	print '\nError: --pe argument takes 2 space-separated arguments.\n'
	exit(1)


"""************************************************
****             LOAD INPUT MODELS
************************************************"""


#	mutation models
#
MUT_MODEL = parseInputMutationModel(MUT_MODEL,1)
if CANCER:
	CANCER_MODEL = parseInputMutationModel(CANCER_MODEL,2)
if MUT_RATE < 0.:
	MUT_RATE = None

#	sequencing error model
#
if SE_RATE < 0.:
	SE_RATE = None
if SE_MODEL == None:
	print 'Using default sequencing error model.'
	SE_MODEL = SIM_PATH+'/models/errorModel_toy.p'
	SE_CLASS = ReadContainer(READLEN,SE_MODEL,SE_RATE)
else:
	# probably need to do some sanity checking
	SE_CLASS = ReadContainer(READLEN,SE_MODEL,SE_RATE)

#	GC-bias model
#
if GC_BIAS_MODEL == None:
	print 'Using default gc-bias model.'
	GC_BIAS_MODEL = SIM_PATH+'/models/gcBias_toy.p'
	[GC_SCALE_COUNT, GC_SCALE_VAL] = pickle.load(open(GC_BIAS_MODEL,'rb'))
	GC_WINDOW_SIZE = GC_SCALE_COUNT[-1]
else:
	[GC_SCALE_COUNT, GC_SCALE_VAL] = pickle.load(open(GC_BIAS_MODEL,'rb'))
	GC_WINDOW_SIZE = GC_SCALE_COUNT[-1]

#	fragment length distribution
#
if PAIRED_END and not(PAIRED_END_ARTIFICIAL):
	print 'Using empirical fragment length distribution.'
	[potential_vals, potential_prob] = pickle.load(open(FRAGLEN_MODEL,'rb'))
	FRAGLEN_VALS = []
	FRAGLEN_PROB = []
	for i in xrange(len(potential_vals)):
		if potential_vals[i] > READLEN:
			FRAGLEN_VALS.append(potential_vals[i])
			FRAGLEN_PROB.append(potential_prob[i])
	# should probably add some validation and sanity-checking code here...
	FRAGLEN_DISTRIBUTION = DiscreteDistribution(FRAGLEN_PROB,FRAGLEN_VALS)
	FRAGMENT_SIZE = FRAGLEN_VALS[mean_ind_of_weighted_list(FRAGLEN_PROB)]

#	Indicate not writing FASTQ reads
#
if NO_FASTQ:
	print 'Bypassing FASTQ generation...'

"""************************************************
****            HARD-CODED CONSTANTS
************************************************"""


# target window size for read sampling. how many times bigger than read/frag length
WINDOW_TARGET_SCALE = 100
# sub-window size for read sampling windows. this is basically the finest resolution
# that can be obtained for targeted region boundaries and GC% bias
SMALL_WINDOW        = 20
# is the mutation model constant throughout the simulation? If so we can save a lot of time
CONSTANT_MUT_MODEL  = True


"""************************************************
****               DEFAULT MODELS
************************************************"""


# fragment length distribution: normal distribution that goes out to +- 6 standard deviations
if PAIRED_END and PAIRED_END_ARTIFICIAL:
	print 'Using artificial fragment length distribution. mean='+str(FRAGMENT_SIZE)+', std='+str(FRAGMENT_STD)
	if FRAGMENT_STD == 0:
		FRAGLEN_DISTRIBUTION = DiscreteDistribution([1],[FRAGMENT_SIZE],degenerateVal=FRAGMENT_SIZE)
	else:
		potential_vals = range(FRAGMENT_SIZE-6*FRAGMENT_STD,FRAGMENT_SIZE+6*FRAGMENT_STD+1)
		FRAGLEN_VALS   = []
		for i in xrange(len(potential_vals)):
			if potential_vals[i] > READLEN:
				FRAGLEN_VALS.append(potential_vals[i])
		FRAGLEN_PROB = [np.exp(-(((n-float(FRAGMENT_SIZE))**2)/(2*(FRAGMENT_STD**2)))) for n in FRAGLEN_VALS]
		FRAGLEN_DISTRIBUTION = DiscreteDistribution(FRAGLEN_PROB,FRAGLEN_VALS)


"""************************************************
****          MORE INPUT ERROR CHECKING
************************************************"""


isInRange(READLEN,           10,1000000, 'Error: -R must be between 10 and 1,000,000')
isInRange(COVERAGE,          0,1000000,  'Error: -c must be between 0 and 1,000,000')
isInRange(PLOIDS,            1,100,      'Error: -p must be between 1 and 100')
isInRange(OFFTARGET_SCALAR,  0,1,        'Error: -to must be between 0 and 1')
if MUT_RATE != -1 and MUT_RATE != None:
	isInRange(MUT_RATE,      0,0.3,      'Error: -M must be between 0 and 0.3')
if SE_RATE != -1 and SE_RATE != None:
	isInRange(SE_RATE,       0,0.3,      'Error: -E must be between 0 and 0.3')
if NJOBS != 1:
	isInRange(NJOBS,         1,1000,     'Error: --job must be between 1 and 1,000')
	isInRange(MYJOB,         1,1000,     'Error: --job must be between 1 and 1,000')
	isInRange(MYJOB,         1,NJOBS,    'Error: job id must be less than or equal to number of jobs')
if N_MAX_QUAL != -1:
	isInRange(N_MAX_QUAL,    1,40,       'Error: -N must be between 1 and 40')


"""************************************************
****                   MAIN()
************************************************"""


def main():

	# index reference
	refIndex = indexRef(REFERENCE)
	if PAIRED_END:
		N_HANDLING = ('random',FRAGMENT_SIZE)
	else:
		N_HANDLING = ('ignore',READLEN)
	indices_by_refName = {refIndex[n][0]:n for n in xrange(len(refIndex))}

	# parse input variants, if present
	inputVariants = []
	if INPUT_VCF != None:
		if CANCER:
			(sampNames, inputVariants) = parseVCF(INPUT_VCF,tumorNormal=True,ploidy=PLOIDS)
			tumorInd  = sampNames.index('TUMOR')
			normalInd = sampNames.index('NORMAL')
		else:
			(sampNames, inputVariants) = parseVCF(INPUT_VCF,ploidy=PLOIDS)
		for k in sorted(inputVariants.keys()):
			inputVariants[k].sort()

	# parse input targeted regions, if present
	inputRegions = {}
	refList      = [n[0] for n in refIndex]
	if INPUT_BED != None:
		with open(INPUT_BED,'r') as f:
			for line in f:
				[myChr,pos1,pos2] = line.strip().split('\t')[:3]
				if myChr not in inputRegions:
					inputRegions[myChr] = [-1]
				inputRegions[myChr].extend([int(pos1),int(pos2)])
		# some validation
		nInBedOnly = 0
		nInRefOnly = 0
		for k in refList:
			if k not in inputRegions:
				nInRefOnly += 1
		for k in inputRegions.keys():
			if not k in refList:
				nInBedOnly += 1
				del inputRegions[k]
		if nInRefOnly > 0:
			print 'Warning: Reference contains sequences not found in targeted regions BED file.'
		if nInBedOnly > 0:
			print 'Warning: Targeted regions BED file contains sequence names not found in reference (regions ignored).'


	# parse input mutation rate rescaling regions, if present
	mutRateRegions = {}
	mutRateValues  = {}
	if MUT_BED != None:
		with open(MUT_BED,'r') as f:
			for line in f:
				[myChr,pos1,pos2,metaData] = line.strip().split('\t')[:4]
				mutStr = re.findall(r"MUT_RATE=.*?(?=;)",metaData+';')
				(pos1,pos2) = (int(pos1),int(pos2))
				if len(mutStr) and (pos2-pos1) > 1:
					# mutRate = #_mutations / length_of_region, let's bound it by a reasonable amount
					mutRate = max([0.0,min([float(mutStr[0][9:]),0.3])])
					if myChr not in mutRateRegions:
						mutRateRegions[myChr] = [-1]
						mutRateValues[myChr]  = [0.0]
					mutRateRegions[myChr].extend([pos1,pos2])
					mutRateValues.extend([mutRate*(pos2-pos1)]*2)

	# initialize output files (part I)
	bamHeader = None
	if SAVE_BAM:
		bamHeader = [copy.deepcopy(refIndex)]
	vcfHeader = None
	if SAVE_VCF:
		vcfHeader = [REFERENCE]
	
	# If processing jobs in parallel, precompute the independent regions that can be process separately
	if NJOBS > 1:
		parallelRegionList  = getAllRefRegions(REFERENCE,refIndex,N_HANDLING,saveOutput=SAVE_NON_N)
		(myRefs, myRegions) = partitionRefRegions(parallelRegionList,refIndex,MYJOB,NJOBS)
		if not len(myRegions):
			print 'This job id has no regions to process, exiting...'
			exit(1)
		for i in xrange(len(refIndex)-1,-1,-1):	# delete reference not used in our job
			if not refIndex[i][0] in myRefs:
				del refIndex[i]
		# if value of NJOBS is too high, let's change it to the maximum possible, to avoid output filename confusion
		corrected_nJobs = min([NJOBS,sum([len(n) for n in parallelRegionList.values()])])
	else:
		corrected_nJobs = 1

	# initialize output files (part II)
	if CANCER:
		OFW = OutputFileWriter(OUT_PREFIX+'_normal',paired=PAIRED_END,BAM_header=bamHeader,VCF_header=vcfHeader,gzipped=GZIPPED_OUT,noFASTQ=NO_FASTQ)
		OFW_CANCER = OutputFileWriter(OUT_PREFIX+'_tumor',paired=PAIRED_END,BAM_header=bamHeader,VCF_header=vcfHeader,gzipped=GZIPPED_OUT,jobTuple=(MYJOB,corrected_nJobs),noFASTQ=NO_FASTQ)
	else:
		OFW = OutputFileWriter(OUT_PREFIX,paired=PAIRED_END,BAM_header=bamHeader,VCF_header=vcfHeader,gzipped=GZIPPED_OUT,jobTuple=(MYJOB,corrected_nJobs),noFASTQ=NO_FASTQ)
	OUT_PREFIX_NAME = OUT_PREFIX.split('/')[-1]


	"""************************************************
	****        LET'S GET THIS PARTY STARTED...
	************************************************"""


	readNameCount = 1	# keep track of the number of reads we've sampled, for read-names
	unmapped_records = []

	for RI in xrange(len(refIndex)):

		# read in reference sequence and notate blocks of Ns
		(refSequence,N_regions) = readRef(REFERENCE,refIndex[RI],N_HANDLING)

		# if we're processing jobs in parallel only take the regions relevant for the current job
		if NJOBS > 1:
			for i in xrange(len(N_regions['non_N'])-1,-1,-1):
				if not (refIndex[RI][0],N_regions['non_N'][i][0],N_regions['non_N'][i][1]) in myRegions:
					del N_regions['non_N'][i]

		# count total bp we'll be spanning so we can get an idea of how far along we are (for printing progress indicators)
		total_bp_span   = sum([n[1]-n[0] for n in N_regions['non_N']])
		currentProgress = 0
		currentPercent  = 0

		# prune invalid input variants, e.g variants that:
		#		- try to delete or alter any N characters
		#		- don't match the reference base at their specified position
		#		- any alt allele contains anything other than allowed characters
		validVariants = []
		nSkipped = [0,0,0]
		if refIndex[RI][0] in inputVariants:
			for n in inputVariants[refIndex[RI][0]]:
				span = (n[0],n[0]+len(n[1]))
				rseq = str(refSequence[span[0]-1:span[1]-1])	# -1 because going from VCF coords to array coords
				anyBadChr = any((nn not in ALLOWED_NUCL) for nn in [item for sublist in n[2] for item in sublist])
				if rseq != n[1]:
					nSkipped[0] += 1
					continue
				elif 'N' in rseq:
					nSkipped[1] += 1
					continue
				elif anyBadChr:
					nSkipped[2] += 1
					continue
				#if bisect.bisect(N_regions['big'],span[0])%2 or bisect.bisect(N_regions['big'],span[1])%2:
				#	continue
				validVariants.append(n)
			print 'found',len(validVariants),'valid variants for '+refIndex[RI][0]+' in input VCF...'
			if any(nSkipped):
				print sum(nSkipped),'variants skipped...'
				print ' - ['+str(nSkipped[0])+'] ref allele does not match reference'
				print ' - ['+str(nSkipped[1])+'] attempting to insert into N-region'
				print ' - ['+str(nSkipped[2])+'] alt allele contains non-ACGT characters'


		# add large random structural variants
		#
		#	TBD!!!


		# determine sampling windows based on read length, large N regions, and structural mutations.
		# in order to obtain uniform coverage, windows should overlap by:
		#		- READLEN, if single-end reads
		#		- FRAGMENT_SIZE (mean), if paired-end reads
		# ploidy is fixed per large sampling window,
		# coverage distributions due to GC% and targeted regions are specified within these windows
		samplingWindows  = []
		ALL_VARIANTS_OUT = {}
		sequences        = None
		if PAIRED_END:
			targSize = WINDOW_TARGET_SCALE*FRAGMENT_SIZE
			overlap  = FRAGMENT_SIZE
		else:
			targSize = WINDOW_TARGET_SCALE*READLEN
			overlap  = READLEN

		print '--------------------------------'
		if ONLY_VCF:
			print 'generating vcf...'
		else:
			print 'sampling reads...'
		tt = time.time()

		for i in xrange(len(N_regions['non_N'])):
			(pi,pf) = N_regions['non_N'][i]
			nTargWindows = max([1,(pf-pi)/targSize])
			bpd = int((pf-pi)/float(nTargWindows))
			bpd += GC_WINDOW_SIZE - bpd%GC_WINDOW_SIZE

			#print len(refSequence), (pi,pf), nTargWindows
			#print structuralVars

			# if for some reason our region is too small to process, skip it! (sorry)
			if nTargWindows == 1 and (pf-pi) < overlap-1:
				#print 'Does this ever happen?'
				continue

			start = pi
			end   = min([start+bpd,pf])
			####print '------------------RAWR:', (pi,pf), bpd
			varsFromPrevOverlap = []
			varsCancerFromPrevOverlap = []
			vindFromPrev = 0
			isLastTime = False
			havePrinted100 = False

			while True:
				
				# which inserted variants are in this window?
				varsInWindow = []
				updated = False
				for j in xrange(vindFromPrev,len(validVariants)):
					vPos = validVariants[j][0]
					if vPos > start and vPos < end:	# update: changed >= to >, so variant cannot be inserted in first position
						varsInWindow.append(tuple([vPos-1]+list(validVariants[j][1:])))	# vcf --> array coords
					if vPos >= end-overlap-1 and updated == False:
						updated = True
						vindFromPrev = j
					if vPos >= end:
						break

				# determine which structural variants will affect our sampling window positions
				structuralVars = []
				for n in varsInWindow:
					bufferNeeded = max([max([abs(len(n[1])-len(alt_allele)),1]) for alt_allele in n[2]]) # change: added abs() so that insertions are also buffered.
					structuralVars.append((n[0]-1,bufferNeeded))	# -1 because going from VCF coords to array coords

				# adjust end-position of window based on inserted structural mutations
				buffer_added = 0
				keepGoing = True
				while keepGoing:
					keepGoing = False
					for n in structuralVars:
						# adding "overlap" here to prevent SVs from being introduced in overlap regions
						# (which can cause problems if random mutations from the previous window land on top of them)
						delta = (end-1) - (n[0] + n[1]) - 2 - overlap
						if delta < 0:
							#print 'DELTA:', delta, 'END:', end, '-->',
							buffer_added = -delta
							end += buffer_added
							####print end
							keepGoing = True
							break
				next_start = end-overlap
				next_end   = min([next_start+bpd,pf])
				if next_end-next_start < bpd:
					end = next_end
					isLastTime = True

				# print progress indicator
				####print 'PROCESSING WINDOW:',(start,end), [buffer_added], 'next:', (next_start,next_end)
				currentProgress += end-start
				newPercent = int((currentProgress*100)/float(total_bp_span))
				if newPercent > currentPercent:
					if newPercent <= 99 or (newPercent == 100 and not havePrinted100):
						sys.stdout.write(str(newPercent)+'% ')
						sys.stdout.flush()
					currentPercent = newPercent
					if currentPercent == 100:
						havePrinted100 = True

				##### if computing only VCF, we can skip this...
				####if ONLY_VCF:
				####	coverage_dat = None
				####	coverage_avg = None
				####else:
				####	# pre-compute gc-bias and targeted sequencing coverage modifiers
				####	nSubWindows  = (end-start)/GC_WINDOW_SIZE
				####	coverage_dat = (GC_WINDOW_SIZE,[])
				####	for j in xrange(nSubWindows):
				####		rInd = start + j*GC_WINDOW_SIZE
				####		if INPUT_BED == None:
				####			tCov = True
				####		elif refIndex[RI][0] in inputRegions:
				####			tCov = not(bisect.bisect(inputRegions[refIndex[RI][0]],rInd)%2) or not(bisect.bisect(inputRegions[refIndex[RI][0]],rInd+GC_WINDOW_SIZE)%2)
				####		else:
				####			tCov = False
				####		if tCov:
				####			tScl = 1.0
				####		else:
				####			tScl = OFFTARGET_SCALAR
				####		gc_v = refSequence[rInd:rInd+GC_WINDOW_SIZE].count('G') + refSequence[rInd:rInd+GC_WINDOW_SIZE].count('C')
				####		gScl = GC_SCALE_VAL[gc_v]
				####		coverage_dat[1].append(1.0*tScl*gScl)
				####	coverage_avg = np.mean(coverage_dat[1])
				####
				##### pre-compute mutation rate tracks
				##### PROVIDED MUTATION RATES OVERRIDE AVERAGE VALUE

				# construct sequence data that we will sample reads from
				if sequences == None:
					sequences = SequenceContainer(start,refSequence[start:end],PLOIDS,overlap,READLEN,[MUT_MODEL]*PLOIDS,MUT_RATE,onlyVCF=ONLY_VCF)
				else:
					sequences.update(start,refSequence[start:end],PLOIDS,overlap,READLEN,[MUT_MODEL]*PLOIDS,MUT_RATE)

				# some inserted variant debugging...
				####print '\n',start,end,end-overlap,'\n',varsFromPrevOverlap,'\n',varsInWindow,'\n'

				# insert variants
				sequences.insert_mutations(varsFromPrevOverlap + varsInWindow)
				all_inserted_variants = sequences.random_mutations()
				#print all_inserted_variants

				# compute coverage modifiers
				coverage_avg = None
				if ONLY_VCF == False:
					coverage_dat = [GC_WINDOW_SIZE,GC_SCALE_VAL,[]]
					if INPUT_BED == None:
						coverage_dat[2] = [1.0]*(end-start)
					else:
						for j in xrange(start,end):
							if not(bisect.bisect(inputRegions[refIndex[RI][0]],j)%2):
								coverage_dat[2].append(1.0)
							else:
								coverage_dat[2].append(OFFTARGET_SCALAR)
					ASDF_TT = time.time()
					if sum(coverage_dat[2]) < LOW_COV_THRESH:
						coverage_avg = 0.0
					else:
						if PAIRED_END:
							coverage_avg = sequences.init_coverage(tuple(coverage_dat),fragDist=FRAGLEN_DISTRIBUTION)
						else:
							coverage_avg = sequences.init_coverage(tuple(coverage_dat))
					#print 'COV:',time.time()-ASDF_TT

				# unused cancer stuff
				if CANCER:
					tumor_sequences = SequenceContainer(start,refSequence[start:end],PLOIDS,overlap,READLEN,[CANCER_MODEL]*PLOIDS,MUT_RATE,coverage_dat)
					tumor_sequences.insert_mutations(varsCancerFromPrevOverlap + all_inserted_variants)
					all_cancer_variants = tumor_sequences.random_mutations()

				# which variants do we need to keep for next time (because of window overlap)?
				varsFromPrevOverlap       = []
				varsCancerFromPrevOverlap = []
				for n in all_inserted_variants:
					if n[0] >= end-overlap-1:
						varsFromPrevOverlap.append(n)
				if CANCER:
					for n in all_cancer_variants:
						if n[0] >= end-overlap-1:
							varsCancerFromPrevOverlap.append(n)

				# if we're only producing VCF, no need to go through the hassle of generating reads
				if ONLY_VCF:
					pass
				else:
					if PAIRED_END:
						readsToSample = int(((end-start)*float(COVERAGE)*coverage_avg)/(2*READLEN))+1
					else:
						readsToSample = int(((end-start)*float(COVERAGE)*coverage_avg)/(READLEN))+1

					# if coverage is so low such that no reads are to be sampled, skip region
					#      (i.e., remove buffer of +1 reads we add to every window)
					if readsToSample == 1 and sum(coverage_dat[2]) < LOW_COV_THRESH:
						readsToSample = 0

					# sample reads
					ASDF2_TT = time.time()
					for i in xrange(readsToSample):

						isUnmapped = []
						if PAIRED_END:
							myFraglen = FRAGLEN_DISTRIBUTION.sample()
							myReadData = sequences.sample_read(SE_CLASS,myFraglen)
							if myReadData == None:	# skip if we failed to find a valid position to sample read
								continue
							if myReadData[0][0] == None:
								isUnmapped.append(True)
							else:
								isUnmapped.append(False)
								myReadData[0][0] += start	# adjust mapping position based on window start
							if myReadData[1][0] == None:
								isUnmapped.append(True)
							else:
								isUnmapped.append(False)
								myReadData[1][0] += start
						else:
							myReadData = sequences.sample_read(SE_CLASS)
							if myReadData == None:	# skip if we failed to find a valid position to sample read
								continue
							if myReadData[0][0] == None:	# unmapped read (lives in large insertion)
								isUnmapped = [True]
							else:
								isUnmapped = [False]
								myReadData[0][0] += start	# adjust mapping position based on window start

						if NJOBS > 1:
							myReadName = OUT_PREFIX_NAME+'-j'+str(MYJOB)+'-'+refIndex[RI][0]+'-r'+str(readNameCount)
						else:
							myReadName = OUT_PREFIX_NAME+'-'+refIndex[RI][0]+'-'+str(readNameCount)
						readNameCount += len(myReadData)

						# if desired, replace all low-quality bases with Ns
						if N_MAX_QUAL > -1:
							for j in xrange(len(myReadData)):
								myReadString = [n for n in myReadData[j][2]]
								for k in xrange(len(myReadData[j][3])):
									adjusted_qual = ord(myReadData[j][3][k])-SE_CLASS.offQ
									if adjusted_qual <= N_MAX_QUAL:
										myReadString[k] = 'N'
								myReadData[j][2] = ''.join(myReadString)

						# if read (or read + mate for PE) are unmapped, put them at end of bam file
						if all(isUnmapped):
							if PAIRED_END:
								unmapped_records.append((myReadName+'/1',myReadData[0],109))
								unmapped_records.append((myReadName+'/2',myReadData[1],157))
							else:
								unmapped_records.append((myReadName+'/1',myReadData[0],4))

						# write read data out to FASTQ and BAM files, bypass FASTQ if option specified
						myRefIndex = indices_by_refName[refIndex[RI][0]]
						if len(myReadData) == 1:
							if NO_FASTQ != True:
								OFW.writeFASTQRecord(myReadName,myReadData[0][2],myReadData[0][3])
							if SAVE_BAM:
								if isUnmapped[0] == False:
									OFW.writeBAMRecord(myRefIndex, myReadName+'/1', myReadData[0][0], myReadData[0][1], myReadData[0][2], myReadData[0][3], samFlag=0)
						elif len(myReadData) == 2:
							if NO_FASTQ != True:
								OFW.writeFASTQRecord(myReadName,myReadData[0][2],myReadData[0][3],read2=myReadData[1][2],qual2=myReadData[1][3])
							if SAVE_BAM:
								if isUnmapped[0] == False and isUnmapped[1] == False:
									OFW.writeBAMRecord(myRefIndex, myReadName+'/1', myReadData[0][0], myReadData[0][1], myReadData[0][2], myReadData[0][3], samFlag=99,  matePos=myReadData[1][0])
									OFW.writeBAMRecord(myRefIndex, myReadName+'/2', myReadData[1][0], myReadData[1][1], myReadData[1][2], myReadData[1][3], samFlag=147, matePos=myReadData[0][0])
								elif isUnmapped[0] == False and isUnmapped[1] == True:
									OFW.writeBAMRecord(myRefIndex, myReadName+'/1', myReadData[0][0], myReadData[0][1], myReadData[0][2], myReadData[0][3], samFlag=105,  matePos=myReadData[0][0])
									OFW.writeBAMRecord(myRefIndex, myReadName+'/2', myReadData[0][0], myReadData[1][1], myReadData[1][2], myReadData[1][3], samFlag=149, matePos=myReadData[0][0], alnMapQual=0)
								elif isUnmapped[0] == True and isUnmapped[1] == False:
									OFW.writeBAMRecord(myRefIndex, myReadName+'/1', myReadData[1][0], myReadData[0][1], myReadData[0][2], myReadData[0][3], samFlag=101,  matePos=myReadData[1][0], alnMapQual=0)
									OFW.writeBAMRecord(myRefIndex, myReadName+'/2', myReadData[1][0], myReadData[1][1], myReadData[1][2], myReadData[1][3], samFlag=153, matePos=myReadData[1][0])
						else:
							print '\nError: Unexpected number of reads generated...\n'
							exit(1)
					#print 'READS:',time.time()-ASDF2_TT

					if not isLastTime:
						OFW.flushBuffers(bamMax=next_start)
					else:
						OFW.flushBuffers(bamMax=end+1)

				# tally up all the variants that got successfully introduced
				for n in all_inserted_variants:
					ALL_VARIANTS_OUT[n] = True

				# prepare indices of next window
				start = next_start
				end   = next_end
				if isLastTime:
					break
				if end >= pf:
					isLastTime = True

		if currentPercent != 100 and not havePrinted100:
			print '100%'
		else:
			print ''
		if ONLY_VCF:
			print 'VCF generation completed in',
		else:
			print 'Read sampling completed in',
		print int(time.time()-tt),'(sec)'

		# write all output variants for this reference
		if SAVE_VCF:
			print 'Writing output VCF...'
			for k in sorted(ALL_VARIANTS_OUT.keys()):
				currentRef = refIndex[RI][0]
				myID       = '.'
				myQual     = '.'
				myFilt     = 'PASS'
				# k[0] + 1 because we're going back to 1-based vcf coords
				OFW.writeVCFRecord(currentRef, str(int(k[0])+1), myID, k[1], k[2], myQual, myFilt, k[4])

		#break

	# write unmapped reads to bam file
	if SAVE_BAM and len(unmapped_records):
		print 'writing unmapped reads to bam file...'
		for umr in unmapped_records:
			if PAIRED_END:
				OFW.writeBAMRecord(-1, umr[0], 0, umr[1][1], umr[1][2], umr[1][3], samFlag=umr[2], matePos=0, alnMapQual=0)
			else:
				OFW.writeBAMRecord(-1, umr[0], 0, umr[1][1], umr[1][2], umr[1][3], samFlag=umr[2], alnMapQual=0)

	# close output files
	OFW.closeFiles()
	if CANCER:
		OFW_CANCER.closeFiles()


if __name__ == '__main__':
	main()



