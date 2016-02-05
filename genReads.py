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

from inputChecking		import requiredField, checkFileOpen, checkDir
from refFunc			import indexRef, readRef, ALLOWED_NUCL
from vcfFunc			import parseVCF
from OutputFileWriter	import OutputFileWriter
from probability		import DiscreteDistribution
from SequenceContainer	import SequenceContainer, ReadContainer, parseInputMutationModel


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
parser.add_argument('-m', type=str,   required=False, metavar='<str>',   default=None,  help="mutation model directory")
parser.add_argument('-M', type=float, required=False, metavar='<float>', default=-1,    help="rescale avg mutation rate to this")
parser.add_argument('-s', type=str,   required=False, metavar='<str>',   default=None,  help="input sample model")
parser.add_argument('-v', type=str,   required=False, metavar='<str>',   default=None,  help="input VCF file")

parser.add_argument('--paired-art', nargs=2, type=int,   required=False, metavar=('<int>','<int>'), default=(0,0), help='paired-end fragment length and std')
parser.add_argument('--paired-emp', nargs=1, type=str,   required=False, metavar='<str>',           default=None,  help='empirical fragment length distribution')
parser.add_argument('--cancer',                          required=False, action='store_true',       default=False, help='produce tumor/normal datasets')
parser.add_argument('-cm',                   type=str,   required=False, metavar='<str>',           default=None,  help="cancer mutation model directory")
parser.add_argument('-cp',                   type=float, required=False, metavar='<float>',         default=0.8,   help="tumor sample purity")
parser.add_argument('--gc-model',   nargs=1, type=str,   required=False, metavar='<str>',           default=None,  help='empirical GC coverage bias distribution')
parser.add_argument('--job',        nargs=2, type=int,   required=False, metavar=('<int>','<int>'), default=(0,0), help='jobs IDs for generating reads in parallel')
parser.add_argument('--bam',                             required=False, action='store_true',       default=False, help='output golden BAM file')
parser.add_argument('--vcf',                             required=False, action='store_true',       default=False, help='output golden VCF file')
parser.add_argument('--rng',                 type=int,   required=False, metavar='<int>',           default=-1,    help='rng seed value')
parser.add_argument('--gz',                              required=False, action='store_true',       default=False, help='gzip output FQ and VCF')
args = parser.parse_args()

# required args
(REFERENCE, READLEN, OUT_PREFIX) = (args.r, args.R, args.o)
# various dataset parameters
(COVERAGE, PLOIDS, INPUT_BED, SE_MODEL, SE_RATE, MUT_MODEL, MUT_RATE, SAMP_MODEL, INPUT_VCF) = (args.c, args.p, args.t, args.e, args.E, args.m, args.M, args.s, args.v)
(CANCER_MODEL, CANCER_PURITY) = (args.cm, args.cp)
# important flags
(CANCER, SAVE_BAM, SAVE_VCF, GZIPPED_OUT) = (args.cancer, args.bam, args.vcf, args.gz)

(FRAGMENT_SIZE, FRAGMENT_STD) = args.paired_art
FRAGLEN_MODEL = args.paired_emp

# if user specified mean/std, use artificial fragment length distribution, otherwise use
# the empirical model specified. If neither, then we're doing single-end reads.
PAIRED_END = False
PAIRED_END_ARTIFICIAL = False
if FRAGMENT_SIZE > 0 and FRAGMENT_STD > 0:
	PAIRED_END = True
	PAIRED_END_ARTIFICIAL = True
elif FRAGLEN_MODEL != None:
	PAIRED_END = True
	PAIRED_END_ARTIFICIAL = False

(MYJOB, NJOBS) = args.job
if MYJOB == 0:
	MYJOB = 1
	NJOBS = 1

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
	SE_MODEL = SIM_PATH+'/models/errorModel_toy.p'
	SE_CLASS = ReadContainer(READLEN,SE_MODEL,SE_RATE)
else:
	print '\nFeature not added yet, sorry!\n'
	exit(1)

#	fragment length distribution
#
if PAIRED_END and not(PAIRED_END_ARTIFICIAL):
	print 'Using empirical fragment length distribution.'
	[FRAGLEN_VALS, FRAGLEN_PROB] = pickle.load(open(FRAGLEN_MODEL,'rb'))
	# should probably add some validation and sanity-checking code here...
	FRAGLEN_DISTRIBUTION = DiscreteDistribution(FRAGLEN_PROB,FRAGLEN_VALS)


"""************************************************
****            HARD-CODED CONSTANTS
************************************************"""


# target window size for read sampling. how many times bigger than read/frag length
WINDOW_TARGET_SCALE = 10
# sub-window size for read sampling windows. this is basically the finest resolution
# that can be obtained for targeted region boundaries and GC% bias
SMALL_WINDOW        = 20


"""************************************************
****               DEFAULT MODELS
************************************************"""


# fragment length distribution: normal distribution that goes out to +- 6 standard deviations
if PAIRED_END and PAIRED_END_ARTIFICIAL:
	print 'Using artificial fragment length distribution.'
	if FRAGMENT_STD == 0:
		FRAGLEN_DISTRIBUTION = DiscreteDistribution([1],[FRAGMENT_SIZE],degenerateVal=FRAGMENT_SIZE)
	else:
		FRAGLEN_VALS = range(FRAGMENT_SIZE-6*FRAGMENT_STD,FRAGMENT_SIZE+6*FRAGMENT_STD+1)
		FRAGLEN_PROB = [np.exp(-(((n-float(FRAGMENT_SIZE))**2)/(2*(FRAGMENT_STD**2)))) for n in FRAGLEN_VALS]
		FRAGLEN_DISTRIBUTION = DiscreteDistribution(FRAGLEN_PROB,FRAGLEN_VALS)


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

	# parse input variants, if present
	inputVariants = []
	if INPUT_VCF != None:
		if CANCER:
			(sampNames, inputVariants) = parseVCF(INPUT_VCF,tumorNormal=True)
			tumorInd  = sampNames.index('TUMOR')
			normalInd = sampNames.index('NORMAL')
		else:
			(sampNames, inputVariants) = parseVCF(INPUT_VCF)
		inputVariants.sort()

	#print sampNames
	#for k in sorted(inputVariants.keys()):
	#	for n in inputVariants[k]:
	#		print k, n

	# parse input targeted regions, if present
	inputRegions = {}
	if INPUT_BED != None:
		with open(INPUT_BED,'r') as f:
			for line in f:
				[myChr,pos1,pos2] = line.strip().split('\t')[:3]
				if myChr not in inputRegions:
					inputRegions[myChr] = []
				inputRegions[myChr].extend([int(pos1),int(pos2)])

	# initialize output files
	bamHeader = None
	if SAVE_BAM:
		bamHeader = [refIndex]
	vcfHeader = None
	if SAVE_VCF:
		vcfHeader = [REFERENCE]
	if CANCER:
		OFW = OutputFileWriter(OUT_PREFIX+'_normal',paired=PAIRED_END,BAM_header=bamHeader,VCF_header=vcfHeader,gzipped=GZIPPED_OUT)
		OFW_CANCER = OutputFileWriter(OUT_PREFIX+'_tumor',paired=PAIRED_END,BAM_header=bamHeader,VCF_header=vcfHeader,gzipped=GZIPPED_OUT)
	else:
		OFW = OutputFileWriter(OUT_PREFIX,paired=PAIRED_END,BAM_header=bamHeader,VCF_header=vcfHeader,gzipped=GZIPPED_OUT)


	"""************************************************
	****                   MAIN()
	************************************************"""


	for RI in xrange(len(refIndex)):

		# read in reference sequence and notate blocks of Ns
		(refSequence,N_regions) = readRef(REFERENCE,refIndex[RI],N_HANDLING)

		# prune invalid input variants, e.g variants that:
		#		- try to delete or alter any N characters
		#		- don't match the reference base at their specified position
		#		- any alt allele contains anything other than allowed characters
		validVariants = []
		nSkipped = 0
		if refIndex[RI][0] in inputVariants:
			for n in inputVariants[refIndex[RI][0]]:
				span = (n[0],n[0]+len(n[1]))
				rseq = str(refSequence[span[0]-1:span[1]-1])	# -1 because going from VCF coords to array coords
				anyBadChr = any((nn not in ALLOWED_NUCL) for nn in [item for sublist in n[2] for item in sublist])
				if rseq != n[1] or 'N' in rseq or anyBadChr:
					nSkipped += 1
					continue
				#if bisect.bisect(N_regions['big'],span[0])%2 or bisect.bisect(N_regions['big'],span[1])%2:
				#	continue
				validVariants.append(n)
		print 'found',len(validVariants),'valid variants for '+refIndex[RI][0]+' in input VCF...'
		if nSkipped:
			print nSkipped,'variants skipped (invalid position or alt allele)'

		# add large random structural variants
		#
		#	TBD!!!

		# determine which structural variants will affect our sampling window positions
		structuralVars = []
		for n in validVariants:
			bufferNeeded = max([max([len(n[1])-len(alt_allele),1]) for alt_allele in n[2]])
			structuralVars.append((n[0]-1,bufferNeeded))	# -1 because going from VCF coords to array coords

		# determine sampling windows based on read length, large N regions, and structural mutations.
		# in order to obtain uniform coverage, windows should overlap by:
		#		- READLEN, if single-end reads
		#		- FRAGMENT_SIZE, if paired-end reads
		# ploidy is fixed per large sampling window,
		# coverage distributions due to GC% and targeted regions are specified within these windows
		samplingWindows  = []
		readNameCount    = 1
		ALL_VARIANTS_OUT = {}
		if PAIRED_END:
			targSize = WINDOW_TARGET_SCALE*FRAGMENT_SIZE
			overlap  = FRAGMENT_SIZE
		else:
			targSize = WINDOW_TARGET_SCALE*READLEN
			overlap  = READLEN
		for i in xrange(len(N_regions['non_N'])):
			(pi,pf) = N_regions['non_N'][i]
			nTargWindows = max([1,(pf-pi)/targSize])
			bpd = int((pf-pi)/float(nTargWindows))
			print (pi,pf), nTargWindows
			print structuralVars

			# adjust end-position of window based on inserted structural mutations
			start = pi
			end   = min([start+bpd,pf])
			currentVariantInd = 0
			varsFromPrevOverlap = []
			varsCancerFromPrevOverlap = []
			vindFromPrev = 0
			isLastTime = False
			while True:
				relevantVars = []
				if len(structuralVars) and currentVariantInd < len(structuralVars):
					prevVarInd = currentVariantInd
					while structuralVars[currentVariantInd][0] <= end:
						delta = (end-1) - (structuralVars[currentVariantInd][0] + structuralVars[currentVariantInd][1])
						if delta <= 0:
							end -= (delta-1)
						currentVariantInd += 1
						if currentVariantInd == len(structuralVars):
							break
					relevantVars = structuralVars[prevVarInd:currentVariantInd]
				next_start = end-overlap
				next_end   = min([next_start+bpd,pf])
				if next_end-next_start < bpd:
					end = next_end
					isLastTime = True
				print 'PROCESSING WINDOW:',(start,end)

				# which inserted variants are in this window?
				varsInWindow = []
				updated = False
				for j in xrange(vindFromPrev,len(inputVariants)):
					vPos = inputVariants[j][0]
					if vPos >= start and vPos < end:
						varsInWindow.append([vPos-1]+inputVariants[j][1:])	# vcf --> array coords
					if vPos >= end-overlap-1 and updated == False:
						updated = True
						vindFromPrev = j
					if vPos >= end:
						break

				# construct sequence data that we will sample reads from
				sequences = SequenceContainer(start,refSequence[start:end],PLOIDS,overlap,READLEN,[MUT_MODEL]*PLOIDS,MUT_RATE)
				sequences.insert_mutations(varsFromPrevOverlap + varsInWindow)
				all_inserted_variants = sequences.random_mutations()

				if CANCER:
					tumor_sequences = SequenceContainer(start,refSequence[start:end],PLOIDS,overlap,READLEN,[CANCER_MODEL]*PLOIDS,MUT_RATE)
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
				
				# for each sampling window, construct sub-windows with coverage information
				covWindows = [COVERAGE for n in xrange((end-start)/SMALL_WINDOW)]
				if (end-start)%SMALL_WINDOW:
					covWindows.append(COVERAGE)
				meanCov = sum(covWindows)/float(len(covWindows))
				if PAIRED_END:
					readsToSample = int(((end-start)*meanCov)/(2*READLEN))+1
				else:
					readsToSample = int(((end-start)*meanCov)/(READLEN))+1

				# sample reads from altered reference
				for i in xrange(readsToSample):
					if PAIRED_END:
						myFraglen = FRAGLEN_DISTRIBUTION.sample()
						myReadData = sequences.sample_read(SE_CLASS,myFraglen)
					else:
						myReadData = sequences.sample_read(SE_CLASS)
				
					myReadName = OUT_PREFIX+'_'+str(readNameCount)
					readNameCount += len(myReadData)

					# write read data out to FASTQ and BAM files
					if len(myReadData) == 1:
						OFW.writeFASTQRecord(myReadName,myReadData[0][2],myReadData[0][3])
						if SAVE_BAM:
							OFW.writeBAMRecord(RI, myReadName+'/1', myReadData[0][0], myReadData[0][1], myReadData[0][2], myReadData[0][3], samFlag=0)
					elif len(myReadData) == 2:
						OFW.writeFASTQRecord(myReadName,myReadData[0][2],myReadData[0][3],read2=myReadData[1][2],qual2=myReadData[1][3])
						if SAVE_BAM:
							OFW.writeBAMRecord(RI, myReadName+'/1', myReadData[0][0], myReadData[0][1], myReadData[0][2], myReadData[0][3], samFlag=99,  matePos=myReadData[1][0])
							OFW.writeBAMRecord(RI, myReadName+'/2', myReadData[1][0], myReadData[1][1], myReadData[1][2], myReadData[1][3], samFlag=147, matePos=myReadData[0][0])
					else:
						print '\nError: Unexpected number of reads generated...\n'
						exit(1)

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

		# write all output variants for this reference
		if SAVE_VCF:
			for k in sorted(ALL_VARIANTS_OUT.keys()):
				currentRef = refIndex[RI][0]
				myID       = '.'
				myQual     = '.'
				myFilt     = 'PASS'
				OFW.writeVCFRecord(currentRef, k[0], myID, k[1], k[2], myQual, myFilt, k[4])

		break

	# close output files
	OFW.closeFiles()
	if CANCER:
		OFW_CANCER.closeFiles()


if __name__ == '__main__':
	main()



