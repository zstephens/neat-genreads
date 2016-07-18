import sys
import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as mpl

# absolute path to the directory above this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-2])
sys.path.append(SIM_PATH+'/py/')

from refFunc import indexRef

parser = argparse.ArgumentParser(description='genMutModel.py')
parser.add_argument('-r', type=str, required=True, metavar='<str>',                    help="* ref.fa")
parser.add_argument('-m', type=str, required=True, metavar='<str>',                    help="* mutations.tsv")
parser.add_argument('-o', type=str, required=True, metavar='<str>',                    help="* output.p")
parser.add_argument('--save-trinuc',required=False,action='store_true', default=False, help='save trinuc counts for ref')
args = parser.parse_args()
(REF, TSV, OUT_PICKLE, SAVE_TRINUC) = (args.r, args.m, args.o, args.save_trinuc)

REF_WHITELIST =  [str(n) for n in xrange(1,30)] + ['x','y','X','Y','mt','Mt','MT']
REF_WHITELIST += ['chr'+n for n in REF_WHITELIST]
VALID_NUCL    =  ['A','C','G','T']
VALID_TRINUC  =  [VALID_NUCL[i]+VALID_NUCL[j]+VALID_NUCL[k] for i in xrange(len(VALID_NUCL)) for j in xrange(len(VALID_NUCL)) for k in xrange(len(VALID_NUCL))]

# given a reference index, grab the sequence string of a specified reference
def getChrFromFasta(refPath,ref_inds,chrName):

	for i in xrange(len(ref_inds)):
		if ref_inds[i][0] == chrName:
			ref_inds_i = ref_inds[i]
			break

	refFile = open(refPath,'r')
	refFile.seek(ref_inds_i[1])
	myDat = ''.join(refFile.read(ref_inds_i[2]-ref_inds_i[1]).split('\n'))
	return myDat

# cluster a sorted list
def clusterList(l,delta):
	outList    = [[l[0]]]
	prevVal    = l[0]
	currentInd = 0
	for n in l[1:]:
		if n-prevVal <= delta:
			outList[currentInd].append(n)
		else:
			currentInd += 1
			outList.append([])
			outList[currentInd].append(n)
		prevVal = n
	return outList

def list_2_countDict(l):
	cDict = {}
	for n in l:
		if n not in cDict:
			cDict[n] = 0
		cDict[n] += 1
	return cDict

# return the mean distance to the median of a cluster
def mean_dist_from_median(c):
	centroid = np.median([n for n in c])
	dists    = []
	for n in c:
		dists.append(abs(n-centroid))
	return np.mean(dists)

# get median value from counting dictionary
def quick_median(countDict):
	midPoint = sum(countDict.values())/2
	mySum    = 0
	myInd    = 0
	sk       = sorted(countDict.keys())
	while mySum < midPoint:
		mySum += countDict[sk[myInd]]
		if mySum >= midPoint:
			break
		myInd += 1
	return myInd

# get median deviation from median of counting dictionary
def median_deviation_from_median(countDict):
	myMedian = quick_median(countDict)
	deviations = {}
	for k in sorted(countDict.keys()):
		d = abs(k-myMedian)
		deviations[d] = countDict[k]
	return quick_median(deviations)


#####################################
#				main()				#
#####################################


def main():

	ref_inds = indexRef(REF)
	refList  = [n[0] for n in ref_inds]

	# how many times do we observe each trinucleotide in the reference?
	TRINUC_REF_COUNT = {}
	# [(trinuc_a, trinuc_b)] = # of times we observed a mutation from trinuc_a into trinuc_b
	TRINUC_TRANSITION_COUNT = {}
	# total count of SNPs
	SNP_COUNT = 0
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

	# load and process variants in each reference sequence individually, for memory reasons...
	for refName in refList:

		if refName not in REF_WHITELIST:
			print refName,'is not in our whitelist, skipping...'
			continue

		print 'reading reference "'+refName+'"...'
		refSequence = getChrFromFasta(REF,ref_inds,refName).upper()
		TOTAL_REFLEN += len(refSequence) - refSequence.count('N')

		# list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
		VDAT_COMMON = []


		""" ##########################################################################
		###						COUNT TRINUCLEOTIDES IN REF						   ###
		########################################################################## """


		if not os.path.isfile(REF+'.trinucCounts'):
			print 'counting trinucleotides in reference...'
			for i in xrange(len(refSequence)-2):
				if i%1000000 == 0 and i > 0:
					print i,'/',len(refSequence)
					#break
				trinuc = refSequence[i:i+3]
				if not trinuc in VALID_TRINUC:
					continue	# skip if trinuc contains invalid characters
				if trinuc not in TRINUC_REF_COUNT:
					TRINUC_REF_COUNT[trinuc] = 0
				TRINUC_REF_COUNT[trinuc] += 1
		else:
			print 'skipping trinuc counts because we found a file...'


		""" ##########################################################################
		###							READ INPUT VARIANTS							   ###
		########################################################################## """


		print 'reading input variants...'
		f = open(TSV,'r')
		isFirst = True
		for line in f:

			if isFirst:
				splt = line.strip().split('\t')
				(c1,c2,c3) = (splt.index('chromosome'),splt.index('chromosome_start'),splt.index('chromosome_end'))
				(m1,m2,m3) = (splt.index('reference_genome_allele'),splt.index('mutated_from_allele'),splt.index('mutated_to_allele'))
				(d_id) = (splt.index('icgc_donor_id'))
				isFirst = False
				continue

			splt = line.strip().split('\t')
			# we have -1 because tsv coords are 1-based, and our reference string index is 0-based
			[chrName,chrStart,chrEnd] = [splt[c1],int(splt[c2])-1,int(splt[c3])-1]
			[allele_ref,allele_normal,allele_tumor] = [splt[m1].upper(),splt[m2].upper(),splt[m3].upper()]
			[donor_id] = [splt[d_id]]

			# hacky, bad.
			if 'chr' not in chrName:
				chrName = 'chr'+chrName
			# hacky, bad.
			if 'chr' not in refName:
				refName = 'chr'+refName

			if chrName != refName:
				continue

			# we want only snps
			# so, no '-' characters allowed, and chrStart must be same as chrEnd
			if '-' not in allele_normal and '-' not in allele_tumor and chrStart == chrEnd:
				trinuc_ref = refSequence[chrStart-1:chrStart+2]
				if not trinuc_ref in VALID_TRINUC:
					continue	# skip ref trinuc with invalid characters
				# only consider positions where ref allele in tsv matches the nucleotide in our reference
				if allele_ref == trinuc_ref[1]:
					trinuc_normal    = refSequence[chrStart-1] + allele_normal + refSequence[chrStart+1]
					trinuc_tumor     = refSequence[chrStart-1] + allele_tumor + refSequence[chrStart+1]
					if not trinuc_normal in VALID_TRINUC or not trinuc_tumor in VALID_TRINUC:
						continue	# skip if mutation contains invalid char
					key = (trinuc_normal,trinuc_tumor)
					if key not in TRINUC_TRANSITION_COUNT:
						TRINUC_TRANSITION_COUNT[key] = 0
					TRINUC_TRANSITION_COUNT[key] += 1
					SNP_COUNT += 1
					VDAT_COMMON.append((chrStart,allele_ref,allele_normal,allele_tumor))
					TOTAL_DONORS[donor_id] = True
				else:
					print '\nError: ref allele in variant call does not match reference.\n'
					print trinuc, allele_ref, allele_normal, allele_tumor
					exit(1)

			# now let's look for indels...
			if '-' in allele_normal: len_normal = 0
			else: len_normal = len(allele_normal)
			if '-' in allele_tumor: len_tumor = 0
			else: len_tumor = len(allele_tumor)
			if len_normal != len_tumor:
				indel_len = len_tumor - len_normal
				if indel_len not in INDEL_COUNT:
					INDEL_COUNT[indel_len] = 0
				INDEL_COUNT[indel_len] += 1
				VDAT_COMMON.append((chrStart,allele_ref,allele_normal,allele_tumor))
				TOTAL_DONORS[donor_id] = True
		f.close()

		#
		# identify common mutations
		#
		percentile_var = 95
		N_DONORS = len(TOTAL_DONORS)
		VDAT_COMMON = list_2_countDict(VDAT_COMMON)
		minVal = int(np.percentile(VDAT_COMMON.values(),percentile_var))
		for k in sorted(VDAT_COMMON.keys()):
			if VDAT_COMMON[k] >= minVal:
				COMMON_VARIANTS.append((refName,k[0],k[1],k[3],VDAT_COMMON[k]/float(N_DONORS)))

		#
		# identify areas that have contained significantly higher random mutation rates
		#
		dist_thresh      = 2000
		percentile_clust = 97
		qptn             = 1000
		# identify regions with disproportionately more variants in them
		VARIANT_POS = sorted([n[0] for n in VDAT_COMMON.keys()])
		clustered_pos = clusterList(VARIANT_POS,dist_thresh)
		byLen  = [(len(clustered_pos[i]),min(clustered_pos[i]),max(clustered_pos[i]),i) for i in xrange(len(clustered_pos))]
		#byLen  = sorted(byLen,reverse=True)
		#minLen = int(np.percentile([n[0] for n in byLen],percentile_clust))
		#byLen  = [n for n in byLen if n[0] >= minLen]
		candidate_regions = []
		for n in byLen:
			bi = int((n[1]-dist_thresh)/float(qptn))*qptn
			bf = int((n[2]+dist_thresh)/float(qptn))*qptn
			candidate_regions.append((n[0]/float(bf-bi),max([0,bi]),min([len(refSequence),bf])))
		minVal = np.percentile([n[0] for n in candidate_regions],percentile_clust)
		for n in candidate_regions:
			if n[0] >= minVal:
				HIGH_MUT_REGIONS.append((refName,n[1],n[2],n[0]))
		# collapse overlapping regions
		for i in xrange(len(HIGH_MUT_REGIONS)-1,0,-1):
			if HIGH_MUT_REGIONS[i-1][2] >= HIGH_MUT_REGIONS[i][1] and HIGH_MUT_REGIONS[i-1][0] == HIGH_MUT_REGIONS[i][0]:
				avgMutRate = 0.5*HIGH_MUT_REGIONS[i-1][3]+0.5*HIGH_MUT_REGIONS[i][3]	# not accurate, but I'm lazy
				HIGH_MUT_REGIONS[i-1] = (HIGH_MUT_REGIONS[i-1][0], HIGH_MUT_REGIONS[i-1][1], HIGH_MUT_REGIONS[i][2], avgMutRate)
				del HIGH_MUT_REGIONS[i]

	#
	# if we didn't count ref trinucs because we found file, read in ref counts from file now
	#
	if os.path.isfile(REF+'.trinucCounts'):
		print 'reading pre-computed trinuc counts...'
		f = open(REF+'.trinucCounts','r')
		for line in f:
			splt = line.strip().split('\t')
			TRINUC_REF_COUNT[splt[0]] = int(splt[1])
		f.close()
	# otherwise, save trinuc counts to file, if desired
	elif SAVE_TRINUC:
		print 'saving trinuc counts to file...'
		f = open(REF+'.trinucCounts','w')
		for trinuc in sorted(TRINUC_REF_COUNT.keys()):
			f.write(trinuc+'\t'+str(TRINUC_REF_COUNT[trinuc])+'\n')
		f.close()


	""" ##########################################################################
	###							COMPUTE PROBABILITIES						   ###
	########################################################################## """


	#for k in sorted(TRINUC_REF_COUNT.keys()):
	#		print k, TRINUC_REF_COUNT[k]
	#
	#for k in sorted(TRINUC_TRANSITION_COUNT.keys()):
	#	print k, TRINUC_TRANSITION_COUNT[k]

	# frequency that each trinuc mutated into anything else
	TRINUC_MUT_PROB = {}
	# frequency that a trinuc mutates into another trinuc, given that it mutated
	TRINUC_TRANS_PROBS = {}

	for trinuc in sorted(TRINUC_REF_COUNT.keys()):
		myCount = 0
		for k in sorted(TRINUC_TRANSITION_COUNT.keys()):
			if k[0] == trinuc:
				myCount += TRINUC_TRANSITION_COUNT[k]
		TRINUC_MUT_PROB[trinuc] = myCount / float(TRINUC_REF_COUNT[trinuc])
		for k in sorted(TRINUC_TRANSITION_COUNT.keys()):
			if k[0] == trinuc:
				TRINUC_TRANS_PROBS[k] = TRINUC_TRANSITION_COUNT[k] / float(myCount)

	# compute average snp and indel frequencies
	totalVar       = SNP_COUNT + sum(INDEL_COUNT.values())
	SNP_FREQ       = SNP_COUNT/float(totalVar)
	AVG_INDEL_FREQ = 1.-SNP_FREQ
	INDEL_FREQ     = {k:(INDEL_COUNT[k]/float(totalVar))/AVG_INDEL_FREQ for k in INDEL_COUNT.keys()}
	AVG_MUT_RATE   = totalVar/float(TOTAL_REFLEN)

	#
	#	print some stuff
	#
	for k in sorted(TRINUC_MUT_PROB.keys()):
		print 'p('+k+' mutates) =',TRINUC_MUT_PROB[k]

	for k in sorted(TRINUC_TRANS_PROBS.keys()):
		print 'p('+k[0]+' --> '+k[1]+' | '+k[0]+' mutates) =',TRINUC_TRANS_PROBS[k]

	for k in sorted(INDEL_FREQ.keys()):
		if k > 0:
			print 'p(ins length = '+str(abs(k))+' | indel occurs) =',INDEL_FREQ[k]
		else:
			print 'p(del length = '+str(abs(k))+' | indel occurs) =',INDEL_FREQ[k]

	for n in COMMON_VARIANTS:
		print n

	for n in HIGH_MUT_REGIONS:
		print n

	print 'p(snp)   =',SNP_FREQ
	print 'p(indel) =',AVG_INDEL_FREQ
	print 'overall average mut rate:',AVG_MUT_RATE
	print 'total variants processed:',totalVar

	#
	# save variables to file
	#
	OUT_DICT = {'AVG_MUT_RATE':AVG_MUT_RATE,
	            'SNP_FREQ':SNP_FREQ,
	            'INDEL_FREQ':INDEL_FREQ,
	            'TRINUC_MUT_PROB':TRINUC_MUT_PROB,
	            'TRINUC_TRANS_PROBS':TRINUC_TRANS_PROBS,
	            'COMMON_VARIANTS':COMMON_VARIANTS,
	            'HIGH_MUT_REGIONS':HIGH_MUT_REGIONS}
	pickle.dump( OUT_DICT, open( OUT_PICKLE, "wb" ) )


if __name__ == "__main__":
	main()




