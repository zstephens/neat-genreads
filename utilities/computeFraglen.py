#
#
#      Compute Fragment Length Model for genReads.py
#                  computeFraglen.py
#
#
#      Usage: samtools view normal.bam | python computeFraglen.py
#
#

import sys
import fileinput
import cPickle as pickle
import numpy as np

FILTER_MAPQUAL  = 10	# only consider reads that are mapped with at least this mapping quality
FILTER_MINREADS = 100	# only consider fragment lengths that have at least this many read pairs supporting it
FILTER_MEDDEV_M = 10	# only consider fragment lengths this many median deviations above the median

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

def median_deviation_from_median(countDict):
	myMedian = quick_median(countDict)
	deviations = {}
	for k in sorted(countDict.keys()):
		d = abs(k-myMedian)
		deviations[d] = countDict[k]
	return quick_median(deviations)

if len(sys.argv) != 1:
	print "Usage: samtools view normal.bam | python computeFraglen.py"
	exit(1)

all_tlens = {}
PRINT_EVERY = 100000
BREAK_AFTER = 1000000
i = 0
for line in fileinput.input():
	splt = line.strip().split('\t')
	samFlag = int(splt[1])
	myRef   = splt[2]
	mapQual = int(splt[4])
	mateRef = splt[6]
	myTlen  = abs(int(splt[8]))

	if samFlag&1 and samFlag&64 and mapQual > FILTER_MAPQUAL:	# if read is paired, and is first in pair, and is confidently mapped...
		if mateRef == '=' or mateRef == myRef:					# and mate is mapped to same reference
			if myTlen not in all_tlens:
				all_tlens[myTlen] = 0
			all_tlens[myTlen] += 1
			i += 1
			if i%PRINT_EVERY == 0:
				print '---',i, quick_median(all_tlens), median_deviation_from_median(all_tlens)
				#for k in sorted(all_tlens.keys()):
				#	print k, all_tlens[k]

			#if i > BREAK_AFTER:
			#	break


med = quick_median(all_tlens)
mdm = median_deviation_from_median(all_tlens)

outVals  = []
outProbs = []
for k in sorted(all_tlens.keys()):
	if k > 0 and k < med + FILTER_MEDDEV_M * mdm:
		if all_tlens[k] >= FILTER_MINREADS:
			print k, all_tlens[k]
			outVals.append(k)
			outProbs.append(all_tlens[k])
countSum = float(sum(outProbs))
outProbs = [n/countSum for n in outProbs]

print '\nsaving model...'
pickle.dump([outVals, outProbs],open('fraglen.p','wb'))


