import fileinput
import cPickle as pickle
import numpy as np

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

all_tlens = {}

PRINT_EVERY = 100000
BREAK_AFTER = 1000000
i = 0
for line in fileinput.input():
	splt = line.strip().split('\t')
	samFlag = int(splt[1])
	myRef   = splt[2]
	mateRef = splt[6]
	myTlen  = abs(int(splt[8]))

	if samFlag&1 and samFlag&64:					# if read is paired, and is first in pair...
		if mateRef == '=' or mateRef == myRef:		# and mate is mapped to same reference
			if myTlen not in all_tlens:
				all_tlens[myTlen] = 0
			all_tlens[myTlen] += 1
			i += 1
			if i%PRINT_EVERY == 0:
				print '---',i, quick_median(all_tlens), median_deviation_from_median(all_tlens)
				#for k in sorted(all_tlens.keys()):
				#	print k, all_tlens[k]

			if i > BREAK_AFTER:
				break



for k in sorted(all_tlens.keys()):
	print k, all_tlens[k]
