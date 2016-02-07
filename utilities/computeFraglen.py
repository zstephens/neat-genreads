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

all_tlens = {}

PRINT_EVERY = 100000
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
				print '---',i, quick_median(all_tlens)
				#for k in sorted(all_tlens.keys()):
				#	print k, all_tlens[k]

for k in sorted(all_tlens.keys()):
	print k, all_tlens[k]
