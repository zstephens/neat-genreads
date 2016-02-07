import fileinput
import cPickle as pickle
import numpy as np

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
				print '---',i
				for k in sorted(all_tlens.keys()):
					print k, all_tlens[k]

for k in sorted(all_tlens.keys()):
	print k, all_tlens[k]
