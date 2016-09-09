#
#
#            computeGC.py
#            Compute GC and coverage model for genReads.py
#
#            Takes output file from bedtools genomecov to generate GC/coverage model
#
#            Usage: bedtools genomecov -d -ibam normal.bam -g reference.fa 
#                   python computeGC.py -r reference.fa -i genomecovfile -W [sliding window length] -o path/to/output_name.p
#
#


import time
import sys
import argparse
import numpy as np
import cPickle as pickle

parser = argparse.ArgumentParser(description='computeGC.py')
parser.add_argument('-i', type=str, required=True, metavar='<str>', help="input.genomecov")
parser.add_argument('-r', type=str, required=True, metavar='<str>', help="reference.fa")
parser.add_argument('-w', type=int, required=True, metavar='<int>', help="sliding window length")
parser.add_argument('-o', type=str, required=True, metavar='<str>', help="output.p")
args = parser.parse_args()

(IN_GCB, REF_FILE, WINDOW_SIZE, OUT_P) = (args.i, args.r, args.w, args.o)

GC_BINS = {n:[] for n in range(WINDOW_SIZE+1)}

print 'reading ref...'
allRefs = {}
f = open(REF_FILE,'r')
for line in f:
	if line[0] == '>':
		refName = line.strip()[1:]
		allRefs[refName] = []
		print refName
		#if refName == 'chr2':
		#	break
	else:
		allRefs[refName].append(line.strip())
f.close()

print 'capitalizing ref...'
for k in sorted(allRefs.keys()):
	print k
	allRefs[k] = ''.join(allRefs[k])
	allRefs[k] = allRefs[k].upper()

print 'reading genomecov file...'
tt = time.time()
f = open(IN_GCB,'r')
currentLine = 0
currentRef  = None
currentCov  = 0
linesProcessed = 0
PRINT_EVERY    = 1000000
STOP_AFTER     = 1000000
for line in f:
	splt = line.strip().split('\t')
	if linesProcessed%PRINT_EVERY == 0:
		print linesProcessed
	linesProcessed += 1

	if currentLine == 0:
		currentRef = splt[0]
		sPos       = int(splt[1])-1

	if currentRef not in allRefs:
		continue

	currentLine += 1
	currentCov  += int(splt[2])

	if currentLine == WINDOW_SIZE:
		currentLine = 0
		seq         = allRefs[currentRef][sPos:sPos+WINDOW_SIZE]
		if 'N' not in seq:
			gc_count = seq.count('G') + seq.count('C')
			GC_BINS[gc_count].append(currentCov)
		currentCov = 0

	#if linesProcessed >= STOP_AFTER:
	#	break

f.close()

runningTot = 0
allMean    = 0.0
for k in sorted(GC_BINS.keys()):
	if len(GC_BINS[k]) == 0:
		print '{0:0.2%}'.format(k/float(WINDOW_SIZE)), 0.0, 0
		GC_BINS[k] = 0
	else:
		myMean = np.mean(GC_BINS[k])
		myLen  = len(GC_BINS[k])
		print '{0:0.2%}'.format(k/float(WINDOW_SIZE)), myMean, myLen
		allMean += myMean * myLen
		runningTot += myLen
		GC_BINS[k] = myMean

avgCov = allMean/float(runningTot)
print 'AVERAGE COVERAGE =',avgCov

y_out = []
for k in sorted(GC_BINS.keys()):
	GC_BINS[k] /= avgCov
	y_out.append(GC_BINS[k])

print 'saving model...'
pickle.dump([range(WINDOW_SIZE+1),y_out],open(OUT_P,'wb'))

print time.time()-tt,'(sec)'

