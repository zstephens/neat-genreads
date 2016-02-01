import time
import sys
import numpy as np

REF_FILE = '/home/n-z/zstephe2/plant_refs/zzz_cleaned/hg19_clean.fa'
REF_FILE = sys.argv[1]
IN_GCB   = '/home/n-z/zstephe2/OICR/Pr_P_PE_635_WG_150224.bam.genomecov'
IN_GCB   = sys.argv[2]

WINDOW_SIZE = 50
WINDOW_SIZE = int(sys.argv[3])

GC_BINS = {n:[] for n in range(WINDOW_SIZE+1)}

print 'reading ref...'
allRefs = {}
f = open(REF_FILE,'r')
for line in f:
	if line[0] == '>':
		refName = line.strip()[1:]
		allRefs[refName] = []
		print refName
		if refName == 'chr2':
			break
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

		if currentRef == 'chr2':
			break

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
	else:
		myMean = np.mean(GC_BINS[k])
		myLen  = len(GC_BINS[k])
		print '{0:0.2%}'.format(k/float(WINDOW_SIZE)), myMean, myLen
		allMean += myMean * myLen
		runningTot += myLen

avgCov = allMean/float(runningTot)
print 'AVERAGE COVERAGE =',avgCov

print time.time()-tt,'(sec)'

