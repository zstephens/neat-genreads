import sys
import time
import os
import random

OK_CHR_ORD   = {ord('A'):True,ord('C'):True,ord('G'):True,ord('T'):True,ord('U'):True}
ALLOWED_NUCL = ['A','C','G','T']

#
#	Index reference fasta
#
def indexRef(refPath):

	tt = time.time()

	fn = None
	if os.path.isfile(refPath+'i'):
		print 'found index '+refPath+'i'
		fn = refPath+'i'
	if os.path.isfile(refPath+'.fai'):
		print 'found index '+refPath+'.fai'
		fn = refPath+'.fai'

	ref_inds = []
	if fn != None:
		fai = open(fn,'r')
		for line in fai:
			splt = line[:-1].split('\t')
			seqLen = int(splt[1])
			offset = int(splt[2])
			lineLn = int(splt[3])
			nLines = seqLen/lineLn
			if seqLen%lineLn != 0:
				nLines += 1
			ref_inds.append((splt[0],offset,offset+seqLen+nLines,seqLen))
		fai.close()
		return ref_inds

	sys.stdout.write('index not found, creating one... ')
	sys.stdout.flush()
	refFile = open(refPath,'r')
	prevR   = None
	prevP   = None
	seqLen  = 0
	while 1:
		data = refFile.readline()
		if not data:
			ref_inds.append( (prevR, prevP, refFile.tell()-len(data), seqLen) )
			break
		if data[0] == '>':
			if prevP != None:
				ref_inds.append( (prevR, prevP, refFile.tell()-len(data), seqLen) )
			seqLen = 0
			prevP  = refFile.tell()
			prevR  = data[1:-1]
		else:
			seqLen += len(data)-1
	refFile.close()

	print '{0:.3f} (sec)'.format(time.time()-tt)
	return ref_inds


#
#	Read in sequence data from reference fasta
#
#	N_unknowns  = True --> all ambiguous characters will be treated as Ns
#	N_handling  = (mode,params)
#		- ('random',read/frag len)      --> all regions of Ns smaller than read or fragment
#                                           length (whichever is bigger) will be replaced
#                                           with uniformly random nucleotides
#		- ('allChr',read/frag len, chr) --> same as above, but replaced instead with a string
#                                           of 'chr's
#		- ('ignore')                    --> do not alter nucleotides in N regions
#
def readRef(refPath,ref_inds_i,N_handling,N_unknowns=True,quiet=False):

	tt = time.time()
	if not quiet:
		sys.stdout.write('reading '+ref_inds_i[0]+'... ')
		sys.stdout.flush()

	refFile = open(refPath,'r')
	refFile.seek(ref_inds_i[1])
	myDat = ''.join(refFile.read(ref_inds_i[2]-ref_inds_i[1]).split('\n'))
	myDat = bytearray(myDat.upper())

	# find N regions
	# data explanation: myDat[N_atlas[0][0]:N_atlas[0][1]] = solid block of Ns
	prevNI = 0
	nCount = 0
	N_atlas = []
	for i in xrange(len(myDat)):
		if myDat[i] == ord('N') or (N_unknowns and myDat[i] not in OK_CHR_ORD):
			if nCount == 0:
				prevNI = i
			nCount += 1
			if i == len(myDat)-1:
				N_atlas.append((prevNI,prevNI+nCount))
		else:
			if nCount > 0:
				N_atlas.append((prevNI,prevNI+nCount))
			nCount = 0

	# handle N base-calls as desired
	N_info = {}
	N_info['all']   = []
	N_info['big']   = []
	N_info['non_N'] = []
	if N_handling[0] == 'random':
		for region in N_atlas:
			N_info['all'].extend(region)
			if region[1]-region[0] <= N_handling[1]:
				for i in xrange(region[0],region[1]):
					myDat[i] = random.choice(ALLOWED_NUCL)
			else:
				N_info['big'].extend(region)
	elif N_handling[0] == 'allChr' and N_handling[2] in OK_CHR_ORD:
		for region in N_atlas:
			N_info['all'].extend(region)
			if region[1]-region[0] <= N_handling[1]:
				for i in xrange(region[0],region[1]):
					myDat[i] = N_handling[2]
			else:
				N_info['big'].extend(region)
	elif N_handling[0] == 'ignore':
		for region in N_atlas:
			N_info['all'].extend(region)
			N_info['big'].extend(region)
	else:
		print '\nERROR: UNKNOWN N_HANDLING MODE\n'
		exit(1)

	habitableRegions = []
	if N_info['big'] == []:
		N_info['non_N'] = [(0,len(myDat))]
	else:
		for i in xrange(0,len(N_info['big']),2):
			if i == 0:
				habitableRegions.append((0,N_info['big'][0]))
			else:
				habitableRegions.append((N_info['big'][i-1],N_info['big'][i]))
		habitableRegions.append((N_info['big'][-1],len(myDat)))
	for n in habitableRegions:
		if n[0] != n[1]:
			N_info['non_N'].append(n)

	if not quiet:
		print '{0:.3f} (sec)'.format(time.time()-tt)
	return (myDat,N_info)

#
#	find all non-N regions in reference sequence ahead of time, for computing jobs in parallel
#
def getAllRefRegions(refPath,ref_inds,N_handling,saveOutput=False):
	outRegions = {}
	fn = refPath+'.nnr'
	if os.path.isfile(fn) and not(saveOutput):
		print 'found list of preidentified non-N regions...'
		f = open(fn,'r')
		for line in f:
			splt = line.strip().split('\t')
			if splt[0] not in outRegions:
				outRegions[splt[0]] = []
			outRegions[splt[0]].append((int(splt[1]),int(splt[2])))
		f.close()
		return outRegions
	else:
		print 'enumerating all non-N regions in reference sequence...'
		for RI in xrange(len(ref_inds)):
			(refSequence,N_regions) = readRef(refPath,ref_inds[RI],N_handling,quiet=True)
			refName = ref_inds[RI][0]
			outRegions[refName] = [n for n in N_regions['non_N']]
		if saveOutput:
			f = open(fn,'w')
			for k in outRegions.keys():
				for n in outRegions[k]:
					f.write(k+'\t'+str(n[0])+'\t'+str(n[1])+'\n')
			f.close()
		return outRegions

#
#	find which of the non-N regions are going to be used for this job
#
def partitionRefRegions(inRegions,ref_inds,myjob,njobs):

	totSize = 0
	for RI in xrange(len(ref_inds)):
		refName = ref_inds[RI][0]
		for region in inRegions[refName]:
			totSize += region[1] - region[0]
	sizePerJob = int(totSize/float(njobs)-0.5)

	regionsPerJob = [[] for n in xrange(njobs)]
	refsPerJob    = [{} for n in xrange(njobs)]
	currentInd    = 0
	currentCount  = 0
	for RI in xrange(len(ref_inds)):
		refName = ref_inds[RI][0]
		for region in inRegions[refName]:
			regionsPerJob[currentInd].append((refName,region[0],region[1]))
			refsPerJob[currentInd][refName] = True
			currentCount += region[1] - region[0]
			if currentCount >= sizePerJob:
				currentCount = 0
				currentInd   = min([currentInd+1,njobs-1])

	relevantRefs = refsPerJob[myjob-1].keys()
	relevantRegs = regionsPerJob[myjob-1]
	return (relevantRefs,relevantRegs)
