import random
import copy
import re
import os
import cPickle as pickle
import numpy as np

from probability import DiscreteDistribution
from cigar import CigarString

MAX_ATTEMPTS = 100

NUCL    = ['A','C','G','T']
TRI_IND = {'AA':0,  'AC':1,  'AG':2,   'AT':3,  'CA':4,  'CC':5,  'CG':6,  'CT':7,
           'GA':8,  'GC':9,  'GG':10,  'GT':11, 'TA':12, 'TC':13, 'TG':14, 'TT':15}
NUC_IND = {'A':0, 'C':1, 'G':2, 'T':3}

#
#	Container for reference sequences, applies mutations
#
class SequenceContainer:
	def __init__(self, xOffset, sequence, ploidy, windowOverlap, readLen, mutationModels=[], mutRate=None, coverageDat=None):

		self.x         = xOffset
		self.ploidy    = ploidy
		self.readLen   = readLen
		self.sequences = [bytearray(sequence) for n in xrange(self.ploidy)]
		self.seqLen    = len(sequence)
		self.cigars    = [[(self.seqLen,'M')] for n in xrange(self.ploidy)]
		self.indelList = [[] for n in xrange(self.ploidy)]
		self.snpList   = [[] for n in xrange(self.ploidy)]
		self.allCigar  = [[] for n in xrange(self.ploidy)]
		self.adj       = [None for n in xrange(self.ploidy)]
		self.blackList = [np.zeros(self.seqLen,dtype='b') for n in xrange(self.ploidy)]
		
		(self.windowSize, coverage_vals) = coverageDat
		self.win_per_read = int(self.readLen/float(self.windowSize)+0.5)
		self.which_bucket = DiscreteDistribution(coverage_vals,range(len(coverage_vals)))

		self.winBuffer = windowOverlap
		for p in xrange(self.ploidy):
			self.blackList[p][-self.winBuffer] = True
			self.blackList[p][-self.winBuffer-1] = True

		if mutationModels == []:
			ml = [DEFAULT_MODEL_1]*self.ploidy
			self.modelData = ml[:self.ploidy]
		else:
			if len(mutationModels) != self.ploidy:
				print '\nError: Number of mutation models recieved is not equal to specified ploidy\n'
				exit(1)
			self.modelData = mutationModels

		# do we need to rescale mutation frequencies?
		if mutRate == None:
			self.mutScalar = 1.0
		else:
			self.mutScalar = float(mutRate)/sum([n[0] for n in self.modelData])

		# init mutation models
		self.models = []
		for n in self.modelData:
			self.models.append([self.mutScalar*n[0],n[1],n[2],n[3],DiscreteDistribution(n[5],n[4]),DiscreteDistribution(n[7],n[6]),[]])
			for m in n[8]:
				self.models[-1][6].append([DiscreteDistribution(m[0],NUCL),
					                       DiscreteDistribution(m[1],NUCL),
					                       DiscreteDistribution(m[2],NUCL),
					                       DiscreteDistribution(m[3],NUCL)])

	def insert_mutations(self, inputList):
		#
		#	TODO!!!!!! user-input variants, determine which ploid to put it on, etc..
		#
		for inpV in inputList:
			whichPloid = []
			wps = inpV[4][0]
			if wps == None:	# if no genotype given, assume heterozygous and choose a single ploid at random
				whichPloid.append(random.randint(0,self.ploidy-1))
				whichAlt = [0]
			else:
				if 'WP=' in wps:
					whichPloid = [int(n) for n in inpV[-1][3:].split(',') if n == '1']
					whichAlt   = [0]*len(whichPloid)
				elif '/' in wps or '|' in wps:
					if '/' in wps:
						splt = wps.split('/')
					else:
						splt = wps.split('|')
					whichPloid = []
					whichAlt   = []
					for i in xrange(len(splt)):
						if splt[i] == '1':
							whichPloid.append(i)
						whichAlt.append(int(splt[i])-1)
						
			for i in xrange(len(whichPloid)):
				p = whichPloid[i]
				myAlt = inpV[2][whichAlt[i]]
				myVar = (inpV[0]-self.x,inpV[1],myAlt)
				inLen = max([len(inpV[1]),len(myAlt)])
				#print myVar, chr(self.sequences[p][myVar[0]])
				if len(inpV[1]) == 1 and len(myAlt) == 1:
					self.snpList[p].append(myVar)
					self.blackList[p][myVar[0]] = True
				else:
					self.indelList[p].append(myVar)
					for k in xrange(myVar[0],myVar[0]+inLen+1):
						self.blackList[p][k] = True

	def random_mutations(self):
		
		#	add random indels
		indelsToAdd = [int(self.seqLen*n[0]*n[2]) for n in self.models]
		all_indels  = [[] for n in self.sequences]
		for i in xrange(self.ploidy):
			myModel  = self.models[i]
			for j in xrange(indelsToAdd[i]):
				if random.random() <= myModel[1]:	# insert homozygous indel
					whichPloid = range(self.ploidy)
				else:								# insert heterozygous indel
					whichPloid = [random.randint(0,self.ploidy-1)]

				# try to find suitable places to insert indels
				eventPos = -1
				for attempt in xrange(MAX_ATTEMPTS):
					eventPos = random.randint(self.winBuffer,self.seqLen-1)
					for p in whichPloid:
						if self.blackList[p][eventPos]:
							eventPos = -1
					if eventPos != -1:
						break
				if eventPos == -1:
					continue

				if random.random() <= myModel[3]:	# insertion
					inLen   = myModel[4].sample()
					# sequence content of random insertions is uniformly random (change this later)
					inSeq   = ''.join([random.choice(NUCL) for n in xrange(inLen)])
					refNucl = chr(self.sequences[i][eventPos])
					myIndel = (eventPos,refNucl,refNucl+inSeq)
				else:								# deletion
					inLen   = myModel[5].sample()
					if eventPos+inLen+1 >= len(self.sequences[i]):	# skip if deletion too close to boundary
						continue
					if inLen == 1:
						inSeq = chr(self.sequences[i][eventPos+1])
					else:
						inSeq = str(self.sequences[i][eventPos+1:eventPos+inLen+1])
					refNucl = chr(self.sequences[i][eventPos])
					myIndel = (eventPos,refNucl+inSeq,refNucl)

				# if event too close to boundary, skip. if event conflicts with other indel, skip.
				skipEvent = False
				if eventPos+len(myIndel[1]) >= self.seqLen-self.winBuffer-1:
					skipEvent = True
				if skipEvent:
					continue
				for p in whichPloid:
					for k in xrange(eventPos,eventPos+inLen+1):
						if self.blackList[p][k]:
							skipEvent = True
				if skipEvent:
					continue

				for p in whichPloid:
					for k in xrange(eventPos,eventPos+inLen+1):
						self.blackList[p][k] = True
					all_indels[p].append(myIndel)

		for i in xrange(len(all_indels)):
			all_indels[i].extend(self.indelList[i])
		all_indels = [sorted(n,reverse=True) for n in all_indels]
		#print all_indels

		#	add random snps
		snpsToAdd = [int(self.seqLen*n[0]*(1.-n[2])) for n in self.models]
		all_snps  = [[] for n in self.sequences]
		for i in xrange(self.ploidy):
			myModel = self.models[i]
			for j in xrange(snpsToAdd[i]):
				if random.random() <= myModel[1]:	# insert homozygous SNP
					whichPloid = range(self.ploidy)
				else:								# insert heterozygous SNP
					whichPloid = [random.randint(0,self.ploidy-1)]

				# try to find suitable places to insert snps
				eventPos = -1
				for attempt in xrange(MAX_ATTEMPTS):
					eventPos = random.randint(self.winBuffer+1,self.seqLen-2)
					for p in whichPloid:
						if self.blackList[p][eventPos]:
							eventPos = -1
					if eventPos != -1:
						break
				if eventPos == -1:
					continue

				refNucl = chr(self.sequences[i][eventPos])
				context = str(chr(self.sequences[i][eventPos-1])+chr(self.sequences[i][eventPos+1]))
				# sample from tri-nucleotide substitution matrices to get SNP alt allele
				newNucl = myModel[6][TRI_IND[context]][NUC_IND[refNucl]].sample()
				mySNP   = (eventPos,refNucl,newNucl)

				for p in whichPloid:
					all_snps[p].append(mySNP)
					self.blackList[p][mySNP[0]] = True

		for i in xrange(len(all_snps)):
			all_snps[i].extend(self.snpList[i])

		# modify reference sequences
		for i in xrange(len(all_snps)):
			for j in xrange(len(all_snps[i])):
				# sanity checking (for debugging purposes)
				vPos = all_snps[i][j][0]
				if all_snps[i][j][1] != chr(self.sequences[i][vPos]):
					print '\nError: Something went wrong!\n', all_snps[i][j], chr(self.sequences[i][vPos]),'\n'
					exit(1)
				else:
					self.sequences[i][vPos] = all_snps[i][j][2]

		adjToAdd = [[] for n in xrange(self.ploidy)]
		for i in xrange(len(all_indels)):
			for j in xrange(len(all_indels[i])):
				# sanity checking (for debugging purposes)
				vPos  = all_indels[i][j][0]
				vPos2 = vPos + len(all_indels[i][j][1])
				#print all_indels[i][j], str(self.sequences[i][vPos:vPos2])
				#print len(self.sequences[i]),'-->',
				if all_indels[i][j][1] != str(self.sequences[i][vPos:vPos2]):
					print '\nError: Something went wrong!\n', all_indels[i][j], str(self.sequences[i][vPos:vPos2]),'\n'
					exit(1)
				else:
					self.sequences[i] = self.sequences[i][:vPos] + bytearray(all_indels[i][j][2]) + self.sequences[i][vPos2:]
					#print 'rawr:', all_indels[i][j]
					adjToAdd[i].append((all_indels[i][j][0],len(all_indels[i][j][2])-len(all_indels[i][j][1])))
					#adjToAdd.append((vPos,))
				#print len(self.sequences[i])
			adjToAdd[i].sort()
			#print adjToAdd[i]

			self.adj[i] = np.zeros(len(self.sequences[i]),dtype='<i4')
			indSoFar = 0
			valSoFar = 0
			for j in xrange(len(self.adj[i])):
				if indSoFar < len(adjToAdd[i]) and j >= adjToAdd[i][indSoFar][0]+1:
					valSoFar += adjToAdd[i][indSoFar][1]
					indSoFar += 1
				self.adj[i][j] = valSoFar

			# precompute cigar strings

			tempSymbolString = ['M']
			prevVal = self.adj[i][0]
			j = 1
			while j < len(self.adj[i]):
				diff = self.adj[i][j] - prevVal
				prevVal = self.adj[i][j]
				if diff > 0:	# insertion
					tempSymbolString.extend(['I']*abs(diff))
					j += abs(diff)
				elif diff < 0:	# deletion
					tempSymbolString.append('D'*abs(diff)+'M')
					j += 1
				else:
					tempSymbolString.append('M')
					j += 1

			#print len(self.sequences[i]), len(tempSymbolString)
			#for j in xrange(len(tempSymbolString)):
			#	print j, tempSymbolString[j]

			for j in xrange(len(tempSymbolString)-self.readLen):
				self.allCigar[i].append(CigarString(listIn=tempSymbolString[j:j+self.readLen]).getString())


		# tally up variants implemented
		countDict = {}
		all_variants = [sorted(all_snps[i]+all_indels[i]) for i in xrange(self.ploidy)]
		for i in xrange(len(all_variants)):
			for j in xrange(len(all_variants[i])):
				all_variants[i][j] = tuple([all_variants[i][j][0]+self.x])+all_variants[i][j][1:]
				t = tuple(all_variants[i][j])
				if t not in countDict:
					countDict[t] = []
				countDict[t].append(i)

		#
		#	TODO: combine multiple variants that happened to occur at same position into single vcf entry
		#

		output_variants = []
		for k in sorted(countDict.keys()):
			output_variants.append(k+tuple([len(countDict[k])/float(self.ploidy)]))
			output_variants[-1] += tuple(['WP='+','.join([str(n) for n in countDict[k]])])
		#for n in output_variants:
		#	print n
		return output_variants


	def sample_read(self, sequencingModel, fragLen=None):
		
		# choose a ploid
		myPloid = random.randint(0,self.ploidy-1)

		# choose a random position within the ploid, and generate quality scores / sequencing errors
		readsToSample = []
		if fragLen == None:
			#rPos = random.randint(0,len(self.sequences[myPloid])-self.readLen-1)	# uniform random

			# decide which subsection of the sequence to sample from using coverage probabilities
			coords_bad = True
			while coords_bad:
				myBucket = max([self.which_bucket.sample() - self.win_per_read, 0])
				coords_to_select_from = [myBucket*self.windowSize,(myBucket+1)*self.windowSize]
				coords_to_select_from[0] += self.adj[myPloid][coords_to_select_from[0]]
				coords_to_select_from[1] += self.adj[myPloid][coords_to_select_from[1]]
				if coords_to_select_from[1] < len(self.sequences[myPloid])-self.readLen:
					coords_bad = False
			rPos = random.randint(coords_to_select_from[0],coords_to_select_from[1]-1)

			# sample read position and call function to compute quality scores / sequencing errors
			rDat = self.sequences[myPloid][rPos:rPos+self.readLen]
			(myQual, myErrors) = sequencingModel.getSequencingErrors(rDat)
			readsToSample.append([rPos,myQual,myErrors,rDat])

		else:
			#rPos1 = random.randint(0,len(self.sequences[myPloid])-fragLen-1)		# uniform random

			# decide which subsection of the sequence to sample from using coverage probabilities
			coords_bad = True
			while coords_bad:
				myBucket = max([self.which_bucket.sample() - self.win_per_read, 0])
				coords_to_select_from = [myBucket*self.windowSize,(myBucket+1)*self.windowSize]
				coords_to_select_from[0] += self.adj[myPloid][coords_to_select_from[0]]
				coords_to_select_from[1] += self.adj[myPloid][coords_to_select_from[0]]	# both ends use index of starting position to avoid issues with reads spanning breakpoints of large events
				rPos1 = random.randint(coords_to_select_from[0],coords_to_select_from[1]-1)
				# for PE-reads, flip a coin to decide if R1 or R2 will be the "covering" read
				if random.randint(1,2) == 1 and rPos1 > fragLen - self.readLen:
					rPos1 -= fragLen - self.readLen
				if rPos1 < len(self.sequences[myPloid])-fragLen:
					coords_bad = False

			rPos2 = rPos1 + fragLen - self.readLen
			rDat1 = self.sequences[myPloid][rPos1:rPos1+self.readLen]
			rDat2 = self.sequences[myPloid][rPos2:rPos2+self.readLen]
			(myQual1, myErrors1) = sequencingModel.getSequencingErrors(rDat1)
			(myQual2, myErrors2) = sequencingModel.getSequencingErrors(rDat2,isReverseStrand=True)
			readsToSample.append([rPos1,myQual1,myErrors1,rDat1])
			readsToSample.append([rPos2,myQual2,myErrors2,rDat2])

		# error format:
		# myError[i] = (type, len, pos, ref, alt)

		# examine sequencing errors to-be-inserted.
		#	- remove deletions that don't have enough bordering sequence content to "fill in"
		# if error is valid, make the changes to the read data
		rOut = []
		for r in readsToSample:
			myCigar = self.allCigar[myPloid][r[0]]
			totalD  = sum([error[1] for error in r[2] if error[0] == 'D'])
			totalI  = sum([error[1] for error in r[2] if error[0] == 'I'])
			availB  = len(self.sequences[myPloid]) - r[0] - self.readLen - 1
			# add buffer sequence to fill in positions that get deleted
			r[3] += self.sequences[myPloid][r[0]+self.readLen:r[0]+self.readLen+totalD]
			expandedCigar = []
			extraCigar    = []
			adj           = 0
			sse_adj       = [0 for n in xrange(self.readLen)]
			anyIndelErr   = False

			# sort by letter (D > I > S) such that we introduce all indel errors before substitution errors
			# secondarily, sort by index
			arrangedErrors = {'D':[],'I':[],'S':[]}
			for error in r[2]:
				arrangedErrors[error[0]].append((error[2],error))
			sortedErrors = []
			for k in sorted(arrangedErrors.keys()):
				sortedErrors.extend([n[1] for n in sorted(arrangedErrors[k])])

			for error in sortedErrors:
				#print r[0], error
				eLen = error[1]
				ePos = error[2]
				if error[0] == 'D' or error[0] == 'I':
					anyIndelErr   = True
					extraCigarVal = []
					if totalD > availB:	# if not enough bases to fill-in deletions, skip all indel erors
						continue
					if expandedCigar == []:
						expandedCigar = CigarString(stringIn=myCigar).getList()

						fillToGo = totalD - totalI
						if fillToGo > 0:
							extraCigarVal = CigarString(stringIn=self.allCigar[myPloid][r[0]+fillToGo]).getList()[-fillToGo:]

					# insert deletion error into read and update cigar string accordingly
					if error[0] == 'D':
						pi = ePos+adj
						pf = ePos+adj+eLen+1
						if str(r[3][pi:pf]) == str(error[3]):
							r[3] = r[3][:pi+1] + r[3][pf:]
							expandedCigar = expandedCigar[:pi+1] + expandedCigar[pf:]
							expandedCigar[pi+1] = 'D'*eLen + expandedCigar[pi+1]
						else:
							print '\nError, ref does not match alt while attempting to insert deletion error!\n'
							exit(1)
						adj -= eLen
						for i in xrange(ePos,len(sse_adj)):
							sse_adj[i] -= eLen

					# insert insertion error into read and update cigar string accordingly
					else:
						if chr(r[3][ePos+adj]) == error[3]:
							r[3] = r[3][:ePos+adj] + error[4] + r[3][ePos+adj+1:]
							expandedCigar = expandedCigar[:ePos+adj] + ['I']*eLen + expandedCigar[ePos+adj+1:]
						else:
							print '\nError, ref does not match alt while attempting to insert insertion error!\n'
							exit(1)
						adj += eLen
						for i in xrange(ePos,len(sse_adj)):
							sse_adj[i] += eLen

				else:	# substitution errors, much easier by comparison...
					if chr(r[3][ePos+sse_adj[ePos]]) == error[3]:
						r[3][ePos+sse_adj[ePos]] = error[4]
					else:
						print '\nError, ref does not match alt while attempting to insert substitution error!\n'
						exit(1)

			if anyIndelErr:
				if len(expandedCigar):
					#print myCigar,'-->',
					relevantCigar = (expandedCigar+extraCigarVal)[:self.readLen]
					myCigar = CigarString(listIn=relevantCigar).getString()
					#print myCigar

				r[3] = r[3][:self.readLen]
				#if len(r[3]) != self.readLen:
				#	print 'AHHHHHH_1'
				#	exit(1)
				#if len(expandedCigar+extraCigarVal) != self.readLen:
				#	print 'AHHHHHH_2'
				#	exit(1)

			rOut.append([r[0]-self.adj[myPloid][r[0]],myCigar,str(r[3]),str(r[1])])

		# rOut[i] = (pos, cigar, read_string, qual_string)
		return rOut


#
#	Container for read data, computes quality scores and positions to insert errors
#
class ReadContainer:
	def __init__(self, readLen, errorModel, reScaledError):

		self.readLen = readLen

		errorDat = pickle.load(open(errorModel,'rb'))
		if len(errorDat) == 6:		# only 1 q-score model present, use same model for both strands
			[initQ1,probQ1,Qscores,offQ,avgError,errorParams] = errorDat
			self.PE_MODELS = False
		elif len(errorDat) == 8:	# found a q-score model for both forward and reverse strands
			#print 'Using paired-read quality score profiles...'
			[initQ1,probQ1,initQ2,probQ2,Qscores,offQ,avgError,errorParams] = errorDat
			self.PE_MODELS = True
			if len(initQ1) != len(initQ2) or len(probQ1) != len(probQ2):
				print '\nError: R1 and R2 quality score models are of different length.\n'
				exit(1)

		self.qErrRate = [0.]*(max(Qscores)+1)
		for q in Qscores:
			self.qErrRate[q] = 10.**(-q/10.)
		self.offQ = offQ

		# errorParams = [SSE_PROB, SIE_RATE, SIE_PROB, SIE_VAL, SIE_INS_FREQ, SIE_INS_NUCL]
		self.errP   = errorParams
		self.errSSE = [DiscreteDistribution(n,NUCL) for n in self.errP[0]]
		self.errSIE = DiscreteDistribution(self.errP[2],self.errP[3])
		self.errSIN = DiscreteDistribution(self.errP[5],NUCL)

		# adjust length to match desired read length
		if self.readLen == len(initQ1):
			self.qIndRemap = range(self.readLen)
		else:
			print 'Warning: Read length of error model ('+str(len(initQ1))+') does not match -R value ('+str(self.readLen)+'), rescaling model...'
			self.qIndRemap = [max([1,len(initQ1)*n/readLen]) for n in xrange(readLen)]

		# adjust sequencing error frequency to match desired rate
		if reScaledError == None:
			self.errorScale = 1.0
		else:
			self.errorScale = reScaledError/avgError
			print 'Warning: Quality scores no longer exactly representative of error probability. Error model scaled by {0:.3f} to match desired rate...'.format(self.errorScale)

		# initialize probability distributions
		self.initDistByPos1        = [DiscreteDistribution(initQ1[i],Qscores) for i in xrange(len(initQ1))]
		self.probDistByPosByPrevQ1 = [None]
		for i in xrange(1,len(initQ1)):
			self.probDistByPosByPrevQ1.append([])
			for j in xrange(len(initQ1[0])):
				if np.sum(probQ1[i][j]) <= 0.:	# if we don't have sufficient data for a transition, use the previous qscore
					self.probDistByPosByPrevQ1[-1].append(DiscreteDistribution([1],[Qscores[j]],degenerateVal=Qscores[j]))
				else:
					self.probDistByPosByPrevQ1[-1].append(DiscreteDistribution(probQ1[i][j],Qscores))

		if self.PE_MODELS:
			self.initDistByPos2        = [DiscreteDistribution(initQ2[i],Qscores) for i in xrange(len(initQ2))]
			self.probDistByPosByPrevQ2 = [None]
			for i in xrange(1,len(initQ2)):
				self.probDistByPosByPrevQ2.append([])
				for j in xrange(len(initQ2[0])):
					if np.sum(probQ2[i][j]) <= 0.:	# if we don't have sufficient data for a transition, use the previous qscore
						self.probDistByPosByPrevQ2[-1].append(DiscreteDistribution([1],[Qscores[j]],degenerateVal=Qscores[j]))
					else:
						self.probDistByPosByPrevQ2[-1].append(DiscreteDistribution(probQ2[i][j],Qscores))

	def getSequencingErrors(self, readData, isReverseStrand=False):

		qOut = [0]*self.readLen
		sErr = []

		if self.PE_MODELS and isReverseStrand:
			myQ = self.initDistByPos2[0].sample()
		else:
			myQ = self.initDistByPos1[0].sample()

		if random.random() < self.qErrRate[myQ]:
			sErr.append(0)
		qOut[0] = myQ + self.offQ
		for i in xrange(1,self.readLen):

			if self.PE_MODELS and isReverseStrand:
				myQ = self.probDistByPosByPrevQ2[self.qIndRemap[i]][myQ].sample()
			else:
				myQ = self.probDistByPosByPrevQ1[self.qIndRemap[i]][myQ].sample()

			if random.random() < self.errorScale*self.qErrRate[myQ]:
				sErr.append(i)
			qOut[i] = myQ + self.offQ
		qOut = ''.join([chr(n) for n in qOut])

		sOut = []
		nDelSoFar = 0
		# don't allow indel errors to occur on subsequent positions
		prevIndel = -2
		# don't allow other sequencing errors to occur on bases removed by deletion errors
		delBlacklist = []

		for ind in sErr[::-1]:	# for each error that we're going to insert...

			# determine error type
			isSub = True
			if ind != 0 and ind != self.readLen-1-max(self.errP[3]) and ind > prevIndel+1:
				if random.random() < self.errP[1]:
					isSub = False

			# errorOut = (type, len, pos, ref, alt)

			if isSub:								# insert substitution error
				myNucl  = chr(readData[ind])
				newNucl = self.errSSE[NUC_IND[myNucl]].sample()
				sOut.append(('S',1,ind,myNucl,newNucl))
			else:									# insert indel error
				indelLen = self.errSIE.sample()
				if random.random() < self.errP[4]:		# insertion error
					myNucl  = chr(readData[ind])
					newNucl = myNucl + ''.join([self.errSIN.sample() for n in xrange(indelLen)])
					sOut.append(('I',len(newNucl)-1,ind,myNucl,newNucl))
				elif ind < self.readLen-2-nDelSoFar:	# deletion error (prevent too many of them from stacking up)
					myNucl  = str(readData[ind:ind+indelLen+1])
					newNucl = chr(readData[ind])
					nDelSoFar += len(myNucl)-1
					sOut.append(('D',len(myNucl)-1,ind,myNucl,newNucl))
					for i in xrange(ind+1,ind+indelLen+1):
						delBlacklist.append(i)
				prevIndel = ind

		# remove blacklisted errors
		for i in xrange(len(sOut)-1,-1,-1):
			if sOut[i][2] in delBlacklist:
				del sOut[i]

		return (qOut,sOut)



"""************************************************
****          DEFAULT MUTATION MODELS
************************************************"""


# parse mutation model pickle file
def parseInputMutationModel(model=None, whichDefault=1):
	if whichDefault == 1:
		outModel = [copy.deepcopy(n) for n in DEFAULT_MODEL_1]
	elif whichDefault == 2:
		outModel = [copy.deepcopy(n) for n in DEFAULT_MODEL_2]
	else:
		print '\nError: Unknown default mutation model specified\n'
		exit(1)

	if model != None:
		pickle_dict = pickle.load(open(model,"rb"))
		outModel[0] = pickle_dict['AVG_MUT_RATE']
		outModel[2] = 1. - pickle_dict['SNP_FREQ']

		insList     = pickle_dict['INDEL_FREQ']
		if len(insList):
			insCount = sum([insList[k] for k in insList.keys() if k >= 1])
			delCount = sum([insList[k] for k in insList.keys() if k <= -1])
			insVals  = [k for k in sorted(insList.keys()) if k >= 1]
			insWght  = [insList[k]/float(insCount) for k in insVals]
			delVals  = [k for k in sorted([abs(k) for k in insList.keys() if k <= -1])]
			delWght  = [insList[-k]/float(delCount) for k in delVals]
		else:	# degenerate case where no indel stats are provided
			insCount = 1
			delCount = 1
			insVals  = [1]
			insWght  = [1.0]
			delVals  = [1]
			delWght  = [1.0]
		outModel[3] = insCount/float(insCount + delCount)
		outModel[4] = insVals
		outModel[5] = insWght
		outModel[6] = delVals
		outModel[7] = delWght

		trinuc_trans_prob = pickle_dict['TRINUC_TRANS_PROBS']
		for k in sorted(trinuc_trans_prob.keys()):
			myInd   = TRI_IND[k[0][0]+k[0][2]]
			(k1,k2) = (NUC_IND[k[0][1]],NUC_IND[k[1][1]])
			outModel[8][myInd][k1][k2] = trinuc_trans_prob[k]
		for i in xrange(len(outModel[8])):
			for j in xrange(len(outModel[8][i])):
				for l in xrange(len(outModel[8][i][j])):
					outModel[8][i][j][l] /= float(sum(outModel[8][i][j]))

	return outModel


# parse mutation model files, returns default model if no model directory is specified
#
# OLD FUNCTION THAT PROCESSED OUTDATED TEXTFILE MUTATION MODELS
def parseInputMutationModel_deprecated(prefix=None, whichDefault=1):
	if whichDefault == 1:
		outModel = [copy.deepcopy(n) for n in DEFAULT_MODEL_1]
	elif whichDefault == 2:
		outModel = [copy.deepcopy(n) for n in DEFAULT_MODEL_2]
	else:
		print '\nError: Unknown default mutation model specified\n'
		exit(1)

	if prefix != None:
		if prefix[-1] != '/':
			prefix += '/'
		if not os.path.isdir(prefix):
			'\nError: Input mutation model directory not found:',prefix,'\n'
			exit(1)

		print 'Reading in mutation model...'
		listing1 = [n for n in os.listdir(prefix) if n[-5:] == '.prob']
		listing2 = [n for n in os.listdir(prefix) if n[-7:] == '.trinuc']
		listing  = sorted(listing1) + sorted(listing2)
		for l in listing:
			f = open(prefix+l,'r')
			fr = [n.split('\t') for n in f.read().split('\n')]
			f.close()

			if '_overall.prob' in l:
				myIns = None
				myDel = None
				for dat in fr[1:]:
					if len(dat) == 2:
						if dat[0] == 'insertion':
							myIns = float(dat[1])
						elif dat[0] == 'deletion':
							myDel = float(dat[1])
				if myIns != None and myDel != None:
					outModel[2] = myIns + myDel
					outModel[3] = myIns / (myIns + myDel)
					print '-',l

			if '_insLength.prob' in l:
				insVals = {}
				for dat in fr[1:]:
					if len(dat) == 2:
						insVals[int(dat[0])] = float(dat[1])
				if len(insVals):
					outModel[4] = sorted(insVals.keys())
					outModel[5] = [insVals[n] for n in outModel[4]]
					print '-',l

			if '_delLength.prob' in l:
				delVals = {}
				for dat in fr[1:]:
					if len(dat) == 2:
						delVals[int(dat[0])] = float(dat[1])
				if len(delVals):
					outModel[6] = sorted(delVals.keys())
					outModel[7] = [delVals[n] for n in outModel[6]]
					print '-',l

			if '.trinuc' == l[-7:]:
				context_ind = TRI_IND[l[-10]+l[-8]]
				p_matrix    = [[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]]
				for i in xrange(len(p_matrix)):
					for j in xrange(len(fr[i])):
						p_matrix[i][j] = float(fr[i][j])
				anyNone = False
				for i in xrange(len(p_matrix)):
					for j in xrange(len(p_matrix[i])):
						if p_matrix[i][j] == -1:
							anyNone = True
				if not anyNone:
					outModel[8][context_ind] = copy.deepcopy(p_matrix)
					print '-',l

	return outModel



DEFAULT_1_OVERALL_MUT_RATE   = 0.001
DEFAULT_1_HOMOZYGOUS_FREQ    = 0.010
DEFAULT_1_INDEL_FRACTION     = 0.005
DEFAULT_1_INS_VS_DEL         = 0.6

DEFAULT_1_INS_LENGTH_VALUES  = [1,2,3,4,5,6,7,8,9,10]
DEFAULT_1_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033]
DEFAULT_1_DEL_LENGTH_VALUES  = [1,2,3,4,5]
DEFAULT_1_DEL_LENGTH_WEIGHTS = [0.3,0.2,0.2,0.2,0.1]

example_matrix_1   = [[0.0, 0.15, 0.7, 0.15],
				      [0.15, 0.0, 0.15, 0.7],
				      [0.7, 0.15, 0.0, 0.15],
				      [0.15, 0.7, 0.15, 0.0]]

DEFAULT_1_TRI_FREQS  = [example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1,
						example_matrix_1]

DEFAULT_MODEL_1 = [DEFAULT_1_OVERALL_MUT_RATE,
				   DEFAULT_1_HOMOZYGOUS_FREQ,
				   DEFAULT_1_INDEL_FRACTION,
				   DEFAULT_1_INS_VS_DEL,
				   DEFAULT_1_INS_LENGTH_VALUES,
				   DEFAULT_1_INS_LENGTH_WEIGHTS,
				   DEFAULT_1_DEL_LENGTH_VALUES,
				   DEFAULT_1_DEL_LENGTH_WEIGHTS,
				   DEFAULT_1_TRI_FREQS]


DEFAULT_2_OVERALL_MUT_RATE   = 0.002
DEFAULT_2_HOMOZYGOUS_FREQ    = 0.200
DEFAULT_2_INDEL_FRACTION     = 0.1
DEFAULT_2_INS_VS_DEL         = 0.3

DEFAULT_2_INS_LENGTH_VALUES  = [1,2,3,4,5,6,7,8,9,10]
DEFAULT_2_INS_LENGTH_WEIGHTS = [0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
DEFAULT_2_DEL_LENGTH_VALUES  = [1,2,3,4,5]
DEFAULT_2_DEL_LENGTH_WEIGHTS = [0.3,0.2,0.2,0.2,0.1]

example_matrix_2   = [[0.0, 0.15, 0.7, 0.15],
				      [0.15, 0.0, 0.15, 0.7],
				      [0.7, 0.15, 0.0, 0.15],
				      [0.15, 0.7, 0.15, 0.0]]

DEFAULT_2_TRI_FREQS  = [example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2,
						example_matrix_2]

DEFAULT_MODEL_2 = [DEFAULT_2_OVERALL_MUT_RATE,
				   DEFAULT_2_HOMOZYGOUS_FREQ,
				   DEFAULT_2_INDEL_FRACTION,
				   DEFAULT_2_INS_VS_DEL,
				   DEFAULT_2_INS_LENGTH_VALUES,
				   DEFAULT_2_INS_LENGTH_WEIGHTS,
				   DEFAULT_2_DEL_LENGTH_VALUES,
				   DEFAULT_2_DEL_LENGTH_WEIGHTS,
				   DEFAULT_2_TRI_FREQS]


