import os
import pickle
import argparse

parser = argparse.ArgumentParser(description='NEAT-genReads V2.0')
parser.add_argument('-r', type=str, required=True, metavar='<str>', help="* ref.fa")
parser.add_argument('-m', type=str, required=True, metavar='<str>', help="* mutations.tsv")
parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output.p")
args = parser.parse_args()
(REF, TSV, OUT_PICKLE) = (args.r, args.m, args.o)

REF_WHITELIST =  [str(n) for n in xrange(1,30)] + ['x','y','X','Y','mt','Mt','MT']
REF_WHITELIST += ['chr'+n for n in REF_WHITELIST]

def getRefList(refPath):

	fn = None
	if os.path.isfile(refPath+'i'):
		#print 'found index '+refPath+'i'
		fn = refPath+'i'
	if os.path.isfile(refPath+'.fai'):
		#print 'found index '+refPath+'.fai'
		fn = refPath+'.fai'

	refList = []
	if fn != None:
		fai = open(fn,'r')
		for line in fai:
			splt = line[:-1].split('\t')
			refList.append(splt[0])
	return refList

def getChrFromFasta(refPath,chrName):

	fn = None
	if os.path.isfile(refPath+'i'):
		#print 'found index '+refPath+'i'
		fn = refPath+'i'
	if os.path.isfile(refPath+'.fai'):
		#print 'found index '+refPath+'.fai'
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
	else:
		print 'go index your reference, you slob!'
		exit(1)

	for i in xrange(len(ref_inds)):
		if ref_inds[i][0] == chrName:
			ref_inds_i = ref_inds[i]
			break

	refFile = open(refPath,'r')
	refFile.seek(ref_inds_i[1])
	myDat = ''.join(refFile.read(ref_inds_i[2]-ref_inds_i[1]).split('\n'))
	return myDat


#####################################
#				main()				#
#####################################


def main():

	refList = getRefList(REF)

	TRINUC_REF_COUNT = {}
	# [(trinuc_a, trinuc_b)] = # of times we observed a mutation from trinuc_a into trinuc_b
	TRINUC_TRANSITION_COUNT = {}

	# load and process variants in each reference sequence individually, for memory reasons...
	for refName in refList:

		if refName not in REF_WHITELIST:
			print refName,'is not in our whitelist, skipping...'
			continue

		print 'reading reference "'+refName+'"...'
		refSequence = getChrFromFasta(REF,refName).upper()


		""" ##########################################################################
		###						COUNT TRINUCLEOTIDES IN REF						   ###
		########################################################################## """


		print 'counting trinucleotides in reference...'
		for i in xrange(len(refSequence)-2):
			if i%1000000 == 0:
				print i,'/',len(refSequence)
			#if i == 1000000:
			#	break
			trinuc = refSequence[i:i+3]
			if 'N' in trinuc:
				continue
			if trinuc not in TRINUC_REF_COUNT:
				TRINUC_REF_COUNT[trinuc] = 0
			TRINUC_REF_COUNT[trinuc] += 1


		""" ##########################################################################
		###						COUNT TRINUCLEOTIDE TRANSITIONS					   ###
		########################################################################## """


		print 'counting trinucleotide transitions...'
		f = open(TSV,'r')
		isFirst = True
		for line in f:

			if isFirst:
				splt = line.strip().split('\t')
				(c1,c2,c3) = (splt.index('chromosome'),splt.index('chromosome_start'),splt.index('chromosome_end'))
				(m1,m2,m3) = (splt.index('reference_genome_allele'),splt.index('mutated_from_allele'),splt.index('mutated_to_allele'))
				isFirst = False
				continue

			splt = line.strip().split('\t')
			# we have -1 because tsv coords are 1-based, and our reference string index is 0-based
			[chrName,chrStart,chrEnd] = [splt[c1],int(splt[c2])-1,int(splt[c3])-1]
			[allele_ref,allele_normal,allele_tumor] = [splt[m1].upper(),splt[m2].upper(),splt[m3].upper()]

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
				if 'N' in trinuc_ref:
					continue
				# only consider positions where ref allele in tsv matches the nucleotide in our reference
				if allele_ref == trinuc_ref[1]:
					trinuc_normal    = refSequence[chrStart-1] + allele_normal + refSequence[chrStart+1]
					trinuc_tumor     = refSequence[chrStart-1] + allele_tumor + refSequence[chrStart+1]
					key = (trinuc_normal,trinuc_tumor)
					if key not in TRINUC_TRANSITION_COUNT:
						TRINUC_TRANSITION_COUNT[key] = 0
					TRINUC_TRANSITION_COUNT[key] += 1
				else:
					print trinuc, allele_ref, allele_normal, allele_tumor

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
	
	#
	#	print some stuff
	#
	for k in sorted(TRINUC_MUT_PROB.keys()):
		print 'p('+k+' mutates) =',TRINUC_MUT_PROB[k]

	for k in sorted(TRINUC_TRANS_PROBS.keys()):
		print 'p('+k[0]+' mutates into '+k[1]+') =',TRINUC_TRANS_PROBS[k]

	#
	# save variables to file
	#
	pickle.dump( [TRINUC_MUT_PROB, TRINUC_TRANS_PROBS], open( OUT_PICKLE, "wb" ) )


if __name__ == "__main__":
	main()




