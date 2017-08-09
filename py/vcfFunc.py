import sys
import time
import os
import re
import random

INCLUDE_HOMS = False
INCLUDE_FAIL = False
CHOOSE_RANDOM_PLOID_IF_NO_GT_FOUND = True

def parseLine(splt,colDict,colSamp):
	
	#	check if we want to proceed..
	ra = splt[colDict['REF']]
	aa = splt[colDict['ALT']]
	# enough columns?
	if len(splt) != len(colDict):
		return None
	# exclude homs / filtered?
	if not(INCLUDE_HOMS) and (aa == '.' or aa == '' or aa == ra):
		return None
	if not(INCLUDE_FAIL) and (splt[colDict['FILTER']] != 'PASS' and splt[colDict['FILTER']] != '.'):
		return None

	#	default vals
	alt_alleles = [aa]
	alt_freqs   = []
	
	gt_perSamp  = []

	#	any alt alleles?
	alt_split = aa.split(',')
	if len(alt_split) > 1:
		alt_alleles = alt_split

	#	check INFO for AF
	af = None
	if 'INFO' in colDict and ';AF=' in ';'+splt[colDict['INFO']]:
		info = splt[colDict['INFO']]+';'
		af   = re.findall(r"AF=.*?(?=;)",info)[0][3:]
	if af != None:
		af_splt = af.split(',')
		while(len(af_splt) < len(alt_alleles)):	# are we lacking enough AF values for some reason?
			af_splt.append(af_splt[-1])			# phone it in.
		if len(af_splt) != 0 and af_splt[0] != '.' and af_splt[0] != '':		# missing data, yay
			alt_freqs = [float(n) for n in af_splt]
	else:
		alt_freqs = [None]*max([len(alt_alleles),1])

	gt_perSamp = None
	#	if available (i.e. we simulated it) look for WP in info
	if len(colSamp) == 0 and 'INFO' in colDict and 'WP=' in splt[colDict['INFO']]:
		info       = splt[colDict['INFO']]+';'
		gt_perSamp = [re.findall(r"WP=.*?(?=;)",info)[0][3:]]
	else:
		#	if no sample columns, check info for GT
		if len(colSamp) == 0 and 'INFO' in colDict and 'GT=' in splt[colDict['INFO']]:
			info       = splt[colDict['INFO']]+';'
			gt_perSamp = [re.findall(r"GT=.*?(?=;)",info)[0][3:]]
		elif len(colSamp):
			fmt = ':'+splt[colDict['FORMAT']]+':'
			if ':GT:' in fmt:
				gtInd = fmt.split(':').index('GT')
				gt_perSamp = [splt[colSamp[iii]].split(':')[gtInd-1] for iii in xrange(len(colSamp))]
				for i in xrange(len(gt_perSamp)):
					gt_perSamp[i] = gt_perSamp[i].replace('.','0')
		if gt_perSamp == None:
			gt_perSamp = [None]*max([len(colSamp),1])

	return (alt_alleles, alt_freqs, gt_perSamp)



def parseVCF(vcfPath,tumorNormal=False,ploidy=2):

	tt = time.time()
	print '--------------------------------'
	sys.stdout.write('reading input VCF...\n')
	sys.stdout.flush()
	
	colDict   = {}
	colSamp   = []
	nSkipped  = 0
	nSkipped_becauseHash = 0
	allVars   = {}	# [ref][pos]
	sampNames = []
	alreadyPrintedWarning = False
	f = open(vcfPath,'r')
	for line in f:

		if line[0] != '#':
			if len(colDict) == 0:
				print '\n\nERROR: VCF has no header?\n'+VCF_FILENAME+'\n\n'
				f.close()
				exit(1)
			splt = line[:-1].split('\t')
			plOut = parseLine(splt,colDict,colSamp)
			if plOut == None:
				nSkipped += 1
			else:
				(aa, af, gt) = plOut

				# make sure at least one allele somewhere contains the variant
				if tumorNormal:
					gtEval = gt[:2]
				else:
					gtEval = gt[:1]
				if None in gtEval:
					if CHOOSE_RANDOM_PLOID_IF_NO_GT_FOUND:
						if not alreadyPrintedWarning:
							print 'Warning: Found variants without a GT field, assuming heterozygous...'
							alreadyPrintedWarning = True
						for i in xrange(len(gtEval)):
							tmp = ['0']*ploidy
							tmp[random.randint(0,ploidy-1)] = '1'
							gtEval[i] = '/'.join(tmp)
					else:
						# skip because no GT field was found
						nSkipped += 1
						continue
				isNonReference = False
				for gtVal in gtEval:
					if gtVal != None:
						if '1' in gtVal:
							isNonReference = True
				if not isNonReference:
					# skip if no genotype actually contains this variant
					nSkipped += 1
					continue

				chrom = splt[0]
				pos   = int(splt[1])
				ref   = splt[3]
				# skip if position is <= 0
				if pos <= 0:
					nSkipped += 1
					continue

				# hash variants to avoid inserting duplicates (there are some messy VCFs out there...)
				if chrom not in allVars:
					allVars[chrom] = {}
				if pos not in allVars[chrom]:
					allVars[chrom][pos] = (pos,ref,aa,af,gtEval)
				else:
					nSkipped_becauseHash += 1
			
		else:
			if line[1] != '#':
				cols = line[1:-1].split('\t')
				for i in xrange(len(cols)):
					if 'FORMAT' in colDict:
						colSamp.append(i)
					colDict[cols[i]] = i
				if len(colSamp):
					sampNames = cols[-len(colSamp):]
					if len(colSamp) == 1:
						pass
					elif len(colSamp) == 2 and tumorNormal:
						print 'Detected 2 sample columns in input VCF, assuming tumor/normal.'
					else:
						print 'Warning: Multiple sample columns present in input VCF. By default genReads uses only the first column.'
				else:
					sampNames = ['Unknown']
				if tumorNormal:
					#tumorInd  = sampNames.index('TUMOR')
					#normalInd = sampNames.index('NORMAL')
					if 'NORMAL' not in sampNames or 'TUMOR' not in sampNames:
						print '\n\nERROR: Input VCF must have a "NORMAL" and "TUMOR" column.\n'
	f.close()

	varsOut = {}
	for r in allVars.keys():
		varsOut[r] = [allVars[r][k] for k in sorted(allVars[r].keys())]
	
	print 'found',sum([len(n) for n in allVars.values()]),'valid variants in input vcf.'
	print ' *',nSkipped,'variants skipped: (qual filtered / ref genotypes / invalid syntax)'
	print ' *',nSkipped_becauseHash,'variants skipped due to multiple variants found per position'
	print '--------------------------------'
	return (sampNames, varsOut)



