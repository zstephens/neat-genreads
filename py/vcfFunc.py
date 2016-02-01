import sys
import time
import os
import re

INCLUDE_HOMS = False
INCLUDE_FAIL = False

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

	#	if no sample columns, check info for GT
	gt_perSamp = None
	if len(colSamp) == 0 and 'INFO' in colDict and 'GT=' in splt[colDict['INFO']]:
		info = splt[colDict['INFO']]+';'
		gt_perSamp = [re.findall(r"GT=.*?(?=;)",info)[0][3:]]
	elif len(colSamp):
		fmt = ':'+splt[colDict['FORMAT']]+':'
		if ':GT:' in fmt:
			gtInd = fmt.split(':').index('GT')
			gt_perSamp = [splt[colSamp[iii]].split(':')[gtInd-1] for iii in xrange(len(colSamp))]
	if gt_perSamp == None:
		gt_perSamp = [None]*max([len(colSamp),1])

	return (alt_alleles, alt_freqs, gt_perSamp)



def parseVCF(vcfPath,tumorNormal=False):

	tt = time.time()
	sys.stdout.write('reading input VCF... ')
	sys.stdout.flush()
	
	colDict   = {}
	colSamp   = []
	nSkipped  = 0
	allVars   = {}	# [ref][pos]
	sampNames = []
	for line in open(vcfPath,'r'):

		if line[0] != '#':
			if len(colDict) == 0:
				print '\n\nERROR: VCF has no header?\n'+VCF_FILENAME+'\n\n'
				exit(1)
			splt = line[:-1].split('\t')
			plOut = parseLine(splt,colDict,colSamp)
			if plOut == None:
				nSkipped += 1
			else:
				(aa, af, gt) = plOut
				chrom = splt[0]
				pos   = int(splt[1])
				ref   = splt[3]
				if chrom not in allVars:
					allVars[chrom] = {}
				if pos not in allVars[chrom]:
					allVars[chrom][pos] = (pos,ref,aa,af,gt)
			
		else:
			if line[1] != '#':
				cols = line[1:-1].split('\t')
				for i in xrange(len(cols)):
					if 'FORMAT' in colDict:
						colSamp.append(i)
					colDict[cols[i]] = i
				if len(colSamp):
					sampNames = cols[-len(colSamp):]
				else:
					sampNames = ['Unknown']
				if tumorNormal:
					#tumorInd  = sampNames.index('TUMOR')
					#normalInd = sampNames.index('NORMAL')
					if 'NORMAL' not in sampNames or 'TUMOR' not in sampNames:
						print '\n\nERROR: Input VCF must have a "NORMAL" and "TUMOR" column.\n'


	varsOut = {}
	for r in allVars.keys():
		varsOut[r] = [allVars[r][k] for k in sorted(allVars[r].keys())]
	
	print '{0:.3f} (sec)'.format(time.time()-tt)
	print nSkipped,'variants skipped (quality filtered or invalid syntax).'
	return (sampNames, varsOut)



