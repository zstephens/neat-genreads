#!/usr/bin/env python
import os
import argparse

def getListOfFiles(inDir,pattern):
	return [inDir+n for n in os.listdir(inDir) if (pattern in n and os.path.getsize(inDir+n))]

TEMP_IND = 0
def stripVCF_header(fn):
	global TEMP_IND
	f = open(fn,'r')
	ftn = fn+'_temp'+str(TEMP_IND)
	f_t = open(ftn,'w')
	hasHeader = False
	for line in f:
		if line[0] == '#':
			if not hasHeader:
				TEMP_IND += 1
			hasHeader = True
		elif hasHeader:
			f_t.write(line)
		else:
			break
	f_t.close()
	f.close()
	if hasHeader:
		return ftn
	else:
		os.system('rm '+ftn)
		return fn

def catListOfFiles(l,outName,gzipped=False):
	for n in l:
		if n[-3:] == '.gz' or n[-5:] == '.gzip':
			gzipped = True
	if gzipped:
		for n in l:
			if not n[-3:] == '.gz' and not n[-5:] == '.gzip':
				print '\nError: Found a mixture of compressed and decompressed files with the specified prefix. Abandoning ship...\n'
				for m in l:
					print m
				print ''
				exit(1)
		cmd = 'cat '+' '.join(sorted(l))+' > '+outName+'.gz'
	else:
		cmd = 'cat '+' '.join(sorted(l))+' > '+outName
	print cmd
	os.system(cmd)

def catBams(l,outName,samtools_exe):
	l_sort = sorted(l)
	tmp = outName+'.tempHeader.sam'
	os.system(samtools_exe+' view -H '+l_sort[0]+' > '+tmp)
	cmd = samtools_exe+' cat -h '+tmp+' '+' '.join(l_sort)+' > '+outName
	print cmd
	os.system(cmd)
	os.system('rm '+tmp)


#####################################
#				main()				#
#####################################

def main():

	parser = argparse.ArgumentParser(description='mergeJobs.py')
	parser.add_argument('-i', type=str, required=True, metavar='<str>', nargs='+', help="* input prefix: [prefix_1] [prefix_2] ...")
	parser.add_argument('-o', type=str, required=True, metavar='<str>',            help="* output prefix")
	parser.add_argument('-s', type=str, required=True, metavar='<str>',            help="* /path/to/samtools")

	args = parser.parse_args()
	(INP,OUP,SAMTOOLS) = (args.i,args.o,args.s)

	inDir = '/'.join(INP[0].split('/')[:-1])+'/'
	if inDir == '/':
		inDir = './'
	#print inDir

	INP_LIST = []
	for n in INP:
		if n[-1] == '/':
			n = n[:-1]
		INP_LIST.append(n.split('/')[-1])
	listing_r1 = []
	listing_r2 = []
	listing_b  = []
	listing_v  = []
	for n in INP_LIST:
		listing_r1 += getListOfFiles(inDir,n+'_read1.fq.job')
		listing_r2 += getListOfFiles(inDir,n+'_read2.fq.job')
		listing_b  += getListOfFiles(inDir,n+'_golden.bam.job')
		if len(listing_v):	# remove headers from vcf files that aren't the first being processed
			initList   = getListOfFiles(inDir,n+'_golden.vcf.job')
			listing_v += [stripVCF_header(n) for n in initList]
		else:
			listing_v  += getListOfFiles(inDir,n+'_golden.vcf.job')
	
	#
	#	merge fq files
	#
	if len(listing_r1):
		catListOfFiles(listing_r1,OUP+'_read1.fq')
	if len(listing_r2):
		catListOfFiles(listing_r2,OUP+'_read2.fq')

	#
	#	merge golden alignments, if present
	#
	if len(listing_b):
		catBams(listing_b,OUP+'_golden.bam',SAMTOOLS)

	#
	#	merge golden vcfs, if present
	#
	if len(listing_v):
		catListOfFiles(listing_v,OUP+'_golden.vcf')


if __name__ == "__main__":
	main()

