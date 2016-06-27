import os
import argparse


parser = argparse.ArgumentParser(description='mergeJobs.py')
parser.add_argument('-i', type=str, required=True, metavar='<str>', help="* input prefix")
parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output prefix")
parser.add_argument('-s', type=str, required=True, metavar='<str>', help="* /path/to/samtools")
args = parser.parse_args()
(INP,OUP,SAMTOOLS) = (args.i,args.o,args.s)


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

	inDir = '/'.join(INP.split('/')[:-1])+'/'
	if inDir == '/':
		inDir = './'
	print inDir
	
	#
	#	merge fq files
	#
	listing_r1 = [inDir+n for n in os.listdir(inDir) if ('_read1.fq.job' in n and os.path.getsize(inDir+n))]
	if len(listing_r1):
		catListOfFiles(listing_r1,OUP+'_read1.fq')
	listing_r2 = [inDir+n for n in os.listdir(inDir) if ('_read2.fq.job' in n and os.path.getsize(inDir+n))]
	if len(listing_r2):
		catListOfFiles(listing_r2,OUP+'_read2.fq')

	#
	#	merge golden alignments, if present
	#
	listing = [inDir+n for n in os.listdir(inDir) if ('_golden.bam.job' in n and os.path.getsize(inDir+n))]
	if len(listing):
		catBams(listing,OUP+'_golden.bam',SAMTOOLS)

	#
	#	merge golden vcfs, if present
	#
	listing = [inDir+n for n in os.listdir(inDir) if ('_golden.vcf.job' in n and os.path.getsize(inDir+n))]
	if len(listing):
		catListOfFiles(listing,OUP+'_golden.vcf')


if __name__ == "__main__":
	main()


