import sys
import os
import re
import gzip
from struct import pack

from biopython_modified_bgzf import BgzfWriter

BAM_COMPRESSION_LEVEL = 6

# return the reverse complement of a string
RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join(RC_DICT[n] for n in s[::-1])

# SAMtools reg2bin function
def reg2bin(a,b):
	b -= 1
	if (a>>14 == b>>14): return ((1<<15)-1)/7 + (a>>14)
	if (a>>17 == b>>17): return ((1<<12)-1)/7 + (a>>17)
	if (a>>20 == b>>20): return  ((1<<9)-1)/7 + (a>>20)
	if (a>>23 == b>>23): return  ((1<<6)-1)/7 + (a>>23)
	if (a>>26 == b>>26): return  ((1<<3)-1)/7 + (a>>26)
	return 0

CIGAR_PACKED = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8}
SEQ_PACKED   = {'=':0, 'A':1, 'C':2, 'M':3, 'G':4, 'R':5, 'S':6, 'V':7,
                'T':8, 'W':9, 'Y':10,'H':11,'K':12,'D':13,'B':14,'N':15}

BUFFER_BATCH_SIZE = 1000		# write out to file after this many reads

#
#	outFQ      = path to output FASTQ prefix
#	paired     = True for PE reads, False for SE
#	BAM_header = [refIndex]
#	VCF_header = [path_to_ref]
#	gzipped    = True for compressed FASTQ/VCF, False for uncompressed
#
class OutputFileWriter:
	def __init__(self, outPrefix, paired=False, BAM_header=None, VCF_header=None, gzipped=False, jobTuple=(1,1), noFASTQ=False):
		
		jobSuffix = ''
		if jobTuple[1] > 1:
			jsl = len(str(jobTuple[1]))
			jsb = '0'*(jsl-len(str(jobTuple[0])))
			jobSuffix = '.job'+jsb+str(jobTuple[0])+'of'+str(jobTuple[1])

		fq1 = outPrefix+'_read1.fq'+jobSuffix
		fq2 = outPrefix+'_read2.fq'+jobSuffix
		bam = outPrefix+'_golden.bam'+jobSuffix
		vcf = outPrefix+'_golden.vcf'+jobSuffix

		self.noFASTQ = noFASTQ
		if not self.noFASTQ:
			if gzipped:
				self.fq1_file = gzip.open(fq1+'.gz', 'wb')
			else:
				self.fq1_file = open(fq1,'w')

			self.fq2_file = None
			if paired:
				if gzipped:
					self.fq2_file = gzip.open(fq2+'.gz', 'wb')
				else:
					self.fq2_file = open(fq2,'w')

		#
		#	VCF OUTPUT
		#
		self.vcf_file = None
		if VCF_header != None:
			if gzipped:
				self.vcf_file = gzip.open(vcf+'.gz', 'wb')
			else:
				self.vcf_file = open(vcf, 'wb')

			# WRITE VCF HEADER (if parallel: only for first job)
			if jobTuple[0] == 1:
				self.vcf_file.write('##fileformat=VCFv4.1\n')
				self.vcf_file.write('##reference='+VCF_header[0]+'\n')
				self.vcf_file.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
				self.vcf_file.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
				#self.vcf_file.write('##INFO=<ID=READS,Number=1,Type=String,Description="Names of Reads Covering this Variant">\n')
				self.vcf_file.write('##INFO=<ID=VMX,Number=1,Type=String,Description="SNP is Missense in these Read Frames">\n')
				self.vcf_file.write('##INFO=<ID=VNX,Number=1,Type=String,Description="SNP is Nonsense in these Read Frames">\n')
				self.vcf_file.write('##INFO=<ID=VFX,Number=1,Type=String,Description="Indel Causes Frameshift">\n')
				self.vcf_file.write('##INFO=<ID=WP,Number=A,Type=Integer,Description="NEAT-GenReads ploidy indicator">\n')
				self.vcf_file.write('##ALT=<ID=DEL,Description="Deletion">\n')
				self.vcf_file.write('##ALT=<ID=DUP,Description="Duplication">\n')
				self.vcf_file.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n')
				self.vcf_file.write('##ALT=<ID=INV,Description="Inversion">\n')
				self.vcf_file.write('##ALT=<ID=CNV,Description="Copy number variable region">\n')
				self.vcf_file.write('##ALT=<ID=TRANS,Description="Translocation">\n')
				self.vcf_file.write('##ALT=<ID=INV-TRANS,Description="Inverted translocation">\n')
				self.vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

		#
		#	BAM OUTPUT
		#
		self.bam_file = None
		if BAM_header != None:
			self.bam_file = BgzfWriter(bam, 'w', compresslevel=BAM_COMPRESSION_LEVEL)

			# WRITE BAM HEADER (if parallel: only for first job)
			if True or jobTuple[0] == 1:
				self.bam_file.write("BAM\1")
				header = '@HD\tVN:1.5\tSO:coordinate\n'
				for n in BAM_header[0]:
					header += '@SQ\tSN:'+n[0]+'\tLN:'+str(n[3])+'\n'
				header += '@RG\tID:NEAT\n'
				headerBytes = len(header)
				numRefs     = len(BAM_header[0])
				self.bam_file.write(pack('<i',headerBytes))
				self.bam_file.write(header)
				self.bam_file.write(pack('<i',numRefs))

				for n in BAM_header[0]:
					l_name = len(n[0])+1
					self.bam_file.write(pack('<i',l_name))
					self.bam_file.write(n[0]+'\0')
					self.bam_file.write(pack('<i',n[3]))

		# buffers for more efficient writing
		self.fq1_buffer = []
		self.fq2_buffer = []
		self.bam_buffer = []

	def writeFASTQRecord(self,readName,read1,qual1,read2=None,qual2=None):
		###self.fq1_file.write('@'+readName+'/1\n'+read1+'\n+\n'+qual1+'\n')
		self.fq1_buffer.append('@'+readName+'/1\n'+read1+'\n+\n'+qual1+'\n')
		if read2 != None:
			####self.fq2_file.write('@'+readName+'/2\n'+RC(read2)+'\n+\n'+qual2[::-1]+'\n')
			self.fq2_buffer.append('@'+readName+'/2\n'+RC(read2)+'\n+\n'+qual2[::-1]+'\n')

	def writeVCFRecord(self, chrom, pos, idStr, ref, alt, qual, filt, info):
		self.vcf_file.write(str(chrom)+'\t'+str(pos)+'\t'+str(idStr)+'\t'+str(ref)+'\t'+str(alt)+'\t'+str(qual)+'\t'+str(filt)+'\t'+str(info)+'\n')

	def writeBAMRecord(self, refID, readName, pos_0, cigar, seq, qual, samFlag, matePos=None, alnMapQual=70):

		myBin     = reg2bin(pos_0,pos_0+len(seq))
		#myBin     = 0	# or just use a dummy value, does this actually matter?

		myMapQual = alnMapQual
		cig_letters = re.split(r"\d+",cigar)[1:]
		cig_numbers = [int(n) for n in re.findall(r"\d+",cigar)]
		cig_ops     = len(cig_letters)
		next_refID = refID
		if matePos == None:
			next_pos = 0
			my_tlen  = 0
		else:
			next_pos = matePos
			if pos_0 < next_pos:
				my_tlen = next_pos + len(seq) - pos_0
			else:
				my_tlen = -pos_0 - len(seq) + next_pos

		encodedCig = ''
		for i in xrange(cig_ops):
			encodedCig += pack('<I',(cig_numbers[i]<<4) + CIGAR_PACKED[cig_letters[i]])
		encodedSeq = ''
		encodedLen = (len(seq)+1)/2
		seqLen     = len(seq)
		if seqLen&1:
			seq += '='
		for i in xrange(encodedLen):
			encodedSeq += pack('<B',(SEQ_PACKED[seq[2*i]]<<4) + SEQ_PACKED[seq[2*i+1]])

		# apparently samtools automatically adds 33 to the quality score string...
		encodedQual = ''.join([chr(ord(n)-33) for n in qual])

		#blockSize = 4 +		# refID 		int32
		#            4 +		# pos			int32
		#            4 +		# bin_mq_nl		uint32
		#            4 +		# flag_nc		uint32
		#            4 +		# l_seq			int32
		#            4 +		# next_refID	int32
		#            4 +		# next_pos		int32
		#            4 +		# tlen			int32
		#            len(readName)+1 +
		#            4*cig_ops +
		#            encodedLen +
		#            len(seq)

		#blockSize = 32 + len(readName)+1 + 4*cig_ops + encodedLen + len(seq)
		blockSize = 32 + len(readName)+1 + len(encodedCig) + len(encodedSeq) + len(encodedQual)

		####self.bam_file.write(pack('<i',blockSize))
		####self.bam_file.write(pack('<i',refID))
		####self.bam_file.write(pack('<i',pos_0))
		####self.bam_file.write(pack('<I',(myBin<<16) + (myMapQual<<8) + len(readName)+1))
		####self.bam_file.write(pack('<I',(samFlag<<16) + cig_ops))
		####self.bam_file.write(pack('<i',seqLen))
		####self.bam_file.write(pack('<i',next_refID))
		####self.bam_file.write(pack('<i',next_pos))
		####self.bam_file.write(pack('<i',my_tlen))
		####self.bam_file.write(readName+'\0')
		####self.bam_file.write(encodedCig)
		####self.bam_file.write(encodedSeq)
		####self.bam_file.write(encodedQual)

		# a horribly compressed line, I'm sorry.
		# (ref_index, position, data)
		self.bam_buffer.append((refID, pos_0, pack('<i',blockSize) + pack('<i',refID) + pack('<i',pos_0) + pack('<I',(myBin<<16) + (myMapQual<<8) + len(readName)+1) + pack('<I',(samFlag<<16) + cig_ops) + pack('<i',seqLen) + pack('<i',next_refID) + pack('<i',next_pos) + pack('<i',my_tlen) + readName+'\0' + encodedCig + encodedSeq + encodedQual))


	def flushBuffers(self,bamMax=None,lastTime=False):
		if (len(self.fq1_buffer) >= BUFFER_BATCH_SIZE or len(self.bam_buffer) >= BUFFER_BATCH_SIZE) or (len(self.fq1_buffer) and lastTime) or (len(self.bam_buffer) and lastTime):
			# fq
			if not self.noFASTQ:
				self.fq1_file.write(''.join(self.fq1_buffer))
				if len(self.fq2_buffer):
					self.fq2_file.write(''.join(self.fq2_buffer))
			# bam
			if len(self.bam_buffer):
				bam_data = sorted(self.bam_buffer)
				if lastTime:
					self.bam_file.write(''.join([n[2] for n in bam_data]))
					self.bam_buffer = []
				else:
					ind_to_stop_at = 0
					for i in xrange(0,len(bam_data)):
						# if we are from previous reference, or have coordinates lower than next window position, it's safe to write out to file
						if bam_data[i][0] != bam_data[-1][0] or bam_data[i][1] < bamMax:
							ind_to_stop_at = i+1
						else:
							break
					self.bam_file.write(''.join([n[2] for n in bam_data[:ind_to_stop_at]]))
					####print 'BAM WRITING:',ind_to_stop_at,'/',len(bam_data)
					if ind_to_stop_at >= len(bam_data):
						self.bam_buffer = []
					else:
						self.bam_buffer = bam_data[ind_to_stop_at:]
			self.fq1_buffer = []
			self.fq2_buffer = []


	def closeFiles(self):
		self.flushBuffers(lastTime=True)
		if not self.noFASTQ:
			self.fq1_file.close()
			if self.fq2_file != None:
				self.fq2_file.close()
		if self.vcf_file != None:
			self.vcf_file.close()
		if self.bam_file != None:
			self.bam_file.close()




