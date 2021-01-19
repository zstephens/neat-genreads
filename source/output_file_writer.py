import sys
import os
import re
import gzip
from struct import pack
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
from Bio.bgzf import *
import pathlib
from source.neat_cigar_rework import CigarString

# from source.biopython_modified_bgzf import BgzfWriter

BAM_COMPRESSION_LEVEL = 6


# TODO figure out why these functions are in this file in the first place
def reverse_complement(dna_string: str) -> str:
    """
    Return the reverse complement of a string from a DNA strand
    :param dna_string: string of DNA
    :return: the reverse compliment of the above string
    """
    rc_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(rc_dict[n] for n in dna_string[::-1])


# SAMtools reg2bin function
def reg2bin(beg: int, end: int):
    """
    Finds the largest superset bin of region. Numeric values taken from hts-specs
    Note: description of this function taken from source code for bamnostic.bai
        (https://bamnostic.readthedocs.io/en/latest/_modules/bamnostic/bai.html)
    :param beg: inclusive beginning position of region
    :param end: exclusive end position of region
    :return: distinct bin ID or largest superset bin of region
    """
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


# takes list of strings, returns numerical flag
def sam_flag(string_list: list) -> int:
    out_val = 0
    string_list = list(set(string_list))
    for n in string_list:
        if n == 'paired':
            out_val += 1
        elif n == 'proper':
            out_val += 2
        elif n == 'unmapped':
            out_val += 4
        elif n == 'mate_unmapped':
            out_val += 8
        elif n == 'reverse':
            out_val += 16
        elif n == 'mate_reverse':
            out_val += 32
        elif n == 'first':
            out_val += 64
        elif n == 'second':
            out_val += 128
        elif n == 'not_primary':
            out_val += 256
        elif n == 'low_quality':
            out_val += 512
        elif n == 'duplicate':
            out_val += 1024
        elif n == 'supplementary':
            out_val += 2048
    return out_val


CIGAR_PACKED = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
SEQ_PACKED = {'=': 0, 'A': 1, 'C': 2, 'M': 3, 'G': 4, 'R': 5, 'S': 6, 'V': 7,
              'T': 8, 'W': 9, 'Y': 10, 'H': 11, 'K': 12, 'D': 13, 'B': 14, 'N': 15}

# TODO figure out an optimum batch size
BUFFER_BATCH_SIZE = 8000  # write out to file after this many reads


# TODO find a better way to write output files
class OutputFileWriter:
    def __init__(self, out_prefix, paired=False, bam_header=None, vcf_header=None, gzipped=False,
                 no_fastq=False, fasta_instead=False):

        self.fasta_instead = fasta_instead
        # TODO Eliminate paired end as an option for fastas
        if self.fasta_instead:
            fq1 = pathlib.Path(out_prefix + '.fasta.gz')
            fq2 = None
        else:
            fq1 = pathlib.Path(out_prefix + '_read1.fq.gz')
            fq2 = pathlib.Path(out_prefix + '_read2.fq.gz')
        bam = pathlib.Path(out_prefix + '_golden.bam')
        vcf = pathlib.Path(out_prefix + '_golden.vcf.gz')

        # TODO Make a fasta-specific method
        self.no_fastq = no_fastq
        if not self.no_fastq:
            if gzipped:
                self.fq1_file = gzip.open(fq1.with_suffix(fq1.suffix + '.gz'), 'wb')
            else:
                self.fq1_file = open(fq1, 'w')

            self.fq2_file = None
            if paired:
                if gzipped:
                    self.fq2_file = gzip.open(fq2.with_suffix(fq2.suffix + '.gz'), 'wb')
                else:
                    self.fq2_file = open(fq2, 'w')

        # VCF OUTPUT
        self.vcf_file = None
        if vcf_header is not None:
            if gzipped:
                self.vcf_file = gzip.open(vcf.with_suffix(vcf.suffix + '.gz'), 'wb')
            else:
                self.vcf_file = open(vcf, 'wb')

            # WRITE VCF HEADER
            self.vcf_file.write('##fileformat=VCFv4.1\n'.encode('utf-8'))
            reference = '##reference=' + vcf_header[0] + '\n'
            self.vcf_file.write(reference.encode('utf-8'))
            self.vcf_file.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'.encode('utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'.encode('utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=VMX,Number=1,Type=String,Description="SNP is Missense in these Read Frames">\n'.encode(
                    'utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=VNX,Number=1,Type=String,Description="SNP is Nonsense in these Read Frames">\n'.encode(
                    'utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=VFX,Number=1,Type=String,Description="Indel Causes Frameshift">\n'.encode('utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=WP,Number=A,Type=Integer,Description="NEAT-GenReads ploidy indicator">\n'.encode(
                    'utf-8'))
            self.vcf_file.write('##ALT=<ID=DEL,Description="Deletion">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=DUP,Description="Duplication">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=INV,Description="Inversion">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=CNV,Description="Copy number variable region">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=TRANS,Description="Translocation">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=INV-TRANS,Description="Inverted translocation">\n'.encode('utf-8'))
            # TODO add sample to vcf output
            self.vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'.encode('utf-8'))

        # BAM OUTPUT
        self.bam_file = None
        if bam_header is not None:
            self.bam_file = BgzfWriter(bam, 'w', compresslevel=BAM_COMPRESSION_LEVEL)

            # WRITE BAM HEADER
            self.bam_file.write("BAM\1")
            header = '@HD\tVN:1.5\tSO:coordinate\n'
            for n in bam_header[0]:
                header += '@SQ\tSN:' + n[0] + '\tLN:' + str(n[3]) + '\n'
            header += '@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n'
            header_bytes = len(header)
            num_refs = len(bam_header[0])
            self.bam_file.write(pack('<i', header_bytes))
            self.bam_file.write(header)
            self.bam_file.write(pack('<i', num_refs))

            for n in bam_header[0]:
                l_name = len(n[0]) + 1
                self.bam_file.write(pack('<i', l_name))
                self.bam_file.write(n[0] + '\0')
                self.bam_file.write(pack('<i', n[3]))

        # buffers for more efficient writing
        self.fq1_buffer = []
        self.fq2_buffer = []
        self.bam_buffer = []

    # TODO add write_fasta_record

    def write_fastq_record(self, read_name, read1, qual1, read2=None, qual2=None, orientation=None):
        (read1, quality1) = (str(read1), qual1)
        if read2 is not None and orientation is True:
            (read2, quality2) = (str(reverse_complement(read2)), qual2[::-1])
        elif read2 is not None and orientation is False:
            (read1, quality1) = (str(reverse_complement(read2)), qual2[::-1])
            (read2, quality2) = (str(read1), qual1)

        if self.fasta_instead:
            self.fq1_buffer.append('>' + read_name + '/1\n' + read1 + '\n')
            if read2 is not None:
                self.fq2_buffer.append('>' + read_name + '/2\n' + read2 + '\n')
        else:
            self.fq1_buffer.append('@' + read_name + '/1\n' + read1 + '\n+\n' + quality1 + '\n')
            if read2 is not None:
                self.fq2_buffer.append('@' + read_name + '/2\n' + read2 + '\n+\n' + quality2 + '\n')

    def write_vcf_record(self, chrom, pos, id_str, ref, alt, qual, filt, info):
        self.vcf_file.write(
            str(chrom) + '\t' + str(pos) + '\t' + str(id_str) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(
                qual) + '\t' + str(filt) + '\t' + str(info) + '\n')

    def write_bam_record(self, ref_id, read_name, pos_0, cigar, seq, qual, sam_flag, mate_pos=None, aln_map_quality=70):

        my_bin = reg2bin(pos_0, pos_0 + len(seq))
        # my_bin     = 0	# or just use a dummy value, does this actually matter?

        my_map_quality = aln_map_quality
        cig_letters = []
        cig_numbers = []
        for item in cigar.items():
            cig_numbers.append(item[0])
            cig_letters.append(item[1])
        # TODO delete the following two lines once the new cigar rework is 100% working
        # cig_letters = re.split(r"\d+", cigar)[1:]
        # cig_numbers = [int(n) for n in re.findall(r"\d+", cigar)]
        cig_ops = len(cig_letters)
        next_ref_id = ref_id
        if mate_pos is None:
            next_pos = 0
            my_t_len = 0
        else:
            next_pos = mate_pos
            if pos_0 < next_pos:
                my_t_len = next_pos + len(seq) - pos_0
            else:
                my_t_len = -pos_0 - len(seq) + next_pos

        encoded_cig = bytearray()
        for i in range(cig_ops):
            encoded_cig.extend(pack('<I', (cig_numbers[i] << 4) + CIGAR_PACKED[cig_letters[i]]))
        encoded_seq = bytearray()
        encoded_len = (len(seq) + 1) // 2
        seq_len = len(seq)
        if seq_len & 1:
            seq += '='
        for i in range(encoded_len):
            # print(seq[2*i], seq[2*i+1])
            encoded_seq.extend(
                pack('<B', (SEQ_PACKED[seq[2 * i].capitalize()] << 4) + SEQ_PACKED[seq[2 * i + 1].capitalize()]))

        # apparently samtools automatically adds 33 to the quality score string...
        encoded_qual = ''.join([chr(ord(n) - 33) for n in qual])

        # block_size = 4 +		# refID 		int32
        #            4 +		# pos			int32
        #            4 +		# bin_mq_nl		uint32
        #            4 +		# flag_nc		uint32
        #            4 +		# l_seq			int32
        #            4 +		# next_ref_id	int32
        #            4 +		# next_pos		int32
        #            4 +		# tlen			int32
        #            len(readName)+1 +
        #            4*cig_ops +
        #            encoded_len +
        #            len(seq)

        # block_size = 32 + len(readName)+1 + 4*cig_ops + encoded_len + len(seq)
        block_size = 32 + len(read_name) + 1 + len(encoded_cig) + len(encoded_seq) + len(encoded_qual)

        ####self.bam_file.write(pack('<i',block_size))
        ####self.bam_file.write(pack('<i',refID))
        ####self.bam_file.write(pack('<i',pos_0))
        ####self.bam_file.write(pack('<I',(my_bin<<16) + (my_map_quality<<8) + len(readName)+1))
        ####self.bam_file.write(pack('<I',(samFlag<<16) + cig_ops))
        ####self.bam_file.write(pack('<i',seq_len))
        ####self.bam_file.write(pack('<i',next_ref_id))
        ####self.bam_file.write(pack('<i',next_pos))
        ####self.bam_file.write(pack('<i',my_tlen))
        ####self.bam_file.write(readName+'\0')
        ####self.bam_file.write(encoded_cig)
        ####self.bam_file.write(encoded_seq)
        ####self.bam_file.write(encoded_qual)

        # a horribly compressed line, I'm sorry.
        # (ref_index, position, data)
        self.bam_buffer.append((ref_id, pos_0, pack('<i', block_size) + pack('<i', ref_id) + pack('<i', pos_0) +
                                pack('<I', (my_bin << 16) + (my_map_quality << 8) + len(read_name) + 1) +
                                pack('<I', (sam_flag << 16) + cig_ops) + pack('<i', seq_len) + pack('<i', next_ref_id) +
                                pack('<i', next_pos) + pack('<i', my_t_len) + read_name.encode('utf-8') +
                                b'\0' + encoded_cig + encoded_seq + encoded_qual.encode('utf-8')))

    def flush_buffers(self, bam_max=None, last_time=False):
        if (len(self.fq1_buffer) >= BUFFER_BATCH_SIZE or len(self.bam_buffer) >= BUFFER_BATCH_SIZE) or (
                len(self.fq1_buffer) and last_time) or (len(self.bam_buffer) and last_time):
            # fq
            if not self.no_fastq:
                self.fq1_file.write(''.join(self.fq1_buffer))
                if len(self.fq2_buffer):
                    self.fq2_file.write(''.join(self.fq2_buffer))
            # bam
            if len(self.bam_buffer):
                bam_data = sorted(self.bam_buffer)
                if last_time:
                    self.bam_file.write(b''.join([n[2] for n in bam_data]))
                    self.bam_buffer = []
                else:
                    ind_to_stop_at = 0
                    for i in range(0, len(bam_data)):
                        # if we are from previous reference, or have coordinates lower
                        # than next window position, it's safe to write out to file
                        if bam_data[i][0] != bam_data[-1][0] or bam_data[i][1] < bam_max:
                            ind_to_stop_at = i + 1
                        else:
                            break
                    self.bam_file.write(b''.join([n[2] for n in bam_data[:ind_to_stop_at]]))
                    ####print 'BAM WRITING:',ind_to_stop_at,'/',len(bam_data)
                    if ind_to_stop_at >= len(bam_data):
                        self.bam_buffer = []
                    else:
                        self.bam_buffer = bam_data[ind_to_stop_at:]
            self.fq1_buffer = []
            self.fq2_buffer = []

    def close_files(self):
        self.flush_buffers(last_time=True)
        if not self.no_fastq:
            self.fq1_file.close()
            if self.fq2_file is not None:
                self.fq2_file.close()
        if self.vcf_file is not None:
            self.vcf_file.close()
        if self.bam_file is not None:
            self.bam_file.close()
