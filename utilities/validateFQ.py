#!/usr/bin/env python

#
#	A quickie tool for validating the correctness of a FASTQ file
#
#	python validateFQ.py read1.fq [read2.fq]
#

import sys

def get4lines(fn):
	l1 = fn.readline().strip()
	l2 = fn.readline().strip()
	l3 = fn.readline().strip()
	l4 = fn.readline().strip()
	if any([l1,l2,l3,l4]) and not all([l1,l2,l3,l4]):
		print '\nError: missing lines:\n'
		print l1+'\n'+l2+'\n'+l3+'\n'+l4+'\n'
		exit(1)
	return (l1,l2,l3,l4)

ALLOWED_QUAL = '!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ'
ALLOWED_NUCL = 'ACGTN'

def validate4lines(l1,l2,l3,l4):
	failed = 0
	# make sure lines contain correct delimiters
	if l1[0] != '@' or l1[-2] != '/' or l3[0] != '+':
		failed = 1
	# make sure seq len == qual length
	if len(l2) != len(l4):
		failed = 2
	# make sure seq string contains only valid characters
	for n in l2:
		if n not in ALLOWED_NUCL:
			failed = 3
	# make sure qual string contains only valid characters
	for n in l4:
		if n not in ALLOWED_QUAL:
			failed = 4
	if failed:
		print '\nError: malformed lines:'
		if failed == 1: print ' ---- invalid delimiters\n'
		elif failed == 2: print ' ---- seq len != qual len\n'
		elif failed == 3: print ' ---- seq contains invalid characters\n'
		elif failed == 4: print ' ---- qual contains invalid characters\n'
		print l1+'\n'+l2+'\n'+l3+'\n'+l4+'\n'
		exit(1)

f1 = open(sys.argv[1],'r')
(l1_r1, l2_r1, l3_r1, l4_r1) = get4lines(f1)
f2 = None
if len(sys.argv) == 3:
	f2 = open(sys.argv[2],'r')
	(l1_r2, l2_r2, l3_r2, l4_r2) = get4lines(f2)

while l1_r1:
	# check line syntax
	validate4lines(l1_r1,l2_r1,l3_r1,l4_r1)
	if f2 != None:
		validate4lines(l1_r2,l2_r2,l3_r2,l4_r2)
		# make sure seq id is same for r1/r2
		if l1_r1[:-1] != l1_r2[:-1]:
			print '\nError: mismatched r1/r2 name:\n'
			print l1_r1+'\n'+l1_r2+'\n'
			exit(1)

	# grab next 4 lines...
	(l1_r1, l2_r1, l3_r1, l4_r1) = get4lines(f1)
	if f2 != None:
		(l1_r2, l2_r2, l3_r2, l4_r2) = get4lines(f2)

if f2 != None:
	f2.close()
f1.close()

print '\nPASSED WITH FLYING COLORS. GOOD DAY.\n'

