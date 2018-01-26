#!/usr/bin/env python

import sys
import os
import gzip
from struct import unpack

BAM_EOF = ['1f', '8b', '08', '04', '00', '00', '00', '00', '00', 'ff', '06', '00', '42', '43', '02', '00', '1b', '00', '03', '00', '00', '00', '00', '00', '00', '00', '00', '00']

def getBytes(fmt,amt):
	if fmt == '<i' or fmt == '<I':
		mySize = 4
	elif fmt == '<c' or fmt == '<b' or fmt == '<B':
		mySize = 1
	else:
		print '\nError, unknown format:',fmt,'\n'
		exit(1)
	if amt == 1:
		fread = f.read(mySize)
		if not fread:
			return None
		return unpack(fmt,fread)[0]
	else:
		fread = f.read(mySize*amt)
		if not fread:
			return None
		return unpack(fmt,fread)

# check eof
IN_BAM = sys.argv[1]
f = open(IN_BAM,'rb')
f.seek(os.path.getsize(IN_BAM)-28)
EOF = [format(ord(n),'02x') for n in f.read()]
print 'EOF_MARKER:  ', ' '.join(EOF)
if EOF != BAM_EOF:
	print '\nWARNING: BAM EOF DOES NOT MATCH EXPECTED STRING.\n'
f.close()

# check other stuff
f = gzip.open(IN_BAM,'rb')

print 'MAGIC STRING:', f.read(4)
l_text = getBytes('<i',1)
print 'l_text:      ', l_text
print 'text:      \n', f.read(l_text)
n_ref = getBytes('<i',1)
print 'n_ref:       ', n_ref

for i in xrange(n_ref):
	l_name = getBytes('<i',1)
	print 'ref'+str(i)+' - l_name:', l_name
	print 'ref'+str(i)+' - name:  ', f.read(l_name)
	print 'ref'+str(i)+' - l_ref: ', getBytes('<i',1)

print '\nEXAMINING ALIGNMENT DATA:\n'
aln_N = 0
while True:
	aln_N += 1
	blockSize = getBytes('<i',1)
	if blockSize == None:
		break
	print '['+str(aln_N)+']:', 'blockSize:', blockSize
	print '-- refID:', getBytes('<i',1)
	print '-- pos:  ', getBytes('<i',1)
	bmqnl = getBytes('<I',1)
	binv  = (bmqnl>>16)&65535
	mapq  = (bmqnl>>8)&255
	lrn   = bmqnl&255
	print '-- bmqnl:', bmqnl, '(bin='+str(binv)+', mapq='+str(mapq)+', l_readname+1='+str(lrn)+')'
	flgnc = getBytes('<I',1)
	flag  = (flgnc>>16)&65535
	ncig  = flgnc&65535
	print '-- flgnc:', flgnc, '(flag='+str(flag)+', ncig='+str(ncig)+')'
	print '-- l_seq:', getBytes('<i',1)
	print '-- nxtID:', getBytes('<i',1)
	print '-- nxtPo:', getBytes('<i',1)
	print '-- tlen: ', getBytes('<i',1)
	print '-- rname:', str([f.read(lrn)])[1:-1]

	f.read(blockSize-32-lrn)
	#print [block]


f.close()
