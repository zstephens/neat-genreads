#!/usr/bin/env python

# Python 3 ready

import sys
import os
import gzip
from struct import unpack

BAM_EOF = ['1f', '8b', '08', '04', '00', '00', '00', '00', '00', 'ff', '06', '00', '42', '43', '02', '00', '1b', '00',
           '03', '00', '00', '00', '00', '00', '00', '00', '00', '00']


def get_bytes(fmt, amt):
    if fmt == '<i' or fmt == '<I':
        my_size = 4
    elif fmt == '<c' or fmt == '<b' or fmt == '<B':
        my_size = 1
    else:
        print('\nError, unknown format:', fmt, '\n')
        exit(1)
    if amt == 1:
        f_read = f.read(my_size)
        if not f_read:
            return None
        return unpack(fmt, f_read)[0]
    else:
        f_read = f.read(my_size * amt)
        if not f_read:
            return None
        return unpack(fmt, f_read)


# check eof
IN_BAM = sys.argv[1]
f = open(IN_BAM, 'rb')
f.seek(os.path.getsize(IN_BAM) - 28)
EOF = [format(n, '02x') for n in f.read()]
print('EOF_MARKER:  ', ' '.join(EOF))
if EOF != BAM_EOF:
    print('\nWARNING: BAM EOF DOES NOT MATCH EXPECTED STRING.\n')
f.close()

# check other stuff
f = gzip.open(IN_BAM, 'rb')

print('MAGIC STRING:', f.read(4))
l_text = get_bytes('<i', 1)
print('l_text:      ', l_text)
print('text:      \n', f.read(l_text))
n_ref = get_bytes('<i', 1)
print('n_ref:       ', n_ref)

for i in range(n_ref):
    l_name = get_bytes('<i', 1)
    print('ref' + str(i) + ' - l_name:', l_name)
    print('ref' + str(i) + ' - name:  ', f.read(l_name))
    print('ref' + str(i) + ' - l_ref: ', get_bytes('<i', 1))

print('\nEXAMINING ALIGNMENT DATA:\n')
aln_N = 0
while True:
    aln_N += 1
    block_size = get_bytes('<i', 1)
    if block_size == None:
        break
    print('[' + str(aln_N) + ']:', 'block_size:', block_size)
    print('-- refID:', get_bytes('<i', 1))
    print('-- pos:  ', get_bytes('<i', 1))
    bmqnl = get_bytes('<I', 1)
    binv = (bmqnl >> 16) & 65535
    mapq = (bmqnl >> 8) & 255
    lrn = bmqnl & 255
    print('-- bmqnl:', bmqnl, '(bin=' + str(binv) + ', mapq=' + str(mapq) + ', l_readname+1=' + str(lrn) + ')')
    flgnc = get_bytes('<I', 1)
    flag = (flgnc >> 16) & 65535
    ncig = flgnc & 65535
    print('-- flgnc:', flgnc, '(flag=' + str(flag) + ', ncig=' + str(ncig) + ')')
    print('-- l_seq:', get_bytes('<i', 1))
    print('-- nxtID:', get_bytes('<i', 1))
    print('-- nxtPo:', get_bytes('<i', 1))
    print('-- tlen: ', get_bytes('<i', 1))
    print('-- rname:', str([f.read(lrn)])[1:-1])

    f.read(block_size - 32 - lrn)
# print [block]

f.close()
