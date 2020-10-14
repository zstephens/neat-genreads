import sys
import time
import os
import random
from Bio.Seq import Seq

#	Index reference fasta
def index_ref(reference_path):
    tt = time.time()

    filename = None
    if os.path.isfile(reference_path + 'i'):
        print('found index ' + reference_path + 'i')
        filename = reference_path + 'i'
    if os.path.isfile(reference_path + '.fai'):
        print('found index ' + reference_path + '.fai')
        filename = reference_path + '.fai'

    ref_inds = []
    if filename is not None:
        fai = open(filename, 'r')
        for line in fai:
            splt = line[:-1].split('\t')
            seq_len = int(splt[1])
            offset = int(splt[2])
            line_ln = int(splt[3])
            n_lines = seq_len // line_ln
            if seq_len % line_ln != 0:
                n_lines += 1
            ref_inds.append((splt[0], offset, offset + seq_len + n_lines, seq_len))
        fai.close()
        return ref_inds

    sys.stdout.write('index not found, creating one... ')
    sys.stdout.flush()
    ref_file = open(reference_path, 'r')
    prev_r = None
    prev_p = None
    seq_len = 0
    while 1:
        data = ref_file.readline()
        if not data:
            ref_inds.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            break
        if data[0] == '>':
            if prev_p is not None:
                ref_inds.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            seq_len = 0
            prev_p = ref_file.tell()
            prev_r = data[1:-1]
        else:
            seq_len += len(data) - 1
    ref_file.close()

    print('{0:.3f} (sec)'.format(time.time() - tt))
    return ref_inds


def read_ref(ref_path, ref_inds_i, n_handling, n_unknowns=True, quiet=False):
    OK_CHR_ORD = {'A': True, 'C': True, 'G': True, 'T': True, 'U': True}
    ALLOWED_NUCL = ['A', 'C', 'G', 'T']
    tt = time.time()
    if not quiet:
        sys.stdout.write('reading ' + ref_inds_i[0] + '... ')
        sys.stdout.flush()

    ref_file = open(ref_path, 'r')
    ref_file.seek(ref_inds_i[1])
    my_dat = ''.join(ref_file.read(ref_inds_i[2] - ref_inds_i[1]).split('\n'))
    my_dat = Seq(my_dat.upper())
    my_dat = my_dat.tomutable()

    # find N regions
    # data explanation: my_dat[n_atlas[0][0]:n_atlas[0][1]] = solid block of Ns
    prev_ni = 0
    n_count = 0
    n_atlas = []
    for i in range(len(my_dat)):
        if my_dat[i] == 'N' or (n_unknowns and my_dat[i] not in OK_CHR_ORD):
            if n_count == 0:
                prev_ni = i
            n_count += 1
            if i == len(my_dat) - 1:
                n_atlas.append((prev_ni, prev_ni + n_count))
        else:
            if n_count > 0:
                n_atlas.append((prev_ni, prev_ni + n_count))
            n_count = 0

    # handle N base-calls as desired
    n_info = {'all': [], 'big': [], 'non_N': []}
    if n_handling[0] == 'random':
        for region in n_atlas:
            n_info['all'].extend(region)
            if region[1] - region[0] <= n_handling[1]:
                for i in range(region[0], region[1]):
                    my_dat[i] = random.choice(ALLOWED_NUCL)
            else:
                n_info['big'].extend(region)
    elif n_handling[0] == 'allChr' and n_handling[2] in OK_CHR_ORD:
        for region in n_atlas:
            n_info['all'].extend(region)
            if region[1] - region[0] <= n_handling[1]:
                for i in range(region[0], region[1]):
                    my_dat[i] = n_handling[2]
            else:
                n_info['big'].extend(region)
    elif n_handling[0] == 'ignore':
        for region in n_atlas:
            n_info['all'].extend(region)
            n_info['big'].extend(region)
    else:
        print('\nERROR: UNKNOWN N_HANDLING MODE\n')
        exit(1)

    habitable_regions = []
    if not n_info['big']:
        n_info['non_N'] = [(0, len(my_dat))]
    else:
        for i in range(0, len(n_info['big']), 2):
            if i == 0:
                habitable_regions.append((0, n_info['big'][0]))
            else:
                habitable_regions.append((n_info['big'][i - 1], n_info['big'][i]))
        habitable_regions.append((n_info['big'][-1], len(my_dat)))
    for n in habitable_regions:
        if n[0] != n[1]:
            n_info['non_N'].append(n)

    if not quiet:
        print('{0:.3f} (sec)'.format(time.time() - tt))
    return my_dat, n_info


#	find all non-N regions in reference sequence ahead of time, for computing jobs in parallel
def get_all_ref_regions(ref_path, ref_inds, n_handling, save_output=False):
    out_regions = {}
    fn = ref_path + '.nnr'
    if os.path.isfile(fn) and not (save_output):
        print('found list of preidentified non-N regions...')
        f = open(fn, 'r')
        for line in f:
            splt = line.strip().split('\t')
            if splt[0] not in out_regions:
                out_regions[splt[0]] = []
            out_regions[splt[0]].append((int(splt[1]), int(splt[2])))
        f.close()
        return out_regions
    else:
        print('enumerating all non-N regions in reference sequence...')
        for RI in range(len(ref_inds)):
            (refSequence, N_regions) = read_ref(ref_path, ref_inds[RI], n_handling, quiet=True)
            ref_name = ref_inds[RI][0]
            out_regions[ref_name] = [n for n in N_regions['non_N']]
        if save_output:
            f = open(fn, 'w')
            for k in out_regions.keys():
                for n in out_regions[k]:
                    f.write(k + '\t' + str(n[0]) + '\t' + str(n[1]) + '\n')
            f.close()
        return out_regions


#	find which of the non-N regions are going to be used for this job
def partition_ref_regions(in_regions, ref_inds, my_job, n_jobs):
    tot_size = 0
    for RI in range(len(ref_inds)):
        ref_name = ref_inds[RI][0]
        for region in in_regions[ref_name]:
            tot_size += region[1] - region[0]
    size_per_job = int(tot_size / float(n_jobs) - 0.5)

    regions_per_job = [[] for n in range(n_jobs)]
    refs_per_job = [{} for n in range(n_jobs)]
    current_ind = 0
    current_count = 0
    for RI in range(len(ref_inds)):
        ref_name = ref_inds[RI][0]
        for region in in_regions[ref_name]:
            regions_per_job[current_ind].append((ref_name, region[0], region[1]))
            refs_per_job[current_ind][ref_name] = True
            current_count += region[1] - region[0]
            if current_count >= size_per_job:
                current_count = 0
                current_ind = min([current_ind + 1, n_jobs - 1])

    relevant_refs = refs_per_job[my_job - 1].keys()
    relevant_regs = regions_per_job[my_job - 1]
    return relevant_refs, relevant_regs
