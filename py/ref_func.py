import sys
import time
import os
import pathlib
import random
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

OK_CHR_ORD = {'A': True, 'C': True, 'G': True, 'T': True, 'U': True}
ALLOWED_NUCL = ['A', 'C', 'G', 'T']


class Dog:
    pass


def index_ref(reference_path: str) -> list:
    """
    Index reference fasta
    :param reference_path: string path to the reference
    :return: reference index in list from
    """
    tt = time.time()

    absolute_reference_location = pathlib.Path(reference_path)

    # sanity check
    if not absolute_reference_location.is_file():
        print("\nProblem reading the reference fasta file.\n")
        sys.exit(1)

    index_filename = None

    # check if the reference file already exists
    if absolute_reference_location.with_suffix('.fai').is_file():
        print('found index ' + absolute_reference_location.with_suffix('.fai'))
        index_filename = absolute_reference_location.with_suffix('.fai')
    elif absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai').is_file():
        print('found index ' + absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai'))
        index_filename = absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai')
    else:
        pass

    ref_indices = []
    if index_filename is not None:
        fai = open(index_filename, 'r')
        for line in fai:
            splt = line[:-1].split('\t')
            # Defined as the number of bases in the contig
            seq_len = int(splt[1])
            # Defined as the byte index where the contig sequence begins
            offset = int(splt[2])
            # Defined as bases per line in the Fasta file
            line_ln = int(splt[3])
            n_lines = seq_len // line_ln
            if seq_len % line_ln != 0:
                n_lines += 1
            # Item 3 in this gives you the byte position of the next contig, I believe
            ref_indices.append((splt[0], offset, offset + seq_len + n_lines, seq_len))
        fai.close()
        return ref_indices

    print('Index not found, creating one... ')
    ref_file = open(absolute_reference_location, 'r')
    prev_r = None
    prev_p = None
    seq_len = 0

    while True:
        data = ref_file.readline()
        if not data:
            ref_indices.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            break
        elif data[0] == '>':
            if prev_p is not None:
                ref_indices.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            seq_len = 0
            prev_p = ref_file.tell()
            prev_r = data[1:-1]
        else:
            seq_len += len(data) - 1
    ref_file.close()

    print('{0:.3f} (sec)'.format(time.time() - tt))
    return ref_indices


def read_ref(ref_path, ref_inds_i, n_handling, n_unknowns=True, quiet=False):
    tt = time.time()
    if not quiet:
        print('reading ' + ref_inds_i[0] + '... ')

    absolute_reference_path = pathlib.Path(ref_path)
    try:
        ref_file = open(absolute_reference_path, 'r')
    except IOError:
        print('\nProblem reading reference file.\n')
        sys.exit(1)

    # TODO convert to SeqIO containers
    ref_file.seek(ref_inds_i[1])
    my_dat = ''.join(ref_file.read(ref_inds_i[2] - ref_inds_i[1]).split('\n'))
    my_dat = Seq(my_dat.upper(), IUPAC.unambiguous_dna)
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
        sys.exit(1)

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

    ref_file.close()

    if not quiet:
        print('{0:.3f} (sec)'.format(time.time() - tt))

    return my_dat, n_info


def get_all_ref_regions(ref_path, ref_inds, n_handling, save_output=False):
    """
    Find all non-N regions in reference sequence ahead of time, for computing jobs in parallel

    :param ref_path:
    :param ref_inds:
    :param n_handling:
    :param save_output:
    :return:
    """
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
            (ref_sequence, N_regions) = read_ref(ref_path, ref_inds[RI], n_handling, quiet=True)
            ref_name = ref_inds[RI][0]
            out_regions[ref_name] = [n for n in N_regions['non_N']]
        if save_output:
            f = open(fn, 'w')
            for k in out_regions.keys():
                for n in out_regions[k]:
                    f.write(k + '\t' + str(n[0]) + '\t' + str(n[1]) + '\n')
            f.close()
        return out_regions


def partition_ref_regions(in_regions, ref_inds, my_job, n_jobs):
    """
    Find which of the non-N regions are going to be used for this job

    :param in_regions:
    :param ref_inds:
    :param my_job:
    :param n_jobs:
    :return:
    """
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
