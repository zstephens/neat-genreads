#!/usr/bin/env python

#
#
#          genSeqErrorModel.py
#          Computes sequencing error model for genReads.py
#
#         
#          Usage: python genSeqErrorModel.py -i input_reads.fq -o path/to/output_name.p
#
#
# Python 3 ready


import os
import gzip
import numpy as np
import argparse
import sys
import pickle

# absolute path to this script
sim_path = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/py/'
sys.path.append(sim_path)

from probability import DiscreteDistribution


def parse_fq(inf):
    #Takes a gzip or sam file and returns the simulation's average error rate, 
    print('reading ' + inf + '...')
    if inf[-3:] == '.gz':
        print('detected gzip suffix...')
        f = gzip.open(inf, 'r')
    else:
        f = open(inf, 'r')

    is_sam = False
    if inf[-4:] == '.sam':
        print('detected sam input...')
        is_sam = True

    r_read = 0
    actual_readlen = 0
    q_dict = {}
    while True:

        if is_sam:
            data4 = f.readline()
            if not len(data4):
                break
            try:
                data4 = data4.split('\t')[10]
            except IndexError:
                break
        # need to add some input checking here? Yup, probably.
        else:
            data1 = f.readline()
            data2 = f.readline()
            data3 = f.readline()
            data4 = f.readline()
            if not all([data1, data2, data3, data4]):
                break

        if actual_readlen == 0:
            if inf[-3:] != '.gz' and not is_sam:
                total_size = os.path.getsize(inf)
                entry_size = sum([len(n) for n in [data1, data2, data3, data4]])
                print('estimated number of reads in file:', int(float(total_size) / entry_size))
            actual_readlen = len(data4) - 1
            print('assuming read length is uniform...')
            print('detected read length (from first read found):', actual_readlen)
            prior_q = np.zeros([actual_readlen, RQ])
            total_q = [None] + [np.zeros([RQ, RQ]) for n in range(actual_readlen - 1)]

        # sanity-check readlengths
        if len(data4) - 1 != actual_readlen:
            print('skipping read with unexpected length...')
            continue

        for i in range(len(data4) - 1):
            q = data4[i] - off_q
            q_dict[q] = True
            prev_q = q
            if i == 0:
                prior_q[i][q] += 1
            else:
                total_q[i][prev_q, q] += 1
                prior_q[i][q] += 1
            

        r_read += 1
        if r_read % print_every == 0:
            print(r_read)
        if max_reads > 0 and r_read >= max_reads:
            break
    f.close()

    # some sanity checking again...
    q_range = [min(q_dict.keys()), max(q_dict.keys())]
    if q_range[0] < 0:
        print('\nError: Read in Q-scores below 0\n')
        exit(1)
    if q_range[1] > RQ:
        print('\nError: Read in Q-scores above specified maximum:', q_range[1], '>', RQ, '\n')
        exit(1)

    print('computing probabilities...')
    prob_q = [None] + [[[0. for m in range(RQ)] for n in range(RQ)] for p in range(actual_readlen - 1)]
    for p in range(1, actual_readlen):
        for i in range(RQ):
            row_sum = float(np.sum(total_q[p][i, :])) + prob_smooth * RQ
            if row_sum <= 0.:
                continue
            for j in range(RQ):
                prob_q[p][i][j] = (total_q[p][i][j] + prob_smooth) / row_sum

    init_q = [[0. for m in range(RQ)] for n in range(actual_readlen)]
    for i in range(actual_readlen):
        row_sum = float(np.sum(prior_q[i, :])) + init_smooth * RQ
        if row_sum <= 0.:
            continue
        for j in range(RQ):
            init_q[i][j] = (prior_q[i][j] + init_smooth) / row_sum

    if plot_stuff:
        mpl.rcParams.update({'font.size': 14, 'font.weight': 'bold', 'lines.linewidth': 3})

        mpl.figure(1)
        Z = np.array(init_q).T
        X, Y = np.meshgrid(range(0, len(Z[0]) + 1), range(0, len(Z) + 1))
        mpl.pcolormesh(X, Y, Z, vmin=0., vmax=0.25)
        mpl.axis([0, len(Z[0]), 0, len(Z)])
        mpl.yticks(range(0, len(Z), 10), range(0, len(Z), 10))
        mpl.xticks(range(0, len(Z[0]), 10), range(0, len(Z[0]), 10))
        mpl.xlabel('Read Position')
        mpl.ylabel('Quality Score')
        mpl.title('Q-Score Prior Probabilities')
        mpl.colorbar()

        mpl.show()

        v_min_log = [-4, 0]
        min_val = 10 ** v_min_log[0]
        q_labels = [str(n) for n in range(q_range[0], q_range[1] + 1) if n % 5 == 0]
        print(q_labels)
        q_ticks_x = [int(n) + 0.5 for n in q_labels]
        q_ticks_y = [(RQ - int(n)) - 0.5 for n in q_labels]

        for p in range(1, actual_readlen, 10):
            current_data = np.array(prob_q[p])
            for i in range(len(current_data)):
                for j in range(len(current_data[i])):
                    current_data[i][j] = max(min_val, current_data[i][j])

            # matrix indices:		pcolormesh plotting:	plot labels and axes:
            #
            #      y				   ^					   ^
            #	   -->				 x |					 y |
            #  x |					    -->					    -->
            #    v 					    y					    x
            #
            # to plot a MxN matrix 'Z' with rowNames and colNames we need to:
            #
            # pcolormesh(X,Y,Z[::-1,:])		# invert x-axis
            # # swap x/y axis parameters and labels, remember x is still inverted:
            # xlim([yMin,yMax])
            # ylim([M-xMax,M-xMin])
            # xticks()
            #

            mpl.figure(p + 1)
            z = np.log10(current_data)
            x, y = np.meshgrid(range(0, len(Z[0]) + 1), range(0, len(Z) + 1))
            mpl.pcolormesh(x, y, z[::-1, :], vmin=v_min_log[0], vmax=v_min_log[1], cmap='jet')
            mpl.xlim([q_range[0], q_range[1] + 1])
            mpl.ylim([RQ - q_range[1] - 1, RQ - q_range[0]])
            mpl.yticks(q_ticks_y, q_labels)
            mpl.xticks(q_ticks_x, q_labels)
            mpl.xlabel('\n' + r'$Q_{i+1}$')
            mpl.ylabel(r'$Q_i$')
            mpl.title('Q-Score Transition Frequencies [Read Pos:' + str(p) + ']')
            cb = mpl.colorbar()
            cb.set_ticks([-4, -3, -2, -1, 0])
            cb.set_ticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'])

        # mpl.tight_layout()
        mpl.show()

    print('estimating average error rate via simulation...')
    q_scores = range(RQ)
    # print (len(init_q), len(init_q[0]))
    # print (len(prob_q), len(prob_q[1]), len(prob_q[1][0]))

    init_dist_by_pos = [DiscreteDistribution(init_q[i], q_scores) for i in range(len(init_q))]
    prob_dist_by_pos_by_prev_q = [None]
    for i in range(1, len(init_q)):
        prob_dist_by_pos_by_prev_q.append([])
        for j in range(len(init_q[0])):
            if np.sum(prob_q[i][j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                prob_dist_by_pos_by_prev_q[-1].append(DiscreteDistribution([1], [q_scores[j]], degenerateVal=q_scores[j]))
            else:
                prob_dist_by_pos_by_prev_q[-1].append(DiscreteDistribution(prob_q[i][j], q_scores))

    count_dict = {}
    for q in q_scores:
        count_dict[q] = 0
    for samp in range(1, n_samp + 1):
        if samp % print_every == 0:
            print(samp)
        my_q = init_dist_by_pos[0].sample()
        count_dict[my_q] += 1
        for i in range(1, len(init_q)):
            my_q = prob_dist_by_pos_by_prev_q[i][my_q].sample()
            count_dict[my_q] += 1

    tot_bases = float(sum(count_dict.values()))
    avg_err = 0.
    for k in sorted(count_dict.keys()):
        eVal = 10. ** (-k / 10.)
        # print k, eVal, count_dict[k]
        avg_err += eVal * (count_dict[k] / tot_bases)
    print('AVG ERROR RATE:', avg_err)

    return init_q, prob_q, avg_err


parser = argparse.ArgumentParser(description='genSeqErrorModel.py')
parser.add_argument('-i', type=str, required=True, metavar='<str>', help="* input_read1.fq (.gz) / input_read1.sam")
parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output.p")
parser.add_argument('-i2', type=str, required=False, metavar='<str>', default=None,
                    help="input_read2.fq (.gz) / input_read2.sam")
parser.add_argument('-p', type=str, required=False, metavar='<str>', default=None, help="input_alignment.pileup")
parser.add_argument('-q', type=int, required=False, metavar='<int>', default=33, help="quality score offset [33]")
parser.add_argument('-Q', type=int, required=False, metavar='<int>', default=41, help="maximum quality score [41]")
parser.add_argument('-n', type=int, required=False, metavar='<int>', default=-1,
                    help="maximum number of reads to process [all]")
parser.add_argument('-s', type=int, required=False, metavar='<int>', default=1000000,
                    help="number of simulation iterations [1000000]")
parser.add_argument('--plot', required=False, action='store_true', default=False, help='perform some optional plotting')
args = parser.parse_args()

(inf, ouf, off_q, max_q, max_reads, n_samp) = (args.i, args.o, args.q, args.Q, args.n, args.s)
(inf2, pile_up) = (args.i2, args.p)

RQ = max_q + 1

init_smooth = 0.
prob_smooth = 0.
print_every = 10000
plot_stuff = args.plot
if plot_stuff:
    print('plotting is desired, lets import matplotlib...')
    import matplotlib.pyplot as mpl


def main():
    q_scores = range(RQ)
    if inf2 == None:
        (init_q, prob_q, avg_err) = parse_fq(inf)
    else:
        (init_q, prob_q, avg_err1) = parse_fq(inf)
        (init_q2, prob_q2, avg_err2) = parse_fq(inf2)
        avg_err = (avg_err1 + avg_err2) / 2.

    #
    #	embed some default sequencing error parameters if no pileup is provided
    #
    if pile_up == None:

        print('Using default sequencing error parameters...')

        # sequencing substitution transition probabilities
        sse_prob = [[0., 0.4918, 0.3377, 0.1705],
                    [0.5238, 0., 0.2661, 0.2101],
                    [0.3754, 0.2355, 0., 0.3890],
                    [0.2505, 0.2552, 0.4942, 0.]]
        # if a sequencing error occurs, what are the odds it's an indel?
        sie_rate = 0.01
        # sequencing indel error length distribution
        sie_prob = [0.999, 0.001]
        sie_val = [1, 2]
        # if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
        sie_ins_freq = 0.4
        # if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
        sie_ins_nucl = [0.25, 0.25, 0.25, 0.25]

    #
    #	otherwise we need to parse a pileup and compute statistics!
    #
    else:
        print('\nPileup parsing coming soon!\n')
        exit(1)

    err_params = [sse_prob, sie_rate, sie_prob, sie_val, sie_ins_freq, sie_ins_nucl]

    #
    #	finally, let's save our output model
    #
    print('saving model...')
    if inf2 == None:
        pickle.dump([init_q, prob_q, q_scores, off_q, avg_err, err_params], open(ouf, 'wb'))
    else:
        pickle.dump([init_q, prob_q, init_q2, prob_q2, q_scores, off_q, avg_err, err_params], open(ouf, 'wb'))


if __name__ == '__main__':
    main()
