import numpy as np
import argparse
import sys
import pickle
import matplotlib.pyplot as mpl
import pathlib
import pysam
from functools import reduce

print("genSeqErrorModel starts running....", file=sys.stderr)
# enables import from neighboring package
sys.path.append(str(pathlib.Path(__file__).resolve().parents[1]))

from source.probability import DiscreteDistribution


def parse_file(input_file, real_q, off_q, max_reads, n_samp, plot_stuff):
    init_smooth = 0.
    prob_smooth = 0.

    # Takes a gzip or sam file and returns the simulation's average error rate,
    print('reading ' + input_file + '...', file=sys.stderr)
    is_aligned = False
    lines_to_read = 0
    try:
        if input_file[-4:] == '.bam' or input_file[-4:] == '.sam':
            print('detected aligned file....', file=sys.stderr)
            stats = pysam.idxstats(input_file).strip().split('\n')
            lines_to_read = reduce(lambda x, y: x + y, [eval('+'.join(l.rstrip('\n').split('\t')[2:])) for l in stats])
            f = pysam.AlignmentFile(input_file)
            is_aligned = True
        else:
            print('detected fastq file....', file=sys.stderr)
            with pysam.FastxFile(input_file) as f:
                for _ in f:
                    lines_to_read += 1
            f = pysam.FastxFile(input_file)
    except FileNotFoundError:
        print("Check input file. Must be fastq, gzipped fastq, or bam/sam file.", file=sys.stderr)
        sys.exit(1)

    actual_readlen = 0
    q_dict = {}
    current_line = 0
    quarters = lines_to_read // 4

    if is_aligned:
        g = f.fetch()
    else:
        g = f

    for read in g:
        if is_aligned:
            qualities_to_check = read.query_alignment_qualities
        else:
            qualities_to_check = read.get_quality_array()
        if actual_readlen == 0:
            actual_readlen = len(qualities_to_check) - 1
            print('assuming read length is uniform...', file=sys.stderr)
            print('detected read length (from first read found):', actual_readlen, file=sys.stderr)
            prior_q = np.zeros([actual_readlen, real_q])
            total_q = [None] + [np.zeros([real_q, real_q]) for n in range(actual_readlen - 1)]

        # sanity-check readlengths
        if len(qualities_to_check) - 1 != actual_readlen:
            print('skipping read with unexpected length...', file=sys.stderr)
            continue

        for i in range(actual_readlen):
            q = qualities_to_check[i]
            q_dict[q] = True
            prev_q = q
            if i == 0:
                prior_q[i][q] += 1
            else:
                total_q[i][prev_q, q] += 1
                prior_q[i][q] += 1

        current_line += 1
        if current_line % quarters == 0:
            print(f'{(current_line/lines_to_read)*100:.0f}%')
        if 0 < max_reads <= current_line:
            break

    f.close()

    # some sanity checking again...
    q_range = [min(q_dict.keys()), max(q_dict.keys())]
    if q_range[0] < 0:
        print('\nError: Read in Q-scores below 0\n', file=sys.stderr)
        exit(1)
    if q_range[1] > real_q:
        print('\nError: Read in Q-scores above specified maximum:', q_range[1], '>', real_q, '\n', file=sys.stderr)
        exit(1)

    print('computing probabilities...', file=sys.stderr)
    prob_q = [None] + [[[0. for m in range(real_q)] for n in range(real_q)] for p in range(actual_readlen - 1)]
    for p in range(1, actual_readlen):
        for i in range(real_q):
            row_sum = float(np.sum(total_q[p][i, :])) + prob_smooth * real_q
            if row_sum <= 0.:
                continue
            for j in range(real_q):
                prob_q[p][i][j] = (total_q[p][i][j] + prob_smooth) / row_sum

    init_q = [[0. for m in range(real_q)] for n in range(actual_readlen)]
    for i in range(actual_readlen):
        row_sum = float(np.sum(prior_q[i, :])) + init_smooth * real_q
        if row_sum <= 0.:
            continue
        for j in range(real_q):
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
        q_ticks_y = [(real_q - int(n)) - 0.5 for n in q_labels]

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
            mpl.ylim([real_q - q_range[1] - 1, real_q - q_range[0]])
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

    print('estimating average error rate via simulation...', file=sys.stderr)
    q_scores = range(real_q)
    # print (len(init_q), len(init_q[0]))
    # print (len(prob_q), len(prob_q[1]), len(prob_q[1][0]))

    init_dist_by_pos = [DiscreteDistribution(init_q[i], q_scores) for i in range(len(init_q))]
    prob_dist_by_pos_by_prev_q = [None]
    for i in range(1, len(init_q)):
        prob_dist_by_pos_by_prev_q.append([])
        for j in range(len(init_q[0])):
            if np.sum(prob_q[i][j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                prob_dist_by_pos_by_prev_q[-1].append(DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
            else:
                prob_dist_by_pos_by_prev_q[-1].append(DiscreteDistribution(prob_q[i][j], q_scores))

    count_dict = {}
    for q in q_scores:
        count_dict[q] = 0
    lines_to_sample = len(range(1, n_samp + 1))
    samp_quarters = lines_to_sample // 4
    for samp in range(1, n_samp + 1):
        if samp % samp_quarters == 0:
            print(f'{(samp/lines_to_sample)*100:.0f}%')
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
    print('AVG ERROR RATE:', avg_err, file=sys.stderr)

    return init_q, prob_q, avg_err


def main():
    parser = argparse.ArgumentParser(description='genSeqErrorModel.py')
    parser.add_argument('-i', type=str, required=True, metavar='<str>',
                        help="* input_read1.fq[.gz] / input_read1.b/sam")
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
    parser.add_argument('--plot', required=False, action='store_true', default=False,
                        help='perform some optional plotting')
    args = parser.parse_args()

    (infile, outfile, off_q, max_q, max_reads, n_samp) = (args.i, args.o, args.q, args.Q, args.n, args.s)
    (infile2, pile_up) = (args.i2, args.p)
    print("Read1 is {}, Read2 is {}, output model is {}".format(infile, infile2, outfile), file=sys.stderr)

    real_q = max_q + 1

    plot_stuff = args.plot

    q_scores = range(real_q)
    if infile2 is None:
        print('Read2 is not input. Generate Seq Error Model on single end reads', file=sys.stderr)
        (init_q, prob_q, avg_err) = parse_file(infile, real_q, off_q, max_reads, n_samp, plot_stuff)
    else:
        print('Read2 is inputted. Generate Seq Error Model on pair-end reads', file=sys.stderr)
        (init_q, prob_q, avg_err1) = parse_file(infile, real_q, off_q, max_reads, n_samp, plot_stuff)
        (init_q2, prob_q2, avg_err2) = parse_file(infile2, real_q, off_q, max_reads, n_samp, plot_stuff)
        avg_err = (avg_err1 + avg_err2) / 2.

    #
    #	embed some default sequencing error parameters if no pileup is provided
    #
    if pile_up == None:

        print('Using default sequencing error parameters...', file=sys.stderr)

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
        print('\nPileup parsing coming soon!\n', file=sys.stderr)
        exit(1)

    err_params = [sse_prob, sie_rate, sie_prob, sie_val, sie_ins_freq, sie_ins_nucl]

    #
    #	finally, let's save our output model
    #
    outfile = pathlib.Path(outfile).with_suffix(".p")
    print('saving model...', file=sys.stderr)
    if infile2 is None:
        pickle.dump([init_q, prob_q, q_scores, off_q, avg_err, err_params], open(outfile, 'wb'))
    else:
        pickle.dump([init_q, prob_q, init_q2, prob_q2, q_scores, off_q, avg_err, err_params], open(outfile, 'wb'))


if __name__ == '__main__':
    print("genSeqErrorModel starts executing the main function....", file=sys.stderr)
    main()
