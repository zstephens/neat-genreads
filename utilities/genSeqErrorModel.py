import os
import sys
import random
import numpy as np
import argparse
import cPickle as pickle

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-2])+'/py/'
sys.path.append(SIM_PATH)

from probability		import DiscreteDistribution

parser = argparse.ArgumentParser(description='fastq_to_qScoreModel.py')
parser.add_argument('-i',  type=str, required=True,  metavar='<str>',                  help="* input_read1.fq")
parser.add_argument('-o',  type=str, required=True,  metavar='<str>',                  help="* output.p")
parser.add_argument('-i2', type=str, required=False, metavar='<str>', default=None,    help="input_read2.fq")
parser.add_argument('-p',  type=str, required=False, metavar='<str>', default=None,    help="input_alignment.pileup")
parser.add_argument('-q',  type=int, required=False, metavar='<int>', default=33,      help="quality score offset [33]")
parser.add_argument('-Q',  type=int, required=False, metavar='<int>', default=41,      help="maximum quality score [41]")
parser.add_argument('-n',  type=int, required=False, metavar='<int>', default=-1,      help="maximum number of reads to process [all]")
parser.add_argument('-s',  type=int, required=False, metavar='<int>', default=1000000, help="number of simulation iterations [1000000]")
args = parser.parse_args()

(INF, OUF, offQ, maxQ, MAX_READS, N_SAMP) = (args.i, args.o, args.q, args.Q, args.n, args.s)
(INF2, PILEUP) = (args.i2, args.p)

RQ = maxQ+1

INIT_SMOOTH = 0.
PRINT_EVERY = 10000
PLOT_STUFF  = False
if PLOT_STUFF:
	import matplotlib.pyplot as mpl

def parseFQ(inf):
	print 'reading '+INF+'...'
	f = open(INF,'r')

	rRead  = 0
	actual_readlen = 0
	qDict  = {}
	while True:
		data1 = f.readline()
		data2 = f.readline()
		data3 = f.readline()
		data4 = f.readline()
		if not all([data1,data2,data3,data4]):
			break

		if actual_readlen == 0:
			actual_readlen = len(data4)-1
			priorQ = np.zeros([actual_readlen,RQ])
			totalQ = [None] + [np.zeros([RQ,RQ]) for n in xrange(actual_readlen-1)]

		# sanity-check readlengths
		if len(data4)-1 != actual_readlen:
			print 'skipping read with unexpected length...'
			continue

		for i in range(len(data4)-1):
			q = ord(data4[i])-offQ
			qDict[q] = True
			if i == 0:
				priorQ[i][q] += 1
			else:
				totalQ[i][prevQ,q] += 1
				priorQ[i][q] += 1
			prevQ = q

		rRead += 1
		if rRead%PRINT_EVERY == 0:
			print rRead
		if MAX_READS > 0 and rRead >= MAX_READS:
			break
	f.close()

	print 'computing probabilities...'
	probQ  = [None] + [[[0. for m in xrange(RQ)] for n in xrange(RQ)] for p in xrange(actual_readlen-1)]
	for p in xrange(1,actual_readlen):
		for i in xrange(RQ):
			rowSum = float(np.sum(totalQ[p][i,:]))
			if rowSum <= 0.:
				continue
			for j in xrange(RQ):
				probQ[p][i][j] = totalQ[p][i][j]/rowSum

	initQ  = [[INIT_SMOOTH for m in xrange(RQ)] for n in xrange(actual_readlen)]
	for i in xrange(actual_readlen):
		rowSum = float(np.sum(priorQ[i,:]))+INIT_SMOOTH*RQ
		if rowSum <= 0.:
			continue
		for j in xrange(RQ):
			initQ[i][j] = (priorQ[i][j]+INIT_SMOOTH)/rowSum

	if PLOT_STUFF:
		mpl.figure(1)
		Z = np.array(initQ).T
		X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
		mpl.pcolormesh(X,Y,Z,vmin=0.,vmax=0.25)
		mpl.axis([0,len(Z[0]),0,len(Z)])
		mpl.yticks(range(0,len(Z),10),range(0,len(Z),10))
		mpl.xticks(range(0,len(Z[0]),10),range(0,len(Z[0]),10))
		mpl.xlabel('Read Position')
		mpl.ylabel('Quality Score')
		mpl.title('Q-Score Prior Probabilities')
		mpl.colorbar()

		mpl.show()

		for p in xrange(1,actual_readlen):
			mpl.figure(p+1)
			Z = np.array(probQ[p])
			X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
			mpl.pcolormesh(X,Y,Z[::-1],vmin=0.,vmax=0.25)
			mpl.axis([0,len(Z[0]),0,len(Z)])
			mpl.yticks(range(0,len(Z)),range(len(Z)-1,-1,-1))
			mpl.xticks(range(0,len(Z[0])),range(0,len(Z[0])))
			mpl.xlabel('next Q')
			mpl.ylabel('current Q')
			mpl.title('Q-Score Transition Probabilities [pos:'+str(p)+']')
			mpl.colorbar()

			mpl.show()

	print 'estimating average error rate via simulation...'
	Qscores = range(RQ)
	#print (len(initQ), len(initQ[0]))
	#print (len(probQ), len(probQ[1]), len(probQ[1][0]))

	initDistByPos        = [DiscreteDistribution(initQ[i],Qscores) for i in xrange(len(initQ))]
	probDistByPosByPrevQ = [None]
	for i in xrange(1,len(initQ)):
		probDistByPosByPrevQ.append([])
		for j in xrange(len(initQ[0])):
			if np.sum(probQ[i][j]) <= 0.:	# if we don't have sufficient data for a transition, use the previous qscore
				probDistByPosByPrevQ[-1].append(DiscreteDistribution([1],[Qscores[j]],degenerateVal=Qscores[j]))
			else:
				probDistByPosByPrevQ[-1].append(DiscreteDistribution(probQ[i][j],Qscores))

	countDict = {}
	for q in Qscores:
		countDict[q] = 0
	for samp in xrange(1,N_SAMP+1):
		if samp%PRINT_EVERY == 0:
			print samp
		myQ = initDistByPos[0].sample()
		countDict[myQ] += 1
		for i in xrange(1,len(initQ)):
			myQ = probDistByPosByPrevQ[i][myQ].sample()
			countDict[myQ] += 1

	totBases = float(sum(countDict.values()))
	avgError = 0.
	for k in sorted(countDict.keys()):
		eVal = 10.**(-k/10.)
		#print k, eVal, countDict[k]
		avgError += eVal * (countDict[k]/totBases)
	print 'AVG ERROR RATE:',avgError

	return (initQ, probQ, avgError)

def main():

	Qscores = range(RQ)
	if INF2 == None:
		(initQ, probQ, avgError) = parseFQ(INF)
	else:
		(initQ, probQ, avgError1)   = parseFQ(INF)
		(initQ2, probQ2, avgError2) = parseFQ(INF2)
		avgError = (avgError1+avgError2)/2.

	#
	#	embed some default sequencing error parameters if no pileup is provided
	#
	if PILEUP == None:

		print 'Using default sequencing error parameters...'

		# sequencing substitution transition probabilities
		SSE_PROB   = [[0.,     0.4918, 0.3377, 0.1705 ],
					  [0.5238,     0., 0.2661, 0.2101 ],
					  [0.3754, 0.2355,     0., 0.3890 ],
					  [0.2505, 0.2552, 0.4942, 0.     ]]
		# if a sequencing error occurs, what are the odds it's an indel?
		SIE_RATE     = 0.01
		# sequencing indel error length distribution
		SIE_PROB     = [0.999,0.001]
		SIE_VAL      = [1,2]
		# if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
		SIE_INS_FREQ = 0.4
		# if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
		SIE_INS_NUCL = [0.25, 0.25, 0.25, 0.25]

	#
	#	otherwise we need to parse a pileup and compute statistics!
	#
	else:
		print '\nPileup parsing coming soon!\n'
		exit(1)

	errorParams  = [SSE_PROB, SIE_RATE, SIE_PROB, SIE_VAL, SIE_INS_FREQ, SIE_INS_NUCL]

	#
	#	finally, let's save our output model
	#
	print 'saving model...'
	if INF2 == None:
		pickle.dump([initQ,probQ,Qscores,offQ,avgError,errorParams],open(OUF,'wb'))
	else:
		pickle.dump([initQ,probQ,initQ2,probQ2,Qscores,offQ,avgError,errorParams],open(OUF,'wb'))

if __name__ == '__main__':
	main()
