import math
import random
import bisect
import copy
import numpy as np

LOW_PROB_THRESH = 1e-12

def mean_ind_of_weighted_list(l):
	myMid = sum(l)/2.0
	mySum = 0.0
	for i in xrange(len(l)):
		mySum += l[i]
		if mySum >= myMid:
			return i

class DiscreteDistribution:
	def __init__(self, weights, values, degenerateVal=None, method='bisect'):

		# some sanity checking
		if not len(weights) or not len(values):
			print '\nError: weight or value vector given to DiscreteDistribution() are 0-length.\n'
			exit(1)

		self.method  = method
		sumWeight    = float(sum(weights))
		
		# if probability of all input events is 0, consider it degenerate and always return the first value
		if sumWeight < LOW_PROB_THRESH:
			self.degenerate = values[0]
		else:
			self.weights = [n/sumWeight for n in weights]
			self.values  = copy.deepcopy(values)
			if len(self.values) != len(self.weights):
				print '\nError: length and weights and values vectors must be the same.\n'
				exit(1)
			self.degenerate = degenerateVal
			# prune values with probability too low to be worth using [DOESN'T REALLY IMPROVE PERFORMANCE]
			####if self.degenerate != None:
			####	for i in xrange(len(self.weights)-1,-1,-1):
			####		if self.weights[i] < LOW_PROB_THRESH:
			####			del self.weights[i]
			####			del self.values[i]
			####	if len(self.weights) == 0:
			####		print '\nError: probability distribution has no usable values.\n'
			####		exit(1)

			if self.method == 'alias':
				K       = len(self.weights)
				q       = np.zeros(K)
				J       = np.zeros(K, dtype=np.int)
				smaller = []
				larger  = []
				for kk, prob in enumerate(self.weights):
					q[kk] = K*prob
					if q[kk] < 1.0:
						smaller.append(kk)
					else:
						larger.append(kk)
				while len(smaller) > 0 and len(larger) > 0:
					small = smaller.pop()
					large = larger.pop()
					J[small] = large
					q[large] = (q[large] + q[small]) - 1.0
					if q[large] < 1.0:
						smaller.append(large)
					else:
						larger.append(large)

				self.a1 = len(J)-1
				self.a2 = J.tolist()
				self.a3 = q.tolist()

			elif self.method == 'bisect':
				self.cumP = np.cumsum(self.weights).tolist()[:-1]
				self.cumP.insert(0,0.)

	def __str__(self):
		return str(self.weights)+' '+str(self.values)+' '+self.method

	def sample(self):

		if self.degenerate != None:
			return self.degenerate

		else:

			if self.method == 'alias':
				r1 = random.randint(0,self.a1)
				r2 = random.random()
				if r2 < self.a3[r1]:
					return self.values[r1]
				else:
					return self.values[self.a2[r1]]

			elif self.method == 'bisect':
				r = random.random()
				return self.values[bisect.bisect(self.cumP,r)-1]


# takes k_range, lambda, [0,1,2,..], returns a DiscreteDistribution object with the corresponding to a poisson distribution
MIN_WEIGHT = 1e-12
def poisson_list(k_range,l):
	if l < MIN_WEIGHT:
		return DiscreteDistribution([1],[0],degenerateVal=0)
	logFactorial_list = [0.0]
	for k in k_range[1:]:
		logFactorial_list.append(np.log(float(k))+logFactorial_list[k-1])
	w_range = [np.exp(k*np.log(l) - l - logFactorial_list[k]) for k in k_range]
	w_range = [n for n in w_range if n >= MIN_WEIGHT]
	if len(w_range) <= 1:
		return DiscreteDistribution([1],[0],degenerateVal=0)
	return DiscreteDistribution(w_range,k_range[:len(w_range)])

# quantize a list of values into blocks
MIN_PROB = 1e-12
QUANT_BLOCKS = 10
def quantize_list(l):
	suml = float(sum(l))
	ls = sorted([n for n in l if n >= MIN_PROB*suml])
	if len(ls) == 0:
		return None
	qi = []
	for i in xrange(QUANT_BLOCKS):
		#qi.append(ls[int((i)*(len(ls)/float(QUANT_BLOCKS)))])
		qi.append(ls[0]+(i/float(QUANT_BLOCKS))*(ls[-1]-ls[0]))
	qi.append(1e12)
	runningList = []
	prevBi = None
	previ  = None
	for i in xrange(len(l)):
		if l[i] >= MIN_PROB*suml:
			bi = bisect.bisect(qi,l[i])
			#print i, l[i], qi[bi-1]
			if prevBi != None:
				if bi == prevBi and previ == i-1:
					runningList[-1][1] += 1
				else:
					runningList.append([i,i,qi[bi-1]])
			else:
				runningList.append([i,i,qi[bi-1]])
			prevBi = bi
			previ  = i
	return runningList

