import math
import random
import bisect
import numpy as np

class DiscreteDistribution:
	def __init__(self, weights, values, degenerateVal=None, method='bisect'):

		self.weights = [n/float(sum(weights)) for n in weights]
		self.values  = values
		self.method  = method
		if len(self.values) != len(self.weights):
			print '\nError: length and weights and values vectors must be the same.\n'
			exit(1)
		self.degenerate = degenerateVal

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

	def sample(self):

		if self.degenerate != None:
			return degenerateVal

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

