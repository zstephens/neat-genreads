import random
import bisect
import copy
import sys
from typing import Union

import numpy as np

LOW_PROB_THRESH = 1e-12


def mean_ind_of_weighted_list(candidate_list: list) -> int:
    """
    Returns the index of the mean of a weighted list

    :param candidate_list: weighted list
    :return: index of mean
    """
    my_mid = sum(candidate_list) / 2.0
    my_sum = 0.0
    for i in range(len(candidate_list)):
        my_sum += candidate_list[i]
        if my_sum >= my_mid:
            return i


class DiscreteDistribution:
    def __init__(self, weights, values, degenerate_val=None, method='bisect'):

        # some sanity checking
        if not len(weights) or not len(values):
            print('\nError: weight or value vector given to DiscreteDistribution() are 0-length.\n')
            sys.exit(1)

        self.method = method
        sum_weight = float(sum(weights))

        # if probability of all input events is 0, consider it degenerate and always return the first value
        if sum_weight < LOW_PROB_THRESH:
            self.degenerate = values[0]
        else:
            self.weights = [n / sum_weight for n in weights]
            # TODO This line is slowing things down and seems unnecessary. Are these "values
            # possibly some thing from another class?
            self.values = copy.deepcopy(values)
            if len(self.values) != len(self.weights):
                print('\nError: length and weights and values vectors must be the same.\n')
                exit(1)
            self.degenerate = degenerate_val

            if self.method == 'alias':
                len_weights = len(self.weights)
                prob_vector = np.zeros(len_weights)
                count_vector = np.zeros(len_weights, dtype=np.int)
                smaller = []
                larger = []
                for kk, prob in enumerate(self.weights):
                    prob_vector[kk] = len_weights * prob
                    if prob_vector[kk] < 1.0:
                        smaller.append(kk)
                    else:
                        larger.append(kk)
                while len(smaller) > 0 and len(larger) > 0:
                    small = smaller.pop()
                    large = larger.pop()
                    count_vector[small] = large
                    prob_vector[large] = (prob_vector[large] + prob_vector[small]) - 1.0
                    if prob_vector[large] < 1.0:
                        smaller.append(large)
                    else:
                        larger.append(large)

                self.a1 = len(count_vector) - 1
                self.a2 = count_vector.tolist()
                self.a3 = prob_vector.tolist()

            elif self.method == 'bisect':
                self.cum_prob = np.cumsum(self.weights).tolist()[:-1]
                self.cum_prob.insert(0, 0.)

            else:
                print("\nUnknown discreet distribution method.\n")

    def __str__(self):
        return str(self.weights) + ' ' + str(self.values) + ' ' + self.method

    def sample(self) -> Union[int, float]:
        """
        This is one of the slowest parts of the code. Or it just gets hit the most times. Will need
        to investigate at some point.
        :return: Since this function is selecting an item from a list, and the list could theoretically be anything,
        then in a broad sense this function returns a list item or a generic object. But I'm fairly confident that most
        of these uses will be lists of ints or floats, but will investigate further
        """

        if self.degenerate is not None:
            return self.degenerate

        else:

            if self.method == 'alias':
                random1 = random.randint(0, self.a1)
                random2 = random.random()
                if random2 < self.a3[random1]:
                    return self.values[random1]
                else:
                    return self.values[self.a2[random1]]

            elif self.method == 'bisect':
                r = random.random()
                return self.values[bisect.bisect(self.cum_prob, r) - 1]


# takes k_range, lambda, [0,1,2,..], returns a DiscreteDistribution object
# with the corresponding to a poisson distribution

def poisson_list(k_range, input_lambda):
    min_weight = 1e-12
    if input_lambda < min_weight:
        return DiscreteDistribution([1], [0], degenerate_val=0)
    log_factorial_list = [0.0]
    for k in k_range[1:]:
        log_factorial_list.append(np.log(float(k)) + log_factorial_list[k - 1])
    w_range = [np.exp(k * np.log(input_lambda) - input_lambda - log_factorial_list[k]) for k in k_range]
    w_range = [n for n in w_range if n >= min_weight]
    if len(w_range) <= 1:
        return DiscreteDistribution([1], [0], degenerate_val=0)
    return DiscreteDistribution(w_range, k_range[:len(w_range)])


# quantize a list of values into blocks
def quantize_list(list_to_quantize):
    min_prob = 1e-12
    quant_blocks = 10
    sum_list = float(sum(list_to_quantize))
    sorted_list = sorted([n for n in list_to_quantize if n >= min_prob * sum_list])
    if len(sorted_list) == 0:
        return None
    qi = []
    for i in range(quant_blocks):
        # qi.append(sorted_list[int((i)*(len(sorted_list)/float(quant_blocks)))])
        qi.append(sorted_list[0] + (i / float(quant_blocks)) * (sorted_list[-1] - sorted_list[0]))
    qi.append(1e12)
    running_list = []
    prev_bi = None
    prev_i = None
    for i in range(len(list_to_quantize)):
        if list_to_quantize[i] >= min_prob * sum_list:
            bi = bisect.bisect(qi, list_to_quantize[i])
            # print i, l[i], qi[bi-1]
            if prev_bi is not None:
                if bi == prev_bi and prev_i == i - 1:
                    running_list[-1][1] += 1
                else:
                    running_list.append([i, i, qi[bi - 1]])
            else:
                running_list.append([i, i, qi[bi - 1]])
            prev_bi = bi
            prev_i = i
    return running_list
