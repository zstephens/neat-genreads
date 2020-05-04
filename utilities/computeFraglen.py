#
#
#      Compute Fragment Length Model for genReads.py
#                  computeFraglen.py
#
#
#      Usage: samtools view normal.bam | python computeFraglen.py
#
#
# Python 3 ready

import fileinput
import pickle
import argparse

FILTER_MAPQUAL = 10  # only consider reads that are mapped with at least this mapping quality
FILTER_MINREADS = 100  # only consider fragment lengths that have at least this many read pairs supporting it
FILTER_MEDDEV_M = 10  # only consider fragment lengths this many median deviations above the median

def quick_median(count_dict):
    midPoint = sum(count_dict.values()) // 2
    mySum = 0
    myInd = 0
    sk = sorted(count_dict.keys())
    while mySum < midPoint:
        mySum += count_dict[sk[myInd]]
        if mySum >= midPoint:
            break
        myInd += 1
    return myInd


def median_deviation_from_median(count_dict):
    myMedian = quick_median(count_dict)
    deviations = {}
    for k in sorted(count_dict.keys()):
        d = abs(k - myMedian)
        deviations[d] = count_dict[k]
    return quick_median(deviations)


def count_frags(file: str) -> dict:
    count_dict = {}
    PRINT_EVERY = 100000
    i = 0
    for line in fileinput.input():
        # Skip all comments and headers
        if line[0] == '#' or line[0] == '@':
            continue
        splt = line.strip().split('\t')
        samFlag = int(splt[1])
        myRef = splt[2]
        mapQual = int(splt[4])
        mateRef = splt[6]
        myTlen = abs(int(splt[8]))

        # if read is paired, and is first in pair, and is confidently mapped...
        if samFlag & 1 and samFlag & 64 and mapQual > FILTER_MAPQUAL:
            # and mate is mapped to same reference
            if mateRef == '=' or mateRef == myRef:
                if myTlen not in count_dict:
                    count_dict[myTlen] = 0
                count_dict[myTlen] += 1
                i += 1
                if i % PRINT_EVERY == 0:
                    print('---', i, quick_median(count_dict), median_deviation_from_median(count_dict))
    return count_dict


def compute_probs(count_dict: dict) -> (list, list):
    values = []
    probabilities = []
    med = quick_median(count_dict)
    mdm = median_deviation_from_median(count_dict)

    for k in sorted(count_dict.keys()):
        if k > 0 and k < med + FILTER_MEDDEV_M * mdm:
            if count_dict[k] >= FILTER_MINREADS:
                print(k, count_dict[k])
                values.append(k)
                probabilities.append(count_dict[k])
    countSum = float(sum(probabilities))
    probabilities = [n / countSum for n in probabilities]
    return values, probabilities

def main():
    parser = argparse.ArgumentParser(description="computeFraglen.py")
    parser.add_argument('-i', type=str, required=True, default=None, help="Sam file input (samtools view name.bam > name.sam")
    parser.add_argument('-o', type=str, required=True, default=None, help="Prefix for output")

    args = parser.parse_args()
    input_file = args.i
    output_prefix = args.o
    output = output_prefix + '.p'

    all_tlens = count_frags(input_file)
    print('\nsaving model...')
    out_vals, out_probs = compute_probs(all_tlens)
    print(out_probs)
    pickle.dump([out_vals, out_probs], open(output, 'wb'))


if __name__ == "__main()":
    main()
