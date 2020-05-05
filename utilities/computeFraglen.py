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
import pysam


def quick_median(count_dict: dict) -> int:
    """
    Finds the median of a counting dictionary
    :param count_dict: the counting dictionary to find the median of
    :return: integer index of the location of the median
    """
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


def median_deviation_from_median(count_dict: dict) -> int:
    """
    calculates the deviation from the median of each element of counting dictionary,
    then returns the median of that dictionary
    :param count_dict: Counting dictionary to analyze
    :return: index of median of the deviations
    """
    myMedian = quick_median(count_dict)
    deviations = {}
    for k in sorted(count_dict.keys()):
        d = abs(k - myMedian)
        deviations[d] = count_dict[k]
    return quick_median(deviations)


def count_frags(file: str) -> dict:
    """
    Takes a sam file input and creates a counting dictionary of the number of reads that are paired,
    first in the pair, confidently mapped and whose pair is mapped to the same reference
    :param file: A sam input file
    :return: A dictionary of the counts of the above reads
    """
    FILTER_MAPQUAL = 10  # only consider reads that are mapped with at least this mapping quality
    count_dict = {}
    PRINT_EVERY = 100000
    i = 0
    with open(file, 'r') as f:
        for line in f:
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
    """
    Computes the probabilities for fragments with at least 100 pairs supporting it and that are at least 10 median
    deviations from the median.
    :param count_dict: A dictionary of fragments with counts
    :return: A list of values that meet the criteria and a list of -their associated probabilities
    """
    FILTER_MINREADS = 100  # only consider fragment lengths that have at least this many read pairs supporting it
    FILTER_MEDDEV_M = 10  # only consider fragment lengths this many median deviations above the median
    values = []
    probabilities = []
    med = quick_median(count_dict)
    mdm = median_deviation_from_median(count_dict)

    for k in sorted(count_dict.keys()):
        if 0 < k < med + FILTER_MEDDEV_M * mdm:
            if count_dict[k] >= FILTER_MINREADS:
                print(k, count_dict[k])
                values.append(k)
                probabilities.append(count_dict[k])
    countSum = float(sum(probabilities))
    probabilities = [n / countSum for n in probabilities]
    return values, probabilities


def main():
    """
    Main function takes 2 arguments:
        input - a samfile input that can be formed by applying samtools to a bam file
        in the follawing way: samtools view nameof.bam > nameof.sam

        output - the prefix of the output. The actual output will be the prefix plus ".p" at the end
        for pickle file. The list of values and list of probabilities are dumped as a list of lists
        into a pickle file on completion of the analysis

    :return: None
    """
    parser = argparse.ArgumentParser(description="computeFraglen.py")
    parser.add_argument('-i', type=str, metavar="input", required=True, default=None,
                        help="Sam file input (samtools view name.bam > name.sam)")
    parser.add_argument('-o', type=str, metavar="output", required=True, default=None, help="Prefix for output")

    args = parser.parse_args()
    input_file = args.i
    output_prefix = args.o
    output = output_prefix + '.p'

    all_tlens = count_frags(input_file)
    print('\nsaving model...')
    out_vals, out_probs = compute_probs(all_tlens)
    print(out_probs)
    pickle.dump([out_vals, out_probs], open(output, 'wb'))


if __name__ == "__main__":
    main()
