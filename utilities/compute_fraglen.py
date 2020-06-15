#
#
#      Compute Fragment Length Model for genReads.py
#                  compute_fraglen.py
#
#
#      Usage: samtools view normal.bam | python compute_fraglen.py
#
#
# Upgraded 5/6/2020 to match Python 3 standards and refactored for easier reading

import pickle
import argparse
import pysam


def quick_median(count_dict: dict) -> int:
    """
    Finds the median of a counting dictionary. A counting dictionary is a representation of a list of data.
    For example, the data [1,1,2,2,2] could be represented in dictionary form as {1:2, 2:3}, indicating that 1 is
    repeated twice and 2 is repeating three times.
    :param count_dict: the counting dictionary to find the median of
    :return: The median of the set
    """
    # The sum of the values of the counting dictionary tells us how long the expanded list is. To find the midpoint
    # We divide by 2
    if not count_dict:
        raise Exception('Empty counting dictionary. Check that sam input was correctly formatted.')
    else:
        mid_point = sum(count_dict.values()) / 2
        my_sum = 0
        my_ind = 0
        sk = sorted(count_dict.keys())
        # Here we basically count up one group of dictionary values at a time until we're at the midpoint
        while my_sum < mid_point:
            my_sum += count_dict[sk[my_ind]]
            if my_sum >= mid_point:
                break
            my_ind += 1
        # Once we've found the midpoint, we calculate the median, which is just the middle value if there are an
        # odd number of values, or the average of the two middle values if there are an even number
        if sum(count_dict.values()) % 2 == 0:
            median = (sk[my_ind] + sk[my_ind-1])//2
        else:
            median = sk[my_ind]
    return median


def median_absolute_deviation(count_dict: dict) -> int:
    """
    Calculates the deviation from the median of each element of counting dictionary,
    then returns the median of that dictionary. A counting dictionary such as {2: 3, 5:2} expands to [2,2,2,5,5]
    :param count_dict: Counting dictionary to analyze
    :return: index of median of the deviations
    """
    my_median = quick_median(count_dict)
    deviations = {}
    for key in sorted(count_dict.keys()):
        # We take the absolute difference between the value and the median
        X_value = abs(key - my_median)
        # We retain the number of items for each element of the dictionary
        # This is equivalent to subtracting the median from each value of the expanded dataset
        deviations[X_value] = count_dict[key]
    # Taking the median of this dataset gives us the median deviation from the median
    # (also called the Median Absolute Deviation)
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
    # Check if the file is sam or bam and decide how to open based on that
    if file[-4:] == ".sam":
        file_to_parse = open(file, 'r')
    elif file[-4:] == ".bam":
        file_to_parse = pysam.AlignmentFile(file, 'rb')
    else:
        print("Unknown file type, must be bam or sam")
        exit(1)

    for item in file_to_parse:
        # Need to convert bam iterable objects into strings for the next part
        line = str(item)
        # Skip all comments and headers
        if line[0] == '#' or line[0] == '@':
            continue
        splt = line.strip().split('\t')
        samFlag = int(splt[1])
        myRef = splt[2]
        mapQual = int(splt[4])
        mateRef = splt[6]
        my_tlen = abs(int(splt[8]))

        # if read is paired, and is first in pair, and is confidently mapped...
        if samFlag & 1 and samFlag & 64 and mapQual > FILTER_MAPQUAL:
            # and mate is mapped to same reference
            if mateRef == '=' or mateRef == myRef:
                if my_tlen not in count_dict:
                    count_dict[my_tlen] = 0
                count_dict[my_tlen] += 1
                i += 1
                if i % PRINT_EVERY == 0:
                    print('---', i, quick_median(count_dict), median_absolute_deviation(count_dict))
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
    mdm = median_absolute_deviation(count_dict)

    for key in sorted(count_dict.keys()):
        if 0 < key < med + FILTER_MEDDEV_M * mdm:
            if count_dict[key] >= FILTER_MINREADS:
                print(key, count_dict[key])
                values.append(key)
                probabilities.append(count_dict[key])
    count_sum = float(sum(probabilities))
    probabilities = [n / count_sum for n in probabilities]
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
    parser = argparse.ArgumentParser(description="compute_fraglen.py")
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
