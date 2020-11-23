#
#
#      Compute Fragment Length Model for gen_reads.source
#                  compute_fraglen.source
#
#
#      Usage: samtools view normal.bam | source compute_fraglen.source
#
#
# Upgraded 5/6/2020 to match Python 3 standards and refactored for easier reading

import pickle
import argparse
import platform

os = platform.system()
if os !='Windows':
    import pysam


def median(datalist: list) -> float:
    """
    Finds the median of a list of data. For this function, the data are expected to be a list of
    numbers, either float or int.
    :param datalist: the list of data to find the median of. This should be a set of numbers.
    :return: The median of the set
    >>> median([2])
    2
    >>> median([2183, 2292, 4064, 4795, 7471, 12766, 14603, 15182, 16803, 18704, 21504, 21677, 23347, 23586, 24612, 24878, 25310, 25993, 26448, 28018, 28352, 28373, 28786, 30037, 31659, 31786, 33487, 33531, 34442, 39138, 39718, 39815, 41518, 41934, 43301])
    25993
    >>> median([1,2,4,6,8,12,14,15,17,21])
    10.0
    """
    # using integer division here gives the index of the midpoint, due to zero-based indexing.
    midpoint = len(datalist)//2

    # Once we've found the midpoint, we calculate the median, which is just the middle value if there are an
    # odd number of values, or the average of the two middle values if there are an even number
    if len(datalist) % 2 == 0:
        median = (datalist[midpoint] + datalist[midpoint-1])/2
    else:
        median = datalist[midpoint]
    return median



def median_absolute_deviation(datalist: list) -> float:
    """
    Calculates the absolute value of the median deviation from the median for each element of of a datalist. 
    Then returns the median of these values.
    :param datalist: A list of data to find the MAD of
    :return: index of median of the deviations
    >>> median_absolute_deviation([2183, 2292, 4064, 4795, 7471, 12766, 14603, 15182, 16803, 18704, 21504, 21677, 23347, 23586, 24612, 24878, 25310, 25993, 26448, 28018, 28352, 28373, 28786, 30037, 31659, 31786, 33487, 33531, 34442, 39138, 39718, 39815, 41518, 41934, 43301])
    7494
    >>> median_absolute_deviation([1,2,4,6,8,12,14,15,17,21])
    5.5
    >>> median_absolute_deviation([0,2])
    1.0
    """
    my_median = median(datalist)
    deviations = []
    for item in datalist:
        # We take the absolute difference between the value and the median
        X_value = abs(item - my_median)
        # This creates a dataset that is the absolute deviations about the median
        deviations.append(X_value)
    # The median of the absolute deviations is the median absolute deviation
    return median(sorted(deviations))


def count_frags(file: str) -> list:
    """
    Takes a sam or bam file input and creates a list of the number of reads that are paired,
    first in the pair, confidently mapped and whose pair is mapped to the same reference
    :param file: A sam input file
    :return: A list of the tlens from the bam/sam file
    """
    FILTER_MAPQUAL = 10  # only consider reads that are mapped with at least this mapping quality
    count_list = []
    # Check if the file is sam or bam and decide how to open based on that
    if file[-4:] == ".sam":
        file_to_parse = open(file, 'r')
    elif file[-4:] == ".bam":
        print("WARNING: Must have pysam installed to read bam files. Pysam does not work on Windows OS.")
        if os != 'Windows':
            file_to_parse = pysam.AlignmentFile(file, 'rb')
        else:
            raise Exception("Your machine is running Windows. Please convert any BAM files to SAM files using samtools prior to input")
    else:
        print("Unknown file type, file extension must be bam or sam")
        exit(1)

    for item in file_to_parse:
        # Need to convert bam iterable objects into strings for the next part
        line = str(item)
        # Skip all comments and headers
        if line[0] == '#' or line[0] == '@':
            continue
        splt = line.strip().split('\t')
        sam_flag = int(splt[1])
        my_ref = splt[2]
        map_qual = int(splt[4])
        mate_ref = splt[6]
        my_tlen = abs(int(splt[8]))
        # if read is paired, and is first in pair, and is confidently mapped...
        if sam_flag & 1 and sam_flag & 64 and map_qual > FILTER_MAPQUAL:
            # and mate is mapped to same reference
            if mate_ref == '=' or mate_ref == my_ref:
                count_list.append(my_tlen)
    count_list = sorted(count_list)
    file_to_parse.close()
    return count_list


def compute_probs(datalist: list) -> (list, list):
    """
    Computes the probabilities for fragments with at least 100 pairs supporting it and that are at least 10 median
    deviations from the median.
    :param datalist: A list of fragments with counts
    :return: A list of values that meet the criteria and a list of their associated probabilities
    """
    FILTER_MINREADS = 100  # only consider fragment lengths that have at least this many read pairs supporting it
    FILTER_MEDDEV_M = 10  # only consider fragment lengths this many median deviations above the median
    values = []
    probabilities = []
    med = median(datalist)
    mad = median_absolute_deviation(datalist)

    for item in list(set(datalist)):
        if 0 < item <= med + FILTER_MEDDEV_M * mad:
            data_count = datalist.count(item)
            if data_count >= FILTER_MINREADS:
                values.append(item)
                probabilities.append(data_count)
    count_sum = float(sum(probabilities))
    probabilities = [n / count_sum for n in probabilities]
    return values, probabilities


def main():
    """
    Main function takes 2 arguments:
        input - a path to a sam or bam file input. Note that sam files can be formed by applying samtools to a bam file
        in the follawing way: samtools view nameof.bam > nameof.sam
        
        output - the string prefix of the output. The actual output will be the prefix plus ".p" at the end
        for pickle file. The list of values and list of probabilities are dumped as a list of lists
        into a pickle file on completion of the analysis

    :return: None
    """
    parser = argparse.ArgumentParser(description="compute_fraglen.source")
    parser.add_argument('-i', type=str, metavar="input", required=True, default=None,
                        help="Sam file input (samtools view name.bam > name.sam)")
    parser.add_argument('-o', type=str, metavar="output", required=True, default=None, help="Prefix for output")

    args = parser.parse_args()
    input_file = args.i
    output_prefix = args.o
    output = output_prefix + '.p'

    all_tlens = count_frags(input_file)
    print('\nSaving model...')
    out_vals, out_probs = compute_probs(all_tlens)
    pickle.dump([out_vals, out_probs], open(output, 'wb'))
    print('\nModel successfully saved.')

if __name__ == "__main__":
    main()