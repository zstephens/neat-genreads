#
#
#   computeGC.py
#   Compute GC and coverage model for genReads.py
#
#   Takes output file from bedtools genomecov to generate GC/coverage model
#
#   Usage: bedtools genomecov -d -ibam input.bam -g reference.fa > genomeCov.dat
#          python computeGC.py -r reference.fa -i genomeCov.dat -w [sliding window length] -o output_name.p
#
#
# Python 3 ready

import time
import argparse
import numpy as np
import pickle
from Bio import SeqIO


def process_fasta(file: str) -> dict:
    ref_dict = {}

    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
        # ref_dict = {rec.id: rec.seq for rec in SeqIO.parse(file, "fasta")}
    except UnicodeDecodeError:
        print("Input file incorrect: -r should specify the reference fasta")
        exit(1)

    if not ref_dict:
        print("Input file incorrect: -r should specify the reference fasta")
        exit(1)

    return ref_dict


def capitalize_ref(input_dict: dict) -> dict:
    for k in sorted(input_dict.keys()):
        print("Capitalizing " + k)
        input_dict[k] = ''.join(input_dict[k])
        input_dict[k] = input_dict[k].upper()


def main():
    parser = argparse.ArgumentParser(description='computeGC.py')
    parser.add_argument('-i', type=str, required=True, metavar='input', help="input.genomecov")
    parser.add_argument('-r', type=str, required=True, metavar='reference', help="reference.fasta")
    parser.add_argument('-o', type=str, required=True, metavar='output prefix',
                        help="prefix for output (/path/to/output)")
    parser.add_argument('-w', type=int, required=False, metavar='sliding window',
                        help="sliding window length [50]", default=50)
    args = parser.parse_args()

    (IN_GCB, REF_FILE, WINDOW_SIZE, OUT_P) = (args.i, args.r, args.w, args.o)

    GC_BINS = {n: [] for n in range(WINDOW_SIZE + 1)}

    print('reading ref...')
    allrefs = process_fasta(REF_FILE)

    print('capitalizing ref...')
    allrefs = capitalize_ref(allrefs)

    print('reading genomecov file...')
    tt = time.time()
    f = open(IN_GCB, 'r')
    currentLine = 0
    currentRef = None
    currentCov = 0
    linesProcessed = 0
    PRINT_EVERY = 1000000
    STOP_AFTER = 1000000
    for line in f:
        splt = line.strip().split('\t')
        if linesProcessed % PRINT_EVERY == 0:
            print(linesProcessed)
        linesProcessed += 1

        # if linesProcessed > STOP_AFTER:
        #	break

        if currentLine == 0:
            currentRef = splt[0]
            sPos = int(splt[1]) - 1

        if currentRef not in allrefs:
            continue

        currentLine += 1
        currentCov += float(splt[2])

        if currentLine == WINDOW_SIZE:
            currentLine = 0
            seq = allrefs[currentRef][sPos:sPos + WINDOW_SIZE]
            if 'N' not in seq:
                gc_count = seq.count('G') + seq.count('C')
                GC_BINS[gc_count].append(currentCov)
            currentCov = 0

    f.close()

    runningTot = 0
    allMean = 0.0
    for k in sorted(GC_BINS.keys()):
        if len(GC_BINS[k]) == 0:
            print('{0:0.2%}'.format(k / float(WINDOW_SIZE)), 0.0, 0)
            GC_BINS[k] = 0
        else:
            myMean = np.mean(GC_BINS[k])
            myLen = len(GC_BINS[k])
            print('{0:0.2%}'.format(k / float(WINDOW_SIZE)), myMean, myLen)
            allMean += myMean * myLen
            runningTot += myLen
            GC_BINS[k] = myMean

    avgCov = allMean / float(runningTot)
    print('AVERAGE COVERAGE =', avgCov)

    y_out = []
    for k in sorted(GC_BINS.keys()):
        GC_BINS[k] /= avgCov
        y_out.append(GC_BINS[k])

    print('saving model...')
    pickle.dump([range(WINDOW_SIZE + 1), y_out], open(OUT_P, 'wb'))

    print(time.time() - tt, '(sec)')


if __name__ == "__main__":
    main()