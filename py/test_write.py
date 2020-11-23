#!/bin/bash/python

import pathlib
import gzip
from timeit import default_timer as timer


class OutputFileWriter1:
    def __init__(self, out_prefix, gzipped, fasta):
        start = timer()

        if gzipped:
            file = pathlib.Path(out_prefix + '.fasta.gz')
        else:
            file = pathlib.Path(out_prefix + '.fasta')

        if gzipped:
            self.file = gzip.open(file, 'wb')
        else:
            self.file = open(file, 'w')

        fasta_path = pathlib.Path(fasta)

        for i in range(100):
            with open(fasta_path, 'r') as file:
                for line in file:
                    self.file.write(line)

        self.file.close()

        end = timer()

        print("It took {} seconds!".format(end-start))


class OutputFileWriter2:
    def __init__(self, out_prefix, gzipped, fasta):
        start = timer()

        if gzipped:
            self.file = pathlib.Path(out_prefix + '.fasta.gz')
        else:
            self.file = pathlib.Path(out_prefix + '.fasta')

        lines = []
        fasta_path = pathlib.Path(fasta)
        for i in range(100):
            with open(fasta_path, 'r') as file:
                for line in file:
                    lines.append(line)
                if i % 10 == 0:
                    if gzipped:
                        with gzip.open(self.file, 'wb') as f:
                            f.write("\n".join(lines))
                            lines = []
                    else:
                        with open(self.file, 'w') as f:
                            f.write("\n".join(lines))
                            lines = []

        end = timer()

        print("It took {} seconds!".format(end - start))


class OutputFileWriter3:
    def __init__(self, out_prefix, gzipped, fasta):
        start = timer()

        if gzipped:
            self.file = pathlib.Path(out_prefix + '.fasta.gz')
        else:
            self.file = pathlib.Path(out_prefix + '.fasta')

        lines = []
        fasta_path = pathlib.Path(fasta)
        for i in range(100):
            with open(fasta_path, 'r') as file:
                for line in file:
                    lines.append(line)

        if gzipped:
            with gzip.open(self.file, 'wb') as f:
                f.write("\n".join(lines))
        else:
            with open(self.file, 'w') as f:
                f.write("\n".join(lines))

        end = timer()

        print("It took {} seconds!".format(end - start))


fasta = '/home/joshfactorial/Documents/neat_data/chr21.fasta'
OutputFileWriter1('test1', False, fasta)
OutputFileWriter2('test2', False, fasta)
# OutputFileWriter3('test3', False, fasta)
