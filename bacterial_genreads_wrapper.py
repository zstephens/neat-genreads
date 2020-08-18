#!/usr/bin/env python

import genReads
import argparse
import random
import pathlib
import gzip
import shutil
# from Bio import SeqIO


class Bacterium:
    def __init__(self, reference: str, name: str):
        self.reference = pathlib.Path(reference)
        self.name = name
        # Temporarily set the reference as the bacterium's file, until it is analyzed
        self.file = pathlib.Path(reference)
        self.analyze()

    def __repr__(self):
        return str(self.name)

    def __str__(self):
        return str(self.name)

    def get_file(self):
        return self.file

    def analyze(self):
        args = ['-r', str(self.reference), '-R', '101', '-o', self.name, '--fa']
        genReads.main(args)

        # The following workaround is due to the fact that genReads writes out a compressed
        # fasta but does not put the .gz extension on it. Also, genReads cannot handle gzipped
        # fasta files, so we further have to unzip it for it to actually work.
        self.file = pathlib.Path().absolute() / (self.name + "_read1.fa")
        print(self.file)
        new_name = self.name + "_read1.fa.gz"
        self.file.rename(pathlib.Path(pathlib.Path().absolute(), new_name))
        self.file = pathlib.Path().absolute() / (self.name + "_read1.fa.gz")
        print(self.file)
        true_path = pathlib.Path().absolute() / (self.name + "_read1.fa")
        unzip_file(self.file, true_path)
        pathlib.Path.unlink(pathlib.Path().absolute() / (self.name + "_read1.fa.gz")) # deletes unused zip file
        self.file = true_path
        # end workaround

        # Now we further have to fix the fasta file, which outputs in a form that doesn't make much sense,
        # so that it can be properly analyzed in the next generation by genreads.
        chromosome_name = ""
        sequence = ""
        with self.file.open() as f:
            first_line = f.readline()
            if first_line.startswith(">"):
                chromosome_name = first_line
            else:
                print("Something went wrong with the fasta file")
                exit(1)
            for line in f:
                if line.startswith(">"):
                    continue
                else:
                    sequence = sequence + line
            f.close()
        # re-write file with just the chrom name and sequence
        self.file.open('w').write(chromosome_name + sequence)



    def remove(self):
        pathlib.Path.unlink(self.file)


def unzip_file(zipped_file: pathlib, unzipped_file: pathlib):
    with gzip.open(zipped_file, 'rb') as f_in:
        with open(unzipped_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def cull(population: list, percentage: float = 0.5) -> list:
    """
    The purpose of this function will be to cull the bacteria created in the model
    :param percentage: percentage of the population to eliminate
    :param population: the list of members to cull
    :return: The list of remaining members
    """
    cull_amount = round(len(population) * percentage)
    print("Culling {} members from population".format(cull_amount))
    for i in range(cull_amount):
        selection = random.choice(population)
        population.remove(selection)
        selection.remove()
    return population


def crossover(population: list) -> list:
    """
    This function will take a list of individuals and do crossover type mixing
    :param population:
    :return new_population:
    """
    new_population = population
    return new_population


def initialize_population(reference: str, pop_size) -> list:
    """
    The purpose of this function is to evolve the bacteria
    :param reference:
    :param pop_size:
    :return population:
    """
    names = []
    for j in range(pop_size):
        names.append("bacterium_0_{}".format(j+1))
    population = []
    for i in range(pop_size):
        new_member = Bacterium(reference, names[i])
        population.append(new_member)
    return population


def evolve_population(population: list, generation: int) -> list:
    """
    This evolves an existing population by doubling them (binary fission), then introducing random mutation to
    each member of the population.
    :param generation: Helps determine the starting point of the numbering system so the bacteria have unique names
    :param population: A list of fasta files representing the bacteria.
    :return:
    """
    children_population = population + population
    names = []
    new_population = []
    for j in range(len(children_population)):
        names.append("bacterium_{}_{}".format(generation, j+1))
    for i in range(len(children_population)):
        print(children_population[i].get_file())
        child = Bacterium(children_population[i].get_file(), names[i])
        new_population.append(child)
    return new_population


def main():
    parser = argparse.ArgumentParser(description='bacterial_genreads_wrapper.py')
    parser.add_argument('-r', type=str, required=True, metavar='reference fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-g', type=int, required=True, metavar='generations', help="Number of generations to run")
    parser.add_argument('-i', type=int, required=True, metavar='initial population', help="Initial population size")
    parser.add_argument('-c', type=float, required=False, metavar='cull pct',
                        help="Percentage of population to cull each cycle (0.5 will keep population relatively stable)",
                        default=0.5)
    args = parser.parse_args()

    (ref_fasta, init_population_size, generations) = (args.r, args.i, args.g)
    cull_percentage = args.c

    population = initialize_population(ref_fasta, init_population_size)

    for i in range(generations):
        new_population = evolve_population(population, i+1)

        new_population = crossover(new_population)

        new_population = cull(new_population, cull_percentage)

        # If all elements get culled, then break the loop
        if not new_population:
            break


if __name__ == '__main__':
    main()
