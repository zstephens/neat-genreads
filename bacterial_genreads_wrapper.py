#!/usr/bin/env python

import genReads
import argparse
import random
import pathlib
# from Bio import SeqIO


class Bacterium:
    def __init__(self, reference: str, name: str, analyze: bool):
        self.reference = pathlib.Path(reference)
        self.name = name
        self.file = ""

        if analyze:
            self.analyze()

    def __repr__(self):
        return str(self.name)

    def __str__(self):
        return str(self.name)

    def analyze(self):
        args = ['-r', str(self.reference), '-R', '101', '-o', self.name, '--fa']
        genReads.main(args)
        self.file = pathlib.Path().absolute() / (self.name + "_read1.fa")

    def remove(self):
        pathlib.unlink(self.file)


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


def initialize_population(reference: str, pop_size):
    """
    The purpose of this function is to evolve the bacteria
    :param reference:
    :param pop_size:
    :return population:
    """
    names = []
    for j in range(pop_size):
        names.append("bacterium_0_" + str(j+1))
    population = []
    for i in range(pop_size):
        new_member = Bacterium(reference, names[i], True)
        population.append(new_member)
    return population


def evolve(population: list, generation: int) -> list:
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
        names.append("bacterium" + '_' + generation + '_' + str(j+1))
    for i in range(len(children_population)):
        child = Bacterium(children_population[i], names[i], True)
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

    (ref_fasta, init_population_size, cycles) = (args.r, args.i, args.C)
    cull_percentage = args.c

    population = initialize_population(ref_fasta, init_population_size)

    for i in range(cycles):
        new_population = evolve(population, i + 1)

        new_population = crossover(new_population)

        new_population = cull(new_population, cull_percentage)
        # If all elements get culled, then break and quit
        if not new_population:
            break


if __name__ == '__main__':
    main()
