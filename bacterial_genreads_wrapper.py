#!/usr/bin/env python

import genReads
import argparse


def cull(population: list) -> list:
    """
    The purpose of this function will be to cull the bacteria created in the model
    :param population:
    :return:
    """
    return population


def create_fasta(neat_input):
    """
    The purpose of this function will be to convert the output from neat-genreads into something
    that can be an input back into neat genreads
    :param neat_input:
    :return fasta_output:
    """
    fasta_output = None
    return fasta_output


def crossover(population: list) -> list:
    """
    This function will take a list of individuals and do crossover type mixing
    :param population:
    :return new_population:
    """
    new_population = population
    return new_population


def evolve(reference, read_length, output_prefix, pop_size):
    """
    The purpose of this function is to evolve the bacteria
    :param reference:
    :param read_length:
    :param output_prefix:
    :param cycles:
    :return population:
    """
    args = ['-r', reference, '-R', str(read_length), '-o', output_prefix]
    population = []
    for i in range(pop_size):
        new_evolution = genReads.main(args)
        new_member = create_fasta(new_evolution)
        population.append(new_member)
    return population


def main():
    parser = argparse.ArgumentParser(description='bacterial_genreads_wrapper.py')
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-R', type=int, required=True, metavar='<int>', help="read length")
    parser.add_argument('-o', type=str, required=True, metavar='<str>', help="output prefix")
    parser.add_argument('-C', type=int, required=True, metavar='<str>', help="Number of cycles to run")
    parser.add_argument('-p', type=int, required=True, metavar='<str>', help="Population size per cycle")
    args = parser.parse_args()

    (ref_fasta, read_len, out_pref, population_size) = (args.r, args.R, args.o, args.p)
    cycles = args.C

    old_fasta = ref_fasta
    candidate_population = evolve(old_fasta, read_len, out_pref, population_size)
    for i in range(cycles):
        new_population = cull(candidate_population)
        # If all elements get culled, then break and quit
        if not new_population:
            break
        new_population = crossover(new_population)
        candidate_population = evolve(new_population)
