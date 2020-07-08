#!/usr/bin/env python

import genReads
import argparse


def cull(fasta):
    """
    The purpose of this function will be to cull the bacteria created in the model
    :param fasta:
    :return:
    """
    return fasta


def evolve(reference, read_length, output_prefix):
    """
    The purpose of this function is to evolve the bacteria
    :param reference:
    :param read_length:
    :param output_prefix
    :return:
    """
    args = ['-r', reference, '-R', str(read_length), '-o', output_prefix]
    new_fasta = genReads.main(args)
    return new_fasta


def main():
    parser = argparse.ArgumentParser(description='bacterial_genreads_wrapper.py')
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-R', type=int, required=True, metavar='<int>', help="read length")
    parser.add_argument('-o', type=str, required=True, metavar='<str>', help="output prefix")
    parser.add_argument('-C', type=str, required=True, metavar='<str>', help="Number of cycles to run")
    args = parser.parse_args()

    (ref_fasta, read_len, out_pref) = (args.r, args.R, args.o)
    cycles = args.C

    old_fasta = ref_fasta
    for i in range(cycles):
        new_fasta = evolve(old_fasta, read_len, out_pref)
        cull(new_fasta)
        old_fasta = new_fasta
