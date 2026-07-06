#!/usr/bin/env python

# CLI for making a PDB from an IDR sequence

import os
import argparse

from dodo.build import pdb_from_sequence as build_pdb_from_sequence

def pdb_from_sequence():
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description='Generate a PDB of an IDR from the sequence alone.')
    # add args
    parser.add_argument('sequence', help='Amino acid sequence of the IDR to make a structure of.')
    parser.add_argument('-o', '--out_path', help='Path to the output file.', required=True)
    parser.add_argument('-m', '--mode', default='predicted', help='Mode to use for generating the structure.')
    parser.add_argument('-n', '--num_models', default=1, help='Number of models to generate.', type=int)
    parser.add_argument('-c', '--no_CONECT_lines', action='store_false', default=True, help='Include CONECT lines in the output PDB.')
    parser.add_argument('-apr', '--attempts_per_residue', default=1000, help='Number of attempts to make per residue.', type=int)
    parser.add_argument('-api', '--attempts_per_IDR', default=50, help='Number of attempts to make per coordinate.', type=int)
    parser.add_argument('--use_pulchra', action='store_true', default=False,
        help='Run PULCHRA on the generated CA trace to rebuild backbone and side-chain atoms.')
    parser.add_argument('--pulchra_executable', default='pulchra',
        help='PULCHRA executable name or full path (default: pulchra).')

    # parser args
    args = parser.parse_args()

    # build pdb
    build_pdb_from_sequence(args.sequence, out_path=args.out_path, mode=args.mode, 
        attempts_per_res=args.attempts_per_residue, CONECT_lines=args.no_CONECT_lines,
        attempts_per_idr=args.attempts_per_IDR, num_models=args.num_models,
        use_pulchra=args.use_pulchra, pulchra_executable=args.pulchra_executable)
