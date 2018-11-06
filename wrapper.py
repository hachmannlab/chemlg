from libgen import library_generator
import argparse

## Argument parser desription

parser = argparse.ArgumentParser(description='This is a pacakge to generate a combinatorial library of molecules based on the building blocks provided. Please provide the building blocks in the a file in either SMILES form or InChi.')
parser.add_argument('-i',"--config_file", action='store', dest='config_file', default='config.dat', help="provide the config file. Default is config.dat file.")
parser.add_argument('-o',"--output_dir", action='store', dest='output_dir', default='./', help="Path to the output directory.")

## defining arguments
args = parser.parse_args()
BB_file = args.config_file
output_dir = args.output_dir

library_generator(BB_file, output_dir)
