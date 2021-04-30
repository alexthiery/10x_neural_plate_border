#!/usr/bin/env python

import os
import sys
import argparse
import loompy


def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input loom file.", metavar='')
    parser.add_argument('-o', '--output', help="Path to concatenated loom output file.", metavar='')
    return parser.parse_args(args)
    

def check_args(args=None):
    if not os.path.isdir(args.input):
        raise Exception("'--input' path is not a directory.")

def get_file_paths(path):
    files = [path+'/'+f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    return files

    
# Subset loom files from directory and concatenate if necessary
def concatenate_loom(args=None):
    input_files = get_file_paths(os.path.abspath(args.input))
    selected_files = [file for file in input_files if file.lower().endswith('loom')]

    if len(selected_files) < len(input_files):
        if(len(selected_files) == 0):
            raise Exception("'--input' must contain '.loom' files.")
        print("\n'--input' directory contains non .loom files.\nThe following files have been subset for merging:", selected_files)

    output_file = os.path.abspath(args.output)

    print("\nCombining loom files:", selected_files, "\nOutput loom file:", output_file, "\n")

    loompy.combine(selected_files, output_file=output_file, key="Accession")


def main(args=None):
    args = parse_args(args)
    check_args(args)
    concatenate_loom(args)

if __name__ == "__main__":
    sys.exit(main())
