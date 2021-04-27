#!/usr/bin/env python

import os
import sys
import errno
import argparse
import loompy

def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input loom file.", metavar='')
    parser.add_argument('-o', '--output', type=dir_path, help="Output file.", metavar='')
    parser.add_argument('-l', '--loom', help="Full path to concatenated loom output file.", metavar='')
    return parser.parse_args(args)

def check_path(path):
    if not os.path.isdir(path) and not :
        return os.path.abspath(path)
    else:
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid path")

def get_file_paths(path):
    path = dir_path(path)
    files = [path+'/'+f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    return files

def evaluate_input(input_path):
    path = os.path.abspath(input_path)
    if os.path.isdir(path):
        path = get_file_paths(path)
        concatenate_loom(paths)
        path = 

    elif os.path.isfile(path):
        path = dir_path(path)

    else:
        raise argparse.ArgumentTypeError(f"'{input_path}' is not a valid path")

    print(path)




def concatenate_loom(args, input_path):

    

    if len(files) > 1:
        if args.loom is None:
            loompy.combine(files, output_file='test', key="Accession")
        else:
            loom_input = loompy.combine(files, output_file='test', key="Accession")
    elif len(files) == 0:
        exit(0)
    # if os.path.isdir(args.input):
    #     print(os.listdir(file_list))

    # return(loompy.combine(file_list, key="Accession"))

def main(args=None):

    # print(os.getcwd()+'test')
    args = parse_args(args)

    evaluate_input(args.input)

    # if args.input.length > 1:
    #     concatenate_loom(args.input)

    # print(parse_args(args))

if __name__ == "__main__":
    sys.exit(main())