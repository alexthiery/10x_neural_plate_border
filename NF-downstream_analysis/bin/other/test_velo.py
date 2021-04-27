#!/usr/bin/env python

import os
import sys
import errno
import argparse
import loompy
import scvelo as scv


def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input loom file.", metavar='')
    # parser.add_argument('-o', '--output', help="Output file.", metavar='')
    parser.add_argument('-si', '--seuratIntersect', help="Pre-filtered H5 seurat object for intersecting with loom", metavar='')
    parser.add_argument('-lo', '--loomOutput', help="Path to concatenated loom output file.", metavar='')
    return parser.parse_args(args)
    

def get_file_paths(path):
    files = [path+'/'+f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    return files

def evaluate_args(args):
    path = os.path.abspath(args.input)
    if os.path.isdir(path):
        path = get_file_paths(path)
        if args.loomOutput is None:
            raise Exception("'--loomOutput must be specified when '--input' is a list of loom files to concatenate.")
            
    elif not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"'--input' path:'{args.input_path}' is not a valid path")
    return(path)


def concatenate_loom(loom_paths, output_path):
    if not isinstance(loom_paths, str):
        loompy.combine(loom_paths, output_file=output_path, key="Accession")
        return(output_path)
    else:
        return(loom_paths)

def prepare_scVelo(loom):
    adata = scv.read(loom, cache=False)



def main(args=None):
    args = parse_args(args)    
    path = evaluate_args(args)
    
    loom_path = concatenate_loom(path, args.loomOutput)

    ds = loompy.connect(loom_path)

    print(ds.ra.keys())
    print(ds.ca.keys())

    temp = ds.view[:, 10:20]

    print(temp.ra.Accession)

    # print(ds.ra[["Accession"]])

    # scv.read("mouseBM.h5ad")


    # prepare_scVelo(loom_path)

if __name__ == "__main__":
    sys.exit(main())