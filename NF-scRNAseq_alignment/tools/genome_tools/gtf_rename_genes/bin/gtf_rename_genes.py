#!/usr/local/bin/python
    
import re

def parse_args(args=None):
    Description = "Reformat nf-core/viralrecon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=str, help="Input gtf path.", metavar='')
    parser.add_argument('-o', '--output', type=str, help="Prefix for output gtf file.", metavar=''
    return parser.parse_args(args)


def rename_genes(gtf, outfile):
    gtf =  open(gtf, 'rt')
    
    # create output file to write to
    outfile = open(outfile, 'a')

    new_gtf = []

    for line in gtf:
        if 'gene_id' in line:
    #       Check if any gene_ids in line
            if any(gene_id in line for gene_id in goi.keys()):
    #           Get matching gene name from dict
                gene_name = [gene_name for gene_id, gene_name in goi.items() if gene_id in line]
    #           If gene_name is already in line - remove and replace
                if 'gene_name' in line:
                    line = re.sub('gene_name.*?;', '', line, flags=re.DOTALL)
                line = line.rstrip() + ' gene_name "' + gene_name[0] + '";\n'
        outfile.write(line)


def main(args=None):
    args = parse_args(args)

    goi = {'ENSGALG00000030902': 'SNAI2'}

    rename_genes(gtf = args.input, 'rt', outfile = args.output)

if __name__ == '__main__':
    sys.exit(main())