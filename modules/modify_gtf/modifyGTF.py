import pandas as pd
import re
import sys
import os

# access ensembl IDs for genes from chrom W in galgal6 but not necessarily in galgal5
manual_names = sys.argv[2]

if os.path.isfile(manual_names) == False:
    sys.exit("W genes file missing")

W_genes = pd.read_csv(manual_names)
W_genes.dropna(subset = ["gal5_ensembl_id"], inplace=True)
W_genes = list(W_genes.gal5_ensembl_id)
        
# make list of chromosomes to edit
lab = ['MT', 'Z', 'W']

# create output file to write to
outfile = open(sys.argv[3], 'a')

fn = sys.argv[1]
if os.path.isfile(fn) == False:
    sys.exit("GTF file missing")

with open(fn, 'rt') as gtf:
    for line in gtf:

        # only search lines with gene_id in order to skip header lines
        if 'gene_id' in line:

            # change names for genes from chroms of interest (lab)
            chr = line.split()[0]
            if chr in lab:
                match = [s for s in lab if chr in s]
                line = re.sub('gene_id "', 'gene_id "'+''.join(match)+'-', line)
                if 'gene_name' in line:
                    line = re.sub('gene_name "', 'gene_name "'+''.join(match)+'-', line)

            # change gene names and gene ids which are on W chrom in galgal6 but not galgal5 (W_genes)
            gene_id = line.split('gene_id')[1].split(';')[0].replace('"', '').strip()
            if gene_id in W_genes:
                line = re.sub('gene_id "', 'gene_id "W-', line)
                if 'gene_name' in line:
                    line = re.sub('gene_name "', 'gene_name "W-', line)

        outfile.write(line)