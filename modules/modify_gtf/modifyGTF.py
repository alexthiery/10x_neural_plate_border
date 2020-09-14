import pandas as pd
import re
import sys
import os

        
# make list of chromosomes to edit
lab = ['MT', 'Z', 'W']

# create output file to write to
outfile = open(sys.argv[2], 'a')

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

        outfile.write(line)