#!/usr/bin/env python

# You can refer to the help manual by `python parse_gtf.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse
import re

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()

# we are asking argparse to require a -i or --input flag on the command line when this
# script is invoked. It will store it in the "filenames" attribute of the object. Here
# we are only asking to provide this script one file: the GTF file we are parsing
# We also ask it to require a value for the -o or --output flag, which will specify
# the name of the file we produce

parser.add_argument("-i", "--input", help='The input file specified will be the GTF file provided by snakemake',dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake',dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()

# if you try running this on the command line and supply it a value for -i or --input
# it will show up here, stored in this object

# try just running this script and supply it a random string for the -i and -o argument
# example: `python parse_gtf.py -i <your_string> -o <another_string>`
# You have now passed command line arguments directly into this python script

# replace the code that comes after this with the code necessary to parse the GTF

#find the gene name and IDs from the specified col
def extract_gene_info(column):
    gene_id_match = re.search(r'gene_id "(.+?)"', column)
    gene_name_match = re.search(r'gene_name "(.+?)"', column)
    if gene_id_match and gene_name_match:
        gene_id = gene_id_match.group(1)
        gene_name = gene_name_match.group(1)
        return gene_id, gene_name
    else:
        return None, None

#open GTF file
with open(args.input, "r") as gtf_file:
    #new file to write the gene ID and gene name pairs
    with open(args.output, "w") as output_file:
        #loop through each line in the GTF file
        for line in gtf_file:
            #skip header lines
            if line.startswith('##'):
                continue
            #split the line into cols
            columns = line.strip().split('\t')
            #if the feature type is "gene"
            if columns[2] == "gene":
                #extract gene ID and gene name from the ninth column
                gene_id, gene_name = extract_gene_info(columns[8])
                if gene_id and gene_name:
                    #write gene ID and gene name pair to the output file
                    output_file.write(f"{gene_id},{gene_name}\n")