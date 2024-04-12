#!/usr/bin/env python

# You can refer to the help manual by `python concat_df.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse
import pandas as pd
import os

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()

# we are asking argparse to require a -i or --input flag on the command line when this
# script is invoked. It will store it in the "filenames" attribute of the object
# we will be passing it via snakemake, a list of all the outputs of verse so we can
# concatenate them into a single matrix using pandas 

parser.add_argument("-i", "--input", help='A list of the VERSE output filenames provided by snakemake', dest="input", required=True, nargs='+')
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake', dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()


# if you try running this on the command line and supply it a value for -i or --input
# it will show up here, stored in this object

# try just running this script and supply it a random string for the -i and -o argument
# example: `python concat_df.py -i <list of files/strings> -o <list of output file>`
# try testing

#print(args.input)
#print(args.output)


#read in each verse counts file and extract counts
def read_file(file):
    df = pd.read_csv(file, sep='\t')
    return df.set_index('gene')['count']

#instatiate empty lists for e/a sample and sample names
dfs = []
sample_names = []

#add the col header from the file name (sample name)
for file in args.input:
    sample_name = os.path.basename(file).split('.')[0]
    sample_names.append(sample_name)
    df = read_file(file)
    dfs.append(df.rename(sample_name))
#merge the files on the gene col
result = pd.concat(dfs, axis=1)

#write the concatenated dataframe to a CSV file
result.to_csv(args.output)