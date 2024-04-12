#!/usr/bin/env python

##TESTING:
#  python filter_cts_mat.py -i results/VERSE/test_verse_concat.csv -o results/VERSE/test_filter_out.csv
# python filter_cts_mat.py -i results/VERSE/verse_concat_all.csv -o results/VERSE/filtered_concat_all_2.csv

#FIX SO THIS INLCUDES ANY SAMPLES THAT HAVE AT LEAST 1 COUNT- NOT SAMPLES THAT HAVE A COUNT FOR E/A SAMPLE


# You can refer to the help manual by `python parse_gtf.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse
import pandas as pd

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()

# we are asking argparse to require a -i or --input flag on the command line when this
# script is invoked. It will store it in the "filenames" attribute of the object. Here
# we are only asking to provide this script one file: the GTF file we are parsing
# We also ask it to require a value for the -o or --output flag, which will specify
# the name of the file we produce

parser.add_argument("-i", "--input", help='The input file specified will be the counts matrix generated from the VERSE files',dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake',dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()


#read data into table
df = pd.read_csv(args.input, index_col = 0)  

#check for non zero value for e/a sample
df_filtered = df[(df != 0).all(axis=1)] #gives incorrect val for whole samples filtered_concat_all.csv , 8007 rows

#write filtered dataframe to output file
df_filtered.to_csv(args.output)