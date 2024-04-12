# Snakemake-RNAseq
Snakemake workflow for RNAseq analysis

Deliverables:
- Run FastQC and MultiQC on all 16 original data files (8 paired reads) for quality control
- Using the primary assembly FASTA of the m39 genome and its matching GTF file, generate a STAR index using default parameters
- Use the STAR index to align all 8 samples to the full m39 genome
- Generated a snakemake rule that runs VERSE on each of the BAM files
- Generated a snakemake rule that calls the provided script, concat_df.py, and concatenates the 8 output files from VERSE into a single counts matrix.
- Generated a snakemake rule that calls the provided script, parse_gtf.py, that produces a delimited .txt file containing the correct mapping of gene IDs to their corresponding gene symbol from the GTF annotation file.
- Generated a snakemake rule and a python script named ’filter_matrix.py` that filters the counts matrix to only retain genes where every sample has a non-zero count.
- Perform differential expression comparing p0 and AD timepoints and display the results in a tibble.  
