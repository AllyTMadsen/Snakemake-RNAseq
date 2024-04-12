library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')


#helper function to create metadata
timepoint_from_sample <- function(x) {
  find <- str_locate(x, "r") 
  timepoint <- str_sub(x, 1, find -1)
  timepoint <- unique(timepoint)
  return(factor(timepoint))
}

#generate summerized experiment object for deseq2
make_se <- function(csv_path, metafile, selected_times) {
  counts_data <- read.csv(csv_path, row.names = 1)
  
  colData <- read.csv(metafile)
  
  colData <- colData[colData$Timepoint %in% selected_times, ]
  
  selected_samples <- colData$Sample
  counts_data <- counts_data[, selected_samples]

  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts_data)),
    colData = colData,
  )
  
  metadata(se) <- list(model = "model")
  
  return(se)
}


#run deseq2 on selected timepoints 
return_deseq_res <- function(se, design) {
  unique_genes <- rownames(se)[!duplicated(rownames(se))]
  se_unique <- se[unique_genes, ]
  
  dds <- DESeqDataSet(se_unique, design = design)
  dds <- DESeq(dds)
  results <- results(dds, contrast = c("Timepoint", "AD", "P0"))
  results_df <- as.data.frame(
    results
  )
  result_list <- list(results_df = results_df, dds= dds)
  return(result_list)
}


#call the functions to run deseq2 on the filtered_concat_all.csv file:
se <- make_se('results/VERSE/filtered_concat_all.csv', 'results/VERSE/metadata_week4.csv', c('P0', 'AD'))
results <- return_deseq_res(se, ~ Timepoint)

gene_column <- rownames(results$results_df)
results_df <- data.frame(gene_id = gene_column, results$results_df)
results_df <- as_tibble(results_df)

#used these lines to send result to tsv file
#results_df <- results$results_df
#results_df <- results_df %>% rownames_to_column(var = "gene") %>% write_tsv("results/VERSE/DESEQ_week4.tsv")

#show results of DESEQ2
print(head(results_df))

#add gene names to results_df:
gene_names = read.csv("results/txn_map_all.csv", stringsAsFactors = FALSE)
colnames(gene_names) <- c("gene_id", "gene_name")

merged_df <- merge(results_df, gene_names, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
merged_df <- merged_df[, c("gene_id", "gene_name", setdiff(colnames(merged_df), c("gene_id", "gene_name")))]
merged_df <- as_tibble(merged_df)

#show results of DESEQ2 with labels
print(head(merged_df))

#used this lines to send labeled results to tsv file
#write_tsv(merged_df, "results/VERSE/DESEQ_Named_week4.tsv")