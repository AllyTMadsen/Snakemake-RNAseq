##
## WEEK 4 R CODE
##
library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('fgsea')
library('ggplot2')

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
se <- make_se('results/VERSE/filtered_concat_all.csv', 'results/VERSE/metadata_week4.csv', c('P0', 'AD','P4', 'P7'))
results <- return_deseq_res(se, ~ Timepoint)

gene_column <- rownames(results$results_df)
results_df <- data.frame(gene_id = gene_column, results$results_df)
results_df <- as_tibble(results_df)

#show results of DESEQ2
print("these are the DESEQ2 results: ")
print(head(results_df))


#add gene names to results_df:
label_results_df <- function(label_file, deseq_results) {
  gene_names = read.csv(label_file, stringsAsFactors = FALSE)
  colnames(gene_names) <- c("gene_id", "gene_name")

  merged_df <- merge(deseq_results, gene_names, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
  merged_df <- merged_df[, c("gene_id", "gene_name", setdiff(colnames(merged_df), c("gene_id", "gene_name")))]
  merged_df <- as_tibble(merged_df)

  return(merged_df)
}

print("here we are checking the labeled deseq results: ")
labeled_df <- label_results_df("results/txn_map_all.csv", results_df)
print(head(labeled_df))


##
##  WEEK 5 STARTS HERE
##
make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene_map <- read.csv(id2gene_path)
  colnames(id2gene_map) <- c("EnsemblID", "MGI_Symbol")

  #print(head(id2gene_map))
  #print(head(labeled_results))
  
  merged_data <- merge(results_df, id2gene_map, by.x = "gene_id", by.y = "EnsemblID", all.x = TRUE)
  filtered_data <- filter(merged_data, !is.na(log2FoldChange))
  sorted_data <- filtered_data[order(-filtered_data$log2FoldChange), ]
  
  ranked <- setNames(sorted_data$log2FoldChange, sorted_data$MGI_Symbol)
  return(ranked)
}

ranked <- make_ranked_log2fc(merged_df, 'results/txn_map_all.csv')
print("this is the head and tail of the ranked list: ")
print(head(ranked))
print(tail(ranked))

run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  
  rnk_list <- rnk_list[!is.na(rnk_list)]
  rnk_list <- rnk_list[is.finite(rnk_list)]

  pathways <- fgsea::gmtPathways(gmt_file_path)
  fgsea_res <- fgsea(pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  )
  return(tibble(fgsea_res))
}

##run FGSEA
fgsea_results <- run_fgsea('/projectnb/bf528/materials/project_1_rnaseq/m2.all.v2023.2.Mm.symbols.gmt', ranked, 15, 500)
fgsea_results <- fgsea_results %>% arrange(NES)
print("Bottom 10 pathways: ")
print(head(fgsea_results, n = 10))
print("Top 10 pathways: ")
print(tail(fgsea_results, n = 10))

#get FPKM counts -- THIS IS DIFFERENT
dds <- results$dds
dds <- estimateSizeFactors(dds)

##TEST fpkm() ON DDS OBJ
fpm_dds <- fpm(dds, robust = TRUE)
print("these are the FPMK values from the dds object: ")
print(head(fpm_dds))

#Normalized counts below
#norm_counts <- counts(dds, normalized=TRUE)
#show the norm counts, note gene IDs are rownames!
#print("these are the normalized counts from the dds object: ")
#print(head(norm_counts))    


#add gene names to counts from dds object:
label_dds_counts <- function(label_file, dds_counts) {
  gene_names = read.csv(label_file, stringsAsFactors = FALSE)
  colnames(gene_names) <- c("gene_id", "gene_name")

  gene_column <- rownames(dds_counts)
  dds_counts <- data.frame(gene_id = gene_column, dds_counts)
  row.names(dds_counts) <- NULL

  merged_counts <- merge(dds_counts, gene_names, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
  return(merged_counts)
}

print("these are the FPKM LABELED counts: ")
labeled_dds_counts <- label_dds_counts("results/txn_map_all.csv", fpm_dds)
print(head(labeled_dds_counts))


##START CODE FOR R PLOTS

##Read in metadata for WEEK 5 plot
metadata <- read.csv('results/VERSE/metadata_week4.csv')
print(metadata)

##Merge the counts and meta data frames for WEEK 5 plots
labeled_dds_counts_long <- labeled_dds_counts %>% pivot_longer(cols = -c(gene_id, gene_name), names_to = "Sample", values_to = "Count")
merged_counts_meta <- labeled_dds_counts_long %>% left_join(metadata, by = "Sample")
merged_counts_meta$gene_name <- as.factor(merged_counts_meta$gene_name)
print(head(merged_counts_meta))

##GRAPH ONE:
#SARCOMERE-- FPKM vs Timepoint
plot_sarc_FPKM <- function(merged_df, output_file) {
  genes_sarc <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
  filtered_df <- merged_counts_meta %>% filter(gene_name %in% genes_sarc)
  median_counts <- filtered_df %>% group_by(gene_name, Timepoint) %>% summarize(med_count = median(Count, na.rm = TRUE)) %>% mutate(Timepoint = factor(Timepoint, levels = c("P0", "P4", "P7", "AD")))

  p <- ggplot(median_counts, aes(x = Timepoint, y = med_count, group = gene_name, shape = gene_name)) +
    scale_shape_manual(values=1:nlevels(median_counts$gene_name)) +
    geom_line(aes(color = median_counts$gene_name)) +
    geom_point(color = 'black', size = 2, show.legend = FALSE) +   
    labs(x = "Timepoint", y = "FPKM", color = "Gene Name") +
    ggtitle("Sarcomere") +
    theme_light() 

  ggsave(output_file, plot = p, width = 10, height = 6, units = "in", dpi = 300)  
}

plot_sarc_FPKM(merged_counts_meta, "figure1D_sarc.png")


##GRAPH TWO:
#MITOCHONDRIA-- FPKM vs Timepoint
plot_mito_FPKM <- function(merged_df, output_file) {
  genes_mito <- c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh')
  filtered_df <- merged_counts_meta %>% filter(gene_name %in% genes_mito)
  median_counts <- filtered_df %>% group_by(gene_name, Timepoint) %>% summarize(med_count = median(Count, na.rm = TRUE)) %>% mutate(Timepoint = factor(Timepoint, levels = c("P0", "P4", "P7", "AD")))

  p <- ggplot(median_counts, aes(x = Timepoint, y = med_count, group = gene_name, shape = gene_name)) +
    scale_shape_manual(values=1:nlevels(median_counts$gene_name)) +
    geom_line(aes(color = median_counts$gene_name)) +
    geom_point(color = 'black', size = 2, show.legend = FALSE) +   
    labs(x = "Timepoint", y = "FPKM", color = "Gene Name") +
    ggtitle("Mitochondria") +
    theme_light() 

  ggsave(output_file, plot = p, width = 10, height = 6, units = "in", dpi = 300)
}

plot_mito_FPKM(merged_counts_meta, "figure1D_mito.png")

##GRAPH THREE:
#CELL CYCLE-- FPKM vs Timepoint
plot_cycle_FPKM <- function(merged_df, output_file) {
genes_cycle <- c('Cdc7', 'E2f8', 'Cdk7', 'Cdc26', 'Cdc6', 'E2f1', 'Cdc27', 'Bora', 'Cdc45', 'rad51', 'Aurkb', 'Cdc23')
filtered_df <- merged_counts_meta %>% filter(gene_name %in% genes_cycle)
  median_counts <- filtered_df %>% group_by(gene_name, Timepoint) %>% summarize(med_count = median(Count, na.rm = TRUE)) %>% mutate(Timepoint = factor(Timepoint, levels = c("P0", "P4", "P7", "AD")))

  p <- ggplot(median_counts, aes(x = Timepoint, y = med_count, group = gene_name, shape = gene_name)) +
    scale_shape_manual(values=1:nlevels(median_counts$gene_name)) +
    geom_line(aes(color = median_counts$gene_name)) +
    geom_point(color = 'black', size = 2, show.legend = FALSE) +   
    labs(x = "Timepoint", y = "FPKM", color = "Gene Name") +
    ggtitle("Cell Cycle") +
    theme_light() 

  ggsave(output_file, plot = p, width = 10, height = 6, units = "in", dpi = 300)
}

plot_cycle_FPKM(merged_counts_meta, "figure1D_cellcycle.png")