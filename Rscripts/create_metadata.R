library(tidyverse)

#timepoint_from_sample function
timepoint_from_sample <- function(x) {
  find <- str_locate(x, "r") 
  timepoint <- str_sub(x, 1, find -1)
  timepoint <- unique(timepoint)
  return(factor(timepoint))
}

#read in the CSV file
data <- read.csv("results/VERSE/filtered_concat_all.csv")

#get col names excluding 'gene' col
sample_names <- colnames(data)[-1]  

#apply the timepoint_from_sample function to each sample name
timepoints <- sapply(sample_names, timepoint_from_sample)

#create a df with sample names and corresponding timepoints
metadata <- data.frame(
  Sample = sample_names,
  Timepoint = as.character(timepoints),
  stringsAsFactors = FALSE
)

#print the metadata dataframe
print(metadata)

#write metadata to a new CSV file for viewing
#write.csv(metadata, "results/VERSE/metadata_week4.csv", row.names = FALSE)