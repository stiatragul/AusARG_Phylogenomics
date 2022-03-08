# filename: sample_mapping.R

library(readr)
library(dplyr)
library(tidyr)
library(purrr)


# sample file -------------------------------------------------------------

sample <- readr::read_csv('data/sample_dcf_prg.csv')

# split strings into genus, species, and id
sample$genus <- sapply(sample$sample, function(x) strsplit(x, "_")[[1]][1])
sample$species <- sapply(sample$sample, function(x) strsplit(x, "_")[[1]][2])
sample$id <- sapply(sample$sample, function(x) strsplit(x, "_")[[1]][2])
  

sample$species_name <- paste(sample$genus, '_', sample$species, ':', sep ='')


# aggregate ---------------------------------------------------------------

sample_tab <- aggregate(sample ~ species_name, unique(sample), paste, collapse = ",")


write.table(sample_tab, 
            file = 'data/mapping_file.txt', sep = "", 
            row.names = F, quote = F,
            col.names = F)
