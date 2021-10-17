# ##############################################################################
#
##  Feature table pre-processing
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages
library(readr)
library(stringr)
library(tidyr)
library(yaml)
library(tidyverse)
library(yaml)
library(matrixStats)

# ##############################################################################
# General 

# This script is for cleaning and filtering feature tables for each individual study

# ##############################################################################
# Get data

parameters <- yaml.load_file('../parameters.yaml')
study.tag <- parameters$all.studies
feat.tag <- parameters$functional.feat

memory.limit(56000)
memory.size(max = T)

# Metadata

meta.all <- read_tsv(file = '../data/meta/meta.crc.tsv')

# Load each feature table and filter each study individually

sink("../feature.count.filtering.txt", append = T)

for (tag in study.tag) {
  
  cat('##############################################################################\n')
  cat('# Starting cleaning', tag, 'study...\n')
  
  for (f.tag in feat.tag) {
    
    fn.path <- paste0("../data/", f.tag, "/", f.tag, "_relab.tsv")
    feat <- read.table(fn.path, 
                       sep = "\t", stringsAsFactors = F, 
                       header = T, check.names = F, 
                       row.names = 1, quote ="", fill = F)
    cat(tag, '-', f.tag, 'gene table loaded...\n')
    
    meta <- meta.all %>%
      filter(Study == tag)
    
    feat <- as.matrix(feat)
    feat <- feat[, meta %>% pull(Sample_ID)]
    print(rowMeans(feat[c(1,2),]))
    feat <- feat[-c(1,2),]
    
    feat<- feat[rowSums(feat[])>0,] # removing all zeroes for the actual number of features per study
    
    # Filter low abundant features
    
    temp.max.ab <- sapply(unique(meta$Group), FUN=function(s){
      rowMaxs(feat[,meta %>% filter(Group == s) %>% pull(Sample_ID)])
    })

    f.idx = rowSums(temp.max.ab >= 1e-06) >= 1
    feat.filt <- feat[f.idx,]
    cat('Out of', nrow(feat), 'features,', 'retaining', sum(f.idx), f.tag, 'features for', tag, 'study after low-abundance filtering...\n')
    
    # Save filtered feature table
    
  filtered.feat <- paste0('../data/', f.tag, '/filtered_', f.tag, '_', tag, '.tsv')
  write.table(feat.filt, file=filtered.feat, quote=FALSE, sep='\t',
              row.names=TRUE, col.names=TRUE)
  cat(f.tag, 'has now been processed...\n')
  }
  
}

#################################   DONE   #####################################
