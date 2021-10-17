# ##############################################################################
#
##  GENE TABLE PREPARATION & MERGING OF SAMPLES IN EACH STUDY
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages

library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(yaml)

memory.limit(56000)
memory.size(max = T)

# ##############################################################################
# Get data

# Meta
meta.all <- read_tsv(file = '../data/meta/meta.crc.2.tsv')

# List of all studies and sample information
load("../data/names.RData")

# Library sizes
library_size <- read.delim("../data/library_size/library_size_raw.txt", header = T, sep = ':')
rownames(library_size) <- library_size$Sample
library_size_filtered <- read.delim("../data/library_size/library_size_filtered.txt", header = T, sep = ',')
rownames(library_size_filtered) <- library_size_filtered$Sample

# ##############################################################################
# Merge samples in gene family table

study.tag <- c("at-crc", "jp-crc", "us-crc", "in-crc", "fr-crc", "de-crc", "it-crc")

for (tag in study.tag) {
  
  # Get gene table
  tables <- list.files("../data/genes/", pattern = tag)
  feat <- read.delim(paste0("../data/genes/", tables))
  
  # Clean feture table
  feat <- as.data.frame(feat) %>% 
    rename_all(funs(str_replace(., "_sorted_1_Abundance.RELAB", ""))) %>% 
    rename("gene" = "X..Gene.Family")
  
  cat(tag, 'gene table loaded...\n')
  
  # Associate tags with a study name
  if (tag == "us-crc") {
      study = "Vogtmann"
      cat(tag, 'associated with', study, '\n')
  }
  
  if (tag == "at-crc") {
    study = "Feng"
    cat(tag, 'associated with', study, '\n')
  }
  
  if (tag == "fr-crc") {
    study = "Zeller"
    cat(tag, 'associated with', study, '\n')
  }
  
  if (tag == "jp-crc") {
    study = "Yamada"
    cat(tag, 'associated with', study, '\n')
  }
  
  if (tag == "in-crc") {
    study = "Indian"
    cat(tag, 'associated with', study, '\n')
  }
  
  if (tag == "it-crc") {
    study = "Thomas"
    cat(tag, 'associated with', study, '\n')
  }
  
  if (tag == "de-crc") {
    study = "Wirbel"
    cat(tag, 'associated with', study, '\n')
  }
  
  studies <- c("Thomas", "Wirbel", "Hannigan", "Feng", "Zeller", "Vogtmann", "Yamada", "Yu", "Qin", "Indian")
    
  for (study in studies) {
    sample_names <- names[[study]]
    
    
    if (study == "Thomas") {
      sample_names$run_accession <- paste0(sample_names$run_accession, '_', sample_names$sample_alias)
      
    }
    
    if (study == "Feng") {
      colnames(sample_names) <- c("study_accession", "sample_alias", "run_accession")
    }
    
    if (study == "Yu") {
      sample_names$sample_alias <- sample_names$run_accession
    }
    rownames(sample_names) <- sample_names$run_accession
    
    run_acc <- sample_names$run_accession
    
    n = 0
    for (run in run_acc) {
      if (paste0(run) %in% colnames(feat)) {
        
        reads <- library_size_filtered[paste0(run,'_sorted'), "Reads"]
        sample_alias <- sample_names[run, "sample_alias"]
        
        if (!sample_alias %in% colnames(feat)) {
          paste_run <- paste0(run)
          
          index <- match(paste_run, colnames(feat))
          stopifnot(length(index)==1)
          
          feat[, index] <- round(feat[, index] * reads)
          
          
          colnames(feat)[index] <- sample_alias
          
          
        } else {
          paste_run <- paste0(run)
          
          index <- match(paste_run, colnames(feat))
          index_alias <- match(sample_alias, colnames(feat))
          stopifnot(length(index)==1 & length(index_alias)== 1)
          
          feat[, index] <- round(feat[, index] * reads)
          
          feat[, sample_alias] <- rowSums(feat[, c(index, index_alias)])
          
          
          feat <- feat[, -index]
        }
        
        
        n = n+1
      }
  }
}
  
  rownames(feat) <- make.names(feat[,1])
  
  check <- colnames(feat) %in% meta.all$Sample_ID
  colnames(feat)[check == F]
  
  feat <- feat[, check == T]
  
  stopifnot(all(colnames(feat) %in% meta.all$Sample_ID))

###############################################################################
## Convert to relative abundances

feat.relab <- prop.table(as.matrix(feat), 2)

###############################################################################
## Save files
  
  if (tag == "at-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('AT-CRC has now been processed...\n')
  } 
  
  if (tag == "fr-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('FR-CRC has now been processed...\n')
  }
  
  if (tag == "in-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('IN-CRC has now been processed...\n')
  } 
  
  if (tag == "us-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('US-CRC has now been processed...\n')
  } 
  
  if (tag == "jp-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('JP-CRC has now been processed...\n')
  } 
  
  if (tag == "de-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('DE-CRC has now been processed...\n')
  } 
  
  if (tag == "it-crc") {
    write.table(feat, file=paste0('../data/genes/',tag,"_genes_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    write.table(feat.relab, file=paste0('../data/genes/',tag,"_genes_full_relab.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    cat('IT-CRC has now been processed...\n')
  } 
  
}

#################################   DONE   #####################################
