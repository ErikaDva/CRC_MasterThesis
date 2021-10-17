# ##############################################################################
#
##  TAXONOMIc TABLE PREPARATION & MERGING OF SAMPLES
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

# ##############################################################################
# Get data

# Meta
meta.all <- read_tsv(file = '../data/meta/meta.crc.2.tsv')

# Feature table
fn.path <- paste0('../data/species/mpa.txt')
feat.tax <- read_delim(fn.path, "\t", escape_double = FALSE, trim_ws = TRUE)

# List of all studies and sample information
load("../data/names.RData")

# Library sizes
library_size <- read.delim("../data/library_size/library_size_raw.txt", header = T, sep = ':')
rownames(library_size) <- library_size$Sample
library_size_filtered <- read.delim("../data/library_size/library_size_filtered.txt", header = T, sep = ',')
rownames(library_size_filtered) <- library_size_filtered$Sample

# ##############################################################################
# Clean feature table

feat <- as.data.frame(feat.tax) %>% 
  rename_all(funs(str_replace(., "_mphlan", ""))) #%>% 
  #select(-NCBI_tax_id)

# ##############################################################################
# Merge samples in taxonomic feature table

studies <- c("Thomas", "Wirbel", "Hannigan", "Feng", "Zeller", "Yu", "Vogtmann", "Yamada", "Qin", "Indian")

for (study in studies) {
  sample_names <- names[[study]]
  
  
  if (study == "Thomas") {
    sample_names$run_accession <- paste0(sample_names$run_accession, '_', sample_names$sample_alias)
    
  }
  if (study == "Hannigan") {
    sample_names$run_accession <- paste0(sample_names$run_accession, '_CRC_Virome')
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
      
      reads <- library_size_filtered[paste0(run), "Reads"]
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

write.table(feat, file=paste0("../data/species/species_full_3.tsv"),
            col.names = T, quote = F, row.names = T, sep = '\t')

write.table(feat.relab, file=paste0("../data/species/species_full_relab_3.tsv"),
            col.names = T, quote = F, row.names = T, sep = '\t')

#################################   DONE   #####################################
