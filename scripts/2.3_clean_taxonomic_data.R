# ##############################################################################
#
##  Species feature table pre-processing for all studies
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages
library("readr")
library("stringr")
library("tidyr")
library("yaml")
library("dplyr")
library(tibble)

# ##############################################################################
# Get data

parameters <- yaml.load_file('../parameters.yaml')
all.studies <- parameters$all.studies

# Feature tables
feat.species <- read.table("../data/species/species_full_relab.tsv", 
                        sep = "\t", stringsAsFactors = F, 
                        header = T, check.names = F, 
                        row.names = 1, quote ="", fill = F)

# Metadata
meta.all <- read_delim("../data/meta/meta.crc.tsv", 
                       delim = "\t", 
                       escape_double = FALSE, 
                       trim_ws = TRUE)

# ##############################################################################
# Remove UNKNOWN reads

#feat <- feat.path[c(-1,-2),]
feat <- feat.species
# ##############################################################################
# Filter low abundant features

library(matrixStats)

### Gene families

#feat.matrix <- as.matrix(feat)

temp.max.ab <- sapply(unique(meta.all$Group), FUN=function(s){
  rowMaxs(feat.matrix[,meta.all %>% filter(Group == s) %>% pull(Sample_ID)])
})

f.idx = rowSums(temp.max.ab >= 1e-03) >= 1
species.filtr <- feat.matrix[f.idx,]
cat('Out of', nrow(feat), 'features,', 'retaining', sum(f.idx), 'species after low-abundance filtering...\n')

# ##############################################################################
# Save filtered species table

filtered.species <- '../data/species/filtered.species.NEW.tsv'
write.table(species.filtr, file=filtered.species, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

# #######################
# End of script
# #######################
