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

all.studies <- c("JP-CRC", "DE-CRC", "FR-CRC", "IT-CRC", "IT-CRC-2", "US-CRC", "AT-CRC")

# Feature tables
feat.species <- read.table("../data/species/species_full_relab.tsv", 
                        sep = "\t", stringsAsFactors = F, 
                        header = T, check.names = F, 
                        row.names = 1, quote ="", fill = F)

feat.og <- read.table("../data/species/mpa.txt", 
                           sep = "\t", stringsAsFactors = F, 
                           header = T, check.names = F, 
                           row.names = 1, quote ="", fill = F)

feat.wirbel <- read.table("../data/species/feat_rel_crc.tsv", 
                         sep = "\t", stringsAsFactors = F, 
                         header = T, check.names = F, 
                         row.names = 1, quote ="", fill = F)

feat <- read.table("../data/species/22.2_mpa_species.txt", 
                   sep = "\t", stringsAsFactors = F, 
                   header = T, check.names = F, 
                   row.names = 1, quote ="", fill = F)
feat.matrix <- prop.table(as.matrix(feat), 2)

# Metadata
meta.all <- read_delim("../data/meta/meta.crc.2.tsv", 
                       delim = "\t", 
                       escape_double = FALSE, 
                       trim_ws = TRUE)

meta.all <- meta.all %>% filter(Study != "CN-CRC")

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
# Save filtered gene table

filtered.species <- '../data/species/filtered.species.NEW.tsv'
write.table(species.filtr, file=filtered.species, quote=FALSE, sep='\t',
            row.names=TRUE, col.names=TRUE)

# #######################
# End of script
# #######################

temp.max.ab = t(sapply(row.names(feat),
                       FUN=function(marker){sapply(unique(meta.all$Group_full),
                                                   FUN=function(group_full, marker){
                                                     max.ab = max(feat[marker, which(meta.all$Study == group_full)])
                                                   },
                                                   marker=marker)}))

### low abundance filter
f.idx = rowSums(temp.max.ab >= ab.cutoff) >= 3 &
  row.names(feat.rel.crc) != '-1'

feat.ab.red = feat.ab.crc[f.idx,]
feat.rel.red = feat.rel.crc[f.idx,]
feat.rar.red = feat.rar.crc[f.idx,]
cat('Retaining', sum(f.idx), 'features after low-abundance filtering...\n')





feat.list <- colnames(feat)
meta.list <- meta.all %>% column_to_rownames(var = "Sample_ID") %>% rownames()
setdiff(meta.list, feat.list)
