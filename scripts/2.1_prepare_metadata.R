# ##############################################################################
#
##  METADATA PREPARATION
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages
library(yaml)
library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)

# ##############################################################################
# Get data

parameters <- yaml.load_file('../parameters.yaml')
data.loc <- parameters$setup$data.location
all.studies <- parameters$all.studies

fn.meta <- paste0(data.loc, 'metadata/meta_all.tsv')
meta.all <- read_tsv(fn.meta)

library_size <- read.delim("../data/library_size/library_size_raw.txt", header = T, sep = ':')
rownames(library_size) <- library_size$Sample
library_size_filtered <- read.delim("../data/library_size/library_size_filtered.txt", header = T, sep = ',')
rownames(library_size_filtered) <- library_size_filtered$Sample

# ##############################################################################
# Fix JP-CRC names

JP_table <- read_excel("../data/meta/Yamada_meta.xlsx", sheet = 4, skip = 2)
JP_table <- JP_table %>%
  mutate(BMI = ifelse(BMI == "unknown", NA, BMI))

meta.JP <- JP_table %>% 
  as_tibble() %>%
  mutate(Group = sub("Stage.*", "CRC", Group)) %>%
  mutate(Group = sub("Healthy", "CTR", Group)) %>%
  mutate(Stage = sub("-", NA, Stage)) %>%
  mutate(Localization = sub("-", NA, `Tumor location`)) %>%
  mutate(BMI = as.double(BMI))

# Add rows to meta.all and reserve the ID's of original JP-CRC
JP.ID <- meta.all %>%
  filter(Study == "JP-CRC") %>%
  mutate(External_ID = sub(".*\\.10", "10", External_ID)) %>%
  mutate(External_ID = sub("-[0-9]","", External_ID))

meta.all <- meta.all %>% 
  filter(Study != "JP-CRC") %>%
  add_row(Sample_ID = NA, External_ID = meta.JP$Subject_ID, 
          Age = meta.JP$Age, Gender = meta.JP$Gender, BMI = meta.JP$BMI,
          Country = "JPN", AJCC_stage = meta.JP$Stage, TNM_stage = NA, Localization = meta.JP$Localization,
          Study = ifelse(meta.JP$Subject_ID %in% JP.ID$External_ID, "JP-CRC", "JP-CRC-2"),
          Sampling_rel_to_colonoscopy = "BEFORE", FOBT = NA, Group = meta.JP$Group, Diabetes = NA,
          Vegetarian = NA, Smoking = NA, Library_Size = NA)

# ##############################################################################
# Check if all sample ID's are in meta.all

Yamada_meta <- read.delim("../data/study_names/Yamada_meta.txt", sep = ',')
sample_ID <- as.character(Yamada_meta$Sample_Name)
sample_name <- Yamada_meta$Sample.Name

Yamada <- data.frame(sample_ID, sample_name)
Yamada$index <- 0

rownames(Yamada) <- Yamada$sample_ID

n=1
for (ID in Yamada$sample_ID) {
  ID <- sub("\\.1", "", ID)
  
  index <- grep(ID, meta.all$External_ID)
  
  if (length(index)==1) {
    Yamada[paste0(ID),"indx"] <- index
  } else {
    cat(n, " ",ID,"\n")
    n=n+1
  }
}


# Remove sample that is not there and correct sample_ID
Yamada <- Yamada[-grep("10234", Yamada$sample_ID),]
Yamada[grep("^1$", Yamada$sample_ID),1] <- "00001" 

# Skip samples with ".1" as it is the same individual sequenced twice & a year apart

for (sample in Yamada$sample_ID) {
  if (grepl(".*\\.1", sample) == T) {
    cat("Skipping ",sample,"\n") 
  } else {
    stopifnot(length(grep(sample, meta.all$External_ID))==1)
    index <- grep(paste0("^",sample,"$"), Yamada$sample_ID)
    name <- Yamada$sample_name[index]
    
    meta.all[grep(paste0("^",sample,"$"), meta.all$External_ID),"Sample_ID"] <- name
  }
}

# ##############################################################################
# Add IN-CRC metadata

IN_table <- read_excel("../data/meta/Indian_meta.xlsx", sheet = 1, range = "A2:M62")
meta.IN <- IN_table %>% 
  as_tibble() %>%
  mutate(Group = sub("Healthy", "CTR", Status)) %>%
  mutate(Stage = sub(1, "I", Stage)) %>%
  mutate(Stage = sub(2, "II", Stage)) %>%
  mutate(Stage = sub(3, "III", Stage)) %>%
  mutate(Stage = sub(4, "IV", Stage)) %>%
  mutate(GENDER = ifelse(GENDER %in% c("M", "F"), GENDER, NA)) %>%
  mutate(Study = ifelse(LOCATION == "Bhopal", "IN-CRC", "IN-CRC-2"))

meta.all <- meta.all %>% 
  add_row(Sample_ID = meta.IN$`Sample ID`, External_ID = NA, 
          Age = meta.IN$AGE, Gender = meta.IN$GENDER, BMI = meta.IN$BMI,
          Country = "IND", AJCC_stage = meta.IN$Stage, TNM_stage = meta.IN$`TNM-Staging`, Localization = meta.IN$LOCALIZATION,
          Study = meta.IN$Study, Sampling_rel_to_colonoscopy = "BEFORE", FOBT = NA, Group = meta.IN$Group, Diabetes = NA,
          Vegetarian = NA, Smoking = NA, Library_Size = NA)

# Rename two IN-CRC samples
meta.all[grep("^A11$", meta.all$Sample_ID), "Sample_ID"] <- "A_11"
meta.all[grep("^A15$", meta.all$Sample_ID), "Sample_ID"] <- "A_15"

# ##############################################################################
# Load study names

studies <- c("Thomas", "Wirbel", "Hannigan", "Yu", "Feng", "Zeller", "Vogtmann")

names <- list()
for (study in studies) {
  fn.feat <- paste0("../data/study_names/", study, "_names.txt")
  names[[study]] <- read.delim(fn.feat)
}

## JP-CRC Study
run_accession <- paste0(Yamada_meta$Run, '_of_', Yamada_meta$BioSample)
sample_alias <- Yamada_meta$BioSample
sample_ID <- Yamada_meta$Sample_Name

Yamada <- data.frame(run_accession, sample_alias, sample_ID)
rownames(Yamada) <- sample_ID

# A list of removed samples (n = 28)
Yamada_28 <- Yamada %>% filter(str_detect(sample_ID, "\\.1")) %>% pull(run_accession)

# Removing samples that were sequenced a year after surgery (JP-CRC)

Yamada <- Yamada %>% filter(!str_detect(sample_ID, "\\.1"))
names[["Yamada"]] <- Yamada

# Qin Study
Qin_meta <- read.delim("../data/study_names/Qin_meta.txt", sep = ',')
Qin <- subset(Qin_meta, select = c("Run", "Alias"))
Qin$Alias <- gsub("bgi-", "", Qin$Alias)
Qin$Alias <- gsub("-", "", Qin$Alias)
for (a in Qin$Alias) {
  if (grepl("UC", a)) {
    tmp <- paste0(a, ".0")
    index <- grep(paste0("^", a, "$"), Qin$Alias)
    Qin[index, "Alias"] <- tmp
  }
  if (grepl("CD", a)) {
    tmp <- paste0(a, ".0")
    index <- grep(paste0("^", a, "$"), Qin$Alias)
    Qin[index, "Alias"] <- tmp
  }
}
colnames(Qin) <- c("run_accession", "sample_alias")
names[["Qin"]] <- Qin

# IN-CRC study
Indian_meta1 <- read.delim("../data/study_names/Indian_names.txt", sep = ",")
Indian_meta2 <- read.delim("../data/study_names/Indian_meta2.txt", sep = ",")
run_accession <- c(Indian_meta1$Run, Indian_meta2$Run)
sample_alias <- c(Indian_meta1$Sample.Name, Indian_meta2$Sample.Name)
Indian <- data.frame(run_accession, sample_alias)
names[["Indian"]] <- Indian

# ##############################################################################
# Prepare library size

#if (meta.exists != "meta.all.tsv") {

meta.all <- meta.all %>%
  mutate(Library_Size_raw = as.integer(NA)) %>%
  mutate(Library_Size_filtered = as.integer(NA))

stopifnot(colnames(meta.all)[17] == "Library_Size")
colnames(meta.all)[17] <- "Library_Size_Wirbel"

# Input library concentration
studies <- c("Thomas", "Wirbel", "Hannigan", "Feng", "Zeller", "Vogtmann", "Yamada", "Yu", "Qin", "Indian")
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
  
  alias <- unique(sample_names$sample_alias)
  
  n = 0
  for (a in alias) {
    if (a %in% meta.all$Sample_ID) {
      
      run_acc <- sample_names$run_accession[grep(paste0("^",a,"$"), sample_names$sample_alias)]
      
      for (run in run_acc) {
        index_R <- grep(paste0("^",run,"$"), library_size$Sample)
        index_R.filtered <- grep(paste0("^",run,"_sorted$"), library_size_filtered$Sample)
        if (length(index_R) == 1) {
          
          reads <- library_size$Reads[index_R]
          reads.filtered <- library_size_filtered$Reads[index_R.filtered]
          
          sum <- sum(reads,
                     as.integer(meta.all[grep(paste0("^",a,"$"), meta.all$Sample_ID), "Library_Size_raw"]),
                     na.rm = T)
          
          sum.filtered <- sum(reads.filtered,
                              as.integer(meta.all[grep(paste0("^",a,"$"), meta.all$Sample_ID), "Library_Size_filtered"]),
                              na.rm = T)
          
          meta.all[grep(paste0("^",a,"$"),
                        meta.all$Sample_ID), "Library_Size_raw"] <- sum
          
          meta.all[grep(paste0("^",a,"$"),
                        meta.all$Sample_ID), "Library_Size_filtered"] <- sum.filtered
        }
      }
      n = n+1
    }
    
  }
  
}
#}

# ##############################################################################
# Fix metadata for stage information

library(stringr)
library(rebus)

meta.all.stages <- meta.all %>%
  mutate(Group = sub("MP", "ADA", Group)) %>%
  mutate(Group = sub("NAA", "ADA", Group)) %>%
  mutate(AJCC_stage = sub("Normal", NA, AJCC_stage)) %>%
  mutate(AJCC_stage = sub("Stage_0", 0, AJCC_stage)) %>%
  mutate(AJCC_stage = str_to_upper(AJCC_stage)) %>%
  mutate(AJCC_stage = sub("MP", NA, AJCC_stage)) %>%
  mutate(AJCC_stage = str_replace_all(AJCC_stage, pattern = or("IIA", "IIC"), "II")) %>%
  mutate(AJCC_stage = str_replace_all(AJCC_stage, pattern = "IIIB", "III")) %>%
  mutate(AJCC_stage = str_replace_all(AJCC_stage, pattern = or("IVA", "IVB"), "IV")) %>%
  mutate(Study = str_replace_all(Study, pattern = "JP-CRC-2", "JP-CRC")) %>%
  mutate(Study = str_replace_all(Study, pattern = "IN-CRC-2", "IN-CRC")) %>%
  mutate(Group = ifelse(AJCC_stage %in% "FEW_POLYPS", "ADA", Group)) %>%
  mutate(AJCC_stage = sub("FEW_POLYPS", NA, AJCC_stage)) %>% 
  mutate(Group_full = Group) %>% 
  mutate(Group_full = coalesce(AJCC_stage, Group_full)) %>% 
  mutate(Group_full = str_replace(Group_full, pattern = fixed("0"), "S0")) %>%
  mutate(Group_full = str_replace(Group_full, pattern = fixed("III"), "S3")) %>%
  mutate(Group_full = str_replace(Group_full, pattern = fixed("IV"), "S4")) %>%
  mutate(Group_full = str_replace(Group_full, pattern = fixed("II"), "S2")) %>%
  mutate(Group_full = str_replace(Group_full, pattern = fixed("I"), "S1"))

# Checks
levels(as.factor(meta.all.stages$AJCC_stage))
levels(as.factor(meta.all.stages$Group))
levels(as.factor(meta.all.stages$Group_full))
levels(as.factor(meta.all.stages$Study))

# Subset only CRC studies + exclude HS (healthy that has recovered from CRC) + remove sample with less than 5000 reads (i.e. "MG100183")
meta.all.stages.crc <- meta.all.stages %>% 
  filter(Study != "US-CRC-2") %>%
  filter(str_detect(Study, "CRC")) %>% 
  filter(Group_full != "HS") %>% 
  filter(Sample_ID != "MG100183") %>% 
  filter(Sample_ID != c("SAMEA3136751", "MMRS68910755ST-27-0-0")) # missing stage information

levels(as.factor(meta.all.stages.crc$AJCC_stage))
levels(as.factor(meta.all.stages.crc$Group_full))

# IT-CRC fully missing stage information

  # Checks
levels(as.factor(meta.all.stages.crc$Study))
levels(as.factor(meta.all.stages.crc$Group))

# CRC metadata checks
print(table(meta.all.stages.crc$Group, meta.all.stages.crc$Study))

# Inspecting AT-CRC for missing stage information
#meta.at <- meta.all.stages.crc[meta.all.stages.crc$Study %in% "AT-CRC", ]
#meta.at.crc <- meta.at[meta.at$Group %in% "CRC", ]
#view(meta.at.crc)
# Sample_ID: SAMEA3136751 has no stage for CRC (will need to be removed) & not available on SRA run selector

#meta.all.stages.crc <- meta.all.stages.crc %>% filter(Sample_ID != "SAMEA3136751")

############################### SAVE METADATA ##################################

write_tsv(meta.all, file = '../data/meta/meta.all_17.tsv')
save(meta.all, file = '../data/meta/meta.all_17.RData')
write_tsv(meta.all.stages.crc, file = '../data/meta/meta.crc.2.tsv')
save(meta.all.stages.crc, file = '../data/meta/meta.crc.2.Rdata')

#################################   DONE   #####################################

old_meta <- meta_all_stages_crc[,1] %>% column_to_rownames(var="Sample_ID") %>% rownames()
new_meta <- meta.all.stages.crc[,1] %>% column_to_rownames(var="Sample_ID") %>% rownames()
setdiff(new_meta, old_meta)
