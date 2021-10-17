# ##############################################################################
#
##  Extract information from models & make predictions
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages
library(tidyverse)
library(SIAMCAT)
library(yaml)

# ##############################################################################
# general 

memory.limit(56000)
memory.size(TRUE)

set.seed(2021)

# ##############################################################################
# Get data

meta <- read_tsv(file = '../data/meta/meta.crc.tsv')

# ##############################################################################
# Parameters

parameters <- yaml.load_file('../parameters.yaml')

all.studies <- parameters$all.studies
ada.studies <- parameters$ada.studies

ml.method="lasso"

# ##############################################################################
# Import models for assesment (not used in making predictions)

# ADA
models <- list()
for (study in ada.studies){
  tag="eggnog" # change the feature tag
  stage="ADA"
  # import all models
  mdl.path <- paste0("../models/", tag, "/", study, "_", stage, "_stage_", ml.method, "_model.RData")
  load(mdl.path)
  cat(tag, stage, 'model loaded for', study, '...\n')
  models[[study]] <- siamcat
}

# CRC
models <- list()
for (study in all.studies){
  tag="eggnog" # change the feature tag
  stage = "CRC"
  # import all models
  mdl.path <- paste0("../models/", tag, "/", study, "_", stage, "_stage_", ml.method, "_model2.RData")
  load(mdl.path)
  cat(tag, stage, 'model loaded for', study, '...\n')
  models[[study]] <- siamcat
}

# ##############################################################################
# Make predictions

tag="species"  # change the feature tag
ml.method = "lasso"

fn.path <- paste0("../data/", tag, "/filtered_", tag, "_full.tsv")
feat.all <- read.table(fn.path, 
                       sep = "\t", stringsAsFactors = F, 
                       header = T, check.names = F, 
                       row.names = 1, quote ="", fill = F)
cat(tag, 'feature table loaded...\n')

if (tag == "species") {
  fn.path <- paste0("../data/", tag, "/filtered.", tag, ".tsv")
  feat.all <- read.table(fn.path, 
                         sep = "\t", stringsAsFactors = F, 
                         header = T, check.names = F, 
                         row.names = 1, quote ="", fill = F)
  cat(tag, 'feature table loaded...\n') 
}

models <- list()

for (study in parameters$all.studies){
  stage = "CRC"
  # import all models
  mdl.path <- paste0("../models/", tag, "/", study, "_", stage, "_stage_", ml.method, "_model_NEW.RData")
  load(mdl.path)
  cat(tag, stage, 'model loaded for', study, '...\n')
  models[[study]] <- siamcat
}

pred.matrix <- matrix(NA, nrow=nrow(meta), 
                      ncol=length(all.studies), 
                      dimnames = list(meta$Sample_ID, 
                                      all.studies))

for (study in parameters$all.studies){
  
  # load model
  siamcat <- models[[study]]
  temp <- rowMeans(pred_matrix(siamcat))
  pred.matrix[names(temp), study] <- temp
  
  # predict other studies
  for (study_ext in setdiff(parameters$all.studies, study)){
    
    meta.test <- meta %>%
      filter(Study == study_ext) %>%
      filter(Group != "ADA")
    
    feat.test <- feat.all[,meta.test %>% pull(Sample_ID)]
    
    meta.test <- data.frame(meta.test)
    rownames(meta.test) <- meta.test$Sample_ID
    
    siamcat.test <- siamcat(feat=feat.test,
                            meta = meta.test, #added
                            label = "Group", case = "CRC") #added
    
    siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
    
    temp <- rowMeans(pred_matrix(siamcat.test))
    pred.matrix[names(temp), study] <- temp
    
    cat("Predictions using", study, "model on", study_ext, 'data --> DONE\n')

  }
  write.table(pred.matrix, file=paste0('../files/',tag,'/predictions_', 
                                         ml.method, '.tsv'), 
              quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
  cat("Predictions for ", tag, 'models complete...\n')
}
