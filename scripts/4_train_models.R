# ######################################################################
#
##  Model building
#
# ######################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages
library(tidyverse)
library(SIAMCAT)
library(yaml)
library(readr)

# ##############################################################################
# general 

memory.limit(56000)
memory.size(TRUE)

set.seed(2021)

# Use argument if the scrip is run from the command line

#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("The analysis tag needs to be provided! Exiting...\n")
#}
#tag <- args[1]

# If run on Desktop select an appropriate tag for the functional classification

tag="kegg"
tag="eggnog"
tag="level4ec"
tag="pfam"
tag="go"

parameters <- yaml.load_file('../parameters.yaml')
# extract parameters
feat.tag <- parameters$functional.feat
all.studies <- parameters$all.studies
stages <- parameters$stages # to target ADA as well

stages.full <- parameters$stages.full
norm.method <- parameters$model.building$norm.method
n.p <- list(log.n0=ifelse(tag %in% c('species', 'genus'),
                          as.numeric(parameters$model.building$log.n0),
                          as.numeric(parameters$model.building$log.n0.func)),
            sd.min.q=as.numeric(parameters$model.building$sd.min.q),
            n.p=as.numeric(parameters$model.building$n.p),
            norm.margin=as.numeric(parameters$model.building$norm.margin))
num.folds <- as.numeric(parameters$model.building$num.folds)
num.resample <- as.numeric(parameters$model.building$num.resample)
ml.method <- parameters$model.building$ml.method
min.nonzero.coeff <- as.numeric(parameters$model.building$min.nonzero.coeff)
modsel.crit <- list(parameters$model.building$modsel.crit)
perform.fs <- FALSE
param.fs <- list()
if (!tag %in% c('species', 'genus')){
  perform.fs <- TRUE
  param.fs.ss <- 
    list(thres.fs = as.numeric(
      parameters$model.building$feature.selection$cutoff),
      method.fs = parameters$model.building$feature.selection$type,
      direction = "absolute") #added
  param.fs.loso <- 
    list(thres.fs = 3200,
         method.fs = parameters$model.building$feature.selection$type,
         direction = "absolute") #added
}

# ##############################################################################
# Get Data

meta.all <- read_tsv(file = '../data/meta/meta.crc.tsv')

# ##############################################################################
# Model Building for CRC group

for (tag in feat.tag) {
  
for (study in all.studies) {

    meta <- meta.all %>% filter(Study == study)
    
    if (study == "IT-CRC-2") {
      fn.path <- paste0("../data/", tag, "/filtered_", tag, "_IT-CRC.tsv")
      feat <- read.table(fn.path, 
                         sep = "\t", stringsAsFactors = F, 
                         header = T, check.names = F, 
                         row.names = 1, quote ="", fill = F)
      cat(tag, 'feature table loaded for', study, '...\n')
    }
    
    fn.path <- paste0("../data/", tag, "/filtered_", tag, "_", study, ".tsv")
    feat.all <- read.table(fn.path, 
                       sep = "\t", stringsAsFactors = F, 
                       header = T, check.names = F, 
                       row.names = 1, quote ="", fill = F)
    
    cat(tag, 'feature table loaded for', study, '...\n')
      
      # single stage model (CRC)
      stage = "CRC"

      meta.train <- meta %>% filter(Group %in% c("CTR", stage))
      
      feat.train <- feat.all[,meta.train %>% pull(Sample_ID)]
      
      meta.train <- data.frame(meta.train)
      rownames(meta.train) <- meta.train$Sample_ID
      
      siamcat <- siamcat(feat=feat.train, meta=meta.train,
                         label = 'Group', case= stage)
      siamcat <- normalize.features(siamcat, norm.method = norm.method,
                                    norm.param = n.p, feature.type = 'original',
                                    verbose=3)
      siamcat <- create.data.split(siamcat, num.folds = num.folds,
                                   num.resample = num.resample)
      siamcat <- train.model(siamcat,
                             method = ml.method,
                             modsel.crit=modsel.crit,
                             min.nonzero.coeff = min.nonzero.coeff,
                             #perform.fs = perform.fs,
                             perform.fs = T,
                             param.fs = param.fs.ss)
      siamcat <- make.predictions(siamcat)
      siamcat <- evaluate.predictions(siamcat)
      save(siamcat, file=paste0('../models/',tag,'/', study, '_', stage, '_stage_', 
                                ml.method ,'_model.RData'))
      
      cat("Successfully trained a single stage model for", tag, study, '_', stage, '\n')
  }
}

# ##############################################################################
# Model Building for ADA group

ada.studies <- parameters$ada.studies

models <- list()

for (tag in feat.tag) {
  
  for (study in ada.studies) {
    
    meta <- meta.all %>% filter(Study == study)
    
    if (study == "IT-CRC-2") {
      fn.path <- paste0("../data/", tag, "/filtered_", tag, "_IT-CRC.tsv")
      feat <- read.table(fn.path, 
                         sep = "\t", stringsAsFactors = F, 
                         header = T, check.names = F, 
                         row.names = 1, quote ="", fill = F)
      cat(tag, 'feature table loaded for', study, '...\n')
    }
    
    fn.path <- paste0("../data/", tag, "/filtered_", tag, "_", study, ".tsv")
    feat.all <- read.table(fn.path, 
                           sep = "\t", stringsAsFactors = F, 
                           header = T, check.names = F, 
                           row.names = 1, quote ="", fill = F)
    cat(tag, 'feature table loaded for', study, '...\n')
    
    # single stage model (ADA or CRC) 
    stage = "ADA"

    meta.train <- meta %>% filter(Group %in% c("CTR", stage))
    
    feat.train <- feat.all[,meta.train %>% pull(Sample_ID)]
    
    meta.train <- data.frame(meta.train)
    rownames(meta.train) <- meta.train$Sample_ID
    
    siamcat <- siamcat(feat=feat.train, meta=meta.train,
                       label = 'Group', case= stage)
    siamcat <- normalize.features(siamcat, norm.method = norm.method,
                                  norm.param = n.p, feature.type = 'original',
                                  verbose=3)
    siamcat <- create.data.split(siamcat, num.folds = num.folds,
                                 num.resample = num.resample)
    siamcat <- train.model(siamcat,
                           method = ml.method,
                           modsel.crit=modsel.crit,
                           min.nonzero.coeff = min.nonzero.coeff,
                           #perform.fs = perform.fs,
                           perform.fs = T,
                           param.fs = param.fs.ss)
    siamcat <- make.predictions(siamcat)
    siamcat <- evaluate.predictions(siamcat)
    save(siamcat, file=paste0('../models/',tag,'/', study, '_', stage, '_stage_', 
                              ml.method ,'_model.RData'))
    
    cat("Successfully trained a single stage model for", tag, study, '_', stage, '\n')
  }
}

# #######################
# End of script
# #######################
