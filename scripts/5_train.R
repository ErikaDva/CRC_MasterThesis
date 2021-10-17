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
cat('Starting model building script\n')
start.time <- proc.time()[1]

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

# overwrite ml method with command line argument
#if (length(args) == 2){
#  ml.method <- args[2]
#}

#feat.all <-  as.matrix(read.table("../data/kegg/feat.kegg.tsv", sep='\t', header=TRUE, 
#                                  stringsAsFactors = FALSE, 
#                                  check.names = FALSE, quote=''))

meta.all <- read_tsv(file = '../data/meta/meta.crc.2.tsv')

meta.all <- read_tsv(file = '../data/meta/meta.crc.2.tsv') %>% filter(Sample_ID != "SAMD00164778")

# ##############################################################################
# Get Data
#feat.path <- paste0('../data/', tag, '/filtered.', tag, '.test.tsv')
#feat.all <- as.matrix(read.table(feat.path, sep='\t', header=TRUE, 
#                                 stringsAsFactors = FALSE, 
#                                 check.names = FALSE, quote=''))

#meta.all <- read_tsv('../data/meta/meta.jp.stages.tsv')
#stopifnot(all(meta$Sample_ID %in% colnames(feat.all)))

# ##############################################################################
# Model Building for CRC group

#models <- list()
#study="AT-CRC"
#ml.method="randomForest"

for (tag in feat.tag) {
  
#  for (study in all.studies) {
    study="JP-CRC"
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
    
    outlier <- "SAMD00164778" # JP-CRC outlier
    feat.all <- feat.all %>% select(-outlier) # an outlier for JP-CRC
    cat(tag, 'feature table loaded for', study, '...\n')
  
    #for (stage in parameters$stages){
      
      stage = "CRC"
      # single stage model (ADA or CRC)
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
      #models[[stage]] <- siamcat
      save(siamcat, file=paste0('../models/',tag,'/', study, '_', stage, '_stage_', 
                                ml.method ,'_model_NEW.RData'))
      
      cat("Successfully trained a single stage model for", tag, study, '_', stage, '\n')
    #}
#  }
}

save(models, file =paste0('../models/',tag,'/all_', tag, "_", 
                            ml.method ,'_models.RData'))

# ##############################################################################
# Model Building for ADA group

ada.studies <- parameters$ada.studies
ada.studies <- c("AT-CRC", "IT-CRC", "FR-CRC")
feat.tag <- c("pfam", "level4ec", "go", "eggnog")

models <- list()
ml.method <- "randomForest"
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
    
    #for (stage in parameters$stages){
    
    stage = "ADA"
    # single stage model (ADA or CRC)
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
    #models[[stage]] <- siamcat
    save(siamcat, file=paste0('../models/',tag,'/', study, '_', stage, '_stage_', 
                              ml.method ,'_model.RData'))
    
    cat("Successfully trained a single stage model for", tag, study, '_', stage, '\n')
    #}
  }
}

# #######################
# End of script
# #######################

sc.obj <- check.associations(
  siamcat,
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  feature.type = "normalized",
  panels = c("fc", "prevalence", "auroc"))


model.interpretation.plot(
  siamcat,
  fn.plot = '../interpretation.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)
