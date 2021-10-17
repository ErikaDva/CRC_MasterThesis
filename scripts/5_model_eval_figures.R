# ######################################################################
#
##  4.2.3 Model evaluation 
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

# Load models

# INDIAN
load("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/models/go/IN-CRC_CRC_stage_lasso_model.RData")

# GERMAN
load("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/models/go/US-CRC_CRC_stage_lasso_model.RData")

load("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/models/KEGG/IT-CRC-2_CRC_stage_lasso_model.RData")

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

check.confounders(siamcat)

sc.obj <-  evaluate.predictions(siamcat)
model.evaluation.plot(sc.obj, fn.plot = "../figures/models/evaluation/ROC_JP.pdf")
model.evaluation.plot(sc.obj, fn.plot = "../figures/models/evaluation/ROC_IN.pdf")
model.evaluation.plot(siamcat, fn.plot = "../figures/models/evaluation/ROC_FR_egg_ADA2.pdf")

model.interpretation.plot(
  sc.obj,
  fn.plot = '../figures/models/evaluation/IN_eval.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)

model.interpretation.plot(
  sc.obj,
  fn.plot = '../figures/models/evaluation/JP_eval.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)

model.interpretation.plot(
  sc.obj,
  fn.plot = '../figures/models/evaluation/IT2_eval.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)
