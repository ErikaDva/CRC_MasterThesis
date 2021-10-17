# ######################################################################
#
##  Model evaluation 
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
siamcat_in <- siamcat

# US
load("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/models/go/US-CRC_CRC_stage_lasso_model.RData")
siamcat_us <- siamcat

# Generating ROC curves

model.evaluation.plot(siamcat_in, fn.plot = "../figures/models/evaluation/ROC_IN.pdf")
model.evaluation.plot(siamcat_us, fn.plot = "../figures/models/evaluation/ROC_US.pdf")

# Generating model interpretation plots

model.interpretation.plot(
  siamcat_in,
  fn.plot = '../figures/models/evaluation/IN_eval.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)

model.interpretation.plot(
  siamcat_us,
  fn.plot = '../figures/models/evaluation/JP_eval.pdf',
  consens.thres = 0.5,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
)
