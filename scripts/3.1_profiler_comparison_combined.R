# ##############################################################################
#
##  4.1.2 Explorative analysis - Figures
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# Packages
library(readr)
library(readxl)
library(tidyverse)
library(ggsci)
library(patchwork)

library(stringr)
library(tidyr)
library(yaml)

library(matrixStats)

# ##############################################################################
# General 

# This script is for cleaning and filtering feature tables for each individual study

# ##############################################################################
# Get data

parameters <- yaml.load_file('../parameters.yaml')
study.tag <- parameters$all.studies
feat.tag <- parameters$functional.feat

memory.limit(56000)
memory.size(max = T)

# Metadata

meta.all <- read_tsv(file = '../data/meta/meta.crc.2.tsv')
feat.prof <- read_excel("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/data_profiler_performance.xlsx")

feat.prof.full <- t(read_excel("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/data_profiler_performance_complete.xlsx"))
feat.prof.full <- as.data.frame(feat.prof.full)
names(feat.prof.full) <- feat.prof.full[1,]
feat.prof.full <- feat.prof.full[-1,]
feat.prof.full <- feat.prof.full %>% rownames_to_column("profiler")

library(reshape2)
feat.prof.full <- melt(feat.prof.full, id.vars="profiler")


# Number of mapped features 

# Before
feat.prof %>%
        ggplot(aes(fct_reorder(profiler, before), before)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label=round(before)), color="black", size = 5, vjust = -0.5) +
        labs(x="Profiler", y="Number of mapped features", title="Comparison of the number of features mapepd by different profilers")
# After
feat.prof %>%
  ggplot(aes(fct_reorder(profiler, after), after)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=round(after)), color="black", size = 5, vjust = -0.5) +
  labs(x="Profiler", y="Number of mapped features", title="Comparison of the number of features mapepd by different profilers")

# Grouped barplot of number of features

temp <- feat.prof[, 1:3]
temp <- melt(temp, id.vars = "profiler")

temp %>%
  ggplot(aes(fct_reorder(profiler, value), value, fill=variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label=round(value)), color="black", size = 4, vjust = -0.35, position=position_dodge(0.85)) +
  labs(x="Profiler", y="Number of mapped features", title="Comparison of the number of features mapped by different profilers")+
  scale_fill_discrete(name = "Features",
                      labels=c("Unfiltered", "Filtered")) +
  scale_y_continuous(breaks = seq(0,45000,5000), minor_breaks = 5000)+
  theme_bw()+
  theme(legend.position = c(0.05, .95),
        legend.justification = c("left", "top"),
        #legend.box.background = element_rect(colour = "black"),
        #legend.background = element_blank(),
        #legend.background = element_rect(colour = "white"),
        legend.background = element_rect(colour = 'black', fill = 'white'),
        panel.grid.major.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggsave("../figures/general/features.pdf", plot = last_plot())
ggsave("../figures/general/features.png", plot = last_plot())

# UNGROUPED plot

temp <- feat.prof[1:5, c(1,6)]
temp <- temp %>% mutate(ungrouped = ungrouped*100)

f2 <- temp %>%
  ggplot(aes(fct_reorder(profiler, ungrouped), ungrouped, fill = profiler)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label=paste0(round(ungrouped), "%")), color="black", size = 5, vjust = -0.35, position=position_dodge(0.85)) +
  labs(x="Profiler", y="Proportion of ungrouped gene families (%)", title="Comparison of the proportions of ungrouped gene families by different profilers")+
  scale_y_continuous(limits = c(0, 60), breaks = seq(0,65,10), minor_breaks = 10)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x = element_blank())

f2

f2_npg = f2 + scale_color_npg() + scale_fill_npg()
f2_npg # final
ggsave("../figures/general/ungrouped.pdf", last_plot())
ggsave("../figures/general/ungrouped.png", last_plot())
