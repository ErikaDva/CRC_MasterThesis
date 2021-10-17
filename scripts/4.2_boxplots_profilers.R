# ######################################################################
#
##  4.2 Boxplots for profiler summary
#
# ######################################################################

# ##############################################################################
# Set-up

setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

library(readr)
library(dplyr)
library(ampvis2)
library(ggplot2)
library(stringr)
library(tibble)
library(patchwork)

memory.limit(56000)
memory.size(T)

# ##############################################################################
# Get data

kegg_relab <- read.table("../data/kegg/filtered_kegg_full.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
go_relab <- read.table("../data/go/filtered_go_full.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
eggnog_relab <- read.table("../data/eggnog/filtered_eggnog_full.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
pfam_relab <- read.table("../data/pfam/filtered_pfam_full.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
level4ec_relab <- read.table("../data/level4ec/filtered_level4ec_full.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
species_relab <- read.table('../data/species/filtered.species.NEW.tsv', sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
meta.all <- read_tsv(file = '../data/meta/meta.crc.2.tsv')

# ##############################################################################
# Study information

study <- meta.all %>% select(c(Sample_ID, Study))

# ##############################################################################
# Counting features

temp <- t(species_relab)
temp <- as.data.frame(temp)
temp$MetaPhlAn3 <- rowSums(temp>0)
temp <- temp %>% rownames_to_column("Sample_ID") %>% select(c(Sample_ID, MetaPhlAn3))

# start
all.counts <- inner_join(study, temp, by = "Sample_ID")

# appending other categories
all.counts <- inner_join(all.counts, temp, by = "Sample_ID")
write.table(all.counts, "../all.counts.tsv")
summary(all.counts$KEGG)

# ##############################################################################
# Generating boxplots

library(tidyverse)
library(ggsci)

study.abb <- c("AT", "CN", "DE", "FR", "IN", "IT", "IT-2", "JP", "US")
study.labels <- c("AT-CRC" = "AT", "CN-CRC" = "CN", "DE-CRC" = "DE", "FR-CRC" = "FR", "IN-CRC" = "IN", "IT-CRC" = "IT", 
                  "IT-CRC-2" = "IT-2", "JP-CRC" = "JP", "US-CRC" = "US")

# species
species_count <- all.counts %>%
  ggplot( aes(x=Study, y=MetaPhlAn3, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Number of features")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank()
  ) +
  ggtitle("MetaPhlAn3") +
  xlab("") +
  scale_x_discrete(labels=study.labels)

species_count

ggsave("../figures/general/species_counts.pdf", plot = last_plot())
ggsave("../figures/general/species_counts.png", plot = last_plot())

# KEGG

kegg_count <- all.counts %>%
  ggplot( aes(x=Study, y=KEGG, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Number of features")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank()
  ) +
  ggtitle("KEGG") +
  xlab("") +
  scale_x_discrete(labels=study.labels)

kegg_count

ggsave("../figures/general/kegg_counts.pdf", plot = last_plot())
ggsave("../figures/general/kegg_counts.png", plot = last_plot())

# eggnog
eggnog_count <- all.counts %>%
  ggplot( aes(x=Study, y=eggNOG, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Number of features")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank()
  ) +
  ggtitle("eggNOG") +
  xlab("") +
  scale_x_discrete(labels=study.labels)

eggnog_count

ggsave("../figures/general/eggnog_counts.pdf", plot = last_plot())
ggsave("../figures/general/eggnog_counts.png", plot = last_plot())

# pfam
pfam_count <- all.counts %>%
  ggplot( aes(x=Study, y=Pfam, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Number of features")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank()
  ) +
  ggtitle("Pfam") +
  xlab("") +
  scale_x_discrete(labels=study.labels)

pfam_count

ggsave("../figures/general/pfam_counts.pdf", plot = last_plot())
ggsave("../figures/general/pfam_counts.png", plot = last_plot())

# GO
go_count <- all.counts %>%
  ggplot( aes(x=Study, y=GO, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Number of features")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank()
  ) +
  ggtitle("GO") +
  xlab("") +
  scale_x_discrete(labels=study.labels)

go_count

ggsave("../figures/general/go_counts.pdf", plot = last_plot())
ggsave("../figures/general/go_counts.png", plot = last_plot())

# enzymes
ec_count <- all.counts %>%
  ggplot( aes(x=Study, y=Enzymes, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Number of features")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank()
  ) +
  ggtitle("Enzymes") +
  xlab("") +
  scale_x_discrete(labels=study.labels)

ec_count

ggsave("../figures/general/ec_counts.pdf", plot = last_plot())
ggsave("../figures/general/ec_counts.png", plot = last_plot())

# combined

library(patchwork)

comb <- (species_count+kegg_count)/(eggnog_count+pfam_count)/(go_count+ec_count)
comb

ggsave("../figures/general/all.counts.comb.pdf", plot = last_plot(), height = 14, width = 10)
ggsave("../figures/general/all.counts.comb.png", plot = last_plot(), height = 14, width = 10)

# ##############################################################################

###### UNGROUPED

# ##############################################################################

kegg_relab <- read.table("../data/kegg/kegg_relab.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
go_relab <- read.table("../data/go/go_relab.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
eggnog_relab <- read.table("../data/eggnog/eggnog_relab.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
pfam_relab <- read.table("../data/pfam/pfam_relab.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)
level4ec_relab <- read.table("../data/level4ec/level4ec_relab.tsv", sep = "\t", stringsAsFactors = F, header = T, check.names = F, row.names = 1, quote ="", fill = F)

# ##############################################################################
# Counting UNGROUPED

temp <- t(level4ec_relab)
temp <- as.data.frame(temp)
temp <- temp %>% rownames_to_column("Sample_ID") %>% select(c(Sample_ID, UNGROUPED)) %>% rename("Enzymes" = "UNGROUPED")

# start
all.ungrouped <- inner_join(study, temp, by = "Sample_ID")

# appending other categories
all.ungrouped <- inner_join(all.ungrouped, temp, by = "Sample_ID")

write.table(all.ungrouped, "../all.ungrouped.percent.tsv")

all.ungrouped$KEGG <- all.ungrouped$KEGG*100
all.ungrouped$eggNOG <- all.ungrouped$eggNOG*100 
all.ungrouped$Pfam <- all.ungrouped$Pfam*100
all.ungrouped$GO <- all.ungrouped$GO*100 
all.ungrouped$Enzymes <- all.ungrouped$Enzymes*100
  
# KEGG

kegg_ungrouped <- all.ungrouped %>%
  ggplot( aes(x=Study, y=KEGG, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Proportion of ungrouped reads (%)")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  ggtitle("KEGG") +
  xlab("") +
  scale_x_discrete(labels=study.labels)+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,110, by = 20))

kegg_ungrouped

ggsave("../figures/general/kegg_ungrouped.pdf", plot = last_plot())
ggsave("../figures/general/kegg_ungrouped.png", plot = last_plot())

# eggnog
eggnog_ungrouped <- all.ungrouped %>%
  ggplot( aes(x=Study, y=eggNOG, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Proportion of ungrouped reads (%)")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())+
  ggtitle("eggNOG") +
  xlab("") +
  scale_x_discrete(labels=study.labels)+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,110, by = 20))

eggnog_ungrouped

ggsave("../figures/general/eggnog_ungrouped.pdf", plot = last_plot())
ggsave("../figures/general/eggnog_ungrouped.png", plot = last_plot())

# pfam
pfam_ungrouped <- all.ungrouped %>%
  ggplot( aes(x=Study, y=Pfam, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Proportion of ungrouped reads (%)")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  ggtitle("Pfam") +
  xlab("") +
  scale_x_discrete(labels=study.labels)+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,110, by = 20))

pfam_ungrouped

ggsave("../figures/general/pfam_ungrouped.pdf", plot = last_plot())
ggsave("../figures/general/pfam_ungrouped.png", plot = last_plot())

# GO
go_ungrouped <- all.ungrouped %>%
  ggplot( aes(x=Study, y=GO, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Proportion of ungrouped reads (%)")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  ggtitle("GO") +
  xlab("") +
  scale_x_discrete(labels=study.labels)+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,110, by = 20))

go_ungrouped

ggsave("../figures/general/go_ungrouped.pdf", plot = last_plot())
ggsave("../figures/general/go_ungroupes.png", plot = last_plot())

# enzymes
ec_ungrouped <- all.ungrouped %>%
  ggplot( aes(x=Study, y=Enzymes, fill=Study)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(y="Proportion of ungrouped reads (%)")+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size = 20),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  ggtitle("Enzymes") +
  xlab("") +
  scale_x_discrete(labels=study.labels)+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,110, by = 20))

ec_ungrouped

ggsave("../figures/general/ec_ungrouped.pdf", plot = last_plot())
ggsave("../figures/general/ec_ungrouped.png", plot = last_plot())

# combined

p_blank <- ggplot()

comb <- (kegg_ungrouped+eggnog_ungrouped)/(pfam_ungrouped+go_ungrouped)/(ec_ungrouped+p_blank)
comb

ggsave("../figures/general/all.ungrouped.comb.pdf", plot = last_plot(), height = 14, width = 10)
ggsave("../figures/general/all.ungrouped.comb.png", plot = last_plot(), height = 14, width = 10)
