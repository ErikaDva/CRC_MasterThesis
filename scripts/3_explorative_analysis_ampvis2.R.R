# ######################################################################
#
##  4.3 Ordination figures (CCA)
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
# Table preparation

ampvis_kegg <- kegg_relab %>% rownames_to_column(var="OTU") %>% subset(select=-SAMD00164778)
ampvis_eggnog <- eggnog_relab %>% rownames_to_column(var="OTU") %>% subset(select=-SAMD00164778) %>% subset(select=-HP)
ampvis_pfam <- pfam_relab %>% rownames_to_column(var="OTU") %>% subset(select=-SAMD00164778)
ampvis_go <- go_relab %>% rownames_to_column(var="OTU") %>% subset(select=-SAMD00164778)
ampvis_level4ec <- level4ec_relab %>% rownames_to_column(var="OTU") %>% subset(select=-SAMD00164778)
ampvis_species <- species_relab %>% rownames_to_column(var="OTU") %>% subset(select=-SAMD00164778)

# ##############################################################################
# Creating an ampvis2 object

kegg <- amp_load(otutable = ampvis_kegg,
                 metadata = meta.all)

eggnog <- amp_load(otutable = ampvis_eggnog,
                   metadata = meta.all)

go <- amp_load(otutable = ampvis_go,
               metadata = meta.all)

level4ec <- amp_load(otutable = ampvis_level4ec,
                     metadata = meta.all)

pfam <- amp_load(otutable = ampvis_pfam,
                 metadata = meta.all)

species <- amp_load(otutable = ampvis_species,
                    metadata = meta.all)

# ##############################################################################
# General ordination by group (groups)

species.cca <- amp_ordinate(species,
                            type = "CCA",
                            distmeasure = "none",
                            constrain = "Group",
                            sample_color_by = "Group",
                            transform = "log",
                            filter_species = 0.001,
                            #species_nlabels = 5,
                            sample_colorframe = T,
                            species_label_taxonomy = "OTU")

p1 = species.cca + labs(title="MetaPhlAn") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p1
ggsave("../figures/ordination/cca/group/species.cca.pdf", last_plot())
ggsave("../figures/ordination/cca/group/species.cca.png", last_plot())

# kegg

kegg.cca <- amp_ordinate(kegg,
                         type = "CCA",
                         distmeasure = "none",
                         constrain = "Group", 
                         sample_color_by = "Group",
                         transform = "log",
                         filter_species = 0.001,
                         #species_nlabels = 5,
                         sample_colorframe = T,
                         species_label_taxonomy = "OTU")
p2 = kegg.cca + labs(title="KEGG") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p2
ggsave("../figures/ordination/cca/group/kegg.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/group/kegg.cca.png", plot = last_plot())

# eggnog

eggnog.cca <- amp_ordinate(eggnog,
                           type = "CCA",
                           distmeasure = "none",
                           constrain = "Group", 
                           sample_color_by = "Group",
                           transform = "log",
                           filter_species = 0.001,
                           sample_colorframe = T,
                           species_label_taxonomy = "OTU")

p3 = eggnog.cca + labs(title="eggNOG") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p3
ggsave("../figures/ordination/cca/group/eggnog.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/group/eggnog.cca.png", plot = last_plot())

# pfam

pfam.cca <- amp_ordinate(pfam,
                         type = "CCA",
                         distmeasure = "none",
                         constrain = "Group", 
                         sample_color_by = "Group",
                         transform = "log",
                         filter_species = 0.001,
                         #species_nlabels = 5,
                         sample_colorframe = T,
                         species_label_taxonomy = "OTU")

p4 = pfam.cca + labs(title="Pfam") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p4
ggsave("../figures/ordination/cca/group/pfam.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/group/pfam.cca.png", plot = last_plot())

# go

go.cca <- amp_ordinate(go,
                       type = "CCA",
                       distmeasure = "none",
                       constrain = "Group", 
                       sample_color_by = "Group",
                       transform = "log",
                       filter_species = 0.001,
                       #species_nlabels = 5,
                       sample_colorframe = T,
                       species_label_taxonomy = "OTU")

p5 = go.cca + labs(title="GO") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p5
ggsave("../figures/ordination/cca/group/go.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/group/go.cca.png", plot = last_plot())

# enzymes

enzymes.cca <- amp_ordinate(level4ec,
                            type = "CCA",
                            distmeasure = "none",
                            constrain = "Group", 
                            sample_color_by = "Group",
                            transform = "log",
                            filter_species = 0.001,
                            #species_nlabels = 5,
                            sample_colorframe = T,
                            species_label_taxonomy = "OTU",
                            print_caption = T)

p6 = enzymes.cca + labs(title="Enzymes") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p6
ggsave("../figures/ordination/cca/group/enzymes.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/group/enzymes.cca.png", plot = last_plot())

p <- (p1+p2)/(p3+p4)/(p5+p6) + plot_layout(guides="collect")

p

ggsave("../figures/ordination/cca/group/all.cca.pdf", plot = last_plot(), height = 14, width = 10)
ggsave("../figures/ordination/cca/group/all.cca.png", plot = last_plot(), height = 14, width = 10)


# ##############################################################################
# General ordination by group (all studies combined)

species.cca <- amp_ordinate(species,
                            type = "CCA",
                            distmeasure = "none",
                            constrain = "Study",
                            sample_color_by = "Study",
                            transform = "log",
                            filter_species = 0.001,
                            #species_nlabels = 5,
                            sample_colorframe = T,
                            species_label_taxonomy = "OTU")

p1 = species.cca + labs(title="MetaPhlAn") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p1
ggsave("../figures/ordination/cca/species.cca.pdf", last_plot())
ggsave("../figures/ordination/cca/species.cca.png", last_plot())

# kegg

kegg.cca <- amp_ordinate(kegg,
                         type = "CCA",
                         distmeasure = "none",
                         constrain = "Study", 
                         sample_color_by = "Study",
                         transform = "log",
                         filter_species = 0.001,
                         #species_nlabels = 5,
                         sample_colorframe = T,
                         species_label_taxonomy = "OTU")
p2 = kegg.cca + labs(title="KEGG") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p2
ggsave("../figures/ordination/cca/kegg.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/kegg.cca.png", plot = last_plot())

# eggnog

eggnog.cca <- amp_ordinate(eggnog,
                           type = "CCA",
                           distmeasure = "none",
                           constrain = "Study", 
                           sample_color_by = "Study",
                           transform = "log",
                           filter_species = 0.001,
                           sample_colorframe = T,
                           species_label_taxonomy = "OTU")

p3 = eggnog.cca + labs(title="eggNOG") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p3
ggsave("../figures/ordination/cca/eggnog.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/eggnog.cca.png", plot = last_plot())

# pfam

pfam.cca <- amp_ordinate(pfam,
                         type = "CCA",
                         distmeasure = "none",
                         constrain = "Study", 
                         sample_color_by = "Study",
                         transform = "log",
                         filter_species = 0.001,
                         #species_nlabels = 5,
                         sample_colorframe = T,
                         species_label_taxonomy = "OTU")

p4 = pfam.cca + labs(title="Pfam") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p4
ggsave("../figures/ordination/cca/pfam.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/pfam.cca.png", plot = last_plot())

# go

go.cca <- amp_ordinate(go,
                       type = "CCA",
                       distmeasure = "none",
                       constrain = "Study", 
                       sample_color_by = "Study",
                       transform = "log",
                       filter_species = 0.001,
                       #species_nlabels = 5,
                       sample_colorframe = T,
                       species_label_taxonomy = "OTU")

p5 = go.cca + labs(title="GO") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p5
ggsave("../figures/ordination/cca/go.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/go.cca.png", plot = last_plot())

# enzymes

enzymes.cca <- amp_ordinate(level4ec,
                            type = "CCA",
                            distmeasure = "none",
                            constrain = "Study", 
                            sample_color_by = "Study",
                            transform = "log",
                            filter_species = 0.001,
                            #species_nlabels = 5,
                            sample_colorframe = T,
                            species_label_taxonomy = "OTU",
                            print_caption = T)

p6 = enzymes.cca + labs(title="Enzymes") + theme(plot.title = element_text(hjust = 0.5, size = 20))
p6
ggsave("../figures/ordination/cca/enzymes.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/cca/enzymes.cca.png", plot = last_plot())

p <- (p1+p2)/(p3+p4)/(p5+p6) + plot_layout(guides="collect")

p

ggsave("../figures/ordination/cca/all.cca.pdf", plot = last_plot(), height = 14, width = 10)
ggsave("../figures/ordination/cca/all.cca.png", plot = last_plot(), height = 14, width = 10)

# ##############################################################################
# CCA of IN & US studies by groups

goin <- amp_subset_samples(go, Study == "IN-CRC")
gous <- amp_subset_samples(go, Study == "US-CRC")
gojp <- amp_subset_samples(go, Study == "JP-CRC")
goit2 <- amp_subset_samples(go, Study == "IT-CRC-2")
gode <- amp_subset_samples(go, Study == "DE-CRC")
gofr <- amp_subset_samples(go, Study == "FR-CRC")

goin.cca <- amp_ordinate(goin,
                         type = "CCA",
                         distmeasure = "none",
                         constrain = "Group",
                         sample_color_by = "Group",
                         transform = "log",
                         filter_species = 0.001,
                         sample_colorframe = T,
                         species_label_taxonomy = "OTU")

ind = goin.cca + labs(title="Indian cohort") + theme(plot.title = element_text(size = 20))
ind
ggsave("../figures/ordination/goin.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/goin.cca.png", plot = last_plot())

gous.cca <- amp_ordinate(gous,
                         type = "CCA",
                         distmeasure = "none",
                         constrain = "Group",
                         sample_color_by = "Group",
                         transform = "log",
                         filter_species = 0.001,
                         sample_colorframe = T,
                         species_label_taxonomy = "OTU")

us = gous.cca + labs(title="US cohort") + theme(plot.title = element_text(size = 20))
us
ggsave("../figures/ordination/goin.cca.pdf", plot = last_plot())
ggsave("../figures/ordination/goin.cca.png", plot = last_plot())

p <- (ind+us) + plot_layout(guides="collect")

p

ggsave("../figures/ordination/cca/indus.cca.pdf", plot = last_plot(), width = 10)
ggsave("../figures/ordination/cca/indus.cca.png", plot = last_plot(), width = 10)
