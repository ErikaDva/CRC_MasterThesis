# ######################################################################
#
##  5. ML figures
#
# ######################################################################

# ##############################################################################
# Set-up

setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

library(readxl)
library(reshape2)
library(ggplot2)
library(yaml)

library(dplyr)
library(ampvis2)

library(stringr)
library(tibble)
library(patchwork)

# ##############################################################################
# Get & prepare data for CRC models

theme_set(theme_bw())

parameters <- yaml.load_file('../parameters.yaml')
col.scheme.heatmap <- parameters$plotting$peformance.cols
all.studies <- parameters$all.studies

# importing result table
crc_models <- read_excel("../crc_models.xlsx")

# converting to long format
crc_models <- melt(crc_models, id.vars="Study")
crc_models <- round(value, 2)
crc_models$AUC <- crc_models$value

study_order <- c("AT-CRC", "CN-CRC","DE-CRC", "FR-CRC", "IN-CRC", "IT-CRC", "IT-CRC-2", "JP-CRC", "US-CRC", "Mean")
profiler_order <- c("MetaPhlAn", "KEGG", "eggNOG", "Pfam", "GO", "Enzymes")

crc_models2 <- crc_models %>% filter(Study != "Mean") %>% filter(variable != "Study average")
crc_heat <- ggplot(crc_models2, aes(x = factor(variable, level = profiler_order), y = factor(Study, level = rev(study_order)), fill = AUC)) +
            geom_tile(color = "black") +
            geom_text(aes(label = round(value,2)), color = "black", size = 4) +
            #scale_fill_gradient(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
            scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.5, 1)) +
            coord_fixed() + 
            labs(title= "Lasso models for CRC prediction") +
            theme(axis.text = element_text(size=10),
                  plot.title = element_text(hjust = 0.5, size = 16, face="bold"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())

crc_heat
p1 = crc_heat

p1


crc_study <- crc_models %>% filter(variable == "Study average") %>% filter(Study != "Mean")

crc_study_heat <- ggplot(crc_study, aes(x = "Mean", y = factor(Study, level = rev(study_order)), fill = AUC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,2)), color = "black", size = 4) +
  #scale_fill_gradient(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.5, 1)) +
  coord_fixed() + 
  #labs(title= "Study \n average") +
  theme(#axis.text = element_text(size=10),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

crc_study_heat

p2 = crc_study_heat

crc_heatmap <- p1 + p2 + p3 + plot_layout(ncol=2)
crc_heatmap

ggsave("../figures/models/all.crc.pdf", plot = last_plot(), dpi = 300, height = 10, width = 8)
ggsave("../figures/models/all.crc.png", plot = last_plot(), dpi = 300, height = 10, width = 8)

crc_profilers <- crc_models %>% filter(Study == "Mean") %>% filter(variable != "Study average")

crc_profilers_heat <- ggplot(crc_profilers, aes(x = factor(variable, level = profiler_order), y = "Mean", fill = AUC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,2)), color = "black", size = 4) +
  #scale_fill_gradient(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.5, 1)) +
  coord_fixed() + 
  #labs(title= "Study \n average") +
  theme(#axis.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),)

crc_profilers_heat

p3 = crc_profilers_heat

# ##############################################################################
# Get & prepare data for ADA models

# importing result table
ada_models <- read_excel("../ada_models.xlsx")

# converting to long format
ada_models <- melt(ada_models, id.vars="Study")
ada_models <- round(value, 2)
ada_models$AUC <- ada_models$value

ada_study_order <- c("AT-CRC", "FR-CRC", "IT-CRC", "JP-CRC")
profiler_order <- c("MetaPhlAn", "KEGG", "eggNOG", "Pfam", "GO", "Enzymes")

ada_models2 <- ada_models %>% filter(Study != "Mean") %>% filter(variable != "Study average")
ada_heat <- ggplot(ada_models2, aes(x = factor(variable, level = profiler_order), y = factor(Study, level = rev(study_order)), fill = AUC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,2)), color = "black", size = 4) +
  #scale_fill_gradient(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.45, 1)) +
  coord_fixed() + 
  labs(title= "Lasso models for ADA prediction") +
  theme(axis.text = element_text(size=10),
        plot.title = element_text(hjust = 0.5, size = 16, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ada_heat
a1 = ada_heat

a1


ada_study <- ada_models %>% filter(variable == "Study average") %>% filter(Study != "Mean")

ada_study_heat <- ggplot(ada_study, aes(x = "Mean", y = factor(Study, level = rev(study_order)), fill = AUC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,2)), color = "black", size = 4) +
  #scale_fill_gradient(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.5, 1)) +
  coord_fixed() + 
  #labs(title= "Study \n average") +
  theme(#axis.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())

ada_study_heat

a2 = ada_study_heat

ada_heatmap <- a1 + a2 + a3 + plot_layout(ncol=2)
ada_heatmap

ggsave("../figures/models/all.ada.pdf", plot = last_plot(), dpi = 300, width = 8)
ggsave("../figures/models/all.ada.png", plot = last_plot(), dpi = 300, width = 8)

ada_profilers <- ada_models %>% filter(Study == "Mean") %>% filter(variable != "Study average")

ada_profilers_heat <- ggplot(ada_profilers, aes(x = factor(variable, level = profiler_order), y = "Mean", fill = AUC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,2)), color = "black", size = 4) +
  #scale_fill_gradient(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.5, 1)) +
  coord_fixed() + 
  #labs(title= "Study \n average") +
  theme(#axis.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),)

ada_profilers_heat

a3 = ada_profilers_heat
