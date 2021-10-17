# ##############################################################################
#
##  External validation
##  Figure showing how well the classifier can be transferred between different datasets
#
# ##############################################################################

# Setting working directory
setwd("~/Desktop/crc_analysis/scripts") #macbook
setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/crc_analysis/scripts") #windows

# ##############################################################################
# Packages
library("tidyverse")
library("cowplot")
library("pROC")
library("yaml")

# for cowplot
theme_set(theme_bw())

parameters <- yaml.load_file('../parameters.yaml')

# ##############################################################################
# general 

set.seed(2021)

tag = "go" # change the profiler tag
ml.method = "lasso"

fn.pred <- paste0('../files/',tag,'/predictions_', ml.method, '.tsv')
if (!file.exists(fn.pred)){
  stop('The prediciton file is not available. Exiting...\n')
}

# ##############################################################################
# Get Data
meta <- read_tsv('../data/meta/meta.crc.2.tsv')
meta <- meta %>% filter(Group != "ADA")
studies <- meta %>% pull(Study) %>% unique

pred.matrix <- read.table(fn.pred, 
                          sep='\t', check.names = FALSE)
pred.matrix$Sample_ID <- rownames(pred.matrix)
pred.matrix <- as_tibble(pred.matrix)

df.all <- inner_join(meta, pred.matrix, by='Sample_ID')

# ##############################################################################
# Calculate AUROCs
auroc.all <- tibble()

for (study.train in studies){
  for (study.test in studies){
    predictor <- df.all %>%
      filter(Study == study.test) %>% 
      pull(study.train) 
    response <- df.all %>%
      filter(Study == study.test) %>% 
      pull(Group)                  
    temp <- roc(predictor=predictor, response = response, ci=TRUE)
    
    auroc.all <- bind_rows(auroc.all, 
                           tibble(study.train=study.train, 
                                  study.test=study.test,
                                  AUC=c(temp$auc)))
    
  }
}

# ##############################################################################
# AUROC heatmap

col.scheme.heatmap <- parameters$plotting$peformance.cols
plot.levels <- parameters$all.studies

g <- auroc.all %>% 
  mutate(study.test=factor(study.test, levels=plot.levels)) %>% 
  mutate(study.train=factor(study.train, levels=rev(plot.levels))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC, size=CV, color=CV)) +
  geom_tile(color = "red4") + theme_bw() +
  # test in tiles
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  #scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.42, 1)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)+
  labs(title= "GO")
g
# model average
g2 <- auroc.all %>% 
  filter(study.test != study.train) %>% 
  group_by(study.train) %>% 
  summarise(AUROC=mean(AUC)) %>% 
  mutate(study.train=factor(study.train, levels=rev(plot.levels))) %>% 
  ggplot(aes(y=study.train, x=1, fill=AUROC)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUROC, digits=2)"), col='black', size=7)+
  scale_fill_gradient(low = "white", high = "firebrick3", limits=c(0.42, 1), 
                       guide=FALSE) + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Model Average') + ylab('')
g2

pdf(paste0('../figures/models/performance_heatmap_', ml.method, "_", tag,'.pdf'), 
    width = 12, height = 7.5, useDingbats = FALSE)
plot_grid(g, g2, rel_widths = c(5/6, 2/6), align = 'h')
dev.off()

cat('Successfully plotted performance heatmap in',
    proc.time()[1]-start.time, 'second...\n')


# #######################
# End of script
# #######################
