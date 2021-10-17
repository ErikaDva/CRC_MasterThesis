###################### FINAL FEATURE TABLE PREPARATION #########################

################################################################################
## Set-up

setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/210105_crc_study/src")
#setwd("~/Desktop/210105_crc_study/src")

### Change feature table names, and convert to absolute counts

#meta.all.stages <- read_tsv(file = '../data/metadata/meta.all.stages.tsv')
load_meta.all.stages()
#profiler <- c("kegg_relab", "eggnog_relab", "go_relab", "pfam_relab", "level4ec_relab")
profiler <- c("kegg", "eggnog", "go", "pfam", "level4ec")

feat.kegg <- data.table::fread("../data/feature_tables_1/kegg.tsv", quote = "")

feat.kegg <- data.table::fread("../data/feature_tables_no_name_1/kegg.tsv", quote = "")
feat.kegg <- remove_rownames(as_tibble(feat.kegg))
feat.kegg <- as.data.frame(feat.kegg)

feat.eggnog <- data.table::fread("../data/feature_tables_no_name_1/eggnog.tsv", quote = "")
feat.eggnog <- as.data.frame(remove_rownames(as_tibble(feat.eggnog)))
feat.eggnog <- as.data.frame(feat.eggnog)

feat.go <- data.table::fread("../data/feature_tables_no_name_1/go.tsv", quote = "")
feat.go <- as.data.frame(remove_rownames(as_tibble(feat.go)))
feat.go <- as.data.frame(feat.go)

feat.level4ec <- data.table::fread("../data/feature_tables_no_name_1/level4ec.tsv", quote = "")
feat.level4ec <- as.data.frame(remove_rownames(as_tibble(feat.level4ec)))
feat.level4ec <- as.data.frame(feat.level4ec)

feat.pfam <- data.table::fread("../data/feature_tables_no_name_1/pfam.tsv", quote = "")
feat.pfam <- as.data.frame(remove_rownames(as_tibble(feat.pfam)))
feat.pfam <- as.data.frame(feat.pfam)


library(data.table)
################################################################################
## Feature table

for (tag in profiler) {
  
  #tables <- list.files("../data/feature_tables_1/", pattern = tag)
  tables <- list.files("../data/feature_tables_no_name_1/", pattern = tag)
  #tables <- list.files("../data/feature_tables_clean/", pattern = tag)
  
  #tables <- list.files("../data/feature_tables_no_name/", pattern = tag)
  
  #feat <- data.table::fread((paste0("../data/feature_tables_1/", tables)))
  #feat <- as.data.frame(fread(paste0("../data/feature_tables_1/", tables), quote = ""))
  
  ######feat <- paste0("feat.", tag)
  #feat <- read.delim(paste0("../data/feature_tables_1/", tables))
  
  #feat <- read.delim(paste0("../data/feature_tables_clean/", tables))
  feat <- read.delim(paste0("../data/feature_tables_no_name_1/", tables))
  
  if (tag == "kegg") {
    suffix = "_kegg"
  }
  
  if (tag == "eggnog") {
    suffix = "_eggnog"
  }
  
  if (tag == "go") {
    suffix = "_go"
  }

  if (tag == "pfam") {
    suffix = "_pfam"
  }

  if (tag == "level4ec") {
    suffix = "_lvl4ec"
  }
  
  studies <- c("Thomas", "Wirbel", "Hannigan", "Feng", "Zeller", "Vogtmann", "Yamada", "Yu", "Qin", "Indian")
  
  for (study in studies) {
    sample_names <- names[[study]]
    
    
    if (study == "Thomas") {
      sample_names$run_accession <- paste0(sample_names$run_accession, '_', sample_names$sample_alias)
      
    }
    if (study == "Hannigan") {
      sample_names$run_accession <- paste0(sample_names$run_accession, '_CRC_Virome')
    }
    
    if (study == "Feng") {
      colnames(sample_names) <- c("study_accession", "sample_alias", "run_accession")
    }
    
    if (study == "Yu") {
      sample_names$sample_alias <- sample_names$run_accession
    }
    rownames(sample_names) <- sample_names$run_accession
    
    run_acc <- sample_names$run_accession
    
    n = 0
    for (run in run_acc) {
      if (paste0(run,suffix) %in% colnames(feat)) {
        
        reads <- library_size_filtered[paste0(run,'_sorted'), "Reads"]
        sample_alias <- sample_names[run, "sample_alias"]
        
        if (!sample_alias %in% colnames(feat)) {
          paste_run <- paste0(run, suffix)
          
          index <- match(paste_run, colnames(feat))
          stopifnot(length(index)==1)
          
            feat[, index] <- round(feat[, index] * reads)
            
            
            colnames(feat)[index] <- sample_alias
          
          
          
          
        } else {
          paste_run <- paste0(run, suffix)
          
          index <- match(paste_run, colnames(feat))
          index_alias <- match(sample_alias, colnames(feat))
          stopifnot(length(index)==1 & length(index_alias)== 1)
          
            feat[, index] <- round(feat[, index] * reads)
            
            feat[, sample_alias] <- rowSums(feat[, c(index, index_alias)])
          
          
          feat <- feat[, -index]
        }
        
        
        n = n+1
      }
      
    }
    
  }
  
  rownames(feat) <- make.names(feat[,1])
  
  check <- colnames(feat) %in% meta.all.stages$Sample_ID
  colnames(feat)[check == F]
  
  #[1] "clade_name"            "NCBI_tax_id"           "MMPU99077057ST"        "MMPU84450604ST"        "MMPU72854103ST"       
  #[6] "MMPU68403337ST"        "MMPU29365221ST"        "MMRS42570301ST-27-0-0" "MMRS11664448ST-27-0-0" "MMRS71238091ST-27-0-0"
  #[11] "MMRS48639115ST-27-0-0" "MMRS67690541ST-27-0-0" "MMRS92727331ST-27-0-0" "SAMD00114983"       
  
  
  feat <- feat[, check == T]
  
  stopifnot(all(colnames(feat) %in% meta.all.stages$Sample_ID))
  
  
  ###############################################################################
  ## Write files

  
  if (tag == "eggnog") {
    #write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
    #            col.names = T, quote = F, row.names = T, sep = '\t')
    write.table(feat, file=paste0('../data/feature_tables_merged_1/',tag,"_full.tsv"),
                 col.names = T, quote = F, row.names = T, sep = '\t')
                
    eggnog_merged <- feat
  } 
    
  if (tag == "kegg") {
  #  write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
  #              col.names = T, quote = F, row.names = T, sep = '\t')

    write.table(feat, file=paste0('../data/feature_tables_merged_1/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    kegg_merged <- feat
  }
  
  if (tag == "go") {
  #  write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
  #              col.names = T, quote = F, row.names = T, sep = '\t')
  
  write.table(feat, file=paste0('../data/feature_tables_merged_1/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')  
      
    go_merged <- feat
  }
  
  if (tag == "level4ec") {
  #  write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
  #              col.names = T, quote = F, row.names = T, sep = '\t')
    write.table(feat, file=paste0('../data/feature_tables_merged_1/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    level4ec_merged <- feat
  }
  
  if (tag == "pfam") {
  #  write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
  #              col.names = T, quote = F, row.names = T, sep = '\t')
  
    write.table(feat, file=paste0('../data/feature_tables_merged_1/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')  
      
    pfam_merged <- feat
  }
}

#################################   DONE   #####################################

# Variables went from 1673 to 1563

