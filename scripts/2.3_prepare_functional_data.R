################# FINAL FUNCTIONAL FEATURE TABLE PREPARATION ###################

################################################################################
## Set-up

setwd("C:/Users/Erika Dvarionaite/iCloudDrive/Desktop/210105_crc_study/src")
#setwd("~/Desktop/210105_crc_study/src")

### Change feature table names, and convert to absolute counts

meta.all <- read_tsv(file = '../data/metadata/meta.all.crc.tsv')

profiler <- c("kegg", "eggnog", "go", "pfam", "level4ec")

################################################################################
## Feature table

for (tag in profiler) {
  
  tables <- list.files("../data/feature_tables/", pattern = tag)
  
  feat <- read.delim(paste0("../data/feature_tables/", tables))
  
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
  
  check <- colnames(feat) %in% meta.all$Sample_ID
  colnames(feat)[check == F]
  
  feat <- feat[, check == T]
  
  stopifnot(all(colnames(feat) %in% meta.all$Sample_ID))
  
  ###############################################################################
  ## Write files

  
  if (tag == "eggnog") {
    write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
                
    eggnog_merged <- feat
  } 
    
  if (tag == "kegg") {
    write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
    
    kegg_merged <- feat
  }
  
  if (tag == "go") {
    write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
      
    go_merged <- feat
  }
  
  if (tag == "level4ec") {
    write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')

    
    level4ec_merged <- feat
  }
  
  if (tag == "pfam") {
    write.table(feat, file=paste0('../data/',tag,'/',tag,"_full.tsv"),
                col.names = T, quote = F, row.names = T, sep = '\t')
   
    pfam_merged <- feat
  }
}

#################################   DONE   #####################################
