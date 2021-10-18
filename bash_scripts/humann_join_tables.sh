#!/bin/bash

# 3. Merging outputs produced for each sample (relab/cpm)
# arg1 = units (relab/cpm), arg2 = functional category (i.e. kegg, level4ec, go, pfam, eggnog)

module load HUManN/3.0.0a4-foss-2020b-Python-3.8.6

start=`date +%s`

humann_join_tables --input ../processed_data/humann/renamed_$2 --verbose --output ../processed_data/humann/merged/humann_${2}_$1.tsv --file_name named_${2}_$1

end=`date +%s`
runtime=$((end-start))
printf "\n Total runtime: $runtime \n"
