#!/bin/bash
module load HUManN/3.0.0a4-foss-2020b-Python-3.8.6

# Given a table of feature values and a mapping of groups to component features, produce a new table with group values in place of feature values
#  = genefamilies / pathabundance
# $1 for groups = uniref90_rxn, uniref50_rxn

start=`date +%s`

n=1
for f in ../processed_data/humann/output_all_unstratified/*_sorted_1_genefamilies_unstratified.tsv
    do
bn=$(basename ${f%_sorted_1_genefamilies_unstratified.tsv})

humann_regroup_table --input $f --output ../processed_data/humann/regroup_$1/${bn}_genefamilies_$1.tsv --custom ../humann/utility_mapping/map_${1}_uniref90.txt.gz

done

end=`date +%s`
runtime=$((end-start))
printf "\n Total runtime: $runtime \n"
