#!/bin/bash
module load HUManN/3.0.0a4-foss-2020b-Python-3.8.6

# 1. Run HUMAnN3 with MetaPhlAn output files
start=`date +%s`
n=1

echo "missing_files" > missing_files60.txt

#files=$(ls /srv/MA/users/sdall16/data/fastq_files/processed_fastq/processed/*_1.fastq.gz | cat | head -n 20)
cat ../humann/mpa60x.txt | while read f
#for f in $files
    do
    bn=$(basename ${f%_sorted_1.fastq.gz})
if test -f "/srv/MA/users/edvari19/processed_data/metaphlan/output/${bn}_mphlan.txt"; then

        humann -i $f -o ../processed_data/humann/output --threads 5 --taxonomic-profile ../processed_data/metaphlan/output/${bn}_mphlan.txt \
	 --protein-database ../humann/uniref --bowtie-options "-p 5" --resume --verbose --remove-temp-output --o-log ../processed_data/humann/logs/${bn}.log
printf "\n File $n: $bn has now been processed "
n=$(expr $n + 1)
else 

echo "$bn tax-file does not exist"
echo "$f" >> missing_files.txt

fi
done


end=`date +%s`
runtime=$((end-start))

printf "\n Total runtime: $runtime "

