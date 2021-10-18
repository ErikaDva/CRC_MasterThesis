#!/bin/bash

# arg1 = threads/CPU

module load metaphlan3/3.0.4-foss-2018a

start=`date +%s`

n=1

files=$(ls /srv/MA/users/sdall16/data/fastq_files/processed_fastq/processed/*_1.fastq.gz | cat )
for f in $files

do
input1=$f
bn=$(basename ${f%_sorted_1.fastq.gz})

metaphlan $input1 --input_type fastq --nproc $1 --tmp_dir /scratch --bowtie2db /space/databases/metaphlan/db_2020-10-26 --bowtie2out /srv/MA/users/edvari19/processed_data/metaphlan/bowtie/${bn}.bowtie2.bz2 -o /srv/MA/users/edvari19/processed_data/metaphlan/output/${bn}_mphlan.txt

echo "File $n: $bn has now been processed"
n=$(expr $n + 1 )
done

end=`date +%s`
runtime=$((end-start))
printf "\n Total runtime: $runtime \n"
