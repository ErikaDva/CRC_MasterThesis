#!/bin/bash
module load metaphlan3/3.0.4-foss-2018a

merge_metaphlan_tables.py ../processed_data/metaphlan/output_all/*mphlan.txt > ../processed_data/metaphlan/merged_metaphlan.txt

