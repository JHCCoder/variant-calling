#!/bin/bash

file_list="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_031325_human_samples.txt"
while IFS=$'\t' read -r bam_file; do
        echo "submitting processing for sample $bam_file"
        sbatch hic_alignment_variant_calling_processing_human_deepseek.sh "$bam_file"
done < $file_list
