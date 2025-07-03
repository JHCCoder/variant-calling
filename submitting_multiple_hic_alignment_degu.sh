#!/bin/bash

file_list="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_051125_2DeguSamples.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_042125_degu_sample_1.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_032825_degu_sample_4.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_032725_degu_sample_3.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_032125_degu_sample_2.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_031625_all_degu_ds_samples_except_6991_PFC.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_031425_all_degu_ds_samples_except_6991_PFC.txt" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/file_list_022125_1sample_6991_PFC_ds.txt"
while IFS=$'\t' read -r read_file1 read_file2 sampleid; do
        echo "submitting processing for sample $sampleid"
        sbatch hic_alignment_variant_calling_processing_degu_deepseek.sh "$read_file1" "$read_file2" "$sampleid"
done < $file_list
