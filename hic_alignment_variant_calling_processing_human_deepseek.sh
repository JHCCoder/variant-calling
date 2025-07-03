#!/bin/bash
#SBATCH -J 031325_bHiC_proc_variantCalling_human_subsampled
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=400G
#SBATCH -t 12:00:00
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user jhc103@ucsd.edu

set -eo pipefail

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

# ---------------------------
# Configuration
# ---------------------------
input_dir="/tscc/projects/ps-renlab2/y2xie/tmp/87.FNIH_DHC_IGM_240925"
genome_file="/tscc/projects/ps-renlab2/ps-renlab/share/bwa_indices/hg38.fa"
output_dir="/tscc/projects/ps-renlab2/jhc103/rosetta-stone-proj/outputs/degu-outputs/outputs-from-1st-alignment/outputs-from-hic-alignment-for-variant-calling"

bam_file_name=$1
name1=${bam_file_name%%.*}
tmp_dir="/tscc/lustre/ddn/scratch/jhc103/temp-dir-fast/human-bam-files-$name1"
echo $tmp_dir
output_dir_sam="$output_dir/$name1"

mkdir -p "$output_dir_sam" "$tmp_dir"
# ---------------------------
# Cleanup Function
# ---------------------------
cleanup() {
    echo "Cleaning up temporary files..."
    rm -f "$tmp_dir/${name1}.subsampled.bam" \
          "$tmp_dir/${name1}.stripped.bam" \
          "$tmp_dir/${name1}.sorted.bam"* \
          "$tmp_dir/${name1}.markdup.bam"* \
          "$tmp_dir/${name1}.filtered.bam"*
}
trap cleanup EXIT

# ---------------------------
# Pipeline Steps
# ---------------------------
echo "PROCESSING $name1"

# New Step 0: Subsample to 800 million reads
echo "Subsampling to target reads..."
target_reads=800000000  # Adjust this number as needed

# Get total reads in original BAM
total_reads=$(samtools view -c "${input_dir}/$bam_file_name")
echo "Total reads in input: $total_reads"

# Calculate subsampling fraction
subsample_fraction=$(echo "scale=5; $target_reads / $total_reads" | bc)
echo "Subsampling fraction: $subsample_fraction"

# Generate random seed for reproducibility (remove $RANDOM for fixed seed)
seed=28  # Replace with a fixed number (e.g., 123) for reproducibility

# Perform subsampling
samtools view -@ 16 -s "$seed$subsample_fraction" -b \
  -o "$tmp_dir/${name1}.subsampled.bam" \
  "${input_dir}/$bam_file_name"

# Verify subsampling
subsampled_reads=$(samtools view -c "$tmp_dir/${name1}.subsampled.bam")
echo "Subsampled reads: $subsampled_reads"
[ "$subsampled_reads" -ge $((target_reads * 9 / 10)) ] || { 
  echo "Subsampling failed: Got $subsampled_reads reads (expected ~$target_reads)"; exit 1; 
}

# Step 1: Strip paired-end information (now using subsampled BAM)
echo "Stripping paired-end information"
samtools view -@ 16 -h "$tmp_dir/${name1}.subsampled.bam" | \
  awk '
    BEGIN {OFS="\t"}
    /^@/ {print; next}
    {
      $2 = and($2, compl(0x1 + 0x40 + 0x80));
      $7 = "*"; $8 = 0; $9 = 0;
      print
    }' | \
  samtools view -@ 16 -b -o "$tmp_dir/${name1}.stripped.bam"
[ -s "$tmp_dir/${name1}.stripped.bam" ] || { echo "Stripping failed"; exit 1; }

# Step 2: Sort and index
echo "Sorting BAM"
samtools sort -@ 16 -m 12G -o "$tmp_dir/${name1}.sorted.bam" "$tmp_dir/${name1}.stripped.bam"
[ -s "$tmp_dir/${name1}.sorted.bam" ] || { echo "Sorting failed"; exit 1; }
samtools index -@ 16 "$tmp_dir/${name1}.sorted.bam"

# Step 3: Mark duplicates
echo "Marking duplicates"
samtools markdup -@ 16 -s "$tmp_dir/${name1}.sorted.bam" "$tmp_dir/${name1}.markdup.bam"
[ -s "$tmp_dir/${name1}.markdup.bam" ] || { echo "Markdup failed"; exit 1; }
samtools index -@ 16 "$tmp_dir/${name1}.markdup.bam"

# Step 4: Filter reads
echo "Filtering reads"
samtools view -@ 16 -F 3332 -b "$tmp_dir/${name1}.markdup.bam" > "$tmp_dir/${name1}.filtered.bam"
[ -s "$tmp_dir/${name1}.filtered.bam" ] || { echo "Filtering failed"; exit 1; }
samtools index -@ 16 "$tmp_dir/${name1}.filtered.bam"

# Step 5: Generate coverage
echo "Generating coverage"
bedtools genomecov -ibam "$tmp_dir/${name1}.filtered.bam" -bga > "$output_dir_sam/${name1}.coverage.bedgraph"

# Step 6: Variant calling
echo "Calling variants"
conda activate variant-calling
ln -sf "$genome_file" "$tmp_dir/hg38.fa"
samtools faidx "$tmp_dir/hg38.fa"

bcftools mpileup --threads 16 -O b -o "$output_dir_sam/${name1}.filtered.bcf" -f "$tmp_dir/hg38.fa" "$tmp_dir/${name1}.filtered.bam"
bcftools call --threads 16 --ploidy 2 -m -v -o "$output_dir_sam/${name1}.filtered.vcf" "$output_dir_sam/${name1}.filtered.bcf"

[ -s "$output_dir_sam/${name1}.filtered.vcf" ] || { echo "Variant calling failed"; exit 1; }

echo "SUCCESS: Pipeline completed for $name1"
