#!/bin/bash
#SBATCH -J 051125_bHiC_proc_variantCalling_degu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16                  # Increased threads
#SBATCH --mem=196G             # Increased memory
#SBATCH -t 1-00:00:00         # Increased time limit
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
name1=$3
read1=${1%.fastq.gz}
read2=${2%.fastq.gz}

genome="octDeg1" #"hifi_121614_haphic"
genome_file="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/OctDegus1_genome/OctDeg1/fasta/genome.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hic_ont_121624/04.build/scaffolds.fa"
data_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/sequencing-reads-HiC"
output_dir="/tscc/projects/ps-renlab2/jhc103/rosetta-stone-proj/outputs/degu-outputs/outputs-from-1st-alignment/outputs-from-hic-alignment-for-variant-calling"
map_dir="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/output/mapped_alignments"
tmp_dir="/scratch/$USER/job_${SLURM_JOBID}" # Old temp I set "/tscc/lustre/ddn/scratch/jhc103/temp-dir-fast/hic-processing/$name1"
mkdir -p tmp_dir

# Create directories
output_dir_sam="$output_dir/$name1"
map_dir_sam="$map_dir/$name1"
mkdir -p "$output_dir_sam" "$map_dir_sam" "$tmp_dir"

# ---------------------------
# Cleanup Function
# ---------------------------
cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "$tmp_dir/${name1}"_*
    rm -f "$map_dir_sam/${name1}_${genome}"_{r1,r2}.bam
}
trap cleanup EXIT

# ---------------------------
# Pipeline Steps
# ---------------------------

# Step 0: Run Trim Galore if needed
if [ ! -f "$data_dir/${read1}_val_1.fq.gz" ] || [ ! -f "$data_dir/${read2}_val_2.fq.gz" ]; then
    echo "Running Trim Galore!"
    trim_galore --cores 16 --paired "$data_dir/$1" "$data_dir/$2" -o "$data_dir"
else
    echo "Trimmed files already exist."
fi

# Step 1: Build index if missing
if [[ ! -f "${genome_file}.bwt.2bit.64" ]]; then
    echo "Building bwa-mem2 index..."
    bwa-mem2 index "$genome_file"
fi

# Step 2: Initial alignment
if [ ! -f "$map_dir_sam/${name1}_${genome}_initial.bam" ]; then
    echo "Running initial alignment..."
    bwa-mem2 mem -SP5M -t 16 "$genome_file" \
        "$data_dir/${read1}_val_1.fq.gz" \
        "$data_dir/${read2}_val_2.fq.gz" | \
    samtools view -@ 16 -bhS - > "$tmp_dir/${name1}_${genome}_initial.bam"
    
    mv "$tmp_dir/${name1}_${genome}_initial.bam" "$map_dir_sam/"
else
    echo "Initial BAM already exists."
fi

# Step 3: Convert to single-end properly
echo "Converting to single-end..."
samtools view -@ 16 -h "$map_dir_sam/${name1}_${genome}_initial.bam" | \
  awk '
    BEGIN {OFS="\t"}
    /^@/ {print; next}
    {
        # Clear all pairing-related flags and fields
        $2 = and($2, compl(0x1 + 0x40 + 0x80 + 0x100 + 0x200 + 0x400 + 0x800));
        $7 = "*"; $8 = 0; $9 = 0;
        print
    }' | \
  samtools view -@ 16 -b -o "$tmp_dir/${name1}_${genome}.bam"

# Step 4: Sort and mark duplicates
echo "Sorting and marking duplicates..."
samtools sort -@ 16 -m 4G -o "$tmp_dir/${name1}_${genome}.sorted.bam" "$tmp_dir/${name1}_${genome}.bam"
samtools markdup -@ 16 -s -l 9 "$tmp_dir/${name1}_${genome}.sorted.bam" "$tmp_dir/${name1}_${genome}.markdup.bam"

# Step 5: Filter reads
echo "Filtering..."
samtools view -@ 16 -F 3332 -b "$tmp_dir/${name1}_${genome}.markdup.bam" > "$tmp_dir/${name1}_${genome}.filtered.bam"
samtools index -@ 16 "$tmp_dir/${name1}_${genome}.filtered.bam"

# Step 6: Generate coverage
echo "Generating coverage..."
bedtools genomecov -ibam "$tmp_dir/${name1}_${genome}.filtered.bam" -bga > "$output_dir_sam/${name1}_${genome}.coverage.bedgraph"

# Step 7: Variant calling
echo "Calling variants..."
conda activate variant-calling
samtools faidx "$genome_file"

bcftools mpileup --threads 16 -Ou -f "$genome_file" "$tmp_dir/${name1}_${genome}.filtered.bam" | \
bcftools call --threads 16 -mv -o "$output_dir_sam/${name1}_${genome}.vcf"

# Final validation
[ -s "$output_dir_sam/${name1}_${genome}.vcf" ] || { echo "Variant calling failed"; exit 1; }

echo "Pipeline completed successfully for $name1"
