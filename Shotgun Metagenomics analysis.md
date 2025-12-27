#!/bin/bash

# ==============================================================================
# Shotgun Metagenomics Pipeline: Raw FastQ to Species-Level Abundance Matrix
# Description: This script performs QC, host decontamination, taxonomic classification
#              (Kraken2), and species abundance estimation (Bracken).
# Tools used: fastp, Bowtie2, Kraken2, Bracken, Python/R
# ==============================================================================

# --------------------------
# 1. Configuration & Environment
# --------------------------

# Set working directory
WORKDIR="/path/to/your/project"
mkdir -p "${WORKDIR}/01_clean_data"
mkdir -p "${WORKDIR}/02_host_removed"
mkdir -p "${WORKDIR}/03_taxonomy"
mkdir -p "${WORKDIR}/04_result"

# Database Paths (CRITICAL: These must be pre-built)
# Host Genome Index (e.g., Human hg38) for decontamination
HOST_INDEX="/path/to/bowtie2_index/GRCh38"
# Kraken2 Database (Standard, PlusPF, or custom)
KRAKEN_DB="/path/to/kraken2_db"

# Computational resources
THREADS=16

# Sample list (e.g., Sample_A, Sample_B)
SAMPLES=("Saliva_01" "Saliva_02" "NPC_01" "NPC_02")

# --------------------------
# 2. Analysis Loop per Sample
# --------------------------

echo "[$(date)] Starting Metagenomics Pipeline..."

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"
    
    # Input files
    R1_RAW="${WORKDIR}/raw_data/${SAMPLE}_R1.fq.gz"
    R2_RAW="${WORKDIR}/raw_data/${SAMPLE}_R2.fq.gz"
    
    # Output filenames
    R1_CLEAN="${WORKDIR}/01_clean_data/${SAMPLE}_R1.clean.fq.gz"
    R2_CLEAN="${WORKDIR}/01_clean_data/${SAMPLE}_R2.clean.fq.gz"
    R1_NONHOST="${WORKDIR}/02_host_removed/${SAMPLE}_R1.nonhost.fq.gz"
    R2_NONHOST="${WORKDIR}/02_host_removed/${SAMPLE}_R2.nonhost.fq.gz"
    KRAKEN_REPORT="${WORKDIR}/03_taxonomy/${SAMPLE}.kreport"
    KRAKEN_OUTPUT="${WORKDIR}/03_taxonomy/${SAMPLE}.kraken"
    BRACKEN_OUTPUT="${WORKDIR}/03_taxonomy/${SAMPLE}.bracken"

    # Step 1: Quality Control (fastp)
    echo "--> Step 1: Running fastp (QC)..."
    fastp -i "$R1_RAW" -I "$R2_RAW" \
          -o "$R1_CLEAN" -O "$R2_CLEAN" \
          --detect_adapter_for_pe \
          --thread "$THREADS" \
          -h "${WORKDIR}/01_clean_data/${SAMPLE}_fastp.html" \
          -j "${WORKDIR}/01_clean_data/${SAMPLE}_fastp.json"

    # Step 2: Host Decontamination (Bowtie2)
    # Align to host genome and keep unmapped reads (--un-conc-gz)
    echo "--> Step 2: Removing host reads (Bowtie2)..."
    bowtie2 -p "$THREADS" -x "$HOST_INDEX" \
            -1 "$R1_CLEAN" -2 "$R2_CLEAN" \
            --very-sensitive \
            --un-conc-gz "${WORKDIR}/02_host_removed/${SAMPLE}_%.nonhost.fq.gz" \
            -S /dev/null 2> "${WORKDIR}/02_host_removed/${SAMPLE}_bowtie2.log"
            
    # Rename bowtie2 output to match expected pattern (it outputs _1.nonhost.fq.gz)
    mv "${WORKDIR}/02_host_removed/${SAMPLE}_1.nonhost.fq.gz" "$R1_NONHOST"
    mv "${WORKDIR}/02_host_removed/${SAMPLE}_2.nonhost.fq.gz" "$R2_NONHOST"

    # Step 3: Taxonomic Classification (Kraken2)
    # Classify reads against the database
    echo "--> Step 3: Running Kraken2 classification..."
    kraken2 --db "$KRAKEN_DB" \
            --threads "$THREADS" \
            --paired \
            --report "$KRAKEN_REPORT" \
            --output "$KRAKEN_OUTPUT" \
            --gzip-compressed \
            "$R1_NONHOST" "$R2_NONHOST"

    # Step 4: Species Abundance Estimation (Bracken)
    # Re-estimate abundance specifically at the Species level ('-l S')
    echo "--> Step 4: Running Bracken (Species Level)..."
    bracken -d "$KRAKEN_DB" \
            -i "$KRAKEN_REPORT" \
            -o "$BRACKEN_OUTPUT" \
            -r 150 \
            -l S \
            -t "$THREADS"

    echo "Finished processing ${SAMPLE}"
done

# --------------------------
# 3. Merge Results (Python Script)
# --------------------------

echo "[$(date)] Merging Bracken results into a single matrix..."

# Embedded Python script to merge all sample .bracken files
# Generates a matrix similar to the one you uploaded, but for Species (s__)
cat <<EOF > "${WORKDIR}/merge_bracken.py"
import glob
import pandas as pd
import os

# Find all bracken output files
files = sorted(glob.glob("${WORKDIR}/03_taxonomy/*.bracken"))

merged_df = pd.DataFrame()

for f in files:
    # Get sample name from filename
    sample_name = os.path.basename(f).replace(".bracken", "")
    
    # Read Bracken file
    # Columns: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads
    df = pd.read_csv(f, sep="\t")
    
    # Keep only name (Species) and new_est_reads (Count) or fraction_total_reads (Abundance)
    # Here we use 'new_est_reads' for raw counts (can be TPM normalized later)
    df = df[['name', 'new_est_reads']].set_index('name')
    df.columns = [sample_name]
    
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = merged_df.join(df, how='outer')

# Fill NaNs with 0
merged_df = merged_df.fillna(0)

# Add "s__" prefix to match standard formatting if needed
merged_df.index = "s__" + merged_df.index
merged_df.index.name = "Species"

# Save to CSV
merged_df.to_csv("${WORKDIR}/04_result/species_abundance_matrix.csv", sep=",")
print(f"Successfully merged {len(files)} samples into 04_result/species_abundance_matrix.csv")
EOF

# Run the merge script
python3 "${WORKDIR}/merge_bracken.py"

# Clean up
rm "${WORKDIR}/merge_bracken.py"

echo "[$(date)] Pipeline completed. Output file: 04_result/species_abundance_matrix.csv"