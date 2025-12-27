#!/bin/bash

# ==============================================================================
# RNA-seq Analysis Pipeline: Raw FastQ to TPM Matrix & Log2FC
# Description: This script performs quality control, alignment, quantification,
#              TPM normalization, and basic Log2FC calculation.
# Tools used: fastp, STAR, featureCounts (Subread), R
# ==============================================================================

# --------------------------
# 1. Configuration & Environment
# --------------------------

# Set working directory (Please modify accordingly)
WORKDIR="/path/to/your/project"
mkdir -p "${WORKDIR}/01_clean_data"
mkdir -p "${WORKDIR}/02_alignment"
mkdir -p "${WORKDIR}/03_quantification"
mkdir -p "${WORKDIR}/04_result"

# Reference genome and annotation paths
GENOME_INDEX="/path/to/STAR_index"
GTF_FILE="/path/to/genome_annotation.gtf"

# Computational resources
THREADS=16

# Sample list (Corresponds to file prefixes, e.g., Sample_1_R1.fq.gz)
# NOTE: Ensure your sample names allow distinguishing groups (e.g., "Control", "Treat")
SAMPLES=("Control_1" "Control_2" "Treat_1" "Treat_2")

# --------------------------
# 2. QC and Alignment Loop
# --------------------------

echo "[$(date)] Starting Analysis Pipeline..."

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}"
    
    # Define input/output filenames
    R1_RAW="${WORKDIR}/raw_data/${SAMPLE}_R1.fq.gz"
    R2_RAW="${WORKDIR}/raw_data/${SAMPLE}_R2.fq.gz"
    R1_CLEAN="${WORKDIR}/01_clean_data/${SAMPLE}_R1.clean.fq.gz"
    R2_CLEAN="${WORKDIR}/01_clean_data/${SAMPLE}_R2.clean.fq.gz"
    BAM_PREFIX="${WORKDIR}/02_alignment/${SAMPLE}"
    
    # Step 1: Quality Control (fastp)
    echo "--> Running fastp for QC..."
    fastp -i "$R1_RAW" -I "$R2_RAW" \
          -o "$R1_CLEAN" -O "$R2_CLEAN" \
          -h "${WORKDIR}/01_clean_data/${SAMPLE}_fastp.html" \
          -j "${WORKDIR}/01_clean_data/${SAMPLE}_fastp.json" \
          --thread "$THREADS" \
          --detect_adapter_for_pe

    # Step 2: Alignment (STAR)
    echo "--> Running STAR alignment..."
    STAR --runThreadN "$THREADS" \
         --genomeDir "$GENOME_INDEX" \
         --readFilesIn "$R1_CLEAN" "$R2_CLEAN" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${BAM_PREFIX}" \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingThreadN 4 \
         --quantMode TranscriptomeSAM GeneCounts

    # (Optional) Index BAM file
    samtools index "${BAM_PREFIX}Aligned.sortedByCoord.out.bam"
done

# --------------------------
# 3. Quantification (featureCounts)
# --------------------------

echo "[$(date)] Running featureCounts..."

BAM_LIST=$(ls ${WORKDIR}/02_alignment/*Aligned.sortedByCoord.out.bam)

featureCounts -T "$THREADS" \
              -p \
              -t exon \
              -g gene_id \
              -a "$GTF_FILE" \
              -o "${WORKDIR}/03_quantification/raw_counts.txt" \
              $BAM_LIST

# --------------------------
# 4. TPM & Log2FC Calculation (Embedded R Script)
# --------------------------

echo "[$(date)] Calculating TPM and Log2FC in R..."

cat <<EOF > "${WORKDIR}/calc_results.R"
# Load raw counts
# Note: skip first row (command info) automatically via comment.char
data <- read.table("${WORKDIR}/03_quantification/raw_counts.txt", header=TRUE, comment.char="#", row.names=1)

# Extract count matrix (samples start from column 6) and gene lengths (column 5)
counts_matrix <- data[, 6:ncol(data)]
gene_lengths <- data\$Length

# Clean column names immediately (remove file paths for clarity)
# Adjust regex based on actual file path structure
colnames(counts_matrix) <- gsub(".*02_alignment\\\\.|Aligned.sortedByCoord.out.bam", "", colnames(counts_matrix))
colnames(counts_matrix) <- gsub(".*\\\\/", "", colnames(counts_matrix)) # Fallback

# --- Function: TPM Calculation ---
calculate_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  scaling_factors <- colSums(rpk) / 1000000
  tpm <- t(t(rpk) / scaling_factors)
  return(as.data.frame(tpm))
}

# 1. Generate TPM Matrix
tpm_matrix <- calculate_tpm(counts_matrix, gene_lengths)

# Save TPM Matrix
write.table(tpm_matrix, file="${WORKDIR}/04_result/gene_tpm_matrix.csv", sep=",", quote=FALSE, col.names=NA)
message("TPM matrix saved to: 04_result/gene_tpm_matrix.csv")

# --- Function: Log2 Fold Change Calculation (Treat vs Control) ---
# Note: For rigorous P-values/FDR, please use DESeq2 or edgeR on raw counts.
# This section calculates basic Log2FC based on TPM for data availability purposes.

# Auto-detect groups based on sample names
control_cols <- grep("Control", colnames(tpm_matrix), ignore.case = TRUE)
treat_cols   <- grep("Treat", colnames(tpm_matrix), ignore.case = TRUE)

if(length(control_cols) > 0 && length(treat_cols) > 0) {
  message("Groups detected for Log2FC: Control vs Treat")
  
  # Calculate group means (Adding pseudocount 1 to handle zeros)
  control_mean <- rowMeans(tpm_matrix[, control_cols, drop=FALSE])
  treat_mean   <- rowMeans(tpm_matrix[, treat_cols, drop=FALSE])
  
  # Calculate Log2FC
  # Formula: Log2(Treat_Mean + 1) - Log2(Control_Mean + 1)
  log2fc <- log2(treat_mean + 1) - log2(control_mean + 1)
  
  # Create result dataframe
  fc_results <- data.frame(
    GeneID = rownames(tpm_matrix),
    Control_Mean_TPM = control_mean,
    Treat_Mean_TPM = treat_mean,
    Log2FC = log2fc
  )
  
  # Save Log2FC File
  write.table(fc_results, file="${WORKDIR}/04_result/gene_log2fc.csv", sep=",", row.names=FALSE, quote=FALSE)
  message("Log2FC file saved to: 04_result/gene_log2fc.csv")
  
} else {
  message("Warning: Could not detect 'Control' and 'Treat' in sample names. Skipping Log2FC calculation.")
}
EOF

# Execute R script
Rscript "${WORKDIR}/calc_results.R"

# Clean up
rm "${WORKDIR}/calc_results.R"

echo "[$(date)] Pipeline completed successfully."