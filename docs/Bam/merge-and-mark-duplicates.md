# BAM File Merging and Duplicate Marking Workflow

## Overview

The processing of high-throughput sequencing data often involves multiple sequencing lanes or runs for a single biological sample to achieve adequate coverage depth or to distribute the sequencing load across multiple flow cells. This necessitates a systematic approach to consolidate the resulting alignment files and identify PCR duplicates that may introduce bias in downstream analyses.

## Sequential Processing Strategy

### 1. BAM File Merging

The initial step involves merging multiple BAM files corresponding to different sequencing lanes of the same biological sample. This consolidation is performed using Picard's `MergeSamFiles` tool, which combines aligned reads while preserving read group information and maintaining coordinate-based sorting.

**Key advantages of merging first:**  
- **Computational efficiency**: Processing a single merged file is more resource-efficient than handling multiple separate files in subsequent steps  
- **Read group preservation**: Maintains lane-specific metadata essential for quality control and bias detection  
- **Coordinate optimization**: Ensures proper genomic coordinate sorting across all merged reads  

### 2. Duplicate Marking

Following the merge operation, PCR and optical duplicates are identified and marked using Picard's `MarkDuplicates` tool. This step is crucial for maintaining data quality and preventing artificial inflation of coverage metrics.

**Scientific rationale for post-merge duplicate detection:**  
- **Cross-lane duplicate identification**: Duplicates may occur across different sequencing lanes, which can only be detected after merging  
- **Improved accuracy**: Global duplicate detection across the entire dataset provides more accurate duplicate identification than lane-specific processing  
- **Standardized metrics**: Consolidated duplicate statistics reflect the true duplication rate for the biological sample  

## Implementation Script

The workflow is implemented through the following SLURM-compatible bash script:

```bash
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=08:00:00
#SBATCH --job-name=merge_bams
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-22
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Récupération du nom du fichier BAM à partir du fichier de liste
BAMNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/bam_list)

# Répertoire de sortie pour les fichiers fusionnés
OUTPUT_DIR=results/04_Polishing
mkdir -p $OUTPUT_DIR

# Nom de l'échantillon extrait du nom du fichier BAM
SAMPLE_NAME=$(echo $BAMNAME | awk -F'[_]' '{print $1"_"$2"_"}')

# Fichiers BAM de différentes lanes pour le même échantillon
LANE_BAMS=($(ls results/02_Mapping/${SAMPLE_NAME}*sorted.bam))

# Fichier de sortie fusionné
MERGED_FILE=$OUTPUT_DIR/${SAMPLE_NAME}merged.bam

# QC directory
METRICS_DIR="qc/markdup"
mkdir -p ${METRICS_DIR}

# Load module
module load picard/2.23.5

# Utilisation de picard MergeSamFiles pour fusionner les fichiers BAM
picard MergeSamFiles \
      $(for bam in "${LANE_BAMS[@]}"; do echo "-I $bam"; done) \
      -O $MERGED_FILE \
      --USE_THREADING true \
      --SORT_ORDER coordinate \
      2> results/11_Reports/mergesam/${SAMPLELANE}-SortSam.out

# Suppression des fichiers temporaires si nécessaire
# rm "${LANE_BAMS[@]}"

echo "Fusion des fichiers BAM pour l'échantillon $SAMPLE_NAME terminée."

picard MarkDuplicates -Xmx16000M \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT ${MERGED_FILE} \
  --OUTPUT ${OUTPUT_DIR}/${SAMPLE_NAME}merged_marked.bam \
  --TMP_DIR /tmp/ \
  --METRICS_FILE ${METRICS_DIR}/${SAMPLE_NAME}-merged_marked.metrics \
  2> results/11_Reports/markduplicates/${SAMPLELANE}_merged_marked.out
```

## Technical Specifications

### Merge Parameters
- **Threading**: Enabled for improved performance (`--USE_THREADING true`)
- **Sort order**: Maintains coordinate-based sorting for downstream compatibility
- **Input handling**: Dynamic input file specification using shell array expansion

### Duplicate Marking Parameters
- **Memory allocation**: 16GB heap space (`-Xmx16000M`) for large datasets
- **Duplicate removal**: Physical removal of duplicate reads (`--REMOVE_DUPLICATES true`)
- **Index creation**: Automatic BAI index generation (`--CREATE_INDEX true`)
- **Quality metrics**: Comprehensive duplicate statistics output

## Quality Control Outputs

The workflow generates several quality control files:  
- **Merge logs**: Detailed processing information for the merge operation  
- **Duplicate metrics**: Statistics on duplication rates and patterns  
- **BAM indices**: Binary indices for efficient random access  

## Computational Considerations

This workflow is optimized for high-performance computing environments with SLURM job scheduling. The array-based execution allows parallel processing of multiple samples while maintaining resource efficiency through single-threaded operations with moderate memory requirements (16GB RAM, 8-hour time limit).

The sequential merge-then-mark approach ensures optimal data quality while minimizing computational overhead, making it suitable for large-scale genomic studies requiring robust duplicate detection across multi-lane sequencing datasets.