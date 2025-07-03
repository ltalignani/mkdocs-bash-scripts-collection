# Indel Realignment Workflow (for UnifiedGenotyper Variant-caller)

## Overview

Indel realignment represents a critical preprocessing step in genomic variant calling pipelines, particularly when employing legacy tools such as GATK's `UnifiedGenotyper`. This process addresses systematic alignment artifacts that occur around insertion and deletion (indel) polymorphisms, which are prevalent in mosquito genomes and can significantly impact variant calling accuracy.

## Scientific Rationale

### Alignment Artifacts in Mosquito Genomes

Anopheles genomes, like those of other dipteran species, exhibit substantial structural variation including frequent indel polymorphisms. Standard read alignment algorithms often produce suboptimal alignments in regions flanking true indels, leading to:

- **False positive SNP calls**: Misaligned reads create apparent sequence mismatches
- **Mapping quality degradation**: Ambiguous alignments reduce confidence scores
- **Allelic bias**: Preferential alignment of reference-matching reads over variant-containing reads

### Legacy Tool Compatibility

The UnifiedGenotyper, while superseded by more recent algorithms, remains relevant for specific applications involving:

- **Historical data consistency**: Maintaining compatibility with existing 1000 Genomes Project mosquito datasets
- **Comparative genomics**: Ensuring methodological consistency across temporal datasets
- **Resource-constrained environments**: Lower computational requirements compared to modern alternatives

## Indel Realignment Methodology

### Two-Stage Process

#### Stage 1: Target Identification
The `RealignerTargetCreator` tool identifies genomic intervals requiring realignment by:
- Detecting regions with elevated mismatch rates
- Identifying potential indel sites based on alignment patterns
- Creating target intervals for focused realignment

#### Stage 2: Local Realignment
The `IndelRealigner` performs local realignment within identified target regions by:
- Constructing alternative alignment hypotheses
- Evaluating alignment quality scores
- Selecting optimal alignments that minimize mismatches


```bash linenums="1"
##!/bin/bash

###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=08:00:00
#SBATCH --job-name=indelrealignment
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=1-26
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Récupération du nom du fichier BAM à partir du fichier de liste
BAMNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/bam_list)

# Répertoire de sortie pour les fichiers
OUTPUT_DIR=results/04_Polishing

# Fichier d'entrée
MARKED_BAM=${OUTPUT_DIR}/${BAMNAME}_marked.bam

# Fichier filtré
FILTERED_BAM=${OUTPUT_DIR}/${BAMNAME}_marked_filtered.bam

# Fichier bam taggé
TAGGED_BAM=${OUTPUT_DIR}/${BAMNAME}_marked_filtered_tagged.bam

# Fichiers d'intervalles à réaligner
INTERVALS_FILE=${OUTPUT_DIR}/${BAMNAME}_marked_filtered_tagged.intervals

# Fichier bam réaligné
REALIGNED_BAM=${OUTPUT_DIR}/${BAMNAME}_marked_filtered_realigned.bam

# Report dir
REPORT_DIR=results/11_Reports/indelrealignment
mkdir -p ${REPORT_DIR}

# Load module
module load samtools/1.15.1

samtools view -h -q 30 ${MARKED_BAM} -o ${FILTERED_BAM}

# Load module
module load picard/2.23.5

picard SetNmMdAndUqTags \
      -R resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
      -I ${FILTERED_BAM} \
      -O ${TAGGED_BAM} > {log} 2>&1

samtools index ${TAGGED_BAM}

# Load module
module load gatk/3.8

gatk3 -T RealignerTargetCreator \
      --num_threads 2 \
      -R resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
      -I ${TAGGED_BAM} \
      --defaultBaseQualities 20 --filter_reads_with_N_cigar \
      --out ${INTERVALS_FILE} \
      2> ${REPORT_DIR}/${BAMNAME}_filtered_intervals.out

gatk3 -Xmx64g -T IndelRealigner \
      -I ${TAGGED_BAM} \
      -R resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
      --targetIntervals ${INTERVALS_FILE} \
      --defaultBaseQualities 20 --filter_reads_with_N_cigar \
      --out ${REALIGNED_BAM} \
      2> ${REPORT_DIR}/${BAMNAME}_filtered_realignment.out

awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' ${INTERVALS_FILE} \
      1> ${OUTPUT_DIR}/${BAMNAME}_realignertargetcreator.bed \
      2> ${REPORT_DIR}/${BAMNAME}_awk_filtered.out
```

## Conclusion

Indel realignment represents an essential preprocessing step for variant calling in mosquito genomics, particularly when employing the legacy tool UnifiedGenotyper. The systematic correction of alignment artifacts around indel sites significantly improves variant calling accuracy and enables robust population genomic analyses of Anopheles species. Integration with 1000 Genomes mosquito datasets ensures methodological consistency and facilitates comparative genomic studies across different research groups and time periods.