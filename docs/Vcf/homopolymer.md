# Homopolymer Run Filtering in VCF Files

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Context: Post-Variant Calling Quality Control

This filtering step is particularly relevant in the context of post-variant calling quality control, especially when variants have been called using `legacy tools such as GATK's UnifiedGenotyper`. While modern variant callers like HaplotypeCaller incorporate sophisticated algorithms that better handle problematic genomic regions including homopolymer runs during the variant calling process, UnifiedGenotyper produces inferior results for indel calling compared to HaplotypeCaller, necessitating additional post-processing steps for homopolymer-associated variants.

## Understanding Homopolymer Runs

### Definition and Characteristics

A homopolymer run (also called homopolymer tract) is a genomic sequence consisting of consecutive identical nucleotides. These sequences are characterized by:

- **Repetitive nucleotide composition**: Stretches of single nucleotides (e.g., AAAAA, TTTTTT, GGGGGG, CCCCCC)
- **Variable length**: Typically ranging from 3 to over 20 consecutive identical bases
- **Genome-wide distribution**: Present throughout all genomic regions but with varying frequencies
- **Sequence context dependency**: Often flanked by complex secondary structures

### Examples of Homopolymer Runs
```
Reference: ATCG[AAAAA]TGCAT
Variant:   ATCG[AAAAAA]TGCAT  (insertion of one A)
Variant:   ATCG[AAAA]TGCAT    (deletion of one A)
```

## Technical Challenges in Homopolymer Regions

### Sequencing Technology Limitations

**Illumina Sequencing Issues:**  
- **Phasing errors**: Synchronization loss during sequencing-by-synthesis  
- **Signal intensity variations**: Difficulty in accurately measuring fluorescent signals  
- **Cluster density effects**: Overlapping signals from adjacent clusters  

**Ion Torrent Specific Problems:**  
- **Homopolymer length miscalls**: Inherent limitation of pH-based detection  
- **Systematic length bias**: Tendency to under-call or over-call homopolymer lengths  
- **Flow signal saturation**: Signal plateau effects in long homopolymer runs  

### Bioinformatics Challenges

**Alignment Artifacts:**  
- **Mapping ambiguity**: Multiple valid alignments for reads spanning homopolymer regions  
- **Indel calling errors**: Misinterpretation of sequencing errors as true variants  
- **Strand bias**: Asymmetric representation of variants on forward and reverse strands  

**Variant Calling Complications:**  
- **False positive indels**: Sequencing errors interpreted as genuine insertions/deletions  
- **Allelic imbalance**: Preferential amplification of certain alleles  
- **Genotyping errors**: Incorrect assignment of homozygous vs. heterozygous states  

## Rationale for Homopolymer Run Filtering

### Quality Control Necessity

The removal of variants in homopolymer runs is justified by several factors:

1. **High false positive rate**: Homopolymer regions exhibit elevated rates of spurious variant calls
2. **Systematic sequencing bias**: Technology-specific artifacts that cannot be easily corrected
3. **Reduced genotyping confidence**: Lower quality scores and unreliable allele frequency estimates
4. **Downstream analysis interference**: Contamination of population genetic and association studies

### Impact on Variant Interpretation

**Clinical Implications:**  
- **Diagnostic accuracy**: Reduced false positive rates in medical genetic testing  
- **Treatment decisions**: More reliable variant classification for therapeutic targeting  
- **Genetic counseling**: Improved confidence in variant pathogenicity assessment  

**Research Applications:**  
- **Population genetics**: Cleaner datasets for demographic inference  
- **Genome-wide association studies**: Reduced noise in association testing  
- **Evolutionary genomics**: More accurate phylogenetic reconstructions  

## Implementation Workflow

The following SLURM script demonstrates the annotation and identification of homopolymer runs using GATK tools:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=11-10:00:00
#SBATCH --job-name=hmplmr
#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=2
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

start_time=$(date +%s)

# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chrom.list)
INPUT_DIR="vcfs_filtered/"

# Répertoire de sortie pour les fichiers
OUTPUT_DIR="vcfs_flt_HRun"
mkdir -p ${OUTPUT_DIR}

# Load module
module load gatk4/4.2.6.1
module load gatk/3.8

gatk IndexFeatureFile -I ${INPUT_DIR}/combined.chr${CHROM}.acc.vcf.gz

gatk3 -T VariantAnnotator \
    -V ${INPUT_DIR}/combined.chr${CHROM}.acc.vcf.gz \
    -R /shared/projects/invalbo/bwambae/resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
    -A HomopolymerRun \
    -o ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.vcf.gz

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

## Technical Implementation Details

### GATK VariantAnnotator Parameters

**Core Arguments:**  
- `-T VariantAnnotator`: Specifies the GATK tool for adding annotations to variants  
- `-V`: Input VCF file path (already filtered by accessibility)  
- `-R`: Reference genome file in FASTA format  
- `-A HomopolymerRun`: Annotation module that identifies homopolymer contexts  
- `-o`: Output VCF file with homopolymer annotations  

**Resource Requirements:**  
- **Memory allocation**: 16GB to handle large VCF files and reference genome  
- **CPU utilization**: 8 cores for parallel processing  
- **Runtime**: Extended time allocation (11 days) for large-scale genomic data  

### HomopolymerRun Annotation

The `HomopolymerRun` annotation adds the following information to the VCF INFO field:

- **HRun**: Length of the longest homopolymer run in the variant context  
- **Context window**: Typically examines ±10 bases around the variant position  
- **Strand consideration**: Evaluates both forward and reverse strand contexts  

### Subsequent Filtering Steps

After annotation, variants can be filtered based on homopolymer length thresholds:

```bash linenums="1"
# Example filtering command (not in original script)
gatk VariantFiltration \
    -V combined.chr${CHROM}.acc.HRun.vcf.gz \
    -filter "HRun > 4" \
    --filter-name "HomopolymerRun" \
    -O combined.chr${CHROM}.acc.HRun.filtered.vcf.gz
```

## Best Practices and Considerations

### Threshold Selection

**Conservative Approach (HRun > 3):**  
- Removes most homopolymer-associated artifacts  
- Maintains high variant quality  
- May remove some true variants  

**Moderate Approach (HRun > 5):**  
- Balances quality and sensitivity  
- Suitable for most research applications  
- Commonly used in population studies  

**Lenient Approach (HRun > 7):**  
- Retains more variants for analysis  
- Higher false positive rate  
- Requires additional quality control  

### Validation Strategies

1. **Sanger sequencing confirmation**: Validate a subset of filtered variants  
2. **Comparison with high-quality datasets**: Benchmark against reference populations  
3. **Technology-specific validation**: Use orthogonal sequencing platforms  
4. **Functional annotation**: Assess impact on protein-coding regions  

## Limitations and Considerations

### Potential Data Loss

**Functional Impact:**  
- Removal of true variants in homopolymer regions  
- Loss of clinically relevant mutations  
- Reduced power for association studies  

**Genomic Context:**  
- Regulatory elements within homopolymer tracts  
- Evolutionary important variations  
- Species-specific adaptation signals  

### Alternative Approaches

**Improved Sequencing Technologies:**  
- Long-read sequencing (PacBio, Oxford Nanopore)  
- Linked-read technologies (10x Genomics)  
- Hybrid assembly approaches  

**Advanced Bioinformatics Methods:**  
- Machine learning-based variant calling  
- Consensus calling from multiple algorithms  
- Probabilistic models for homopolymer regions  

## Conclusion

Homopolymer run filtering represents a critical quality control step in variant calling pipelines, particularly for data generated using legacy tools or technologies prone to homopolymer-associated artifacts. While this filtering approach may result in the loss of some true variants, it significantly improves the overall quality and reliability of variant calls by removing systematic biases and false positives. The implementation should be tailored to the specific sequencing technology, analysis objectives, and downstream applications, with careful consideration of the trade-offs between sensitivity and specificity.

Future developments in sequencing technologies and bioinformatics algorithms continue to improve the accuracy of variant calling in challenging genomic regions, potentially reducing the need for such aggressive filtering approaches. However, for current standard practices, homopolymer run filtering remains an essential component of robust variant calling pipelines.