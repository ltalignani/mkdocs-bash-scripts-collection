# Variant Calling in Genomics

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Introduction

Variant calling is a fundamental step in genomics that identifies differences between sequenced genomes and a reference genome. These differences, or variants, include single nucleotide polymorphisms (SNPs), insertions, deletions (indels), and structural variants. Accurate variant calling is crucial for understanding genetic diversity, disease susceptibility, and evolutionary relationships.

The process involves aligning sequencing reads to a reference genome and then identifying positions where the aligned reads differ from the reference. Various computational tools have been developed to perform this task, each with specific strengths and applications.

## GATK HaplotypeCaller

The Genome Analysis Toolkit (GATK) HaplotypeCaller is currently the gold standard for variant calling in human genomics and is widely adopted across various species. This tool uses a sophisticated approach that assembles haplotypes in regions with variation, making it particularly effective at calling variants in complex genomic regions.

### Key Features

HaplotypeCaller employs a local de novo assembly approach that:  
- Identifies active regions where variation is likely to occur  
- Assembles possible haplotypes in these regions  
- Realigns reads to the most likely haplotypes  
- Calls variants based on the assembled haplotypes  

This methodology provides superior accuracy compared to position-based callers, especially for indels and complex variants.

### Implementation Example

The following script demonstrates a typical HaplotypeCaller workflow for generating genomic variant call format (GVCF) files:

```bash
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --job-name=haplotype-gvcf
#SBATCH --time=15-23:00:00
#SBATCH -p long
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # ntasks
#SBATCH --cpus-per-task 4
#SBATCH --mem=48G
#SBATCH -o Cluster_logs/%x-%A-%a.out
#SBATCH -e Cluster_logs/%x-%A-%a.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
#SBATCH --array=1-284 # 71 samples x 4 chromosomes
###################################################################
module load gatk4/4.2.6.1
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/ae_sample_scaffolds)
SAMPLE=$(echo $FILE | awk '{print $1}')
SCAFFOLD=$(echo $FILE | awk '{print $2}')
# Define path to .bam files
BAM_DIR="results/04_Polishing"
# Check if directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Le dossier $BAM_DIR n'existe pas."
    exit 1
fi
OUTPUT_DIR="results/05_Variants"/${SCAFFOLD}
mkdir -p "${OUTPUT_DIR}"
# Create directory for log files
LOG_FOLDER="results/11_Reports/haplotypecaller"
mkdir -p "${LOG_FOLDER}"
# Run HaplotyCaller: 
# -Xmx64g: defines the maximum memory size that the JVM can allocate to the application during execution
gatk --java-options "-Xmx48g" HaplotypeCaller \
  -R resources/genomes/AalbF5_filtered.fasta \
  -I ${BAM_DIR}/${SAMPLE}_marked.bam \
  -O ${OUTPUT_DIR}/${SAMPLE}-${SCAFFOLD}.g.vcf.gz \
  -ERC GVCF \
  -L ${SCAFFOLD} \
  2> ${LOG_FOLDER}/${SAMPLE}-${SCAFFOLD}-hapcall-gvcf.out
```

This script utilizes SLURM job arrays to process multiple samples and chromosomes in parallel, demonstrating the scalability required for large-scale genomic projects.


The output of HaplotypeCaller in `-ERC GVCF` mode is a `genomic VCF` (gVCF) file that records genotype likelihoods and variant calls for every genomic position, including non-variant sites. It enables joint genotyping by capturing reference confidence information. 

Unfortunately, the GVCF format is not always supported by downstream analysis tools, which often require a `standard VCF` file. To convert GVCF files into a conventional VCF format, the GATK pipeline provides the `GenomicsDBImport` tool to aggregate multiple GVCFs into a GenomicsDB datastore, followed by `GenotypeGVCFs`, which performs joint genotyping across samples and produces a final multisample VCF suitable for most variant analysis workflows. 

[See GenomicsDBImport and GenotypeGVCFs in GATK4 Workflow](genomicsdbimport.md#the-gvf-based-workflow)

## GATK UnifiedGenotyper

The UnifiedGenotyper was one of the earlier variant callers in the GATK suite and served as a workhorse for many genomic studies. However, it has been largely superseded by HaplotypeCaller due to several limitations in its algorithm.

### Current Status and Applications

UnifiedGenotyper is **increasingly deprecated** and is now **rarely used** in modern genomics workflows. Its usage has been largely **discontinued** for most applications, with one notable exception: **projects involving the 1000 Genomes Anopheles gambiae dataset**. This specific use case maintains UnifiedGenotyper compatibility due to:

- Historical continuity with existing analyses
- Standardized protocols established for malaria vector genomics
- Specific parameter optimizations for Anopheles gambiae genetic characteristics
- Integration with the Ag1000G (Anopheles gambiae 1000 Genomes) project pipeline

### Implementation for Anopheles gambiae Projects

The following script shows a typical UnifiedGenotyper implementation for Anopheles gambiae genomic data:

```bash
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --job-name=UnifiedGenotyper
#SBATCH --time=2-23:00:00
#SBATCH -p long
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # ntasks
#SBATCH --cpus-per-task 4
#SBATCH --mem=48G
#SBATCH --array=1-4
#SBATCH -o Cluster_logs/%x-%A-%a.out
#SBATCH -e Cluster_logs/%x-%A-%a.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
###################################################################
# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chrom.list)
# Load module
module load gatk/3.8
BAM_LIST="/shared/projects/invalbo/AG3/SN/info_files/bam.list"
REF="/shared/projects/invalbo/AG3/SN/resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
OUTPUT_DIR="results/05_Variants"
mkdir -p "${OUTPUT_DIR}"
# Create directory for log files
LOG="results/11_Reports/unifiedgenotyper"
mkdir -p "${LOG}"
gatk3 -Xmx48G -T UnifiedGenotyper --num_threads 4 --num_cpu_threads_per_data_thread 4 -I ${BAM_LIST} -R ${REF} -L ${CHROM} \
      --out ${OUTPUT_DIR}/SN.ag3.${CHROM}.vcf --genotype_likelihoods_model BOTH --genotyping_mode DISCOVERY --heterozygosity 0.015 \
      --heterozygosity_stdev 0.05 --indel_heterozygosity 0.001 --downsampling_type BY_SAMPLE -dcov 250 --output_mode EMIT_ALL_SITES \
      --min_base_quality_score 17 -stand_call_conf 0.0 -contamination 0.0 -A DepthPerAlleleBySample \
      -A RMSMappingQuality -A Coverage -A FisherStrand -A StrandOddsRatio -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A QualByDepth \
      -A ReadPosRankSumTest -XA ExcessHet -XA InbreedingCoeff -XA MappingQualityZero -XA HaplotypeScore -XA SpanningDeletions -XA ChromosomeCounts \
      &> ${LOG}/SN.ag3.${CHROM}.log
```

Note the specialized parameters optimized for Anopheles gambiae genetics, including heterozygosity rates and quality thresholds tailored to this species.

## BCFtools for Variant Calling

BCFtools, part of the SAMtools suite, provides an alternative approach to variant calling that is particularly valued for its computational efficiency and integration with other bioinformatics tools.

### Algorithm and Approach

BCFtools uses a Bayesian approach to variant calling through its `mpileup` and `call` commands. The process involves:

1. **Pileup generation**: Creating a summary of read alignments at each position
2. **Likelihood calculation**: Computing genotype likelihoods based on base qualities and alignment characteristics
3. **Variant calling**: Applying statistical models to identify the most likely genotypes

### Advantages of BCFtools

- **Computational efficiency**: Generally faster than GATK tools, especially for large datasets
- **Memory efficiency**: Lower memory requirements make it suitable for resource-constrained environments
- **Flexibility**: Extensive command-line options for fine-tuning calling parameters
- **Integration**: Seamless integration with SAMtools and other tools in the suite

### Typical BCFtools Workflow

```bash
# Generate pileup and call variants
bcftools mpileup -f reference.fasta input.bam | bcftools call -mv -O z -o output.vcf.gz

# Alternative approach with separate steps
bcftools mpileup -f reference.fasta input.bam > pileup.vcf
bcftools call -mv pileup.vcf -O z -o variants.vcf.gz
```

### When to Choose BCFtools

BCFtools is particularly suitable for:
- High-throughput projects requiring computational efficiency
- Population genomics studies with large sample sizes
- Projects with limited computational resources
- Analyses requiring tight integration with SAMtools workflows

## Comparative Considerations

### Accuracy vs. Speed Trade-offs

- **GATK HaplotypeCaller**: Highest accuracy, especially for complex variants, but computationally intensive
- **BCFtools**: Good balance of accuracy and speed, particularly efficient for SNP calling
- **UnifiedGenotyper**: Historical tool with limited current applications

### Project-Specific Recommendations

- **Human clinical genomics**: GATK HaplotypeCaller is the standard
- **Population genomics**: BCFtools or HaplotypeCaller depending on accuracy requirements
- **Anopheles gambiae 1000 Genomes projects**: UnifiedGenotyper for continuity with established protocols
- **Resource-limited environments**: BCFtools for efficiency

## Best Practices

### Quality Control

Regardless of the chosen variant caller, implementing robust quality control measures is essential:

- **Base quality filtering**: Remove low-quality bases from consideration
- **Mapping quality assessment**: Filter poorly mapped reads
- **Depth filtering**: Apply appropriate coverage thresholds
- **Allele frequency filtering**: Consider population-specific allele frequencies

### Post-Calling Processing

All variant calling workflows should include:

- **Variant quality score recalibration (VQSR)**: Particularly important for GATK workflows
- **Hard filtering**: Apply quality thresholds when VQSR is not feasible
- **Annotation**: Add functional and population annotations
- **Validation**: Confirm critical variants through independent methods

## Conclusion

The choice of variant calling tool depends on specific project requirements, computational resources, and accuracy needs. While GATK HaplotypeCaller has become the gold standard for most applications, BCFtools remains valuable for efficient processing of large datasets. UnifiedGenotyper, though largely deprecated, continues to serve specific niche applications, particularly in Anopheles gambiae genomics where protocol continuity is essential.

As sequencing technologies and analytical methods continue to evolve, the field of variant calling will likely see further innovations in accuracy, efficiency, and specialized applications for different genomic contexts.