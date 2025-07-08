# Variant and Invariant Site Filtering for Population Genomics

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Overview

In population genomics analyses, the accurate characterization of genetic diversity requires comprehensive assessment of both polymorphic (variant) and monomorphic (invariant) sites across the genome. Traditional variant calling pipelines focus primarily on identifying and filtering variable positions, often discarding invariant sites that are crucial for calculating accurate population genetic parameters such as nucleotide diversity (π), Tajima's D, and demographic inference statistics.

## Rationale for Variant-Invariant Separation

### The Importance of Invariant Sites

Invariant sites represent genomic positions where all sequenced individuals share the same allele. These sites are essential for:

**Population Genetic Statistics:**  
- **Nucleotide diversity calculations**: Accurate π estimates require knowledge of both polymorphic and monomorphic sites  
- **Site frequency spectrum**: Complete allele frequency distributions need invariant site counts  
- **Demographic inference**: Coalescent-based methods require total sequence length information  
- **Selection analysis**: Identifying regions under purifying selection through reduced diversity  

**Comparative Genomics:**  
- **Divergence estimates**: Calculating substitution rates between populations or species  
- **Conservation analysis**: Identifying functionally constrained genomic regions  
- **Mutation rate estimation**: Determining background mutation rates across different genomic contexts  

### Differential Filtering Requirements

Variants and invariants require distinct filtering strategies due to their fundamental differences:  

**Variant Sites:**  
- Higher susceptibility to technical artifacts  
- Complex allelic configurations requiring sophisticated quality metrics  
- Strand bias and mapping quality considerations  
- Allele balance and Hardy-Weinberg equilibrium assessment  

**Invariant Sites:**  
- Simpler quality assessment based primarily on coverage depth  
- Lower risk of systematic bias  
- Focused on sequencing completeness rather than allelic accuracy  
- Reduced computational complexity for quality control  

## Implementation Workflow

### Step 1: Variant and Invariant Separation

The initial step involves separating variant and invariant sites into distinct VCF files for independent processing:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=12:00:00
#SBATCH --job-name=var_invar
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-4
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

start_time=$(date +%s)

# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chrom.list)

# Répertoire de sortie pour les fichiers
INPUT_DIR="vcfs_flt_HRun"

OUTPUT_DIR="vcfs_flt_HRun_slct_vrt"
mkdir -p ${OUTPUT_DIR}

LOG_DIR="logs"
mkdir -p ${LOG_DIR}

# Load module
module load gatk4/4.2.6.1

gatk IndexFeatureFile -I ${INPUT_DIR}/combined.chr${CHROM}.acc.HRun.vcf.gz

gatk SelectVariants \
 -V ${INPUT_DIR}/combined.chr${CHROM}.acc.HRun.vcf.gz \
 -select-type SNP \
 -xl-select-type INDEL \
 -O ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.vcf.gz \
 2>&1 >> ${LOG_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.log

gatk SelectVariants \
 -V ${INPUT_DIR}/combined.chr${CHROM}.acc.HRun.vcf.gz \
 -select-type NO_VARIATION \
 -xl-select-type INDEL \
 -O ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.INVARIANTS.vcf.gz \
 2>&1 >> ${LOG_DIR}/combined.chr${CHROM}.acc.HRun.INVARIANTS.log

end_time=$(date +%s)  # Capture l'heure de fin
duration=$((end_time - start_time))  # Calcule la durée
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

#### Technical Details

**GATK SelectVariants Parameters:**  
- `-select-type SNP`: Selects only single nucleotide polymorphisms for variant file  
- `-select-type NO_VARIATION`: Selects invariant sites for invariant file  
- `-xl-select-type INDEL`: Excludes insertion/deletion variants from both outputs  
- `-V`: Input VCF file path  
- `-O`: Output VCF file path  

**File Organization:**  
- Separate logging for variants and invariants facilitates troubleshooting  
- Chromosome-specific processing enables parallel execution  
- Compressed VCF output reduces storage requirements  


### Step 2: Intermediate Statistics Generation

Quality assessment is performed on both variant and invariant datasets to inform filtering decisions:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=01:30:00
#SBATCH --job-name=stats
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-5
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

start_time=$(date +%s)

# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chrom.list)

# Répertoire de sortie pour les fichiers
OUTPUT_DIR="vcfs_flt_HRun_slct_vrt"

# Load module
module load bcftools/1.15.1

bcftools stats ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.INVARIANTS.vcf.gz > ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.INVARIANTS.stats

bcftools stats ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.vcf.gz > ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.stats

end_time=$(date +%s)  # Capture l'heure de fin
duration=$((end_time - start_time))  # Calcule la durée
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

#### Statistical Analysis

**BCFtools stats output includes:**  
- **Site counts**: Total number of variants and invariants  
- **Quality distributions**: QUAL score histograms for threshold selection  
- **Depth statistics**: Coverage depth distributions across sites  
- **Allele frequency spectra**: Distribution of minor allele frequencies (variants only)  
- **Transition/transversion ratios**: Mutation pattern assessment (variants only)  

### Step 3: Differential Filtering and Flagging

Variants and invariants are subjected to distinct filtering criteria reflecting their different quality requirements.

Before applying any filtering criteria, it is interesting to examine the distribution of DP (depth of coverage) values in each VCF file. This data-driven approach helps determine appropriate thresholds for filtering both variant and invariant sites, ensuring that quality filters are adapted to the characteristics of the dataset.

```python linenums="1" 
pip install pysam numpy
import pysam
import matplotlib.pyplot as plt
import numpy as np

# Load VCF file
vcf_file = "combined.chr2L.acc.HRun.VARIANTS.vcf.gz"

# Define intervals
bins = list(range(0, 1551, 50))
interval_labels = [f"{bins[i]}-{bins[i+1]-1}" for i in range(len(bins)-1)]

# Initialize dict to count frequencies
dp_counts = {label: 0 for label in interval_labels}

# Read VCF file
with pysam.VariantFile(vcf_file, "r") as vcf:
    for record in vcf:
        dp = record.info.get("DP")
        if dp is not None:
            for i, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
                if low <= dp < high:
                    dp_counts[interval_labels[i]] += 1
                    break

# Prepare data for graph
intervals = list(dp_counts.keys())
frequencies = list(dp_counts.values())

# Draw graph
plt.figure(figsize=(12, 6))
plt.bar(intervals, frequencies, color="skyblue", edgecolor="black")
plt.xticks(rotation=45)
plt.xlabel("Intervalle de DP")
plt.ylabel("Fréquence")
plt.title("Distribution des valeurs de DP dans le champ INFO")
plt.tight_layout()

# Save and show graph
plt.savefig("combined.HRun.acc.VARIANTS.2L.distribution_DP.png")
plt.show()
```

Variants and invariants are filtered using different criteria, reflecting their distinct error profiles and interpretation requirements. For example, thresholds for MQ, FS, and DP are applied only to variants, while a simpler DP filter is applied to invariant sites. The script below performs these filtering steps using GATK’s VariantFiltration after indexing the input VCFs. The DP<300 cutoff, among others, was chosen based on the DP distribution plots generated with the Python script:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=01:00:00
#SBATCH --job-name=vfiltABAG
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=2
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Les VCFs ont été filtrés par accessibilité. Il est donc inutile de filtrer les HighDP.

start_time=$(date +%s)

# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chrom.list)

# Répertoire de sortie pour les fichiers
OUTPUT_DIR="vcfs_flt_HRun_slct_vrt"

LOG_DIR="Cluster_logs"
mkdir -p ${LOG_DIR}

# Load module
module load gatk4/4.2.6.1

gatk IndexFeatureFile -I ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.vcf.gz
gatk IndexFeatureFile -I ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.INVARIANTS.vcf.gz

gatk VariantFiltration \
 -R resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
 -O ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.flag.vcf.gz \
 -V ${OUTPUT_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.vcf.gz \
 --filter-expression "MQ<40.0" \
 --filter-name "MQfilter" \
 --filter-expression "DP<300" \
 --filter-name "LowDP" \
 --filter-expression "FS>60.0" \
 --filter-name "FSfilter" \
 --missing-values-evaluate-as-failing \
 2>&1 >> ${LOG_DIR}/combined.chr${CHROM}.acc.HRun.VARIANTS.flag.log

gatk VariantFiltration \
 -R resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
 -O ${OUTPUT_DIR}/AB.AG.HC.acc.INVARIANTS.flag.${CHROM}.vcf.gz \
 -V ${OUTPUT_DIR}/AB.AG.HC.acc.INVARIANTS.${CHROM}.vcf.gz \
 --filter-expression "DP<300" \
 --filter-name "LowDP" \
 --missing-values-evaluate-as-failing \
 2>&1 >> ${LOG_DIR}/AB.AG.HC.acc.INVARIANTS.flag.${CHROM}.log

end_time=$(date +%s)  # Capture l'heure de fin
duration=$((end_time - start_time))  # Calcule la durée
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

#### Filtering Criteria Analysis

**Variant Site Filters:**  
- **MQ < 40.0**: Removes sites with poor mapping quality, indicating potential alignment artifacts  
- **DP < 300**: Excludes undercovered sites with insufficient statistical power for genotype calling  
- **FS > 60.0**: Filters sites with excessive strand bias, suggesting technical artifacts  
- **--missing-values-evaluate-as-failing**: Treats missing annotation values as failing filters  

**Invariant Site Filters:**  
- **DP < 300**: Applies only depth filtering, as other quality metrics are less relevant for monomorphic sites  
- **Simplified criteria**: Reflects the lower complexity of invariant site assessment  

#### Quality Metric Explanations

**Mapping Quality (MQ):**  
- Root mean square of mapping qualities of reads supporting the variant  
- Values < 40 indicate poor uniqueness of genomic position  
- Critical for excluding paralogous sequence alignments  

**Fisher Strand Bias (FS):**  
- Phred-scaled probability of strand bias in variant calls  
- Values > 60 suggest systematic sequencing artifacts  
- Particularly important for distinguishing true variants from technical noise  

**Depth (DP):**  
- Total read depth across all samples at the site  
- Minimum threshold ensures adequate statistical power  
- Balances sensitivity with computational efficiency  

### Step 4: Final Filtering and Site Selection

The final step removes all flagged sites, retaining only high-quality variants and invariants:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=01:00:00
#SBATCH --job-name=snpfltABAG
#SBATCH -p fast
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

# Retrait des SNPs n'étant pas au statut FILTER=PASS

# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chroms_list)

# Répertoire de sortie pour les fichiers
OUTPUT_DIR=results/05_Variants

LOG_DIR="results/11_Reports/snpfiltering"
mkdir -p ${LOG_DIR}

# Load module
module load gatk4/4.2.6.1

gatk SelectVariants \
 -V ${OUTPUT_DIR}/AB.AG.HC.acc.VARIANTS.flag.${CHROM}.vcf.gz \
 --exclude-filtered \
 -O ${OUTPUT_DIR}/AB.AG.HC.acc.VARIANTS.flag.pass.${CHROM}.vcf.gz \
 2>&1 >> ${LOG_DIR}/AB.AG.HC.acc.VARIANTS.flag.pass.${CHROM}.log

gatk SelectVariants \
 -V ${OUTPUT_DIR}//AB.AG.HC.acc.INVARIANTS.flag.${CHROM}.vcf.gz \
 --exclude-filtered \
 -O ${OUTPUT_DIR}/AB.AG.HC.acc.INVARIANTS.flag.pass.${CHROM}.vcf.gz \
 2>&1 >> ${LOG_DIR}/AB.AG.HC.acc.INVARIANTS.flag.pass.${CHROM}.log

end_time=$(date +%s)  # Capture l'heure de fin
duration=$((end_time - start_time))  # Calcule la durée
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

#### Final Selection Parameters

**--exclude-filtered**: Removes all sites that failed any filter criterion, retaining only sites marked as "PASS"  

**Quality Control Verification**: Final datasets should be assessed for:  
- Retained site counts and genomic distribution  
- Quality metric distributions post-filtering  
- Transition/transversion ratios (variants)  
- Coverage uniformity across chromosomes  
 
## Population Genomics Applications

### Downstream Analysis Compatibility

**Nucleotide Diversity Calculations:**
```bash linenums="1"
# Example calculation using both datasets
total_sites=$(bcftools view -H variants.pass.vcf.gz | wc -l) + $(bcftools view -H invariants.pass.vcf.gz | wc -l)
polymorphic_sites=$(bcftools view -H variants.pass.vcf.gz | wc -l)
pi=$(echo "scale=6; $polymorphic_sites / $total_sites" | bc)
```

**Site Frequency Spectrum Analysis:**  
- Variants provide polymorphic site frequencies  
- Invariants contribute to the monomorphic class (frequency = 0 or 1)  
- Combined datasets enable accurate demographic inference  

**Selective Sweep Detection:**  
- Reduced diversity in invariant sites indicates purifying selection  
- Elevated diversity in variant sites suggests balancing selection  
- Comparative analysis identifies selection signatures  

### Computational Considerations

**Storage Requirements:**  
- Invariant sites typically outnumber variants by 10-100 fold  
- Compressed VCF formats essential for storage efficiency  
- Indexed files enable rapid region-specific access  

**Processing Efficiency:**  
- Parallel chromosome processing reduces computational time  
- Separate filtering pipelines optimize resource utilization  
- Modular workflow facilitates parameter optimization  

## Best Practices and Recommendations

### Filter Threshold Selection

**Conservative Approach:**  
- Higher depth requirements (DP > 500) for critical analyses  
- Stricter quality thresholds (MQ > 50) for publication-quality datasets  
- Additional filters (BaseQRankSum, ReadPosRankSum) for variant sites  

**Moderate Approach:**  
- Balanced sensitivity and specificity (current implementation)  
- Suitable for most population genetic analyses   
- Adequate for demographic inference and diversity estimation  

**Lenient Approach:**  
- Lower thresholds (DP > 100) for sample-limited studies  
- Relaxed quality criteria for exploratory analyses    
- Requires additional downstream validation  

### Quality Control Validation

**Statistical Assessment:**  
1. **Ti/Tv ratios**: Expected values ~2.0-2.1 for human, variable for other species  
2. **Depth distributions**: Should follow expected coverage patterns  
3. **Allele frequency spectra**: Should match population genetic expectations  
4. **Genomic distribution**: Uniform distribution across accessible regions  

**Comparative Analysis:**  
- Benchmark against high-quality reference datasets  
- Cross-platform validation using orthogonal technologies  
- Population-specific optimization based on demographic history  

### Integration Strategies

**Downstream Merging:**
```bash linenums="1"
# Combine filtered variants and invariants for analysis
bcftools concat \
  variants.pass.vcf.gz \
  invariants.pass.vcf.gz \
  -a -D -Oz -o combined.filtered.vcf.gz
```

**Analysis-Specific Subsets:**  
- Coding regions: Enhanced filtering for functional impact assessment  
- Intergenic regions: Relaxed filtering for neutral evolution studies  
- Regulatory elements: Intermediate filtering for selection analysis  

## Limitations and Considerations

### Technical Limitations

**Ascertainment Bias:**  
- Variant discovery influenced by reference genome quality  
- Population-specific variants may be systematically missed  
- Affects demographic inference accuracy  

**Coverage Heterogeneity:**  
- Uneven sequencing depth across genomic regions  
- Systematic bias in variant detection sensitivity  
- Requires careful normalization in comparative studies  

### Biological Considerations

**Population Structure:**  
- Filtering thresholds may need population-specific adjustment  
- Admixed populations require modified quality criteria  
- Demographic history affects optimal parameter selection  

**Evolutionary Constraints:**  
- Functional regions may require different filtering approaches  
- Selective pressure influences variant quality distributions  
- Conservation levels affect optimal sensitivity/specificity balance  

## Future Developments

### Methodological Improvements

**Machine Learning Approaches:**  
- Automated filter threshold optimization  
- Context-specific quality assessment  
- Population-adapted filtering strategies  

**Long-Read Integration:**  
- Improved structural variant detection  
- Enhanced repetitive region characterization  
- Reduced reference bias in variant calling  

### Analytical Advances

**Comprehensive Genomic Datasets:**  
- Integration of regulatory annotations  
- Functional impact prediction  
- Evolutionary constraint incorporation  

**Population-Specific Resources:**  
- Ethnically diverse reference panels  
- Population-specific quality metrics  
- Adaptive filtering frameworks  

## Conclusion

The separation and differential filtering of variants and invariants represents a sophisticated approach to population genomic analysis that maximizes both data quality and analytical power. By recognizing the distinct characteristics and requirements of polymorphic and monomorphic sites, this workflow enables accurate estimation of population genetic parameters while maintaining computational efficiency.

The implementation described here provides a robust framework for population genomic studies, with modular components that can be adapted to specific research questions and sample characteristics. The careful balance between sensitivity and specificity, combined with comprehensive quality control measures, ensures the reliability of downstream analyses while preserving the statistical power necessary for demographic inference and selection analysis.

Success in population genomics depends critically on the quality of the underlying data, and the differential treatment of variants and invariants represents a key advancement in achieving this goal. As sequencing technologies and analytical methods continue to evolve, the principles outlined in this workflow will remain essential for extracting meaningful biological insights from genomic variation data.