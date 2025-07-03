# Genotype and Individual-Based Missing Data Filtering in VCF Files

Accurate variant analysis requires not only high-confidence genotypes but also careful handling of **missing data**. Two complementary filtering strategies are commonly applied to Variant Call Format (VCF) files to improve data quality and downstream interpretability: **site-based (horizontal) filtering** and **individual-based (vertical) filtering**.

## Horizontal Filtering: Site-Level Missing Data

This type of filtering focuses on the **genomic positions (sites)** themselves. It removes sites for which a defined by user proportion of samples lack genotype calls. Excessive missingness at a site may result from poor alignment, low coverage, or systematic sequencing issues. Filtering such positions helps reduce noise and avoid biases in population-level statistics. 

In this project, we retained only those positions where **at least 80% of sites** had a valid genotype call.

The script below shows how this filtering step was performed per chromosome using `vcftools`:


```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=02:30:00
#SBATCH --job-name=GTsABAG
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

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chroms_list)
OUTPUT_DIR=results/05_Variants

module load vcftools/0.1.16 htslib/1.14

vcftools \
 --gzvcf ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.sort.${CHROM}.vcf.gz \
 --max-missing 0.8 \
 --recode \
 --recode-INFO-all \
 --stdout \
 | bgzip -c > ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.gt80.${CHROM}.vcf.gz \
 && tabix ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.gt80.${CHROM}.vcf.gz

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

## Vertical Filtering: Individual-Level Missing Data

Once low-quality sites are removed, the next step is to evaluate the quality of each sample (individual) based on the proportion of missing genotype calls across retained sites. Samples with extensive missing data may reflect poor sequencing, contamination, or technical failure. Removing such individuals prevents them from skewing population genetic inferences or association analyses. Here, we excluded individuals with more than 70% missing genotypes.

The script below identifies and filters out these low-quality individuals:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=02:00:00
#SBATCH --job-name=ifltABAG
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

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chroms_list)
OUTPUT_DIR=results/05_Variants
LOG_DIR="results/11_Reports/ind70"
mkdir -p ${LOG_DIR}

module load vcftools/0.1.16 htslib/1.14

SUBSET_VCF=${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.gt80.${CHROM}.vcf.gz
OUT=${OUTPUT_DIR}/AB.AG.HC.vcf.subset.${CHROM}.indv

# Calculate missingness per individual
vcftools \
 --gzvcf ${SUBSET_VCF} \
 --missing-indv \
 --out ${OUT}

# Identify individuals with >70% missing data
mawk '$5 > "0.7"' ${OUT}.imiss | cut -f1 > ${OUTPUT_DIR}/AB.AG.HC.lowDP.indv

VCF_IN=${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.gt80.${CHROM}.vcf.gz
VCF_OUT=${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.gt80.ind70.${CHROM}.vcf.gz

# Remove those individuals from the VCF
vcftools --gzvcf ${VCF_IN} \
        --remove-indels \
        --remove ${OUTPUT_DIR}/AB.AG.HC.lowDP.indv \
        --recode --stdout | bgzip -c > \
        ${VCF_OUT}

end_time=$(date +%s)
duration=$((end_time - start_time))
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

## Conclusion

Together, these two layers of missing data filtering—horizontal (site-based) and vertical (individual-based)—enhance the overall reliability of the genotype dataset. This is especially important in large-scale studies where incomplete or uneven data can lead to biases in genetic diversity estimates, population structure, or association mapping. Filtering thresholds (e.g., 80% completeness per site and 70% per individual) should be selected based on study design, sample size, and sequencing depth.