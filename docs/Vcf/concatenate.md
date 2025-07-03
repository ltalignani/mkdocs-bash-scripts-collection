# Concatenate Variants and Invariants into a Unified VCF File

After separate processing and filtering of **variant** and **invariant** genomic sites, a critical step in the variant calling pipeline is to **recombine** these data into a single VCF file per chromosome. 

This step is referred to as concatenation, where two VCF files—typically from the same set of individuals but representing different sets of genomic positions (e.g., variants and invariants)—are stacked vertically, combining their records into a single file. This contrasts with a merge, which combines VCFs from different samples or individuals, aligning them by genomic position to produce a unified multisample VCF. 

This unified file enables consistent downstream analysis across all genomic positions, regardless of their variant status. 

To achieve this, a **concatenation** step is performed using `bcftools concat`, followed by **sorting** the concatenated VCF by genomic coordinates with `bcftools sort`. This ensures that the resulting file is both complete and properly ordered for compatibility with standard genomic analysis tools.

Below is the SLURM batch script used to perform these operations in parallel across chromosomes on an HPC cluster:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=01:00:00
#SBATCH --job-name=cctABAG
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-5
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Retrieve the chromosome name from the list using the SLURM array index
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chroms_list)

# Define the output directory
OUTPUT_DIR=results/05_Variants
LOG_DIR="results/11_Reports/concat"
mkdir -p ${LOG_DIR}

# Load bcftools module
module load bcftools/1.15.1

# Index both variant and invariant VCFs
bcftools index --threads 8 ${OUTPUT_DIR}/AB.AG.HC.acc.VARIANTS.flag.pass.gq20.${CHROM}.vcf.gz
bcftools index --threads 8 ${OUTPUT_DIR}/AB.AG.HC.acc.INVARIANTS.flag.pass.${CHROM}.vcf.gz

# Concatenate the two VCF files (variants and invariants)
bcftools concat --threads 8 \
        ${OUTPUT_DIR}/AB.AG.HC.acc.VARIANTS.flag.pass.gq20.${CHROM}.vcf.gz \
        ${OUTPUT_DIR}/AB.AG.HC.acc.INVARIANTS.flag.pass.${CHROM}.vcf.gz \
        --allow-overlaps \
        -Oz -o ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.${CHROM}.vcf.gz \
        2>&1 >> ${LOG_DIR}/AB.AG.HC.acc.pass.gq20.${CHROM}.log

# Sort the concatenated VCF by genomic position
bcftools sort -T tmp/ \
        ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.${CHROM}.vcf.gz \
        -Oz -o ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.sort.${CHROM}.vcf.gz \
        2>&1 >> ${LOG_DIR}/AB.AG.HC.acc.pass.gq20.sort.${CHROM}.log

# Clean up the unsorted merged file
rm ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.${CHROM}.vcf.gz

# Index the final sorted VCF
bcftools index --threads 8 ${OUTPUT_DIR}/AB.AG.HC.acc.pass.gq20.sort.${CHROM}.vcf.gz
```

