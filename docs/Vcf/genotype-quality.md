# Filtering VCF Files by Genotype Quality (GQ)

## Genotype Quality (GQ)
Genotype Quality (`GQ`) is a crucial metric in variant calling that estimates the confidence in the genotype assigned to a sample at a given site. GQ values are phred-scaled, meaning that higher values reflect a lower probability of genotyping errors. Filtering variants based on GQ allows for the exclusion of low-confidence genotype calls, which helps improve the reliability of downstream analyses such as population genetics or genotype-phenotype associations.

In a VCF file, the GQ value is found in the sample-specific `FORMAT` fields, typically following the GT (Genotype) field. Here’s an example of a simplified VCF line with the GQ value highlighted:

```
#CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO                    FORMAT            Sample1
2L       123456  .       A       G       99      PASS    DP=520;AF=0.5           GT:GQ:DP          0/1:35:520
```

In this example:  
`The FORMAT field is GT:GQ:DP`
: indicating that each sample entry will provide values for Genotype (GT), Genotype Quality (GQ), and Read Depth (DP).  


`For Sample1, the corresponding values are 0/1:35:520, which means:`  
: GT = 0/1 (heterozygous),  
: GQ = 35 → this is the genotype quality,  
: DP = 520 (number of reads covering the site).  

Thus, the GQ value (35) indicates a reasonably high confidence in the assigned heterozygous genotype at that position.


## Filtering a VCF File by GQ 
In this project, a GQ threshold of 20 was applied to retain only genotypes with a high confidence level (i.e., an error probability of less than 1%). The filtering was performed using `VCFtools`, a widely used utility for processing VCF files. The script below illustrates how this filtering step was implemented for each chromosome on a high-performance computing cluster using SLURM:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=00:45:00
#SBATCH --job-name=GQABAG
#SBATCH -p fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=2
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

start_time=$(date +%s)

# Récupération du nom du fichier BAM à partir du fichier de liste
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/chroms_list)

# Répertoire de sortie pour les fichiers
OUTPUT_DIR=results/05_Variants

LOG_DIR="results/11_Reports/gq20"
mkdir -p ${LOG_DIR}

# Load module
module load vcftools/0.1.16 htslib/1.14

vcftools --gzvcf ${OUTPUT_DIR}/AB.AG.HC.acc.VARIANTS.flag.pass.${CHROM}.vcf.gz --minGQ 20 --recode --stdout | bgzip -c > ${OUTPUT_DIR}/AB.AG.HC.acc.VARIANTS.flag.pass.gq20.${CHROM}.vcf.gz


end_time=$(date +%s)  # Capture l'heure de fin
duration=$((end_time - start_time))  # Calcule la durée
echo "The script completed successfully in $((duration / 60)) minutes and $((duration % 60)) seconds."
```

## Filtering Criteria Analysis

- **--gzvcf**: Specifies the path to the input compressed VCF file (.vcf.gz). Here, it processes a variant file specific to one chromosome (${CHROM}), already filtered by flags (e.g., DP, MQ, FS).
- **--minGQ 20**: Filters out genotypes with GQ < 20, retaining only high-confidence genotypes (error probability < 1%).
- **--recode**: Instructs VCFtools to reconstruct a valid VCF file with only the remaining entries after filtering.
- **--stdout**: Sends the output to standard output, allowing piping to another program.
- **bgzip -c**: Compresses the output using bgzip (a block compression tool compatible with tabix indexing).
- **> ...vcf.gz**: Redirects the compressed output to a new VCF file, named to reflect that it includes only variants with GQ ≥ 20.