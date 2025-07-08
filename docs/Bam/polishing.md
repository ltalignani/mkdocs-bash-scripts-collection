# Polishing of BAM Files

<h2 class="no-toc">Table of Content</h2>

[TOC]

After mapping reads to a reference genome and marking duplicates, polishing the resulting BAM files is essential to ensure the integrity of downstream analyses. This step resolves specific errors often detected by tools such as Picard ValidateSamFile, including:  
- **ERROR:INVALID_TAG_NM**  
- **ERROR:MATE_NOT_FOUND**  

This section describes a two-step polishing pipeline that corrects these errors using samtools calmd and picard FixMateInformation.


## Step 1 — Fixing NM and MD Tags with samtools calmd

samtools calmd recomputes the NM (edit distance) and MD (mismatch positions) tags based on the reference genome. These tags are crucial for tools that perform variant calling or error detection, such as Picard.

We apply this correction to all BAM files with the suffix *marked.bam:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=23:00:00
#SBATCH --job-name=calmd
#SBATCH -p fast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=1G
#SBATCH -o Cluster_logs/%x-%j-%N.out
#SBATCH -e Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Treat and correct ERROR:INVALID_TAG_NM found by Picard ValidateSamFile

# Define path to .bam files
BAM_FOLDER="results/04_Polishing"

# Check if directory exists
if [ ! -d "$BAM_FOLDER" ]; then
    echo "Le dossier $BAM_FOLDER n'existe pas."
    exit 1
fi

# Create directory for log files
LOG_FOLDER="results/11_Reports/samtools"
mkdir -p "$LOG_FOLDER"

# Iterate through all .bam files and submit job
for BAM_FILE in "$BAM_FOLDER"/*marked.bam; do
    # Get file name without extension
    FILENAME_NO_EXT=$(basename "${BAM_FILE%.bam}")

    # Submit job for each file
    sbatch <<EOF
#!/bin/bash
#SBATCH -A invalbo
#SBATCH -p fast
#SBATCH -N 1
#SBATCH --time=00:30:00
#SBATCH --job-name=calmd_${FILENAME_NO_EXT}
#SBATCH --output=Cluster_logs/calmd_${FILENAME_NO_EXT}_%j.out
#SBATCH --error=Cluster_logs/calmd_${FILENAME_NO_EXT}_%j.err
#SBATCH --cpus-per-task 4
#SBATCH --mem=2G
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL

module load samtools/1.15.1

samtools calmd -@ 4 -b \
    "$BAM_FILE" \
    resources/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa \
    > "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".nm.bam \
    2> "${LOG_FOLDER}"/"${FILENAME_NO_EXT}".calmd.log

EOF

done
```

This script:  
- **Iterates** over all *marked.bam files  
- **Runs** samtools calmd to fix NM/MD tags  
- **Writes** a .nm.bam output file and logs errors  



## Step 2 — Fixing Mate Pair Information with Picard

To resolve MATE_NOT_FOUND errors, we use Picard FixMateInformation. However, this tool requires the input BAM file to be sorted by queryname. The process is as follows:

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=00:02:00
#SBATCH --job-name=fixmate
#SBATCH -p fast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=1G
#SBATCH -o Cluster_logs/%x-%j-%N.out
#SBATCH -e Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Treat and correct ERROR:MATE_NOT_FOUND found by Picard ValidateSamFile
# Use Picard ValidateSamFile then samtools calmd first.


# Define path to .bam files
BAM_FOLDER="results/04_Polishing"

# Check if directory exists
if [ ! -d "$BAM_FOLDER" ]; then
    echo "Le dossier $BAM_FOLDER n'existe pas."
    exit 1
fi

# Create directory for log files
LOG_FOLDER="results/11_Reports/samtools"
mkdir -p "$LOG_FOLDER"

# Iterate through all .bam files and submit job
for BAM_FILE in "$BAM_FOLDER"/*marked.nm.bam; do
    # Get file name without extension
    FILENAME_NO_EXT=$(basename "${BAM_FILE%.nm.bam}")

    # Submit job for each file
    sbatch <<EOF
#!/bin/bash
#SBATCH -A invalbo
#SBATCH -p fast
#SBATCH -N 1
#SBATCH --time=01:00:00
#SBATCH --job-name=fixmate_${FILENAME_NO_EXT}
#SBATCH --output=Cluster_logs/fixmate_${FILENAME_NO_EXT}_%j.out
#SBATCH --error=Cluster_logs/fixmate_${FILENAME_NO_EXT}_%j.err
#SBATCH --cpus-per-task 4
#SBATCH --mem=2G
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL

module load picard/2.23.5

# Sort by queryname
picard SortSam \
        -I "${BAM_FILE}" \
        -O "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".sorted_by_queryname.bam \
        --SORT_ORDER queryname --TMP_DIR tmp/ > "${LOG_FOLDER}"/"${FILENAME_NO_EXT}".sorted_by_queryname.log 2>&1

# Apply fixmate correction
picard FixMateInformation \
        -I "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".sorted_by_queryname.bam \
        -O "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".fixed_mates.bam \
        --TMP_DIR tmp/ > "${LOG_FOLDER}"/"${FILENAME_NO_EXT}".fixed_mates.log 2>&1

# Sort by coordinates
picard SortSam \
        -I "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".fixed_mates.bam \
        -O "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".fixed.bam \
        --SORT_ORDER coordinate --CREATE_INDEX true \
        --TMP_DIR tmp/ > "${LOG_FOLDER}"/"${FILENAME_NO_EXT}".sorted_by_coordinates.log 2>&1

rm -rf "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".sorted_by_queryname.bam
rm -rf "${BAM_FOLDER}"/"${FILENAME_NO_EXT}".fixed_mates.bam

EOF

done
```

This script performs:  
	1.	Sorting by queryname  
	2.	Running FixMateInformation to link properly paired reads  
	3.	Sorting again by coordinate (required by GATK, IGV, etc.)  
	4.	Indexing the final BAM  
	5.	Cleaning up intermediate files  


## Final Output

After these two polishing steps, each sample has a clean and valid BAM file:

results/04_Polishing/Sample_1.fixed.bam

These files:  
- **Are free of common structural and tag-based errors**  
- **Contain updated NM and MD fields**  
- **Have correctly paired reads**
- **Are sorted and indexed**

They are now ready for high-quality downstream analysis such as variant calling, coverage analysis, or visualization in genome browsers.

