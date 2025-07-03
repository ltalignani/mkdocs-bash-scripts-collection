# BAM File Quality Control

BAM (Binary Alignment/Map) files are a key output of read mapping processes, and ensuring their integrity is crucial before proceeding to downstream analyses such as variant calling or transcript quantification. Several tools are available to perform different types of quality checks on BAM files. This section covers the most commonly used tools and strategies, along with their associated scripts in an HPC environment using SLURM.

## Quick Validation with samtools quickcheck

samtools quickcheck performs a fast sanity check of BAM or CRAM files. It ensures that each file:  
- Is readable  
- Has the correct format  
- Is not truncated  

It does not check file content consistency or alignment logic, but is ideal for catching incomplete or corrupted files early in a pipeline.

This validation step can be performed just after conversion in Bam format of a sam file, following the mapping/sorting of reads, or after the deduplication step.

A list of Bam files have to be created first:

```bash linenums="1"
#!/bin/bash

# ls -1: display only one file per line, making the output easier to read.
ls -1 results/02_Mapping | grep -E '_sorted\.bam' | sed 's/_sorted\.bam//g' > info_files/bam_list
```

Then, we can run `samtools quickcheck`:
```bash linenums="1"
#!/bin/bash

module load samtools/1.15.1

BAM_LIST="info_files/bam.list"

while IFS= read -r bam_file; do
    echo "Checking $bam_file"
    samtools quickcheck -v "$bam_file"
    if [ $? -ne 0 ]; then
        echo "Error found in $bam_file"
    else
        echo "$bam_file is OK"
    fi
done < "$BAM_LIST"
```

This script loops through a list of BAM files and runs samtools quickcheck -v on each of them, reporting any issues.


## In-depth Quality Reports with qualimap

### Multi-sample Overview: qualimap multi-bamqc

qualimap multi-bamqc is used to summarize multiple BAM files in a single integrated report. It gives a global overview of:  
- Read coverage distribution  
- GC bias  
- Mapping quality  
- Insert size distribution  

It needs an input list in a specific format, which can be created with this script:

```bash linenums="1"
#!/bin/bash

# .bam directory
bam_folder="results/04_Polishing"

# qualimap multi-bamqc output file
output_file="info_files/qualimap_sample_file"

# Pattern to be removed form filename
pattern="_merged_marked.bam"

# Check if directory exists
if [ ! -d "$bam_folder" ]; then
    echo "The directory $bam_folder does not exist."
    exit 1
fi

# Creates the qualimap multi-bamqc input file
for bam_file in "$bam_folder"/*merged_marked.bam; do
    # Remove pattern
    sample_name=$(basename "${bam_file//$pattern}")

    # Add pathway
    echo -e "$sample_name\t$bam_file" >> "$output_file"
done

echo "The file has been created : $output_file"
```

Then Qualimap multi-bamqc can be performed:

```bash linenums="1"
#!/bin/bash
#SBATCH --job-name=qualimap_multibamqc
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
...

OUTPUT_DIR="qc/qualimap_multibamqc/"
mkdir -p ${OUTPUT_DIR}

module load qualimap/2.2.2b

qualimap multi-bamqc -c -d info_files/qualimap_sample_file --run-bamqc \
    -outdir ${OUTPUT_DIR} -nr 10000 \
    --java-mem-size=32G -r
```

This script assumes that the file qualimap_sample_file contains paths and labels of all BAM files to include in the report.

### Per-sample Reports: qualimap bamqc

When detailed metrics are needed for each sample individually, qualimap bamqc is the tool of choice.

```bash linenums="1"
#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
...

BAM_FOLDER="results/04_Polishing"
OUTPUT_FOLDER="qc/qualimap"
LOG_FOLDER="results/11_Reports/qualimap"
mkdir -p "${OUTPUT_FOLDER}" "${LOG_FOLDER}"

for BAM_FILE in "${BAM_FOLDER}"/*_marked.bam; do
    FILENAME_NO_EXT=$(basename "${BAM_FILE%_marked.bam}")

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=qlmp_${FILENAME_NO_EXT}
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

module load qualimap/2.2.2b

qualimap bamqc -bam "${BAM_FILE}" -c -nt 16 \
    --java-mem-size=64G \
    -outdir "${OUTPUT_FOLDER}/${FILENAME_NO_EXT}" \
    &> "${LOG_FOLDER}/qualimap_${FILENAME_NO_EXT}.log"
EOF
done
```

Key Differences:

| Feature  | multi-bamqc                         | bamqc                       |
| :------- | :---------------------------------  | :-------------------------  |
| Input    | Many BAMs via a sample file         | Single BAM file             |
| Output   | Single HTML report with comparisons | One report per BAM          |
| Use case | Project-level summary               | Per-sample diagnostics      |
| Runtime  | Faster overall                      | Can be heavy for large BAMs |


## Read Statistics with samtools stats

samtools stats generates a comprehensive textual summary of alignment metrics such as:
	•	Total reads
	•	Mapped reads
	•	Insert sizes
	•	Base composition
	•	Error rates

```bash
#!/bin/bash
#SBATCH --job-name=stats
#SBATCH --cpus-per-task=16
...

BAM_FOLDER="results/04_Polishing"
OUTPUT_FOLDER="qc/samtools"
LOG_FOLDER="results/11_Reports/samtools"
mkdir -p "$OUTPUT_FOLDER" "$LOG_FOLDER"

for BAM_FILE in "$BAM_FOLDER"/*marked.bam; do
    FILENAME_NO_EXT=$(basename "${BAM_FILE%.bam}")

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=stats_${FILENAME_NO_EXT}
#SBATCH --cpus-per-task=16
...

module load samtools/1.15.1

samtools stats -@ 16 \
    -r resources/genomes/GCF_035046485.1_AalbF5_genomic.fna \
    "$BAM_FILE" \
    1> "${OUTPUT_FOLDER}/${FILENAME_NO_EXT}_stats.txt" \
    2> "${LOG_FOLDER}/${FILENAME_NO_EXT}.stats.log"
EOF
done
```

The result is a plain-text .txt file with hundreds of metrics, suitable for parsing or plotting with tools like plot-bamstats.


## Structural Validation with picard ValidateSamFile

To detect deeper issues in file structure and metadata tags, Picard ValidateSamFile is a reference tool. It checks:
- Invalid mate information  
- NM and MD tag consistency  
- Missing headers  
- Unexpected reference contigs  

```bash linenums="1"
#!/bin/bash
#SBATCH --job-name=validatesam
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
...

BAM_FOLDER="results/04_Polishing"
OUTPUT_FOLDER="qc/validatesam"
LOG_FOLDER="results/11_Reports/validatesamfiles"
mkdir -p "${OUTPUT_FOLDER}" "${LOG_FOLDER}"

for BAM_FILE in "${BAM_FOLDER}"/*marked.bam; do
    FILENAME_NO_EXT=$(basename "${BAM_FILE%.bam}")

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=validatesam_${FILENAME_NO_EXT}
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

module load picard/2.23.5

picard ValidateSamFile -Xmx32000m \
    -I "${BAM_FILE}" \
    -R resources/genomes/GCF_035046485.1_AalbF5_genomic.fna \
    -O "${OUTPUT_FOLDER}/${FILENAME_NO_EXT}_ValidateSam.txt" \
    -M SUMMARY \
    2> ${LOG_FOLDER}/${FILENAME_NO_EXT}_validate_bam.log
EOF
done
```

Output:  
- A text summary report of all found issues  
- A log file containing additional details and error traces  


## Summary Table

|         Tool         |              Purpose              |    Output Type    | Runtime Impact |
| :------------------: | :-------------------------------: | :---------------: | :------------: |
| samtools quickcheck  | Basic format and integrity check  | Terminal messages |   Very fast    |
|    qualimap bamqc    |       Per-sample QC metrics       |     HTML, PDF     |     Medium     |
| qualimap multi-bamqc |  Global QC overview of many BAMs  |       HTML        |      Fast      |
|    samtools stats    |        Detailed read stats        |       Text        |      Fast      |
|   ValidateSamFile    | Deep structure and tag validation |       Text        |     Medium     |
