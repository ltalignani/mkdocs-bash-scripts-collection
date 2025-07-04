# Fastq Contamination Chacking With Fastq-Screen
![Fastq-screen logo](../assets/FQS_logo.png)

## Overview

`Fastq-screen` is a quality control tool designed to screen libraries of short reads in Fastq format against a set of reference databases. The primary purpose of this tool is to identify potential contamination sources and assess the composition of sequencing libraries by aligning reads against multiple reference genomes or databases simultaneously.

## Principle of Operation

`Fastq-screen` operates by taking a subset of reads from input Fastq files and aligning them against user-defined reference databases using fast alignment algorithms. The tool provides quantitative assessment of how many reads align to each reference database, enabling researchers to identify potential contamination, confirm the expected organism composition, and detect unexpected sequences in their datasets.

### Core Methodology

The screening process follows these key steps:

1. **Read Sampling**: A representative subset of reads is extracted from the input FastQ file(s) to reduce computational overhead while maintaining statistical significance  
2. **Multi-Database Alignment**: Reads are aligned against multiple reference databases simultaneously using configurable alignment parameters  
3. **Classification**: Reads are classified based on their alignment patterns:  
   - Uniquely mapping reads (align to only one database)  
   - Multi-mapping reads (align to multiple databases)  
   - Unmapped reads (fail to align to any database)  
4. **Statistical Reporting**: Results are presented as percentages and counts for each reference database  

## Technical Implementation

### Alignment Strategy

FastQ Screen utilizes Bowtie2 as the default aligner, though it supports multiple alignment tools including:  
- Bowtie2 (default)  
- Bowtie  
- BWA  
- HISAT2  

The tool employs a mapping strategy optimized for contamination detection rather than comprehensive genome mapping, using parameters that favor sensitivity over specificity to capture potential contaminants.

### Database Configuration

Reference databases are configured through a configuration file that specifies:  
- Database paths and names  
- Aligner-specific parameters  
- Organism or contamination source labels  
- Custom alignment settings per database  

Common reference databases include:  
- Host organism genomes  
- Common laboratory contaminants (E. coli, yeast, mycoplasma)  
- Adapter sequences  
- Vector sequences  
- PhiX control sequences  

### Performance Optimization

To maintain computational efficiency, FastQ Screen implements several optimization strategies:  
- **Subset Analysis**: By default, only the first 100,000 reads are analyzed  
- **Parallel Processing**: Multi-threading support for simultaneous database screening  
- **Memory Management**: Efficient handling of large reference databases  
- **Configurable Parameters**: Adjustable alignment stringency and sampling parameters  

## Output and Interpretation

### Graphical Output

FastQ Screen generates publication-ready plots showing:  
- Percentage of reads mapping to each database  
- Stacked bar charts displaying unique vs. multi-mapping reads  
- Heat maps for multiple sample comparisons  

### Tabular Results

Detailed tabular output includes:  
- Read counts and percentages for each database  
- Mapping statistics (unique, multi-mapping, unmapped)  
- Quality metrics and alignment parameters used  

### Interpretation Guidelines

Results interpretation focuses on:  
- **Expected Alignments**: High percentage alignment to the target organism database  
- **Contamination Detection**: Unexpected alignments to bacterial, viral, or other databases  
- **Adapter Content**: Presence of sequencing adapters or technical sequences  
- **Cross-contamination**: Alignments to multiple related organisms  

## Applications in Genomics

### Quality Control Workflows

FastQ Screen serves as a critical component in NGS quality control pipelines by:  
- Validating sample identity and purity  
- Detecting cross-contamination between samples  
- Identifying technical artifacts and adapter sequences  
- Assessing library preparation quality  

### Contamination Source Identification

The tool effectively identifies various contamination sources:  
- **Bacterial Contamination**: Detection of common laboratory bacterial strains  
- **Human DNA Contamination**: Identification of human DNA in non-human samples  
- **Environmental Contamination**: Detection of environmental microorganisms  
- **Reagent Contamination**: Identification of contaminating sequences from reagents  

### Multi-Species Analysis

For complex samples containing multiple organisms:  
- Metagenomic sample composition assessment  
- Host-pathogen interaction studies  
- Symbiotic relationship analysis  
- Environmental sample characterization  

## Best Practices

### Database Selection

Optimal database selection should include:  
- Target organism genome(s)  
- Common laboratory contaminants  
- Potential environmental contaminants specific to sample type  
- Technical sequences (adapters, vectors, controls)  

### Parameter Optimization

Key parameters to consider:  
- Number of reads to screen (balance between accuracy and speed)  
- Alignment stringency settings  
- Multi-mapping tolerance thresholds  
- Minimum alignment length requirements  

### Integration with Pipelines

FastQ Screen integrates effectively with:  
- MultiQC for consolidated quality control reporting  
- Automated bioinformatics pipelines  
- High-throughput screening workflows  
- Custom quality control frameworks  

## Limitations and Considerations

### Technical Limitations

- **Database Dependency**: Results quality depends on comprehensive and current reference databases  
- **Short Read Specificity**: Optimized for short-read sequencing technologies  
- **Alignment Sensitivity**: May miss highly divergent sequences or novel organisms  
- **Computational Resources**: Memory requirements scale with database size  

### Biological Considerations

- **Phylogenetic Relationships**: Closely related organisms may show cross-mapping  
- **Repetitive Sequences**: Highly repetitive regions may cause multi-mapping artifacts  
- **Sequence Similarity**: Conservative interpretation needed for evolutionarily related species  
- **Database Completeness**: Absence of alignment doesn't guarantee absence of organism  

## Command Line Usage

### Basic Commands

#### Standard Single-End Screening
```bash linenums="1"
fastq_screen --conf /path/to/fastq_screen.conf sample.fastq.gz
```

#### Paired-End Screening
```bash linenums="1"
fastq_screen --conf /path/to/fastq_screen.conf sample_R1.fastq.gz sample_R2.fastq.gz
```

#### Multiple File Processing
```bash linenums="1"
fastq_screen --conf /path/to/fastq_screen.conf *.fastq.gz
```

### Advanced Usage Examples

#### Custom Subset Size
```bash linenums="1"
# Screen only the first 50,000 reads for faster processing
fastq_screen --conf /path/to/fastq_screen.conf --subset 50000 sample.fastq.gz
```

#### Multi-threading
```bash linenums="1"
# Use 8 threads for parallel processing
fastq_screen --conf /path/to/fastq_screen.conf --threads 8 sample.fastq.gz
```

#### Output Directory Specification
```bash linenums="1"
# Specify custom output directory
fastq_screen --conf /path/to/fastq_screen.conf --outdir /path/to/results/ sample.fastq.gz
```

#### Force Overwrite Existing Results
```bash linenums="1"
# Overwrite existing output files
fastq_screen --conf /path/to/fastq_screen.conf --force sample.fastq.gz
```

### Specialized Screening Options

#### Tag-based Output
```bash linenums="1"
# Add custom tag to output files
fastq_screen --conf /path/to/fastq_screen.conf --tag custom_label sample.fastq.gz
```

#### Aligner Selection
```bash linenums="1"
# Use BWA instead of default Bowtie2
fastq_screen --conf /path/to/fastq_screen.conf --aligner bwa sample.fastq.gz
```

#### Quiet Mode
```bash linenums="1"
# Suppress standard output messages
fastq_screen --conf /path/to/fastq_screen.conf --quiet sample.fastq.gz
```

### Configuration File Examples

#### Basic Configuration Structure
```bash linenums="1"
# Example fastq_screen.conf file content
DATABASE    Human    /path/to/human_genome/human_index    bowtie2
DATABASE    Mouse    /path/to/mouse_genome/mouse_index    bowtie2
DATABASE    E_coli   /path/to/ecoli_genome/ecoli_index    bowtie2
DATABASE    Yeast    /path/to/yeast_genome/yeast_index    bowtie2
DATABASE    PhiX     /path/to/phix_genome/phix_index      bowtie2
DATABASE    Adapters /path/to/adapters/adapter_index     bowtie2
```

#### Advanced Configuration with Custom Parameters
```bash linenums="1"
# Configuration with custom bowtie2 parameters
DATABASE    Human    /path/to/human_genome/human_index    bowtie2    --very-sensitive --end-to-end
DATABASE    Bacteria /path/to/bacteria_db/bacteria_index bowtie2    --local --very-fast
```

### Batch Processing Scripts

#### Shell Script for Multiple Samples

```bash linenums="1"
#!/bin/bash
# Batch processing script
CONF_FILE="/path/to/fastq_screen.conf"
INPUT_DIR="/path/to/fastq_files/"
OUTPUT_DIR="/path/to/results/"

for file in ${INPUT_DIR}*.fastq.gz; do
    echo "Processing $file"
    fastq_screen --conf $CONF_FILE --outdir $OUTPUT_DIR --threads 4 $file
done
```

#### Paired-End Batch Processing
```bash linenums="1"
#!/bin/bash
# Process paired-end files
CONF_FILE="/path/to/fastq_screen.conf"
INPUT_DIR="/path/to/fastq_files/"

for r1_file in ${INPUT_DIR}*_R1.fastq.gz; do
    r2_file=${r1_file/_R1/_R2}
    echo "Processing pair: $r1_file and $r2_file"
    fastq_screen --conf $CONF_FILE --threads 8 $r1_file $r2_file
done
```

### Fastq-screen on a HPC Cluster

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=2-23:00:00
#SBATCH --job-name=fastqscreen  
#SBATCH -p fast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 8
#SBATCH --mem=24G
#SBATCH --array 1-22
#SBATCH -o Cluster_logs/%x-%j-%N.out
#SBATCH -e Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Recover of fastq file name
SAMPLELANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/fastq_files)

# Locate INPUT directory 
INPUT_FOLDER="raw"

LOG_FOLDER="results/11_Reports/fastq-screen"

module load fastq-screen/0.15.3 bwa/0.7.17

fastq_screen -q --threads 8 \
            --conf info_files/fastq-screen.conf \
            --aligner bwa \
            --subset 1000 \
            --outdir qc/fastq-screen/ \
            ${INPUT_FOLDER}/*.fastq.gz \
            &> ${LOG_FOLDER}/${SAMPLELANE}-fastq-screen.log
```


### Integration with Quality Control Pipelines

#### Combined with FastQC
```bash linenums="1"
# Run FastQC and FastQ Screen sequentially
fastqc sample.fastq.gz
fastq_screen --conf /path/to/fastq_screen.conf sample.fastq.gz
```

#### Integration with MultiQC
```bash linenums="1"
# Generate FastQ Screen reports for MultiQC compilation
fastq_screen --conf /path/to/fastq_screen.conf *.fastq.gz
multiqc .
```

### Troubleshooting Commands

#### Test Configuration
```bash linenums="1"
# Test configuration file validity
fastq_screen --conf /path/to/fastq_screen.conf --help
```

#### Verbose Output for Debugging
```bash linenums="1"
# Enable verbose output for troubleshooting
fastq_screen --conf /path/to/fastq_screen.conf --verbose sample.fastq.gz
```

#### Version Information
```bash linenums="1"
# Check FastQ Screen version
fastq_screen --version
```

## Conclusion

FastQ Screen represents an essential tool in the NGS quality control arsenal, providing rapid and comprehensive contamination screening capabilities. Its ability to simultaneously screen against multiple reference databases makes it invaluable for maintaining data quality and integrity in genomics research. When properly configured and interpreted, FastQ Screen enables researchers to confidently assess the composition and quality of their sequencing libraries, ultimately contributing to more reliable and reproducible genomics analyses.