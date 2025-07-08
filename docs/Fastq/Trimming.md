# Trimmomatic: Quality-Based Read Trimming

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Overview

`Trimmomatic` is a comprehensive read trimming tool for Illumina NGS data that implements flexible quality-based trimming algorithms. The software performs sequence quality assessment and trimming operations through a modular pipeline architecture, enabling researchers to remove low-quality bases, adapter sequences, and other artifacts that could compromise downstream analysis accuracy.

## Principle of Operation

`Trimmomatic` operates through a multi-step processing pipeline where each trimming operation is implemented as a discrete module. The tool processes reads sequentially through user-defined trimming steps, applying quality-based filtering and trimming algorithms to maximize both read quality and retention rates.

### Core Methodology

The trimming process is based on several key algorithmic approaches:

1. **Sliding Window Quality Trimming**: Implements a sliding window algorithm that evaluates average quality scores across defined window sizes
2. **Leading/Trailing Quality Trimming**: Removes low-quality bases from read termini based on quality thresholds
3. **Adaptive Adapter Trimming**: Identifies and removes adapter sequences through seed-based alignment algorithms
4. **Length-Based Filtering**: Applies minimum length thresholds to ensure read utility for downstream applications

## Technical Implementation

### Quality Assessment Algorithms

#### Sliding Window Trimming

The sliding window algorithm evaluates sequence quality through:  
- **Window Size Definition**: Configurable window sizes (typically 4-5 bases)  
- **Quality Threshold Application**: Average quality score thresholds within each window  
- **Progressive Evaluation**: Sequential window advancement along read length  
- **Optimal Cutting Point Identification**: Determination of trimming positions that maximize retained sequence quality  

#### Phred Quality Score Integration

Trimmomatic utilizes Phred quality scores for decision-making processes:  
- **Quality Score Interpretation**: Conversion of ASCII-encoded quality scores to probability values  
- **Threshold-Based Decisions**: Application of user-defined quality thresholds for trimming operations  
- **Statistical Quality Assessment**: Probabilistic evaluation of base-calling accuracy  

### Adapter Removal Mechanisms

#### Seed-Based Alignment

Adapter detection employs sophisticated alignment algorithms:  
- **Seed Identification**: Short exact matches between reads and adapter sequences  
- **Extension Algorithms**: Extension of seed matches to identify full adapter sequences  
- **Mismatch Tolerance**: Configurable mismatch thresholds for adapter identification  
- **Partial Adapter Detection**: Identification of adapter sequences at read termini  

#### Palindrome Mode

For paired-end sequencing, Trimmomatic implements palindrome detection:  
- **Insert Size Evaluation**: Assessment of library insert sizes relative to read length  
- **Complementary Sequence Detection**: Identification of adapter-adapter ligation events  
- **Read-Through Artifact Removal**: Elimination of sequences resulting from short insert libraries  

### Processing Modes

#### Single-End Processing

Single-end mode implements:  
- Quality-based trimming algorithms  
- Adapter sequence removal  
- Length-based filtering  
- Output quality assessment  

#### Paired-End Processing

Paired-end mode maintains read pair integrity through:  
- **Synchronous Processing**: Simultaneous processing of forward and reverse reads  
- **Pair Integrity Maintenance**: Preservation of read pairing relationships  
- **Orphan Read Handling**: Management of reads whose pairs fail quality thresholds  
- **Output Stream Segregation**: Separation of paired and unpaired reads post-processing  

## Trimming Operations

### ILLUMINACLIP

Adapter and quality trimming operation with parameters:  
- **Adapter File**: FASTA file containing adapter sequences  
- **Seed Mismatches**: Number of mismatches allowed in seed region  
- **Palindrome Clip Threshold**: Accuracy threshold for palindrome matches  
- **Simple Clip Threshold**: Accuracy threshold for simple adapter matches  

### SLIDINGWINDOW

Quality-based trimming using sliding window approach:  
- **Window Size**: Number of bases in sliding window  
- **Required Quality**: Average quality threshold within window  
- **Trimming Strategy**: 3' end trimming based on quality degradation  

### LEADING/TRAILING

Terminal base removal based on quality thresholds:  
- **LEADING**: Removal of low-quality bases from 5' end  
- **TRAILING**: Removal of low-quality bases from 3' end  
- **Quality Threshold**: Minimum acceptable quality score  

### MINLEN

Length-based filtering:  
- **Minimum Length**: Shortest acceptable read length post-trimming  
- **Read Retention**: Ensures reads maintain sufficient length for alignment  

### CROP/HEADCROP

Fixed-length trimming operations:  
- **CROP**: Trimming reads to specified maximum length  
- **HEADCROP**: Removal of specified number of bases from 5' end  

## Command Line Usage

### Basic Commands

#### Single-End Trimming
```bash
java -jar trimmomatic-0.39.jar SE -phred33 \
    input.fastq.gz output_trimmed.fastq.gz \
    ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

#### Paired-End Trimming
```bash
java -jar trimmomatic-0.39.jar PE -phred33 \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
    output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

### Advanced Usage Examples

#### High-Quality Trimming
```bash
# Stringent quality trimming for sensitive applications
java -jar trimmomatic-0.39.jar PE -phred33 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_paired.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_paired.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
    LEADING:5 TRAILING:5 \
    SLIDINGWINDOW:4:20 \
    MINLEN:50
```

#### RNA-seq Specific Trimming
```bash
# Optimized for RNA-seq with stranded libraries
java -jar trimmomatic-0.39.jar PE -phred33 \
    RNAseq_R1.fastq.gz RNAseq_R2.fastq.gz \
    RNAseq_R1_paired.fastq.gz RNAseq_R1_unpaired.fastq.gz \
    RNAseq_R2_paired.fastq.gz RNAseq_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:25
```

#### Memory-Optimized Processing
```bash
# Large file processing with memory optimization
java -Xmx4g -jar trimmomatic-0.39.jar PE -phred33 \
    -threads 8 \
    large_R1.fastq.gz large_R2.fastq.gz \
    large_R1_paired.fastq.gz large_R1_unpaired.fastq.gz \
    large_R2_paired.fastq.gz large_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

### Specialized Trimming Strategies

#### Nextera Adapter Removal
```bash
# Nextera library preparation adapter trimming
java -jar trimmomatic-0.39.jar PE -phred33 \
    nextera_R1.fastq.gz nextera_R2.fastq.gz \
    nextera_R1_paired.fastq.gz nextera_R1_unpaired.fastq.gz \
    nextera_R2_paired.fastq.gz nextera_R2_unpaired.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

#### Custom Adapter Trimming
```bash
# Custom adapter sequences
java -jar trimmomatic-0.39.jar PE -phred33 \
    custom_R1.fastq.gz custom_R2.fastq.gz \
    custom_R1_paired.fastq.gz custom_R1_unpaired.fastq.gz \
    custom_R2_paired.fastq.gz custom_R2_unpaired.fastq.gz \
    ILLUMINACLIP:custom_adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

#### Fixed-Length Trimming
```bash
# Uniform read length generation
java -jar trimmomatic-0.39.jar PE -phred33 \
    variable_R1.fastq.gz variable_R2.fastq.gz \
    fixed_R1_paired.fastq.gz fixed_R1_unpaired.fastq.gz \
    fixed_R2_paired.fastq.gz fixed_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    HEADCROP:15 \
    CROP:100 \
    LEADING:3 TRAILING:3 \
    MINLEN:75
```

### Batch Processing Scripts

#### Shell Script for Multiple Samples
```bash
#!/bin/bash
# Batch trimming script for multiple paired-end samples

TRIMMOMATIC_JAR="/path/to/trimmomatic-0.39.jar"
ADAPTER_FILE="/path/to/TruSeq3-PE.fa"
INPUT_DIR="/path/to/raw_fastq/"
OUTPUT_DIR="/path/to/trimmed_fastq/"

mkdir -p $OUTPUT_DIR

for r1_file in ${INPUT_DIR}*_R1.fastq.gz; do
    sample_name=$(basename $r1_file _R1.fastq.gz)
    r2_file=${INPUT_DIR}${sample_name}_R2.fastq.gz
    
    echo "Processing sample: $sample_name"
    
    java -jar $TRIMMOMATIC_JAR PE -phred33 -threads 4 \
        $r1_file $r2_file \
        ${OUTPUT_DIR}${sample_name}_R1_paired.fastq.gz \
        ${OUTPUT_DIR}${sample_name}_R1_unpaired.fastq.gz \
        ${OUTPUT_DIR}${sample_name}_R2_paired.fastq.gz \
        ${OUTPUT_DIR}${sample_name}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done
```

#### Quality Control Integration
```bash
#!/bin/bash
# Combined quality assessment and trimming pipeline

SAMPLE_PREFIX="sample"
FASTQC_OUTPUT="/path/to/fastqc_output/"
TRIMMED_OUTPUT="/path/to/trimmed_output/"

# Pre-trimming quality assessment
fastqc -o $FASTQC_OUTPUT ${SAMPLE_PREFIX}_R1.fastq.gz ${SAMPLE_PREFIX}_R2.fastq.gz

# Trimming operation
java -jar trimmomatic-0.39.jar PE -phred33 \
    ${SAMPLE_PREFIX}_R1.fastq.gz ${SAMPLE_PREFIX}_R2.fastq.gz \
    ${TRIMMED_OUTPUT}${SAMPLE_PREFIX}_R1_paired.fastq.gz \
    ${TRIMMED_OUTPUT}${SAMPLE_PREFIX}_R1_unpaired.fastq.gz \
    ${TRIMMED_OUTPUT}${SAMPLE_PREFIX}_R2_paired.fastq.gz \
    ${TRIMMED_OUTPUT}${SAMPLE_PREFIX}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

# Post-trimming quality assessment
fastqc -o $FASTQC_OUTPUT \
    ${TRIMMED_OUTPUT}${SAMPLE_PREFIX}_R1_paired.fastq.gz \
    ${TRIMMED_OUTPUT}${SAMPLE_PREFIX}_R2_paired.fastq.gz
```

### Batch Processing On a HPC Cluster

```bash linenums="1"
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=2-23:00:00
#SBATCH --job-name=trim
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 8
#SBATCH --mem=16G
#SBATCH --array 1-56
#SBATCH -o Cluster_logs/%x-%j-%N.out
#SBATCH -e Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Recover of fastq file name
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/fastq_files)

# Create OUTPUT directory
OUTPUT_FOLDER="results/01_Trimming"
mkdir -p "$OUTPUT_FOLDER"

# Create directory for log files
LOG_FOLDER="results/11_Reports"
mkdir -p "$LOG_FOLDER"

module load trimmomatic/0.39 fastqc/0.12.1

trimmomatic PE -threads 16 \
    -phred33 -trimlog "$LOG_FOLDER/trimmomatic_${SAMPLE}.log" \
    "raw/${SAMPLE}_R1.fastq.gz" \
    "raw/${SAMPLE}_R2.fastq.gz" \
    "$OUTPUT_FOLDER/${SAMPLE}_trimmomatic_R1.fastq.gz" \
    "$OUTPUT_FOLDER/${SAMPLE}_trimmomatic_unpaired_R1.fastq.gz" \
    "$OUTPUT_FOLDER/${SAMPLE}_trimmomatic_R2.fastq.gz" \
    "$OUTPUT_FOLDER/${SAMPLE}_trimmomatic_unpaired_R2.fastq.gz" \
    ILLUMINACLIP:resources/adapters/TruSeq2-PE.fa:2:30:15 LEADING:20 TRAILING:3 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:50

fastqc --threads 16 -o qc/fastqc-post-trim/ "$OUTPUT_FOLDER/${SAMPLE}_trimmomatic_R1.fastq.gz"
fastqc --threads 16 -o qc/fastqc-post-trim/ "$OUTPUT_FOLDER/${SAMPLE}_trimmomatic_R2.fastq.gz"
```

## Performance Considerations

### Computational Requirements

#### Memory Usage
- **Single-End**: Minimal memory requirements (typically <1GB)  
- **Paired-End**: Moderate memory usage scaling with file size  
- **Large Files**: Memory optimization through streaming algorithms  

#### Threading Support
- **Multi-threading**: Parallel processing capability for improved performance  
- **I/O Optimization**: Efficient file handling for compressed formats  
- **Scalability**: Linear performance scaling with thread count  

### Algorithm Efficiency

#### Time Complexity
- **Linear Processing**: O(n) time complexity relative to input size  
- **Sliding Window**: Constant time window evaluation  
- **Adapter Matching**: Efficient seed-based alignment algorithms  

## Output Interpretation

### Quality Metrics

#### Trimming Statistics
- **Input Read Count**: Total reads processed  
- **Output Read Count**: Reads surviving quality filters  
- **Trimming Rates**: Percentage of reads requiring trimming  
- **Length Distributions**: Pre- and post-trimming read length profiles  

#### Quality Improvements
- **Quality Score Distributions**: Comparison of quality profiles  
- **Adapter Contamination**: Reduction in adapter sequence content  
- **Overall Quality Enhancement**: Improvement in dataset quality metrics  

## Best Practices

### Parameter Optimization

#### Quality Thresholds
- **Conservative Approach**: Higher quality thresholds for critical applications  
- **Balanced Strategy**: Moderate thresholds preserving read count and quality  
- **Application-Specific**: Parameter adjustment based on downstream analysis requirements  

#### Adapter Trimming
- **Accurate Adapter Files**: Use of appropriate adapter sequences for library preparation method  
- **Stringency Levels**: Adjustment of mismatch tolerance based on sequence quality  
- **Palindrome Detection**: Optimization for paired-end library characteristics  

### Integration Strategies

#### Pipeline Integration
- **Quality Control Workflows**: Integration with FastQC and MultiQC  
- **Alignment Pipelines**: Preprocessing for genome and transcriptome alignment  
- **Variant Calling**: Quality optimization for SNP and indel detection  

## Limitations and Considerations

### Technical Limitations

#### Algorithm Constraints
- **Quality Score Dependency**: Reliance on accurate quality score calibration  
- **Adapter Database Requirements**: Need for comprehensive adapter sequence databases  
- **Short Read Optimization**: Designed primarily for short-read sequencing technologies  

#### Biological Considerations
- **Over-Trimming Risk**: Potential for excessive sequence removal  
- **Bias Introduction**: Possible introduction of length or quality biases  
- **Context Sensitivity**: Need for application-specific parameter optimization  

### Quality Control Validation

#### Post-Trimming Assessment
- **Quality Metrics Evaluation**: Comprehensive quality assessment post-trimming  
- **Read Length Distributions**: Analysis of length profiles and retention rates  
- **Downstream Impact**: Evaluation of trimming effects on alignment and analysis quality  

## Conclusion

Trimmomatic provides a robust and flexible framework for NGS read quality improvement through sophisticated trimming algorithms. Its modular architecture enables researchers to implement customized quality control strategies while maintaining high processing efficiency. When properly configured and validated, Trimmomatic significantly enhances downstream analysis quality by removing low-quality sequences and technical artifacts, ultimately contributing to more accurate and reliable genomics research outcomes.