# SAM/BAM File Sorting

## Overview

SAM/BAM file sorting is a fundamental preprocessing step in genomic analysis pipelines that organizes alignment records according to specific criteria. Proper sorting is essential for downstream applications including variant calling, coverage analysis, and visualization tools. The sorting process optimizes data access patterns and enables efficient indexing for random access operations.

## Sorting Strategies

### Coordinate-Based Sorting

Coordinate sorting arranges alignment records by their genomic position, following the order:  
1. **Reference sequence**: Sorted by reference sequence name (typically chromosomal order)  
2. **Position**: Sorted by leftmost mapping position (POS field)  
3. **Strand orientation**: Forward strand alignments before reverse strand at the same position  

This sorting order is mandatory for most downstream applications and enables:  
- Efficient variant calling algorithms  
- Proper mate-pair processing  
- Optimized coverage calculations  
- Compatible indexing with BAI/CSI formats  

### Query Name Sorting

Query name sorting organizes records alphabetically by read identifier (QNAME field). This sorting strategy is particularly useful for:  
- Mate-pair analysis and validation  
- Duplicate detection algorithms  
- Quality control assessments  
- Conversion to FASTQ format  

### Unsorted Files

Unsorted BAM files maintain the original alignment order, typically reflecting the input sequence order. While computationally faster to generate, unsorted files have limited utility for most analytical applications.

## Sorting Algorithms and Performance

### Memory-Based Sorting

Modern sorting implementations utilize in-memory algorithms for optimal performance:  
- **Quicksort variants**: Efficient for smaller datasets fitting in available RAM  
- **Merge sort**: Stable sorting with predictable O(n log n) performance  
- **Hybrid approaches**: Combine multiple algorithms based on data characteristics  

### External Sorting

For datasets exceeding available memory, external sorting algorithms are employed:  
- **Multi-way merge**: Divides data into memory-sized chunks  
- **Temporary file management**: Utilizes disk storage for intermediate results  
- **I/O optimization**: Minimizes disk access through efficient buffering  

## Implementation Tools

### SAMtools

SAMtools provides the standard implementation for BAM sorting:

```bash
# Coordinate sorting
samtools sort -o output_sorted.bam input.bam

# Query name sorting
samtools sort -n -o output_qname_sorted.bam input.bam

# Memory optimization
samtools sort -m 2G -@ 8 -o output_sorted.bam input.bam
```

### Picard Tools

Picard SortSam offers advanced sorting capabilities with extensive configuration options:

```bash
picard SortSam \
    INPUT=input.sam \
    OUTPUT=output_sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
```

### GNU Sort

For specialized applications, GNU sort can process SAM files directly:

```bash
# Header-aware coordinate sorting
(samtools view -H input.sam; samtools view input.sam | \
sort -k3,3V -k4,4n) | samtools view -b > output_sorted.bam
```

## Computational Considerations

### Memory Requirements

Memory allocation significantly impacts sorting performance:  
- **Minimum requirements**: 2-4 GB for typical whole-genome datasets  
- **Optimal allocation**: 8-16 GB enables in-memory sorting for most applications  
- **Memory scaling**: Linear relationship between file size and optimal memory allocation  

### Parallelization

Modern sorting implementations support parallel processing:  
- **Thread-level parallelism**: Multiple threads for concurrent sorting operations  
- **Process-level parallelism**: Independent sorting of file segments  
- **Distributed computing**: Cluster-based sorting for extremely large datasets  

### Storage Considerations

Temporary storage requirements during sorting:  
- **Disk space**: 2-3x input file size for intermediate files  
- **I/O bandwidth**: High-speed storage improves sorting performance  
- **Network considerations**: Minimize network I/O during sorting operations  

## Quality Control and Validation

### Sort Order Verification

Validation ensures proper sorting implementation:  
- **Coordinate validation**: Verify monotonic position ordering within chromosomes  
- **Reference sequence order**: Confirm proper chromosomal ordering  
- **Flag consistency**: Validate mate-pair relationships in coordinate-sorted files  

### Performance Monitoring

Key metrics for sorting performance assessment:  
- **Processing time**: Total wall-clock time for sorting operation  
- **Memory utilization**: Peak and average memory usage  
- **I/O statistics**: Read/write operations and throughput  
- **Error rates**: Validation of sorting correctness  

## Practical Implementation Example

The following example demonstrates a production-ready BAM sorting pipeline using SLURM job scheduling:

```bash
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=2-23:00:00
#SBATCH --job-name=sortsam
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH --mem=30G
#SBATCH --array 1-56
#SBATCH -o Cluster_logs/%x-%j-%N.out
#SBATCH -e Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Extract sample identifier from job array
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/fastq_files)
echo ${SAMPLE}

# Load required software module
module load picard/2.23.5

# Execute coordinate-based sorting with indexing
picard SortSam -Xmx40G \
  -I results/02_Mapping/${SAMPLE}.sam \
  -O results/02_Mapping/${SAMPLE}_sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  2> results/11_Reports/sortsam/${SAMPLE}-SortSam.out
```

### Script Components Analysis

**Resource Allocation:**
- **Memory**: 30 GB SLURM allocation with 40 GB Java Virtual Machine (JVM) heap  
- **CPU**: 4 cores per task for parallel processing  
- **Time**: Extended 2-day limit for large datasets  

**Processing Parameters:**
- **Sort Order**: Coordinate-based sorting for downstream compatibility  
- **Index Creation**: Automatic BAI index generation  
- **Memory Management**: 500,000 records in RAM for optimal performance  

**Array Processing:**
- **Parallel Execution**: 56 simultaneous sorting jobs  
- **Sample Management**: Dynamic sample selection from input file  
- **Error Handling**: Comprehensive logging and email notifications  

## Best Practices

### Pre-Sorting Considerations

1. **Input Validation**: Verify SAM/BAM file integrity (with samtools quickcheck) before sorting
2. **Disk Space**: Ensure adequate temporary storage (3x input file size)
3. **Memory Planning**: Allocate sufficient RAM to avoid external sorting
4. **Backup Strategy**: Preserve original files until validation completion

### Optimization Strategies

1. **Memory Tuning**: Adjust heap size based on available system memory
2. **Parallel Processing**: Utilize multiple CPU cores for enhanced performance
3. **I/O Optimization**: Use high-speed storage for temporary files
4. **Compression**: Consider compression level trade-offs between size and speed

### Post-Sorting Validation

1. **Sort Verification**: Validate proper coordinate ordering
2. **Index Generation**: Create BAI/CSI indices for random access
3. **Statistics Collection**: Generate alignment statistics for quality assessment
4. **File Integrity**: Verify BGZF compression and checksums

## Integration with Analysis Pipelines

### Workflow Integration

Sorting typically occurs between alignment and downstream analysis:  
1. **Read Alignment**: Generate initial SAM/BAM files  
2. **Quality Filtering**: Remove low-quality alignments  
3. **Sorting**: Coordinate-based organization  
4. **Indexing**: Generate access indices  
5. **Analysis**: Variant calling, coverage analysis, etc.  

### Automation Considerations

- **Dependency Management**: Ensure proper software module loading
- **Error Handling**: Implement robust error detection and recovery
- **Resource Monitoring**: Track computational resource utilization
- **Scalability**: Design for varying dataset sizes and computational environments

## Troubleshooting Common Issues

### Memory-Related Problems

- **OutOfMemoryError**: Increase JVM heap size or reduce records in RAM
- **System Memory**: Monitor system memory usage during sorting
- **Swap Usage**: Avoid excessive swap file utilization

### Performance Bottlenecks

- **I/O Limitations**: Optimize temporary file location and storage speed
- **CPU Utilization**: Balance thread count with available cores
- **Network Latency**: Minimize network-based file operations

### Data Integrity Issues

- **Truncated Files**: Verify complete file transfer before sorting
- **Corruption Detection**: Implement checksum validation
- **Format Compliance**: Ensure proper SAM/BAM format adherence

## Conclusion

BAM file sorting represents a critical preprocessing step that significantly impacts the efficiency and accuracy of downstream genomic analyses. Proper implementation requires careful consideration of computational resources, algorithm selection, and quality control measures. The integration of sorting operations into automated pipelines enables scalable processing of large-scale genomic datasets while maintaining data integrity and analytical reproducibility.