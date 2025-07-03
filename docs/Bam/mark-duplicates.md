# MarkDuplicates

## Overview

MarkDuplicates is a critical preprocessing tool in genomic analysis pipelines that identifies and manages duplicate reads in aligned sequencing data. Duplicate reads arise from various sources during library preparation and sequencing, potentially introducing bias in downstream analyses. This tool implements sophisticated algorithms to detect, mark, or remove duplicate reads while preserving the highest quality representatives, thereby improving the accuracy of variant calling and other genomic analyses.

## Biological Context of Duplicate Reads

### Sources of Duplication

Duplicate reads in next-generation sequencing data originate from multiple sources:

#### PCR Duplicates
- **Library amplification**: PCR steps during library preparation create identical copies
- **Cluster amplification**: Bridge PCR on flow cells generates optical duplicates
- **Uneven amplification**: Preferential amplification of certain DNA fragments

#### Optical Duplicates
- **Cluster misidentification**: Sequencing instruments incorrectly separate adjacent clusters
- **Signal bleeding**: Fluorescent signals from adjacent clusters interfere
- **Image processing artifacts**: Base calling software misinterprets cluster boundaries

#### Sequencing Artifacts
- **Index hopping**: Incorrect assignment of reads to samples in multiplexed runs
- **Carry-over contamination**: Residual DNA from previous sequencing runs
- **Template switching**: Polymerase jumping between templates during amplification

### Impact on Downstream Analysis

Duplicate reads can significantly affect analytical outcomes:

#### Variant Calling Bias
- **False positive variants**: Duplicates artificially inflate variant allele frequencies
- **Coverage distortion**: Uneven duplicate distribution creates coverage artifacts
- **Quality score inflation**: Multiple observations of the same molecule inflate confidence

#### Population Genetics Analysis
- **Allele frequency skewing**: Duplicates bias population-level statistics
- **Linkage disequilibrium**: Artificial correlation between nearby variants
- **Demographic inference**: Incorrect population size and migration estimates

#### Copy Number Analysis
- **Depth-based calling**: Duplicates create false copy number gains
- **Breakpoint detection**: Duplicates obscure structural variant boundaries
- **Normalization issues**: Uneven duplication affects read depth normalization

## Duplicate Detection Algorithms

### Coordinate-Based Detection

The primary algorithm for duplicate identification relies on alignment coordinates:

#### Single-End Reads
- **Primary criteria**: Identical 5' mapping coordinates on the same strand
- **Secondary criteria**: Identical CIGAR strings for complex alignments
- **Quality assessment**: Retention of highest mapping quality reads

#### Paired-End Reads
- **Template coordinates**: Identical 5' coordinates for both read pairs
- **Orientation consistency**: Proper pair orientation requirements
- **Insert size validation**: Consistent template length measurements

### Sequence-Based Detection

Advanced algorithms incorporate sequence information:  
- **Exact sequence matching**: Identical read sequences indicate duplicates  
- **Fuzzy matching**: Sequence similarity with allowable mismatches  
- **Barcode integration**: Unique molecular identifiers (UMIs) for accurate detection  

### Optical Duplicate Detection

Specialized algorithms for optical duplicates:  
- **Pixel distance calculation**: Physical distance between clusters on flow cells  
- **Tile-based analysis**: Localized duplicate detection within sequencing tiles  
- **Platform-specific parameters**: Instrument-specific distance thresholds  

## Implementation Strategies

### Marking vs. Removal

Two primary approaches for managing duplicates:

#### Marking Duplicates
- **FLAG modification**: Sets bit 0x400 in BAM FLAG field
- **Preservation**: Maintains all reads while identifying duplicates
- **Flexibility**: Allows downstream tools to handle duplicates appropriately
- **Reversibility**: Enables recovery of original data if needed

#### Removing Duplicates
- **Physical removal**: Eliminates duplicate reads from output files
- **Storage efficiency**: Reduces file size and storage requirements
- **Processing speed**: Accelerates downstream analysis by reducing data volume
- **Irreversibility**: Permanent loss of duplicate read information

### Quality-Based Selection

When multiple duplicates exist, selection criteria include:  
- **Mapping quality**: Highest MAPQ score indicates best alignment  
- **Base quality**: Sum of base quality scores across read length  
- **Alignment score**: Aligner-specific scoring metrics  
- **Read completeness**: Preference for reads with fewer soft-clipped bases  

## Tool Parameters and Configuration

### Core Parameters

Essential configuration options for MarkDuplicates:

| Parameter          | Description                | Default     | Recommendation                    |
| ------------------ | -------------------------- | ----------- | --------------------------------- |
| REMOVE_DUPLICATES  | Remove duplicate reads     | false       | true for storage optimization     |
| CREATE_INDEX       | Generate BAI index         | false       | true for downstream compatibility |
| ASSUME_SORTED      | Input is coordinate sorted | false       | true for sorted BAM files         |
| MAX_RECORDS_IN_RAM | Memory buffer size         | 500000      | Adjust based on available RAM     |
| TMP_DIR            | Temporary file directory   | system temp | High-speed storage location       |

### Advanced Parameters

Specialized options for specific use cases:

| Parameter                        | Description                    | Usage                       |
| -------------------------------- | ------------------------------ | --------------------------- |
| OPTICAL_DUPLICATE_PIXEL_DISTANCE | Optical duplicate threshold    | Platform-specific tuning    |
| READ_NAME_REGEX                  | Custom read name parsing       | Non-standard naming schemes |
| DUPLICATE_SCORING_STRATEGY       | Quality scoring method         | SUM_OF_BASE_QUALITIES       |
| CLEAR_DT                         | Remove existing duplicate tags | Reprocessing scenarios      |
| ADD_PG_TAG_TO_READS              | Add program group tags         | Pipeline tracking           |

### Memory Management

Optimal memory configuration strategies:  
- **Heap allocation**: 1-2 GB per million reads as baseline  
- **Buffer sizing**: Balance memory usage with I/O efficiency  
- **Garbage collection**: Tune JVM parameters for sustained performance  
- **Temporary storage**: Ensure adequate disk space for intermediate files  

## Practical Implementation Example

The following example demonstrates production-ready duplicate marking using SLURM job scheduling:

```bash
#!/bin/bash
###################configuration slurm##############################
#SBATCH -A invalbo
#SBATCH --time=2-23:00:00
#SBATCH --job-name=mark_bams
#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --array=1-56
#SBATCH --output=Cluster_logs/%x-%j-%N.out
#SBATCH --error=Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
####################################################################

# Extract BAM filename from job array list
BAMNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" info_files/bam_list)

# Define output directory structure
OUTPUT_DIR=results/04_Polishing
mkdir -p ${OUTPUT_DIR}

# Define output file path
MARKED_BAM=${OUTPUT_DIR}/${BAMNAME}_marked.bam

# Create quality control directory
METRICS_DIR="qc/markdup"
mkdir -p ${METRICS_DIR}

# Create reporting directory
REPORT_DIR=results/11_Reports/markduplicates
mkdir -p ${REPORT_DIR}

# Load required software module
module load picard/2.23.5

# Execute duplicate marking with removal
picard MarkDuplicates -Xmx128000m \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT results/02_Mapping/${BAMNAME}_sorted.bam \
  --OUTPUT ${MARKED_BAM} \
  --TMP_DIR /tmp/ \
  --METRICS_FILE ${METRICS_DIR}/${BAMNAME}.metrics \
  &> ${REPORT_DIR}/${BAMNAME}_marked.out
```

### Script Analysis

**Resource Allocation:**  
- **Memory**: 128 GB allocation for processing large BAM files  
- **CPU**: 2 cores for I/O intensive operations  
- **Time**: Extended processing window for comprehensive datasets  
- **Array processing**: 56 parallel jobs for batch processing  

**File Management:**  
- **Input validation**: Systematic processing of BAM file lists  
- **Directory structure**: Organized output and quality control directories  
- **Temporary storage**: Efficient use of local temporary directories  

**Quality Control Integration:**  
- **Metrics collection**: Comprehensive duplicate statistics generation  
- **Reporting**: Detailed logging for process monitoring  
- **Indexing**: Automatic BAI index creation for downstream compatibility  

**Processing Configuration:**  
- **Duplicate removal**: Physical elimination of duplicate reads  
- **Index creation**: Automatic generation of BAM indices  
- **Memory optimization**: Efficient memory allocation for large datasets  

## Metrics and Quality Assessment

### Duplicate Statistics

MarkDuplicates generates comprehensive metrics including:

#### Library-Level Metrics
- **Total reads**: Complete read count in input file  
- **Duplicate reads**: Number of identified duplicates  
- **Duplication rate**: Percentage of reads identified as duplicates  
- **Optical duplicate rate**: Percentage of optically duplicated reads  

#### Molecular Metrics
- **Unique molecules**: Estimated number of distinct DNA molecules  
- **Molecules with duplicates**: Count of molecules generating duplicates  
- **Mean reads per molecule**: Average amplification rate  

#### Quality Metrics
- **Mapping quality distribution**: Quality scores of retained vs. duplicate reads  
- **Insert size distribution**: Template length statistics for paired-end reads  
- **Duplicate class distribution**: Breakdown by duplicate type  

### Interpretation Guidelines

#### Acceptable Duplication Rates
- **Whole genome sequencing**: 5-15% typical for well-prepared libraries  
- **Exome sequencing**: 10-25% due to target enrichment bias  
- **RNA sequencing**: 20-50% depending on library complexity  
- **Amplicon sequencing**: 30-70% expected for targeted approaches  

#### Quality Indicators
- **Low optical duplicates**: <5% indicates proper cluster density  
- **Consistent across lanes**: Similar rates suggest uniform processing  
- **Insert size distribution**: Proper library size selection  
- **Mapping quality retention**: High-quality reads preferentially retained  

## Alternative Tools and Approaches

### SAMtools Implementations

Basic duplicate marking capabilities:
```bash
# Mark duplicates using SAMtools
samtools markdup input_sorted.bam output_marked.bam

# Remove duplicates using SAMtools
samtools markdup -r input_sorted.bam output_deduped.bam
```

### Specialized Tools

#### UMI-Based Deduplication
- **UMI-tools**: Handles unique molecular identifiers  
- **Picard MarkDuplicates with UMI**: Enhanced UMI support  
- **Custom solutions**: Laboratory-specific UMI strategies  

#### High-Performance Implementations
- **Sambamba**: Parallelized duplicate marking  
- **Biobambam2**: Optimized for large-scale processing  
- **Custom implementations**: Specialized for specific sequencing platforms  

### Cloud-Based Solutions

Modern genomics platforms offer managed duplicate detection:  
- **Google Cloud Genomics**: Integrated duplicate marking  
- **AWS Genomics**: Scalable duplicate detection  
- **Azure Genomics**: Cloud-native processing pipelines  

## Integration with Analysis Pipelines

### Workflow Positioning

Duplicate marking typically occurs after sorting and before variant calling:  

1. **Raw sequencing data**: FASTQ files from sequencing platforms  
2. **Quality control**: Adapter trimming and quality filtering  
3. **Read alignment**: Mapping to reference genome  
4. **Coordinate sorting**: Organization by genomic position  
5. **Read group assignment**: Metadata correction and standardization  
6. **Duplicate marking**: Identification and management ← *Current process*  
7. **Base quality recalibration**: Systematic error correction  
8. **Variant calling**: Genomic variant identification  
9. **Variant filtering**: Quality-based variant selection  

### Downstream Compatibility

Proper duplicate handling ensures compatibility with:  
- **GATK Best Practices**: Standard genomics pipeline requirements  
- **Variant callers**: Accurate allele frequency estimation  
- **Population genetics tools**: Unbiased allele frequency calculations  
- **Structural variant detection**: Accurate breakpoint identification  

## Performance Optimization

### Computational Strategies

#### Memory Optimization
- **Heap tuning**: Optimal JVM memory allocation  
- **Buffer management**: Efficient I/O buffer sizing  
- **Garbage collection**: Appropriate GC strategy selection  

#### I/O Optimization
- **Temporary storage**: High-speed storage for intermediate files  
- **Parallel processing**: Multi-threaded duplicate detection  
- **Streaming algorithms**: Reduced memory footprint for large files  

#### Scalability Considerations
- **Distributed computing**: Cluster-based processing for large datasets  
- **Cloud deployment**: Elastic scaling based on workload  
- **Resource monitoring**: Real-time performance tracking  

### Best Practices

#### Pre-Processing Optimization
1. **Coordinate sorting**: Ensure proper BAM file sorting  
2. **Read group validation**: Verify proper read group assignment  
3. **File integrity**: Validate BAM file completeness  
4. **Resource planning**: Estimate computational requirements  

#### Processing Optimization
1. **Memory allocation**: Size heap based on input file characteristics  
2. **Temporary storage**: Use local, high-speed storage  
3. **Parallel execution**: Utilize job arrays for batch processing  
4. **Error handling**: Implement robust error detection and recovery  

#### Post-Processing Validation
1. **Metrics review**: Analyze duplicate statistics for quality assessment  
2. **File validation**: Verify output file integrity  
3. **Index generation**: Ensure proper BAM indexing  
4. **Downstream testing**: Validate compatibility with analysis tools  

## Troubleshooting Common Issues

### Memory-Related Problems

#### OutOfMemoryError
- **Symptoms**: JVM heap space exhaustion during processing  
- **Solutions**: Increase heap size, reduce records in RAM, optimize GC  
- **Prevention**: Estimate memory requirements based on file size  

#### System Memory Exhaustion
- **Symptoms**: System becomes unresponsive during processing  
- **Solutions**: Reduce concurrent jobs, optimize memory allocation  
- **Prevention**: Monitor system memory usage patterns  

### Performance Issues

#### Slow Processing
- **Causes**: Inadequate memory, slow storage, inefficient parameters  
- **Solutions**: Optimize memory allocation, use faster storage, tune parameters  
- **Monitoring**: Track processing rates and resource utilization  

#### I/O Bottlenecks
- **Causes**: Network storage, insufficient disk bandwidth  
- **Solutions**: Use local storage, optimize I/O patterns  
- **Prevention**: Benchmark storage performance before processing  

### Data Integrity Issues

#### Incomplete Processing
- **Symptoms**: Truncated output files, missing indices  
- **Solutions**: Verify input file integrity, check disk space  
- **Prevention**: Implement comprehensive error checking  

#### Metrics Inconsistencies
- **Symptoms**: Unexpected duplication rates, anomalous statistics  
- **Solutions**: Validate input data quality, review processing parameters  
- **Prevention**: Establish quality control thresholds  

## Quality Control and Validation

### Pre-Processing Checks

Essential validation steps before duplicate marking:  
1. **File integrity**: Verify BAM file completeness and accessibility  
2. **Sort order**: Confirm coordinate-based sorting  
3. **Read groups**: Validate proper read group assignment  
4. **Sample consistency**: Verify sample identification accuracy  

### Processing Monitoring

Real-time monitoring during duplicate marking:  
1. **Resource utilization**: Track memory, CPU, and I/O usage  
2. **Progress tracking**: Monitor processing completion rates  
3. **Error detection**: Identify processing failures early  
4. **Quality metrics**: Review duplicate statistics during processing  

### Post-Processing Validation

Comprehensive validation after duplicate marking:  
1. **Output integrity**: Verify complete file generation  
2. **Metrics analysis**: Review duplicate statistics for anomalies  
3. **Indexing verification**: Confirm proper BAM index generation  
4. **Downstream compatibility**: Test with analysis tools  

## Conclusion

MarkDuplicates represents a fundamental component of modern genomic analysis pipelines, providing essential functionality for managing duplicate reads that can significantly impact analytical accuracy. Proper implementation requires careful consideration of computational resources, algorithm selection, and quality control measures. The tool's integration into automated pipelines enables scalable processing of large-scale genomic datasets while maintaining data integrity and analytical reliability. Understanding the biological sources of duplication, implementing appropriate detection strategies, and maintaining rigorous quality control standards are essential for generating high-quality genomic datasets suitable for downstream analysis and clinical applications.