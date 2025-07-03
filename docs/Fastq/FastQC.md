# FastQC

## Overview

FastQC is a widely-used bioinformatics tool developed by the Babraham Institute for quality control analysis of high-throughput sequencing data. It provides comprehensive statistical assessments and visual representations of sequence quality metrics, enabling researchers to identify potential issues in their sequencing datasets before downstream analysis.

## Key Features

### Platform Compatibility
- **Cross-platform**: Available for Linux, macOS, and Windows
- **Multiple interfaces**: Command-line and graphical user interface (GUI)
- **Batch processing**: Simultaneous analysis of multiple files
- **Format support**: FASTQ, SAM, BAM, and compressed formats

### Analysis Modules
FastQC performs twelve distinct quality control analyses, each generating specific metrics and visualizations.

## Quality Control Modules

### 1. Basic Statistics
Provides fundamental dataset characteristics:   
- Total number of sequences  
- Sequence length distribution  
- GC content percentage  
- Sequence encoding format detection  

### 2. Per Base Sequence Quality
Evaluates Phred quality scores across all base positions:  
- **Pass**: Median quality ≥25 for all positions  
- **Warning**: Lower quartile quality <20 or median <25 for any position  
- **Fail**: Lower quartile quality <10 or median <20 for any position  

### 3. Per Tile Sequence Quality
Assesses quality variation across flow cell tiles (Illumina data):  
- Identifies problematic tiles or regions  
- Detects systematic quality degradation patterns  
- Useful for troubleshooting sequencing run issues  

### 4. Per Sequence Quality Scores
Analyzes the distribution of mean quality scores per read:  
- **Pass**: Peak at quality score >27  
- **Warning**: Peak between 20-27  
- **Fail**: Peak <20  

### 5. Per Base Sequence Content
Examines nucleotide composition across read positions:   
- **Pass**: Difference between A/T and G/C <10%  
- **Warning**: Difference between A/T and G/C 10-20%  
- **Fail**: Difference between A/T and G/C >20%  

### 6. Per Sequence GC Content
Compares observed GC content distribution to theoretical normal distribution:  
- **Pass**: Distribution matches expected normal curve  
- **Warning**: Sum of deviations >15% of reads  
- **Fail**: Sum of deviations >20% of reads  

### 7. Per Base N Content
Quantifies the proportion of uncalled bases (N) at each position:  
- **Pass**: N content <5% at all positions  
- **Warning**: N content 5-20% at any position  
- **Fail**: N content >20% at any position  

### 8. Sequence Length Distribution
Analyzes the distribution of sequence lengths:  
- **Pass**: All sequences same length  
- **Warning**: Variable lengths present  
- **Fail**: Sequences with zero length detected  

### 9. Sequence Duplication Levels
Identifies potential PCR amplification artifacts:  
- **Pass**: Non-unique sequences <20%  
- **Warning**: Non-unique sequences 20-50%  
- **Fail**: Non-unique sequences >50%  

### 10. Overrepresented Sequences
Detects sequences comprising >0.1% of the total:  
- Identifies potential contamination  
- Detects adapter sequences  
- Highlights PCR artifacts  

### 11. Adapter Content
Searches for common sequencing adapters:  
- Illumina Universal Adapter  
- Illumina Small RNA Adapter  
- Nextera Transposase Sequence  
- SOLID Small RNA Adapter  

### 12. Kmer Content
Identifies enriched kmers that may indicate bias:  
- **Pass**: No kmers with positional bias  
- **Warning**: Kmers with positional bias <1%  
- **Fail**: Kmers with positional bias >1%  

## Command Line Usage

### Basic Syntax
```bash
fastqc [options] seqfile1 seqfile2 .. seqfileN 
```

### Common Parameters
```bash
# Basic analysis
fastqc sample.fastq

# Specify output directory
fastqc -o /path/to/output sample.fastq

# Batch processing
fastqc *.fastq

# Specify number of threads
fastqc -t 4 sample.fastq

# Quiet mode (suppress progress output)
fastqc -q sample.fastq

# Extract results to directory
fastqc --extract sample.fastq
```

### Advanced Options
```bash
# Custom adapter sequences
fastqc -a adapter_file.txt sample.fastq

# Specify sequence format
fastqc -f fastq sample.fastq

# Set memory limit (Java heap size)
fastqc -Xmx2g sample.fastq

# Generate report without images
fastqc --nogroup sample.fastq
```

## Output Files

### HTML Report
- Interactive web-based report
- Graphical visualizations
- Pass/Warning/Fail status indicators
- Detailed explanations of each module

### Text Summary
- Tab-delimited summary statistics
- Machine-readable format
- Suitable for automated parsing
- Integration with pipeline workflows

### Data Files
- Raw data underlying each analysis
- Extracted when using `--extract` option
- Enables custom downstream analysis

## Interpretation Guidelines

### Quality Thresholds
| Metric               | Good         | Acceptable    | Poor           |
| -------------------- | ------------ | ------------- | -------------- |
| Per base quality     | >30          | 20-30         | <20            |
| Per sequence quality | Peak >27     | Peak 20-27    | Peak <20       |
| GC content           | Expected ±5% | Expected ±10% | Expected >±10% |
| Duplication level    | <20%         | 20-50%        | >50%           |

### Common Issues and Solutions

#### Low Quality Scores
- **Causes**: Sequencing chemistry degradation, over-clustering
- **Solutions**: Quality trimming, read filtering

#### High Duplication Levels
- **Causes**: PCR amplification bias, low library complexity
- **Solutions**: Duplicate removal, library optimization

#### Adapter Contamination
- **Causes**: Incomplete adapter removal, insert size issues
- **Solutions**: Adapter trimming, library size selection

#### Unusual GC Content
- **Causes**: Contamination, amplification bias, sample composition
- **Solutions**: Contamination screening, bias correction

## Integration with Workflows

### Preprocessing Pipelines
FastQC is typically integrated at multiple stages:  
1. **Raw data assessment**: Initial quality evaluation  
2. **Post-trimming validation**: Verify improvement after processing  
3. **Final validation**: Confirm data suitability for analysis  

### Automation Tools
- **MultiQC**: Aggregate FastQC reports across samples
- **Nextflow**: Workflow management integration
- **Snakemake**: Rule-based workflow incorporation
- **Galaxy**: Web-based platform integration

## Performance Considerations

### Memory Requirements
- Typical usage: 1-2 GB RAM per file
- Large files: Scale memory allocation accordingly
- Java heap size: Adjust with `-Xmx` parameter

### Processing Speed
- Typical throughput: 1-10 million reads per minute
- Factors affecting speed: File size, compression, storage I/O
- Parallelization: Use multiple threads for batch processing

## Best Practices

### Routine Quality Control
- Run FastQC on all sequencing datasets
- Compare results across samples and runs
- Establish baseline quality metrics for your platform

### Troubleshooting
- Correlate quality issues with sequencing run metadata
- Use FastQC results to optimize library preparation
- Document quality control decisions for reproducibility

### Report Interpretation
- Consider experiment-specific quality requirements
- Understand that some "failures" may be expected (e.g., RNA-seq bias)
- Use multiple QC metrics together for comprehensive assessment

## Limitations

### Analysis Scope
- Limited to sequence-level quality assessment
- Does not evaluate biological significance
- Cannot detect all types of experimental artifacts

### Threshold Sensitivity
- Generic thresholds may not suit all applications
- Some warnings may be acceptable for specific analyses
- Requires domain expertise for proper interpretation

## Conclusion

FastQC serves as an essential first step in sequencing data quality assessment, providing standardized metrics and visualizations that guide preprocessing decisions and identify potential experimental issues. Its comprehensive analysis modules, combined with intuitive reporting, make it an indispensable tool in modern genomics workflows. Proper interpretation of FastQC results, combined with understanding of experimental context, enables researchers to make informed decisions about data processing and analysis strategies.