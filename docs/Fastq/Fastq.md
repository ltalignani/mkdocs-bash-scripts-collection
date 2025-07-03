# FASTQ Format

## Overview

The FASTQ format is a text-based file format widely used in bioinformatics to store nucleotide sequences along with their corresponding quality scores. Originally developed for the Sanger sequencing platform, FASTQ has become the de facto standard for representing high-throughput sequencing data from modern platforms including Illumina, Oxford Nanopore, and PacBio systems.

## Format Structure

Each sequence record in a FASTQ file consists of exactly four lines:

1. **Header line**: Begins with `@` followed by the sequence identifier and optional description
2. **Sequence line**: Contains the raw nucleotide sequence (A, T, G, C, N)
3. **Separator line**: Begins with `+` and optionally repeats the sequence identifier
4. **Quality line**: Contains quality scores encoded as ASCII characters

### Example Record

```
@SEQ_ID_001 description text
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

## Quality Score Encoding

Quality scores represent the confidence in each base call, typically expressed as Phred quality scores:

**Phred Quality Score Formula:**
```
Q = -10 × log₁₀(P)
```

Where `P` is the probability that the base call is incorrect.

### Common Encoding Schemes
FASTQ files combine nucleotide sequences and corresponding quality scores. Each quality score is `Phred-scaled`, meaning it’s a logarithmic measure of the probability that a base was called incorrectly. To store these scores efficiently in text format, they are encoded as `ASCII characters`, and different sequencing platforms (and versions) have used different encoding schemes.

| Encoding             | ASCII Range | Phred Score Range | Offset |
| -------------------- | ----------- | ----------------- | ------ |
| Sanger/Illumina 1.8+ | 33-126      | 0-93              | 33     |
| Illumina 1.3-1.7     | 64-126      | 0-62              | 64     |
| Solexa               | 59-126      | -5-62             | 64     |

### Quality Score Interpretation

| Phred Score | Error Probability | Accuracy | ASCII (Sanger) |
| ----------- | ----------------- | -------- | -------------- |
| 10          | 1 in 10           | 90%      | +              |
| 20          | 1 in 100          | 99%      | 5              |
| 30          | 1 in 1000         | 99.9%    | ?              |
| 40          | 1 in 10000        | 99.99%   | I              |

## File Extensions and Compression

### Standard Extensions
- `.fastq` - Uncompressed FASTQ file
- `.fq` - Alternative uncompressed extension
- `.fastq.gz` - Gzip-compressed FASTQ file
- `.fq.gz` - Alternative compressed extension

### Compression Benefits
Compression typically reduces file size by 70-80%, significantly improving storage efficiency and transfer speeds while maintaining data integrity.

## Paired-End Sequencing

For paired-end sequencing data, reads are typically stored in separate files:

- `sample_R1.fastq` - Forward reads (Read 1)
- `sample_R2.fastq` - Reverse reads (Read 2)

Corresponding reads maintain identical identifiers, often with `/1` and `/2` suffixes or `1:` and `2:` flags.

## Best Practices

### Data Integrity
- Always validate FASTQ format before analysis
- Verify quality score encoding scheme
- Check for consistent record structure (4 lines per read)

### Storage Considerations
- Use compression for long-term storage
- Maintain raw FASTQ files as primary data
- Implement proper backup strategies

### Quality Control
- Assess quality score distributions
- Monitor sequence length distributions
- Identify adapter contamination
- Evaluate per-base quality trends

## Common Applications

### Preprocessing
- Quality trimming and filtering
- Adapter removal
- Contamination screening
- Error correction

### Analysis Workflows
- Genome assembly
- Read mapping and alignment
- Variant calling
- RNA-seq expression analysis
- Metagenomics profiling

## Technical Considerations

### Memory Management  
Large FASTQ files require efficient processing strategies:  
- Stream-based processing for memory efficiency  
- Parallel processing for performance optimization  
- Chunked analysis for very large datasets  

### Format Validation  
Key validation checks include:  
- Consistent 4-line record structure  
- Valid nucleotide characters (A, T, G, C, N)  
- Matching sequence and quality string lengths  
- Proper ASCII encoding ranges  

## Related Formats

- **FASTA**: Sequence-only format without quality scores
- **SAM/BAM**: Aligned sequence format with additional mapping information
- **CRAM**: Compressed reference-based sequence format

## Conclusion

The FASTQ format remains fundamental to modern genomics workflows, providing essential quality information alongside sequence data. Understanding its structure, encoding schemes, and best practices is crucial for effective bioinformatics analysis and ensuring data integrity throughout the sequencing-to-analysis pipeline.