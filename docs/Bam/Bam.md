# BAM Format Specification

## Overview

The Binary Alignment/Map (BAM) format is a compressed binary representation of the Sequence Alignment/Map (SAM) format, designed for efficient storage and retrieval of high-throughput sequencing alignment data. BAM files serve as the standard format for storing aligned sequence reads in computational genomics and bioinformatics pipelines.

## Format Architecture

### Binary Structure

BAM files employ a block-compressed format using the BGZF (Blocked GNU Zip Format) compression scheme, which combines the efficiency of gzip compression with random access capabilities. The file structure consists of:

1. **Header section**: Contains metadata including reference sequence information, read group details, and program records
2. **Alignment section**: Stores individual alignment records in binary format
3. **Index files** (optional): Separate `.bai` or `.csi` files enabling rapid random access

### Header Components

The BAM header contains essential metadata organized into record types:

- **@HD**: Header line specifying format version and sort order
- **@SQ**: Reference sequence dictionary with sequence names and lengths
- **@RG**: Read group information including sequencing platform and sample identifiers
- **@PG**: Program records documenting the software pipeline used for alignment

## Alignment Record Structure

Each alignment record in BAM format contains the following fields:

### Mandatory Fields

| Field | Type    | Description                                   |
| ----- | ------- | --------------------------------------------- |
| QNAME | String  | Query template name (read identifier)         |
| FLAG  | Integer | Bitwise flags indicating alignment properties |
| RNAME | String  | Reference sequence name                       |
| POS   | Integer | 1-based leftmost mapping position             |
| MAPQ  | Integer | Mapping quality score (Phred-scaled)          |
| CIGAR | String  | Compact Idiosyncratic Gapped Alignment Report |
| RNEXT | String  | Reference name of the mate/next read          |
| PNEXT | Integer | Position of the mate/next read                |
| TLEN  | Integer | Template length                               |
| SEQ   | String  | Segment sequence                              |
| QUAL  | String  | ASCII-encoded base quality scores             |

### Optional Fields

BAM records support extensible optional fields (tags) following the TAG:TYPE:VALUE format, enabling storage of additional alignment information such as:

- Alignment scores and statistics
- Original quality scores
- Color space information
- Custom annotations

## FLAG Field Interpretation

The FLAG field employs bitwise encoding to represent alignment properties:

| Bit   | Value | Description                              |
| ----- | ----- | ---------------------------------------- |
| 0x1   | 1     | Template having multiple segments        |
| 0x2   | 2     | Each segment properly aligned            |
| 0x4   | 4     | Segment unmapped                         |
| 0x8   | 8     | Next segment unmapped                    |
| 0x10  | 16    | SEQ reverse complemented                 |
| 0x20  | 32    | SEQ of next segment reverse complemented |
| 0x40  | 64    | First segment in template                |
| 0x80  | 128   | Last segment in template                 |
| 0x100 | 256   | Secondary alignment                      |
| 0x200 | 512   | Not passing quality controls             |
| 0x400 | 1024  | PCR or optical duplicate                 |
| 0x800 | 2048  | Supplementary alignment                  |

## CIGAR Operations

The CIGAR string describes the alignment between query and reference sequences using operation codes:

| Operation | Code | Description                                          |
| --------- | ---- | ---------------------------------------------------- |
| M         | 0    | Alignment match (can be sequence match or mismatch)  |
| I         | 1    | Insertion to the reference                           |
| D         | 2    | Deletion from the reference                          |
| N         | 3    | Skipped region from the reference                    |
| S         | 4    | Soft clipping (clipped sequences present in SEQ)     |
| H         | 5    | Hard clipping (clipped sequences NOT present in SEQ) |
| P         | 6    | Padding (silent deletion from padded reference)      |
| =         | 7    | Sequence match                                       |
| X         | 8    | Sequence mismatch                                    |

## Quality Scores

Base quality scores in BAM files follow the Phred scale encoding:  
- ASCII values represent quality scores with an offset of 33  
- Quality score Q relates to error probability P: Q = -10 × log₁₀(P)  
- Higher scores indicate greater confidence in base calls  

## Indexing and Random Access

BAM files support efficient random access through index files:

### BAI Index Format
- Traditional coordinate-based indexing
- Supports files up to 2³¹ bases per reference sequence
- Uses 16 kb linear index intervals

### CSI Index Format
- Enhanced indexing supporting longer reference sequences
- Configurable indexing parameters
- Better performance for highly variable coverage

## Compression and Performance

The BGZF compression scheme provides several advantages:

1. **Block-level compression**: Independent 64 KB blocks enable parallel processing
2. **Random access**: Seek operations without full decompression
3. **Compression efficiency**: Typically 3-5x size reduction compared to SAM
4. **Error detection**: Built-in CRC32 checksums for data integrity

## File Operations and Tools

Standard operations on BAM files include:

### Viewing and Conversion
- `samtools view`: Convert between SAM/BAM formats and filter records
- `samtools flagstat`: Generate alignment statistics

### Sorting and Indexing
- `samtools sort`: Sort alignments by coordinate or read name
- `samtools index`: Generate BAI/CSI index files

### Quality Control
- `samtools stats`: Comprehensive alignment statistics
- `samtools depth`: Calculate per-position coverage depth

## Best Practices

### File Handling
- Always sort BAM files by coordinate for optimal performance
- Generate index files for random access operations
- Validate file integrity using checksum verification

### Storage Considerations
- Use appropriate compression levels based on storage vs. access speed requirements
- Consider reference-based compression for highly redundant datasets
- Implement proper backup strategies for irreplaceable sequencing data

### Performance Optimization
- Utilize parallel processing capabilities of BGZF format
- Optimize memory usage for large-scale genomic analyses
- Implement streaming approaches for memory-constrained environments

## Specifications and Standards

The BAM format specification is maintained by the SAM/BAM Format Specification Working Group and follows semantic versioning. The current specification defines:

- Binary encoding schemes for all data types
- Mandatory and optional field requirements
- Compatibility requirements across different software implementations
- Extension mechanisms for format evolution

## Related Formats

### CRAM Format
- Reference-based compression achieving higher compression ratios
- Maintains compatibility with BAM-based tools
- Suitable for long-term archival storage

### SAM Format
- Human-readable text representation of alignment data
- Useful for debugging and manual inspection
- Less efficient for storage and processing

## Conclusion

The BAM format represents a cornerstone technology in modern genomics, providing efficient storage and manipulation of alignment data while maintaining compatibility across diverse bioinformatics software ecosystems. Its binary structure, combined with robust compression and indexing capabilities, enables scalable analysis of large-scale sequencing datasets essential for contemporary genomic research.