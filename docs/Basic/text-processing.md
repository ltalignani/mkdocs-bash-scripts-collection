# Text Processing in Bash

Text processing is one of the most powerful features of UNIX systems. Whether you're analyzing data files, processing log files, or manipulating configuration files, mastering text processing commands will dramatically improve your productivity. This guide covers essential tools from basic text manipulation to advanced pattern matching and data extraction.

## Table of Content

[TOC]

## Basic Text Processing Commands

### `cat` - Concatenate and Display Files

The `cat` command is the foundation of text processing, allowing you to display and combine text files.

**Basic usage:**
```bash
cat file.txt                        # Display file contents
cat file1.txt file2.txt             # Display multiple files
cat file1.txt file2.txt > combined.txt    # Combine files
```

**Useful options:**
```bash
cat -n file.txt                     # Show line numbers
cat -A file.txt                     # Show all characters (including hidden ones)
cat -s file.txt                     # Squeeze multiple blank lines into one
```

**Creating files with cat:**
```bash
cat > new_file.txt                  # Create file and type content (Ctrl+D to finish)
cat >> existing_file.txt            # Append to existing file
```

### `head` and `tail` - Display File Portions

These commands are essential for examining large files.

**head - First lines:**
```bash
head file.txt                       # First 10 lines (default)
head -n 20 file.txt                 # First 20 lines
head -n 5 *.txt                     # First 5 lines of all .txt files
```

**tail - Last lines:**
```bash
tail file.txt                       # Last 10 lines (default)
tail -n 15 file.txt                 # Last 15 lines
tail -f logfile.txt                 # Follow file changes (useful for logs)
tail -F logfile.txt                 # Follow file, handle rotation
```

**Combining head and tail:**
```bash
head -n 50 file.txt | tail -n 10    # Lines 41-50 of file
```

### `wc` - Word, Line, and Character Count

The `wc` command provides quick statistics about text files.

```bash
wc file.txt                         # Lines, words, characters
wc -l file.txt                      # Count lines only
wc -w file.txt                      # Count words only
wc -c file.txt                      # Count characters only
wc -m file.txt                      # Count characters (multibyte aware)
```

**Practical examples:**
```bash
wc -l *.txt                         # Count lines in all text files
ls | wc -l                          # Count files in directory
cat file.txt | wc -w                # Count words using pipes
```

## Pattern Matching and Searching

### `grep` - Search Text Patterns

The `grep` command searches for patterns in text files and is one of the most important text processing tools.

**Basic syntax:**
```bash
grep "pattern" file.txt
```

**Common options:**
```bash
grep -i "pattern" file.txt          # Case-insensitive search
grep -n "pattern" file.txt          # Show line numbers
grep -v "pattern" file.txt          # Show lines that DON'T match
grep -r "pattern" directory/        # Recursive search in directory
grep -l "pattern" *.txt             # Show only filenames that match
grep -c "pattern" file.txt          # Count matching lines
```

**Advanced pattern matching:**
```bash
grep "^start" file.txt              # Lines starting with "start"
grep "end$" file.txt                # Lines ending with "end"
grep "^$" file.txt                  # Empty lines
grep -E "pattern1|pattern2" file.txt # Multiple patterns (extended regex)
grep -w "word" file.txt             # Match whole words only
```

**Practical examples:**
```bash
# Find all Python files containing a specific function
grep -r "def process_data" --include="*.py" .

# Find error messages in log files
grep -i "error" /var/log/*.log

# Find email addresses in text files
grep -E "[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}" file.txt

# Search for IP addresses
grep -E "([0-9]{1,3}\.){3}[0-9]{1,3}" file.txt
```

### `grep` with Regular Expressions

Regular expressions (regex) make `grep` extremely powerful:

**Basic regex patterns:**
```bash
grep "colou\?r" file.txt            # Match "color" or "colour"
grep "test[123]" file.txt           # Match "test1", "test2", or "test3"
grep "file[0-9]" file.txt           # Match "file" followed by any digit
grep "^[A-Z]" file.txt              # Lines starting with uppercase letter
grep "[0-9]\{3\}" file.txt          # Three consecutive digits
```

**Extended regex with -E:**
```bash
grep -E "test[0-9]+" file.txt       # "test" followed by one or more digits
grep -E "(jpg|png|gif)$" file.txt   # Lines ending with image extensions
grep -E "^[a-zA-Z]+:" file.txt      # Lines starting with word followed by colon
```

## Text Sorting and Uniqueness

### `sort` - Sort Lines of Text

The `sort` command arranges lines in various orders.

**Basic sorting:**
```bash
sort file.txt                       # Sort alphabetically
sort -r file.txt                    # Reverse sort
sort -n file.txt                    # Numeric sort
sort -u file.txt                    # Sort and remove duplicates
```

**Advanced sorting:**
```bash
sort -k2 file.txt                   # Sort by second field
sort -k2,2 file.txt                 # Sort by second field only
sort -t: -k3 /etc/passwd            # Sort by third field, using : as delimiter
sort -k1,1 -k2,2n file.txt          # Sort by first field, then by second numerically
```

**Practical examples:**
```bash
# Sort files by size
ls -la | sort -k5 -n                # Sort by file size (5th column)

# Sort IP addresses correctly
sort -t. -k1,1n -k2,2n -k3,3n -k4,4n ip_addresses.txt

# Sort CSV file by specific column
sort -t, -k3 data.csv               # Sort CSV by third column
```

### `uniq` - Remove or Report Duplicate Lines

The `uniq` command works with sorted input to handle duplicates.

**Basic usage:**
```bash
sort file.txt | uniq                # Remove duplicate lines
sort file.txt | uniq -c             # Count occurrences of each line
sort file.txt | uniq -d             # Show only duplicate lines
sort file.txt | uniq -u             # Show only unique lines
```

**Practical examples:**
```bash
# Find most common words in a file
tr ' ' '\n' < file.txt | sort | uniq -c | sort -nr

# Find duplicate IP addresses in log files
grep -E "([0-9]{1,3}\.){3}[0-9]{1,3}" access.log | sort | uniq -d

# Count unique visitors
cut -d' ' -f1 access.log | sort | uniq | wc -l
```

## Field and Column Processing

### `cut` - Extract Columns from Text

The `cut` command extracts specific columns or fields from text.

**Extract by character position:**
```bash
cut -c1-10 file.txt                 # Characters 1-10
cut -c1,5,10 file.txt               # Characters 1, 5, and 10
cut -c10- file.txt                  # From character 10 to end
```

**Extract by field (delimiter-separated):**
```bash
cut -d: -f1 /etc/passwd             # First field, colon-separated
cut -d, -f2,4 data.csv              # Fields 2 and 4, comma-separated
cut -d$'\t' -f1-3 file.txt          # First 3 fields, tab-separated
```

**Practical examples:**
```bash
# Extract usernames from /etc/passwd
cut -d: -f1 /etc/passwd

# Extract IP addresses from log files (assuming they're in first column)
cut -d' ' -f1 access.log | sort | uniq

# Extract specific columns from CSV
cut -d, -f1,3,5 data.csv > selected_columns.csv
```

### `tr` - Translate Characters

The `tr` command transforms characters in text streams.

**Basic character translation:**
```bash
tr 'a' 'b' < file.txt               # Replace 'a' with 'b'
tr 'a-z' 'A-Z' < file.txt           # Convert lowercase to uppercase
tr 'A-Z' 'a-z' < file.txt           # Convert uppercase to lowercase
```

**Delete characters:**
```bash
tr -d '0-9' < file.txt              # Delete all digits
tr -d ' ' < file.txt                # Delete all spaces
tr -d '\n' < file.txt               # Delete newlines (join lines)
```

**Squeeze characters:**
```bash
tr -s ' ' < file.txt                # Squeeze multiple spaces into one
tr -s '\n' < file.txt               # Remove empty lines
```

**Practical examples:**
```bash
# Clean up text files
tr -d '\r' < dos_file.txt > unix_file.txt    # Remove carriage returns

# Convert spaces to underscores
tr ' ' '_' < file.txt

# Extract words (replace non-letters with newlines)
tr -cs 'A-Za-z' '\n' < file.txt
```

## Introduction to `sed` - Stream Editor

`sed` is a powerful stream editor that can perform complex text transformations. We'll start with simple examples and build up to more advanced usage.

### Basic `sed` Operations

**Substitution (find and replace):**
```bash
sed 's/old/new/' file.txt           # Replace first occurrence per line
sed 's/old/new/g' file.txt          # Replace all occurrences
sed 's/old/new/2' file.txt          # Replace second occurrence per line
```

**Case-insensitive substitution:**
```bash
sed 's/old/new/gi' file.txt         # Global, case-insensitive replacement
```

**In-place editing:**
```bash
sed -i 's/old/new/g' file.txt       # Modify file directly
sed -i.bak 's/old/new/g' file.txt   # Create backup before modifying
```

### `sed` Line Operations

**Delete lines:**
```bash
sed '3d' file.txt                   # Delete line 3
sed '2,5d' file.txt                 # Delete lines 2-5
sed '/pattern/d' file.txt           # Delete lines containing pattern
sed '/^$/d' file.txt                # Delete empty lines
```

**Print specific lines:**
```bash
sed -n '5p' file.txt                # Print only line 5
sed -n '10,20p' file.txt            # Print lines 10-20
sed -n '/pattern/p' file.txt        # Print lines containing pattern
```

**Add/insert lines:**
```bash
sed '3i\New line' file.txt          # Insert line before line 3
sed '3a\New line' file.txt          # Add line after line 3
sed '/pattern/a\New line' file.txt  # Add line after pattern match
```

### Practical `sed` Examples

**Configuration file editing:**
```bash
# Comment out lines containing a pattern
sed 's/^/#/' config.txt             # Comment out all lines
sed '/database/s/^/#/' config.txt   # Comment out lines with "database"

# Uncomment lines
sed 's/^#//' config.txt             # Remove # from beginning of lines
sed '/server/s/^#//' config.txt     # Uncomment lines containing "server"
```

**Log file processing:**
```bash
# Extract specific parts of log entries
sed 's/.*\[\([^]]*\)\].*/\1/' access.log    # Extract timestamp from [timestamp]

# Clean up log files
sed 's/[0-9]\{1,3\}\.[0-9]\{1,3\}\.[0-9]\{1,3\}\.[0-9]\{1,3\}/IP_ADDRESS/g' access.log
```

**Data cleaning:**
```bash
# Remove trailing whitespace
sed 's/[[:space:]]*$//' file.txt

# Convert multiple spaces to single space
sed 's/  */ /g' file.txt

# Remove Windows line endings
sed 's/\r$//' file.txt
```

## Introduction to `awk` - Pattern Scanning and Processing

`awk` is a powerful programming language designed for text processing. It excels at processing structured text files, especially those with columns.

### Basic `awk` Concepts

**Structure of awk:**
```bash
awk 'pattern { action }' file.txt
```

**Built-in variables:**  
- `NR`: Number of records (lines) processed  
- `NF`: Number of fields in current record  
- `$0`: Entire current record  
- `$1, $2, ...`: Individual fields  
- `FS`: Field separator (default: space/tab)  

### Simple `awk` Examples

**Print specific fields:**
```bash
awk '{print $1}' file.txt           # Print first field
awk '{print $1, $3}' file.txt       # Print first and third fields
awk '{print $NF}' file.txt          # Print last field
awk '{print NR, $0}' file.txt       # Print line number and entire line
```

**Field separators:**
```bash
awk -F: '{print $1}' /etc/passwd    # Use colon as field separator
awk -F, '{print $2}' data.csv       # Use comma as field separator
awk 'BEGIN{FS=","} {print $1}' data.csv    # Set field separator in BEGIN block
```

**Pattern matching:**
```bash
awk '/pattern/ {print}' file.txt    # Print lines containing pattern
awk '$1 == "value" {print}' file.txt # Print lines where first field equals "value"
awk '$3 > 100 {print}' file.txt     # Print lines where third field > 100
```

### `awk` Calculations and Processing

**Mathematical operations:**
```bash
awk '{sum += $1} END {print sum}' file.txt       # Sum first column
awk '{print $1 * $2}' file.txt                  # Multiply first two columns
awk '{avg += $1; count++} END {print avg/count}' file.txt  # Average of first column
```

**String operations:**
```bash
awk '{print length($0)}' file.txt               # Print length of each line
awk '{print toupper($1)}' file.txt              # Convert first field to uppercase
awk '{print substr($1, 1, 3)}' file.txt        # Print first 3 characters of first field
```

**Conditional processing:**
```bash
awk '{if ($1 > 50) print "High:", $0; else print "Low:", $0}' file.txt
awk '$3 > 100 {print $1, "is high"}' file.txt
awk 'NR > 1 {print}' file.txt                   # Skip header line
```

### Practical `awk` Examples

**CSV processing:**
```bash
# Calculate total of a column in CSV
awk -F, '{sum += $3} END {print "Total:", sum}' data.csv

# Print rows where a condition is met
awk -F, '$4 > 1000 {print $1, $2}' sales.csv

# Count records by category
awk -F, '{count[$2]++} END {for (i in count) print i, count[i]}' data.csv
```

**Log file analysis:**
```bash
# Count HTTP status codes
awk '{print $9}' access.log | sort | uniq -c

# Calculate average response time (assuming time is in field 10)
awk '{sum += $10; count++} END {print "Average:", sum/count}' access.log

# Extract and count unique IP addresses
awk '{print $1}' access.log | sort | uniq -c | sort -nr
```

**System administration:**
```bash
# Find processes using most CPU (from ps output)
ps aux | awk 'NR > 1 {print $3, $11}' | sort -nr | head -10

# Parse /etc/passwd for user information
awk -F: '{print $1 " has shell " $7}' /etc/passwd

# Monitor disk usage
df -h | awk 'NR > 1 {print $5, $6}' | sort -nr
```

## Advanced Text Processing Techniques

### Combining Commands with Pipes

The real power of UNIX text processing comes from combining commands:

**Complex pipelines:**
```bash
# Find most common words in a file
cat file.txt | tr -s ' ' '\n' | tr '[:upper:]' '[:lower:]' | sort | uniq -c | sort -nr | head -10

# Analyze web server logs
cat access.log | awk '{print $1}' | sort | uniq -c | sort -nr | head -10

# Process CSV data
cat data.csv | sed '1d' | cut -d, -f3 | sort -n | tail -10
```

### Using `xargs` for Batch Processing

`xargs` converts input into arguments for other commands:

```bash
# Find and process multiple files
find . -name "*.txt" | xargs grep "pattern"

# Delete files found by find
find . -name "*.tmp" -print0 | xargs -0 rm

# Process files with custom commands
ls *.txt | xargs -I {} cp {} backup/{}
```

### Text Processing with `paste` and `join`

**paste - Merge lines:**
```bash
paste file1.txt file2.txt           # Merge lines side by side
paste -d, file1.txt file2.txt       # Use comma as delimiter
paste -s file.txt                   # Merge all lines into one
```

**join - Join files on common field:**
```bash
join -t, file1.csv file2.csv        # Join CSV files on first field
join -1 2 -2 1 file1.txt file2.txt  # Join using field 2 of file1 and field 1 of file2
```

## Practical Examples and Use Cases

### Example 1: Processing CSV Data

**Sample CSV file (sales.csv):**
```csv
Date,Product,Quantity,Price
2024-01-01,Widget A,10,25.50
2024-01-02,Widget B,5,15.75
2024-01-03,Widget A,8,25.50
```

**Processing tasks:**
```bash
# Remove header and calculate total revenue
tail -n +2 sales.csv | awk -F, '{revenue += $3 * $4} END {print "Total Revenue:", revenue}'

# Find best-selling product
tail -n +2 sales.csv | awk -F, '{sales[$2] += $3} END {for (p in sales) print p, sales[p]}' | sort -k2 -nr | head -1

# Convert to different format
awk -F, 'NR > 1 {printf "%s: %d units of %s at $%.2f each\n", $1, $3, $2, $4}' sales.csv
```

### Example 2: Log File Analysis

**Sample log entry:**
```
192.168.1.100 - - [01/Jan/2024:12:00:00 +0000] "GET /index.html HTTP/1.1" 200 1024
```

**Analysis tasks:**
```bash
# Extract IP addresses and count occurrences
awk '{print $1}' access.log | sort | uniq -c | sort -nr | head -10

# Find 404 errors
awk '$9 == 404 {print $1, $7}' access.log

# Calculate total bytes transferred
awk '{sum += $10} END {print "Total bytes:", sum}' access.log

# Find peak traffic hours
awk '{gsub(/\[/, "", $4); gsub(/:.*/, "", $4); print $4}' access.log | sort | uniq -c | sort -nr
```

### Example 3: Configuration File Management

**Process configuration files:**
```bash
# Extract uncommented lines
sed '/^#/d; /^$/d' config.conf

# Add comments to specific lines
sed '/database/s/^/# /' config.conf

# Update configuration values
sed 's/^port=.*/port=8080/' config.conf

# Extract configuration values
awk -F= '/^[^#]/ {print $1, $2}' config.conf
```

### Example 4: Data Cleaning and Transformation

**Clean messy data:**
```bash
# Remove extra whitespace
sed 's/^[[:space:]]*//; s/[[:space:]]*$//' file.txt

# Standardize date format
sed 's/\([0-9]\{1,2\}\)\/\([0-9]\{1,2\}\)\/\([0-9]\{4\}\)/\3-\2-\1/g' file.txt

# Convert phone numbers to standard format
sed 's/(\([0-9]\{3\}\)) \([0-9]\{3\}\)-\([0-9]\{4\}\)/\1-\2-\3/g' file.txt

# Extract email addresses
grep -E -o "\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b" file.txt
```

## Performance Tips and Best Practices

### Efficient Text Processing

**Choose the right tool:**
- `grep`: Fast searching and filtering
- `sed`: Simple substitutions and line operations
- `awk`: Field processing and calculations
- `cut`: Simple column extraction
- `tr`: Character transformations

**Performance considerations:**
```bash
# Use grep for simple searches (faster than sed/awk)
grep "pattern" file.txt

# Use cut for simple column extraction
cut -d, -f1 file.csv

# Use awk for complex processing
awk 'complex_logic' file.txt

# Process large files efficiently
# Instead of: cat huge_file.txt | grep pattern
# Use: grep pattern huge_file.txt
```

### Memory and Processing Tips

**Handle large files:**
```bash
# Process files in chunks
split -l 10000 large_file.txt chunk_

# Use streaming where possible
tail -f logfile.txt | grep "ERROR"

# Avoid loading entire file into memory
# Good: grep pattern file.txt
# Avoid: cat file.txt | grep pattern (for large files)
```

### Common Pitfalls and Solutions

**Quoting and escaping:**
```bash
# Problem: Special characters in patterns
grep '$100' file.txt              # Wrong: $ has special meaning
grep '\$100' file.txt             # Correct: Escape the $

# Problem: Spaces in filenames
grep pattern file name.txt        # Wrong: Treated as two files
grep pattern "file name.txt"      # Correct: Quote the filename
```

**Field separators:**
```bash
# Problem: Inconsistent whitespace
awk '{print $2}' file.txt         # May not work with multiple spaces
awk '{print $2}' file.txt         # awk handles multiple spaces correctly

# Problem: Wrong delimiter
cut -d' ' -f2 file.csv            # Wrong: CSV uses commas
cut -d, -f2 file.csv              # Correct: Use comma delimiter
```

## Summary and Next Steps

Text processing is a fundamental skill in UNIX environments. Key takeaways:

**Essential commands to master:**  
- `grep` for searching and filtering  
- `sed` for substitutions and line operations  
- `awk` for field processing and calculations  
- `cut`, `sort`, `uniq` for data manipulation  
- `tr` for character transformations  

**Best practices:**  
- Start with simple commands and build complexity gradually  
- Use pipes to combine commands for powerful workflows  
- Choose the right tool for each task  
- Test commands on small files before processing large datasets  
- Always backup important files before in-place editing  

**Continue learning:**  
- Practice with real data files  
- Learn more advanced regex patterns  
- Explore shell scripting to automate text processing tasks  
- Study more advanced `awk` programming features  
- Learn about other text processing tools like `perl` and `python`  

The combination of these tools provides incredibly powerful text processing capabilities. With practice, you'll be able to handle complex data manipulation tasks efficiently and elegantly.