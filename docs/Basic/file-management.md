# File Management in Bash

Effective file management is crucial for productivity in UNIX environments. This guide will teach you advanced techniques for organizing, searching, and manipulating files using bash commands. These skills will help you handle large datasets, maintain organized project structures, and automate repetitive tasks.

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Understanding File Systems and Paths

### Absolute vs Relative Paths

Understanding paths is fundamental to file management:

**Absolute paths** start from the root directory (`/`) and specify the complete location:
```bash
/home/username/documents/project/data.txt
```

**Relative paths** are relative to your current working directory:
```bash
documents/project/data.txt    # From home directory
../project/data.txt          # From documents directory
./data.txt                   # From project directory
```

**Special path symbols:**  
- `~`: Home directory (`/home/username`)  
- `.`: Current directory  
- `..`: Parent directory  
- `-`: Previous directory (used with `cd`)  

### File Extensions and Types

Unlike Windows, UNIX systems don't rely on file extensions to determine file types. However, extensions are still useful for organization:

```bash
file document.txt        # Determine actual file type
file image.jpg          # Check if it's really an image
file script.py          # Verify file content type
```

## Advanced Directory Operations

### Creating Complex Directory Structures

The `mkdir` command with the `-p` option can create entire directory trees:

```bash
# Create a project structure
mkdir -p project/{src,docs,tests,data/{raw,processed}}

# This creates:
# project/
# ├── src/
# ├── docs/
# ├── tests/
# └── data/
#     ├── raw/
#     └── processed/
```

**Verify the structure:**
```bash
tree project/           # If tree is available
# or
find project/ -type d   # List all directories
```

### Bulk Directory Operations

Create multiple directories at once:
```bash
mkdir week{1..10}                    # Creates week1, week2, ..., week10
mkdir -p courses/{math,physics,chemistry}/{homework,notes,exams}
```

### Directory Navigation Shortcuts

```bash
cd -                    # Switch to previous directory
pushd /path/to/dir     # Save current directory and go to new one
popd                   # Return to saved directory
dirs                   # Show directory stack
```

## Advanced File Operations

### Copying Files and Directories

The `cp` command has many powerful options:

```bash
# Basic copying
cp source.txt destination.txt

# Copy with preservation of attributes
cp -p file.txt backup.txt           # Preserve timestamps and permissions
cp -a directory/ backup_directory/  # Archive mode (preserves everything)

# Interactive and verbose copying
cp -iv source.txt destination.txt   # Ask before overwriting, show progress

# Copy multiple files to directory
cp file1.txt file2.txt file3.txt target_directory/

# Copy with pattern matching
cp *.txt backup_directory/          # Copy all .txt files
cp data_*.csv analysis/             # Copy files matching pattern
```

### Advanced Moving and Renaming

```bash
# Rename multiple files (simple cases)
mv old_name.txt new_name.txt

# Move multiple files
mv *.log logs_directory/

# Rename with pattern (using bash parameter expansion)
for file in *.txt; do
    mv "$file" "${file%.txt}.backup"
done

# Move with backup of existing files
mv -b source.txt destination.txt    # Creates destination.txt~
```

### Safe File Deletion

```bash
# Always use -i for interactive deletion
rm -i unwanted_file.txt

# Remove empty directories
rmdir empty_directory/

# Remove directory and contents (be very careful!)
rm -ri directory_to_delete/         # Interactive recursive deletion

# Remove files older than a certain date
find . -name "*.tmp" -mtime +30 -delete    # Delete .tmp files older than 30 days
```

## File Search and Discovery

### Using `find` for File Search

The `find` command is incredibly powerful for locating files:

**Basic syntax:**  
```bash
find [path] [expression]
```

**Search by name:**  
```bash
find . -name "*.txt"                # Find all .txt files  
find . -name "data*"                # Find files starting with "data"  
find . -iname "*.PDF"               # Case-insensitive search  
find /home -name "config.yml"       # Search in specific directory
```

**Search by type:**
```bash
find . -type f                      # Find files only
find . -type d                      # Find directories only
find . -type l                      # Find symbolic links
```

**Search by size:**
```bash
find . -size +100M                  # Files larger than 100MB
find . -size -1k                    # Files smaller than 1KB
find . -size 50c                    # Files exactly 50 bytes
```

**Search by modification time:**
```bash
find . -mtime -7                    # Modified in last 7 days
find . -mtime +30                   # Modified more than 30 days ago
find . -mmin -60                    # Modified in last 60 minutes
```

**Search by permissions:**
```bash
find . -perm 755                    # Files with exact permissions
find . -perm -644                   # Files with at least these permissions
find . -executable                  # Executable files
```

**Combining search criteria:**
```bash
find . -name "*.log" -size +10M     # Large log files
find . -type f -name "*.tmp" -mtime +7 -delete    # Delete old temp files
```

### Using `locate` for Fast Search

The `locate` command uses a database for fast searching:

```bash
locate filename.txt                 # Fast search (requires updatedb)
locate -i filename                  # Case-insensitive
locate "*.pdf"                      # Search with patterns
```

**Update the locate database:**
```bash
sudo updatedb                       # Update database (run as administrator)
```

### Using `which` and `whereis`

Find executable files and their locations:

```bash
which python                        # Find location of python executable
which -a python                     # Find all python executables in PATH
whereis python                      # Find binary, source, and manual pages
```

## File Content Operations

### Viewing and Analyzing Files

**Get file information:**
```bash
file document.txt                   # Determine file type
stat document.txt                   # Detailed file information
ls -lh document.txt                 # Human-readable file size
```

**Count lines, words, and characters:**
```bash
wc document.txt                     # Lines, words, characters
wc -l document.txt                  # Count lines only
wc -w document.txt                  # Count words only
wc -c document.txt                  # Count characters only
```

### Comparing Files

```bash
diff file1.txt file2.txt            # Show differences between files
diff -u file1.txt file2.txt         # Unified diff format
diff -r dir1/ dir2/                 # Compare directories recursively
```

### Splitting and Joining Files

**Split large files:**
```bash
split -l 1000 large_file.txt        # Split into files of 1000 lines each
split -b 100M large_file.bin        # Split into 100MB chunks
split -d -l 1000 data.txt data_     # Split with numeric suffixes
```

**Join files back together:**
```bash
cat data_00 data_01 data_02 > combined_data.txt
```

## File Permissions and Ownership

### Understanding File Permissions

File permissions are displayed in the format: `drwxrwxrwx`

- First character: file type (`d` = directory, `-` = regular file, `l` = link)
- Next 9 characters: permissions for owner, group, and others
- Each group of 3: read (`r`), write (`w`), execute (`x`)

**Example permission breakdown:**
```bash
ls -l document.txt
# -rw-r--r-- 1 username group 1024 Jan 15 10:30 document.txt
# ││││││││││
# │││││││└└└─ others: read
# │││││││
# ││││└└└──── group: read, execute
# ││││
# │└└└─────── owner: read, write
# │
# └────────── regular file
```

### Changing Permissions

**Numeric method (octal):**
```bash
chmod 755 script.sh                 # rwxr-xr-x (owner: all, group/others: read+execute)
chmod 644 document.txt              # rw-r--r-- (owner: read+write, others: read)
chmod 600 private_file.txt          # rw------- (owner only)
```

**Symbolic method:**
```bash
chmod u+x script.sh                 # Add execute for owner
chmod g-w document.txt              # Remove write for group
chmod o+r public_file.txt           # Add read for others
chmod a+r readme.txt                # Add read for all (owner, group, others)
chmod u=rw,g=r,o=r document.txt     # Set specific permissions
```

**Recursive permission changes:**
```bash
chmod -R 755 directory/             # Apply to directory and all contents
chmod -R u+x scripts/               # Make all files in scripts/ executable for owner
```

### Changing Ownership

```bash
chown username file.txt             # Change owner
chown username:groupname file.txt   # Change owner and group
chown -R username:group directory/  # Change ownership recursively
```

## File Compression and Archives

### Creating Archives with `tar`

The `tar` command creates and extracts archives:

**Create archives:**
```bash
tar -cf archive.tar files/          # Create archive
tar -czf archive.tar.gz files/      # Create compressed archive (gzip)
tar -cjf archive.tar.bz2 files/     # Create archive with bzip2 compression
```

**Extract archives:**
```bash
tar -xf archive.tar                 # Extract archive
tar -xzf archive.tar.gz             # Extract gzip archive
tar -xjf archive.tar.bz2            # Extract bzip2 archive
```

**List archive contents:**
```bash
tar -tf archive.tar                 # List contents without extracting
tar -tzf archive.tar.gz             # List contents of compressed archive
```

**Common tar options:**  
- `-c`: Create archive   
- `-x`: Extract archive  
- `-t`: List contents  
- `-f`: Specify filename  
- `-z`: Use gzip compression  
- `-j`: Use bzip2 compression  
- `-v`: Verbose output  
- `-C`: Change to directory before operation  

### Working with Compressed Files

```bash
# Gzip compression
gzip file.txt                       # Compress file (creates file.txt.gz)
gunzip file.txt.gz                  # Decompress file
zcat file.txt.gz                    # View compressed file without extracting

# Zip archives
zip archive.zip file1.txt file2.txt # Create zip archive
zip -r archive.zip directory/       # Create zip archive recursively
unzip archive.zip                   # Extract zip archive
unzip -l archive.zip                # List zip contents
```

## Symbolic Links and Hard Links

### Creating and Managing Links

**Symbolic links (soft links):**
```bash
ln -s /path/to/original linkname    # Create symbolic link
ln -s ../data/file.txt shortcut.txt # Relative symbolic link
```

**Hard links:**
```bash
ln original.txt hardlink.txt        # Create hard link
```

**Managing links:**
```bash
ls -l linkname                      # Check where symbolic link points
readlink linkname                   # Show target of symbolic link
find . -type l                      # Find all symbolic links
```

### Differences Between Link Types

**Symbolic links:**  
- Point to a path (can be relative or absolute)  
- Can link to files or directories  
- Can link across filesystems  
- Become broken if target is moved or deleted  

**Hard links:**  
- Point to the same inode as the original file  
- Cannot link to directories  
- Cannot link across filesystems  
- Remain valid even if original filename is deleted  

## File Monitoring and Watching

### Real-time File Monitoring

```bash
tail -f logfile.txt                 # Follow file changes in real-time
tail -F logfile.txt                 # Follow file, handle rotation
watch -n 1 "ls -la"                 # Watch command output every second
```

### Monitoring Directory Changes

```bash
# If inotify-tools is available
inotifywait -m -e create,delete,modify directory/
```

## Batch File Operations

### Using Loops for Batch Processing

**Process multiple files:**
```bash
# Rename multiple files
for file in *.txt; do
    mv "$file" "${file%.txt}.backup"
done

# Convert file formats (example)
for file in *.md; do
    pandoc "$file" -o "${file%.md}.pdf"
done

# Process files in subdirectories
find . -name "*.log" -exec gzip {} \;
```

### Using Parameter Expansion

Bash parameter expansion provides powerful string manipulation:

```bash
filename="document.txt"
echo ${filename%.*}         # Remove extension: "document"
echo ${filename##*.}        # Get extension: "txt"
echo ${filename/.txt/.bak}  # Replace extension: "document.bak"

# Useful for batch renaming
for file in *.jpeg; do
    mv "$file" "${file%.jpeg}.jpg"
done
```

## File System Navigation and Organization

### Organizing Project Files

**Recommended directory structure:**
```bash
project/
├── README.md
├── src/                    # Source code
├── docs/                   # Documentation
├── tests/                  # Test files
├── data/                   # Data files
│   ├── raw/               # Original data
│   └── processed/         # Processed data
├── scripts/               # Utility scripts
├── results/               # Output files
└── archive/               # Old versions
```

**Create the structure:**
```bash
mkdir -p project/{src,docs,tests,data/{raw,processed},scripts,results,archive}
```

### File Naming Conventions

**Good practices:**  
- Use descriptive names  
- Avoid spaces (use underscores or hyphens)  
- Include dates in format YYYY-MM-DD  
- Use consistent extensions  

**Examples:**
```bash
# Good
data_2024-01-15_experiment_results.csv
analysis_script_v2.py
meeting_notes_2024-01-15.md

# Avoid
data file.csv
script.py
notes.txt
```

## Practical Examples and Use Cases

### Example 1: Organizing Downloaded Files

```bash
# Create organization structure
mkdir -p ~/Downloads/{documents,images,videos,archives,software}

# Move files by type
mv ~/Downloads/*.{pdf,doc,docx,txt} ~/Downloads/documents/
mv ~/Downloads/*.{jpg,jpeg,png,gif} ~/Downloads/images/
mv ~/Downloads/*.{mp4,avi,mkv} ~/Downloads/videos/
mv ~/Downloads/*.{zip,tar,gz} ~/Downloads/archives/
mv ~/Downloads/*.{deb,rpm,dmg,exe} ~/Downloads/software/
```

### Example 2: Backup Script

```bash
#!/bin/bash
# Simple backup script

backup_date=$(date +%Y-%m-%d)
backup_dir="/backup/daily_$backup_date"

mkdir -p "$backup_dir"
tar -czf "$backup_dir/documents_backup.tar.gz" ~/Documents/
tar -czf "$backup_dir/projects_backup.tar.gz" ~/Projects/

echo "Backup completed: $backup_dir"
```

### Example 3: Finding and Cleaning Old Files

```bash
# Find files older than 30 days
find ~/Downloads -type f -mtime +30

# Find large files (>100MB)
find ~ -size +100M -type f

# Clean temporary files
find /tmp -name "*.tmp" -mtime +7 -delete
find ~ -name "*~" -delete                    # Remove backup files
find ~ -name ".DS_Store" -delete             # Remove macOS metadata
```

### Example 4: Duplicate File Detection

```bash
# Find duplicate files by size
find . -type f -exec ls -la {} \; | sort -k5 -n | uniq -d -f4

# Find files with same name in different directories
find . -type f -printf "%f\n" | sort | uniq -d
```

## Troubleshooting Common Issues

### Permission Denied Errors

```bash
# Check file permissions
ls -la filename

# Fix common permission issues
chmod 644 document.txt              # For regular files
chmod 755 script.sh                 # For executable files
chmod 755 directory/                # For directories
```

### Disk Space Issues

```bash
# Check disk usage
df -h                               # Show disk space usage
du -sh *                            # Show directory sizes
du -sh * | sort -hr                 # Sort by size (largest first)

# Find large files
find . -size +100M -type f
```

### Broken Symbolic Links

```bash
# Find broken symbolic links
find . -type l -exec test ! -e {} \; -print

# Remove broken symbolic links
find . -type l -exec test ! -e {} \; -delete
```

## Best Practices for File Management

1. **Use descriptive names** that indicate the file's purpose and content  
2. **Organize files hierarchically** with logical directory structures  
3. **Regular backups** of important files  
4. **Clean up regularly** to remove unnecessary files  
5. **Use version control** for code and important documents  
6. **Set appropriate permissions** to protect sensitive files  
7. **Document your organization system** so others can understand it  
8. **Use archives** for long-term storage of completed projects  
9. **Test commands** on sample files before applying to important data  
10. **Keep track of file locations** using consistent naming and organization  

## Summary

Effective file management in bash involves mastering a variety of commands and techniques. Key skills include:

- Understanding file system navigation and paths  
- Using `find` for powerful file searching  
- Managing permissions and ownership appropriately  
- Creating and extracting archives  
- Implementing batch operations for efficiency  
- Organizing files systematically  
- Monitoring file changes when needed  
  
Practice these techniques regularly, and you'll develop the skills needed to manage files efficiently in any UNIX environment. Remember to always test commands on sample data before applying them to important files, and maintain regular backups of critical data.