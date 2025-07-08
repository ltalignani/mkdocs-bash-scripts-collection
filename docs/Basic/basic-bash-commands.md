# Introduction to Basic Bash Commands

Welcome to the world of UNIX command-line interfaces! This guide will introduce you to essential bash commands that form the foundation of working in a UNIX environment. As a Master's student, mastering these commands will significantly enhance your productivity and open up powerful possibilities for file management, data processing, and system administration.

## Table of contents

[TOC]


## Getting Help: Your First Tools (--help & man)

Before diving into specific commands, let's learn how to get help when you're stuck or need more information about a command.

### The `--help` Option

Most commands in bash support the `--help` option, which provides a quick overview of the command's usage and available options.

**Syntax:**  
```bash
command --help
```

**Example:**  
```bash
ls --help
```

This will display a concise summary of how to use the `ls` command, including all its options and their meanings. The `--help` option is perfect for quick reference when you remember the command but need to check its options.

### The `man` Command (Manual Pages)

The `man` command provides access to the system's manual pages, which contain comprehensive documentation for commands, system calls, and configuration files.

**Syntax:**  
```bash
man command_name
```

**Example:**  
```bash
man ls
```

This opens the complete manual page for the `ls` command. Manual pages are more detailed than `--help` and include:

- **NAME**: Brief description of the command  
- **SYNOPSIS**: Usage syntax  
- **DESCRIPTION**: Detailed explanation of what the command does  
- **OPTIONS**: Complete list of all available options  
- **EXAMPLES**: Usage examples  
- **SEE ALSO**: Related commands  

**Navigation in man pages:**  
- Use arrow keys or `j`/`k` to scroll up and down  
- Press `Space` to scroll down one page  
- Press `q` to quit the manual page  
- Press `/` followed by a search term to search within the page  


!!! question "What does `tr` command?"
      This exercise is your first introduction to the man command. You will discover
      a command you are not familiar with. Let's take the example of the `tr` command.  
      1. Using the man command, explain what the `tr` command does.   
      2. Must at least one option (short or long) be passed as an argument to the `tr` command?  
      3. Can more than one option (short or long) be passed as an argument to the `tr` command?  
      4. Must at least one argument be passed to the `tr` command?  
      5. How many arguments (other than options) can be passed to the `tr` command?  
      6. How are arguments (other than options) passed to the tr command?  


??? abstract "Answers"
      1. From the NAME part of the command `tr`, we learn that it is used to ‘Convert or delete characters’.  
      2. The SYNOPSIS part of the man of the tr command contains the text below:  
         `SYNOPSIS`
                  `tr`  `[OPTION]`  `...`  `ENSEMBLE1`  `[ENSEMBLE2]`
         The notation with square brackets [ ] around OPTION indicates that the argument marked OPTION is optional. It is therefore not necessary to pass an option as an argument to the `tr` command.
      3. The three dots after `[OPTION]` indicate that you can specify more than one.
      4. The absence of a bracket around ENSEMBLE1 indicates that this argument is necessary, so must pass at least one argument to the tr command.
      5. `[ENSEMBLE2]` indicates that a second set can be specified, and as it is not followed by three dots, this means that a maximum of two sets can be passed as arguments to the tr command.
      6. The `DESCRIPTION` part of the `tr` command man contains the text below:  
         ```
         DESCRIPTION  
               Convert, compress and/or eliminate characters read from the standard input and write them to the standard output.

               -c, -C, --complement
                              use the complement of the ENSEMBLE1
         (...)
               --version
                              Display software name and version and exit
         ```


## File and Directory Operations (pwd, ls, cd, mkdir, rmdir)

### `pwd` - Print Working Directory

The `pwd` command shows you exactly where you are in the filesystem hierarchy.

**Syntax:**  
```bash
pwd
```

**Example:**  
```bash
pwd
# Output: /home/username/documents
```

This command is essential for understanding your current location before performing other operations.

### `ls` - List Directory Contents

The `ls` command displays the contents of directories. It's one of the most frequently used commands.

**Basic syntax:**  
```bash
ls [options] [directory]  
```

**Common options:**
- `-l`: Long format (detailed information including permissions, owner, size, date)  
- `-a`: Show all files, including hidden files (those starting with `.`)  
- `-h`: Human-readable file sizes (when used with `-l`)  
- `-t`: Sort by modification time  
- `-r`: Reverse the order of sorting  

**Examples:**  
```bash
ls                    # List current directory
ls -l                 # Detailed listing
ls -la                # Detailed listing including hidden files
ls -lh                # Detailed listing with human-readable sizes
ls /home/username     # List contents of specific directory
```

### `cd` - Change Directory

The `cd` command allows you to navigate between directories.

**Syntax:**  
```bash
cd [directory]
```

**Special directory references:**  
- `~`: Home directory  
- `.`: Current directory  
- `..`: Parent directory  
- `-`: Previous directory  

**Examples:**  
```bash
cd /home/username     # Navigate to specific directory
cd ~                  # Go to home directory
cd ..                 # Go up one level
cd -                  # Go back to previous directory
cd                    # Go to home directory (same as cd ~)
```

### `mkdir` - Make Directory

The `mkdir` command creates new directories.

**Syntax:**  
```bash
mkdir [options] directory_name
```

**Common options:**  
- `-p`: Create parent directories as needed  
- `-m`: Set permissions for the new directory  

**Examples:**  
```bash
mkdir new_folder                    # Create a single directory
mkdir -p documents/projects/2024    # Create nested directories
mkdir folder1 folder2 folder3       # Create multiple directories
```

### `rmdir` - Remove Directory

The `rmdir` command removes empty directories.

**Syntax:**  
```bash
rmdir directory_name
```

**Example:**  
```bash
rmdir empty_folder
```

**Note:** This command only works on empty directories. For directories with content, use `rm -r` (covered later).  

## File Operations (touch, cp, mv, rm)

### `touch` - Create Empty Files or Update Timestamps

The `touch` command creates new empty files or updates the timestamp of existing files.

**Syntax:**  
```bash
touch filename
```

**Examples:**  
```bash
touch new_file.txt                    # Create a new empty file
touch file1.txt file2.txt file3.txt   # Create multiple files
touch existing_file.txt               # Update timestamp of existing file
```

### `cp` - Copy Files and Directories

The `cp` command copies files and directories from one location to another.

**Syntax:**  
```bash
cp [options] source destination
```

**Common options:**  
- `-r` or `-R`: Copy directories recursively  
- `-i`: Interactive mode (ask before overwriting)  
- `-v`: Verbose mode (show what's being copied)  

**Examples:**  
```bash
cp file.txt backup.txt              # Copy file to new name
cp file.txt /home/username/         # Copy file to different directory
cp -r documents/ backup_documents/  # Copy directory and all contents
cp -i file.txt existing_file.txt    # Copy with confirmation prompt
```

### `mv` - Move/Rename Files and Directories

The `mv` command moves files and directories, and can also rename them.

**Syntax:**  
```bash
mv [options] source destination
```

**Common options:**  
- `-i`: Interactive mode (ask before overwriting)  
- `-v`: Verbose mode (show what's being moved)  

**Examples:**  
```bash
mv old_name.txt new_name.txt        # Rename a file
mv file.txt /home/username/         # Move file to different directory
mv documents/ /backup/              # Move directory
mv -i file.txt existing_file.txt    # Move with confirmation prompt
```

### `rm` - Remove Files and Directories

The `rm` command deletes files and directories. **Use with caution!**

**Syntax:**  
```bash
rm [options] file_or_directory
```

**Common options:**  
- `-r` or `-R`: Remove directories and their contents recursively  
- `-i`: Interactive mode (ask for confirmation)  
- `-f`: Force deletion without prompting  
- `-v`: Verbose mode (show what's being deleted)  

**Examples:**  
```bash
rm file.txt                   # Delete a file
rm -i important_file.txt      # Delete with confirmation
rm -r old_directory/          # Delete directory and all contents
rm -rf temp_folder/           # Force delete directory (be very careful!)
```

**Warning:** The `rm` command permanently deletes files. There's no "trash" or "recycle bin" in the command line. Always double-check before using `rm`, especially with the `-r` and `-f` options.

## Viewing File Contents (cat, less, head, tail)

### `cat` - Display File Contents

The `cat` command displays the entire contents of a file.

**Syntax:**  
```bash
cat filename
```

**Examples:**  
```bash
cat document.txt              # Display file contents
cat file1.txt file2.txt       # Display multiple files
cat > new_file.txt            # Create file and enter content (Ctrl+D to finish)
```

### `less / zless` - View File Contents Page by Page

The `less` command allows you to view file contents one page at a time, which is useful for large files. `zless` is the same command but for compressed files, finishing by `.gz`

**Syntax:**
```bash
less filename
zless -S filename.gz
```

**Navigation in less:**  
- Arrow keys or `j`/`k`: Move up and down  
- `Space`: Move forward one page  
- `b`: Move backward one page  
- `g`: Go to beginning of file  
- `G`: Go to end of file  
- `S`: Chop long lines rather than folding them
- `/pattern`: Search for pattern  
- `q`: Quit  

**Example:**  
```bash
zless -S sample1_R1.fastq.gz
```

### `head` - Display First Lines of a File

The `head` command shows the first few lines of a file (default is 10 lines).

**Syntax:**  
```bash
head [options] filename
```

**Common options:**  
- `-n number`: Show specific number of lines

**Examples:**  
```bash
head document.txt        # Show first 10 lines
head -n 5 document.txt   # Show first 5 lines
```

### `tail` - Display Last Lines of a File

The `tail` command shows the last few lines of a file (default is 10 lines).

**Syntax:**  
```bash
tail [options] filename
```

**Common options:**  
- `-n number`: Show specific number of lines  
- `-f`: Follow the file (useful for log files)  

**Examples:**  
```bash
tail document.txt        # Show last 10 lines
tail -n 20 document.txt  # Show last 20 lines
tail -f logfile.txt      # Follow log file updates (Ctrl+C to stop)
```

## File Permissions and Information (chmod, chown)

### `chmod` - Change File Permissions

The `chmod` command modifies file and directory permissions.

**Syntax:**  
```bash
chmod [options] permissions filename
```

**Permission notation:**  
In Unix, file permissions are represented by a 10-character string like `-rwxr-xr--`, where the first character indicates the file type (-: file; d: directory; l: symbolic link...), and the next nine specify read (r), write (w), and execute (x) permissions for the `owner`, `group`, and `others` (thus, each coded on three characters). 

```bash
ls -l document.txt
# -rwxr-xr-- 1 username group 1024 Jan 15 10:30 document.txt
# ││││││││││
# │││││││└└└─ others: read
# │││││││
# ││││└└└──── group: read, execute
# ││││
# │└└└─────── owner: read, write
# │
# └────────── regular file
```


These permissions can also be represented numerically (e.g., 751) using octal notation: 

- `r` (read) = 4  
- `w` (write) = 2  
- `x` (execute) = 1  

This system controls who can access or modify a file.  
In the example above, the type is a file the owner can read, write and execute (4+2+1=7). The group can read and execute (4+1=5) and others can only read (1). The octal notation in this case is 751.

**How to change permissions on a file:**  
```bash
chmod 755 script.sh      # rwxr-xr-x (owner: all, group/others: read+execute)
chmod 644 document.txt   # rw-r--r-- (owner: read+write, group/others: read only)
chmod +x script.sh       # Add execute permission
chmod -w document.txt    # Remove write permission
```

### `chown` - Change File Ownership

The `chown` command changes the owner and group of files and directories.

**Syntax:**  
```bash
chown [options] owner:group filename
```

**Example:**  
```bash
chown username:usergroup file.txt

# Change owner to 'alice' and group to 'researchers' for a single file
sudo chown alice:researchers my_project/file1.txt

# Change owner only (group remains unchanged)
sudo chown bob my_project/file2.txt

# Recursively change owner and group for the entire directory
sudo chown -R alice:researchers my_project/
```

**Note:** You typically need administrator privileges to change ownership.

## Process Management (ps, kill)

### `ps` - Display Running Processes

The `ps` command shows information about running processes.

**Syntax:**  
```bash
ps [options]
```

**Common options:**  
- `aux`: Show all processes with detailed information  
- `ef`: Show all processes in full format  

**Examples:**  
```bash
ps                # Show processes in current session
ps aux            # Show all processes with details
ps aux | grep python  # Show only Python processes
```

### `kill` - Terminate Processes

The `kill` command sends signals to processes, typically to terminate them.

**Syntax:**  
```bash
kill [options] process_id
```

**Common options:**  
- `-9`: Force kill (SIGKILL)  
- `-15`: Graceful termination (SIGTERM, default)  

**Examples:**  
```bash
kill 1234         # Terminate process with ID 1234
kill -9 1234      # Force terminate process 1234
```

## System Information (whoami, date, df, du)

### `whoami` - Display Current Username

The `whoami` command shows your current username.

**Syntax:**  
```bash
whoami
```

### `date` - Display or Set Date

The `date` command shows the current date and time.

**Syntax:**  
```bash
date [options]
```

**Example:**  
```bash
date                    # Show current date and time
date +"%Y-%m-%d %H:%M"  # Show formatted date
```

### `df` - Display Filesystem Disk Space

The `df` command shows disk space usage for filesystems.

**Syntax:**  
```bash
df [options]
```

**Common options:**  
- `-h`: Human-readable format  

**Example:**  
```bash
df -h              # Show disk usage in human-readable format
```

### `du` - Display Directory Space Usage

The `du` command shows disk space usage for directories.

**Syntax:**  
```bash
du [options] directory
```

**Common options:**  
- `-h`: Human-readable format  
- `-s`: Summary (total size only)  

**Examples:**  
```bash
du -h documents/        # Show space usage of documents directory
du -sh documents/       # Show total size of documents directory
```

## Best Practices and Tips 

1. **Always use `--help` or `man` when uncertain** about a command's options or usage.  

2. **Be careful with destructive commands** like `rm`, especially with options like `-r` and `-f`.  

3. **Use tab completion** to avoid typing errors and speed up your work. Press Tab while typing to auto-complete filenames and commands.

4. **Use relative and absolute paths appropriately:**  
   - Relative paths: `documents/file.txt` (relative to current directory)
   - Absolute paths: `/home/username/documents/file.txt` (full path from root)

5. **Practice in a safe environment** before using commands on important files.  

6. **Use descriptive names** for files and directories to make navigation easier.  

7. **Combine commands with pipes (`|`)** to create powerful workflows:  
   ```bash
   ls -la | grep "document"    # List files and filter for those containing "document"
   ```

## Next Steps

Now that you've learned these basic commands, practice using them regularly. Try creating directories, copying files, and exploring your filesystem. As you become more comfortable, you can learn about:

- Text processing commands (`grep`, `sed`, `awk`)
- File compression (`tar`, `gzip`)
- Network commands (`curl`, `wget`)
- Environment variables and shell scripting
- Advanced file operations and regular expressions

Remember, the command line is a powerful tool that becomes more intuitive with practice. Don't hesitate to use `--help` and `man` pages whenever you need clarification!