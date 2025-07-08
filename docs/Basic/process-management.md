# Process Management in Bash

## Introduction

In UNIX-like systems, a **process** is a running instance of a program. When you execute a command in the terminal, you create a process. Understanding how to manage these processes is crucial for effective system administration and development work.

This guide will teach you the essential commands and concepts for managing processes in bash, even if you've never used a UNIX environment before.

<h2 class="no-toc">Table of Content</h2>

[TOC]

## What is a Process?

A process is simply a program that is currently running on your system. Every time you run a command, open an application, or execute a script, you create a new process. Each process has:

- A unique **Process ID (PID)** - a number that identifies the process
- A **Parent Process ID (PPID)** - the ID of the process that created it
- **Memory allocation** - the amount of RAM the process uses
- **CPU time** - how much processor time the process has consumed

## Basic Process Commands

### `ps` - Display Running Processes

The `ps` command shows you information about currently running processes.

**Basic syntax:**
```bash
ps [options]
```

**Common examples:**

```bash
# Show processes for current user
ps

# Show all processes with detailed information
ps aux

# Show processes in a tree format (shows parent-child relationships)
ps -ef --forest
```

**Example output of `ps aux`:**
```
USER       PID  %CPU %MEM    VSZ   RSS TTY      STAT START   TIME COMMAND
john      1234  0.5  2.1  12345  4567 pts/0    Ss   09:30   0:01 bash
john      5678  0.0  0.8   8901  2345 pts/0    R+   09:35   0:00 ps aux
```

**Column explanations:**  
- **USER**: Who owns the process  
- **PID**: Process ID  
- **%CPU**: Percentage of CPU usage  
- **%MEM**: Percentage of memory usage  
- **COMMAND**: The command that started the process  

### `top` - Real-time Process Monitor

The `top` command displays processes in real-time, updating every few seconds.

```bash
# Start top (press 'q' to quit)
top

# Show processes for a specific user
top -u username
```

**Useful `top` commands while running:**  
- `q`: Quit top  
- `k`: Kill a process (you'll be prompted for the PID)  
- `h`: Show help  
- `M`: Sort by memory usage  
- `P`: Sort by CPU usage  

### `htop` - Enhanced Process Monitor

`htop` is an improved version of `top` with a more user-friendly interface.

```bash
# Start htop (may need to install first)
htop
```

**Navigation in htop:**
- Use arrow keys to navigate  
- `F9`: Kill selected process  
- `F10`: Quit htop  

## Process States

Processes can be in different states:

- **Running (R)**: Currently executing  
- **Sleeping (S)**: Waiting for an event (like user input)  
- **Stopped (T)**: Process has been stopped (paused)  
- **Zombie (Z)**: Process has finished but parent hasn't cleaned it up yet  

## Background and Foreground Processes

### Running Commands in Background

You can run commands in the background using the `&` symbol:

```bash
# Run a command in the background
long_running_command &

# Example: copy a large file in the background
cp large_file.txt backup_file.txt &
```

### `jobs` - List Active Jobs

The `jobs` command shows processes started from your current shell:

```bash
# List all jobs
jobs

# List jobs with process IDs
jobs -l
```

**Example output:**
```
[1]+  Running                 cp large_file.txt backup_file.txt &
[2]-  Stopped                 vim document.txt
```

### `fg` and `bg` - Foreground and Background Control

```bash
# Bring job number 1 to foreground
fg %1

# Send job number 1 to background
bg %1

# Bring the most recent job to foreground
fg
```

### Stopping and Resuming Processes

- **Ctrl+Z**: Stop (pause) the current foreground process  
- **Ctrl+C**: Terminate the current foreground process  

**Example workflow:**
```bash
# Start editing a file
vim document.txt

# Press Ctrl+Z to stop vim and return to shell
# [1]+  Stopped                 vim document.txt

# Continue editing in the background
bg %1

# Bring it back to foreground when needed
fg %1
```

## Process Termination

### `kill` - Terminate Processes

The `kill` command sends signals to processes to terminate them.

**Basic syntax:**
```bash
kill [signal] PID
```

**Common signals:**
- **SIGTERM (15)**: Polite termination request (default)  
- **SIGKILL (9)**: Force termination (cannot be ignored)  
- **SIGSTOP (19)**: Stop process  
- **SIGCONT (18)**: Continue stopped process  

**Examples:**
```bash
# Politely ask process 1234 to terminate
kill 1234

# Force kill process 1234
kill -9 1234

# Kill using signal name
kill -SIGTERM 1234

# Kill all processes with a specific name
killall firefox
```

### `killall` - Kill Processes by Name

```bash
# Kill all processes named "firefox"
killall firefox

# Kill all processes by a specific user
killall -u username
```

### `pkill` - Kill Processes by Pattern

```bash
# Kill processes whose names contain "chrome"
pkill chrome

# Kill processes owned by specific user
pkill -u username
```

## Process Monitoring and Information

### `pgrep` - Find Process IDs

```bash
# Find PIDs of processes named "firefox"
pgrep firefox

# Find PIDs with more details
pgrep -l firefox
```

### `lsof` - List Open Files

Shows which files are being used by which processes:

```bash
# List all open files
lsof

# List files opened by a specific process
lsof -p 1234

# List processes using a specific file
lsof /path/to/file

# List processes using a specific port
lsof -i :8080
```

## Practical Examples

### Example 1: Finding and Killing a Misbehaving Process

```bash
# 1. Find the process
ps aux | grep firefox

# 2. Kill it using PID (replace 1234 with actual PID)
kill 1234

# Alternative: kill by name
killall firefox
```

### Example 2: Running a Long Task in Background

```bash
# Start a backup process in background
tar -czf backup.tar.gz /home/user/documents &

# Check if it's still running
jobs

# Monitor system resources
top
```

### Example 3: Managing a Stopped Process

```bash
# Start a text editor
nano myfile.txt

# Stop it with Ctrl+Z
# [1]+  Stopped                 nano myfile.txt

# List stopped jobs
jobs

# Resume in background
bg %1

# Bring back to foreground when ready
fg %1
```

## Process Priority and Nice Values

Processes have priority levels that determine how much CPU time they get. The `nice` command allows you to start processes with different priorities.

```bash
# Start a process with lower priority (nice value 10)
nice -n 10 long_running_command

# Start a process with higher priority (requires sudo)
sudo nice -n -10 important_command
```

**Nice values range from -20 (highest priority) to 19 (lowest priority).**

### `renice` - Change Priority of Running Process

```bash
# Change priority of process 1234 to nice value 5
renice 5 1234

# Change priority of all processes by user
renice 10 -u username
```

## Best Practices

1. **Always try gentle termination first**: Use `kill` (SIGTERM) before `kill -9` (SIGKILL)  

2. **Monitor resource usage**: Use `top` or `htop` to identify resource-hungry processes  

3. **Use background processes for long tasks**: Don't let long-running commands block your terminal  

4. **Be careful with `kill -9`**: It forces termination without cleanup, which can cause data loss  

5. **Check before killing**: Make sure you're killing the right process by checking the PID  

## Common Scenarios and Solutions

### Scenario 1: Terminal Becomes Unresponsive
```bash
# Open another terminal and find the problematic process
ps aux | grep -i stuck_program

# Kill it
kill -9 PID
```

### Scenario 2: Too Many Background Jobs
```bash
# List all jobs
jobs

# Kill specific job
kill %1

# Kill all jobs
kill $(jobs -p)
```

### Scenario 3: Finding Resource-Heavy Processes
```bash
# Sort by CPU usage
ps aux --sort=-%cpu | head -10

# Sort by memory usage
ps aux --sort=-%mem | head -10
```

## Summary

Process management is a fundamental skill in UNIX systems. The key commands to remember are:

- `ps aux`: List all processes  
- `top`/`htop`: Monitor processes in real-time  
- `jobs`: List jobs started from current shell  
- `kill PID`: Terminate a process  
- `killall name`: Kill processes by name  
- `fg`/`bg`: Move jobs between foreground and background  
- `Ctrl+Z`: Stop current process  
- `Ctrl+C`: Terminate current process  

With these commands, you'll be able to effectively manage processes in any UNIX environment. Practice these commands in a safe environment to become comfortable with process management!