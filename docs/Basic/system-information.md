# System Information in Bash

## Introduction

Understanding your system's configuration, resources, and current state is crucial for system administration, troubleshooting, and optimization. This guide covers essential bash commands for gathering system information in Unix-like environments. These commands will help you monitor hardware, software, processes, and system performance.

## Table of Content

[TOC]

## System Overview Commands

### `uname` - System Information

The `uname` command provides basic system information including kernel name, version, and architecture.

**Basic usage:**
```bash
uname
```
Output: `Linux`

**Show all information:**
```bash
uname -a
```
**Example output:**
```
Linux ubuntu 5.15.0-72-generic #79-Ubuntu SMP Wed Apr 19 08:22:18 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux
```

**Common options:**
```bash
uname -s    # Kernel name (Linux, Darwin, etc.)
uname -r    # Kernel release version
uname -v    # Kernel version
uname -m    # Machine hardware name (x86_64, arm64, etc.)
uname -p    # Processor type
uname -o    # Operating system
```

### `hostname` - System Name

**Show hostname:**
```bash
hostname
```

**Show fully qualified domain name:**
```bash
hostname -f
```

**Show IP address:**
```bash
hostname -I
```

### `whoami` and `id` - User Information

**Current user:**
```bash
whoami
```

**Detailed user information:**
```bash
id
```
**Example output:**
```
uid=1000(john) gid=1000(john) groups=1000(john),4(adm),24(cdrom),27(sudo),30(dip)
```

**Show all logged-in users:**
```bash
who
```

**Show current user and login time:**
```bash
w
```

## Hardware Information

### `lscpu` (Linux) / `sysctl` (macOS) - CPU Information

**Detailed CPU information (Linux):**
```bash
lscpu
```

**Example output:**
```
Architecture:                    x86_64
CPU op-mode(s):                  32-bit, 64-bit
Byte Order:                      Little Endian
CPU(s):                          4
On-line CPU(s) list:             0-3
Thread(s) per core:              2
Core(s) per socket:              2
Socket(s):                       1
NUMA node(s):                    1
Vendor ID:                       GenuineIntel
CPU family:                      6
Model:                           142
Model name:                      Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz
```

On macOS, there is no equivalent to `lscpu`. The command `sysctl` allows you to obtain more detailed information in isolation. It combines the `lscpu` and `lshw` commands under Linux. You can chain commands together to obtain the equivalent :

```bash
#!/bin/bash
echo "Architecture      : $(uname -m)"
echo "CPU(s)            : $(sysctl -n hw.ncpu)"
echo "Core(s) per CPU   : $(sysctl -n machdep.cpu.core_count)"
echo "Thread(s) per CPU : $(sysctl -n machdep.cpu.thread_count)"
echo "Model name        : $(sysctl -n machdep.cpu.brand_string)"
echo "CPU MHz           : $(sysctl -n hw.cpufrequency | awk '{print $1/1000000 " MHz"}')"
echo "L1 Cache          : $(sysctl -n hw.l1dcachesize) bytes"
echo "L2 Cache          : $(sysctl -n hw.l2cachesize) bytes"
echo "L3 Cache          : $(sysctl -n hw.l3cachesize 2>/dev/null || echo "N/A")"
echo "Memory            : $(( $(sysctl -n hw.memsize) / 1024 / 1024 )) MB"
```

### `lshw` - Hardware Information

**Complete hardware information (requires root):**
```bash
sudo lshw
```

**Short format:**
```bash
sudo lshw -short
```

**Specific hardware class:**
```bash
sudo lshw -class processor
sudo lshw -class memory
sudo lshw -class disk
sudo lshw -class network
```

### `lspci` (Linux) / `system_profiler` (macOS) - PCI Devices

**List all PCI devices:**
```bash
lspci                           # Linux
system_profiler SPPciDataType   # macOS
```

**Verbose output:**
```bash
lspci -v
```

**Tree format:**
```bash
lspci -t
```

### `lsusb` (Linux) / `system_profiler` (macOS) - USB Devices

**List USB devices:**
```bash
lsusb                           # Linux
system_profiler SPUSBDataType   # macOS
system_profiler SPUSBDataType | grep -E 'Product ID|Vendor ID|Speed|Location ID|Serial Number|Manufacturer'
```

**Verbose output:**
```bash
lsusb -v
```

### `lsblk` (Linux) / `diskutil` (macOS) - Block Devices

**List block devices:**
```bash
lsblk               # Linux
diskutil list       # macOS
```

**Example output (Linux):**
```
NAME   MAJ:MIN RM   SIZE RO TYPE MOUNTPOINT
sda      8:0    0   500G  0 disk 
├─sda1   8:1    0   512M  0 part /boot/efi
├─sda2   8:2    0     1G  0 part /boot
└─sda3   8:3    0 498.5G  0 part /
```

**Example output (macOS):**
```
/dev/disk0 (internal, physical):
   #:                       TYPE NAME                    SIZE       IDENTIFIER
   0:      GUID_partition_scheme                        *1.0 TB     disk0
   1:             Apple_APFS_ISC Container disk1         524.3 MB   disk0s1
   2:                 Apple_APFS Container disk3        2994.7 GB   disk0s2
   3:        Apple_APFS_Recovery Container disk2         5.4 GB     disk0s3

/dev/disk3 (synthesized):
   #:                       TYPE NAME                    SIZE       IDENTIFIER
   0:      APFS Container Scheme -                      +994.7 GB   disk3
                                 Physical Store disk0s2
   1:                APFS Volume Macintosh HD            11.3 GB    disk3s1
   2:              APFS Snapshot com.apple.os.update-... 11.3 GB    disk3s1s1
   3:                APFS Volume Preboot                 7.1 GB     disk3s2
   4:                APFS Volume Recovery                1.0 GB     disk3s3
   5:                APFS Volume Data                   2901.6 GB   disk3s5
   6:                APFS Volume VM                      2.1 GB     disk3s6
```


**Show filesystem information:**
```bash
lsblk -f
```

## Memory Information

### `free` (Linux) / `vm_stat` (macOS) - Memory Usage

**Show memory usage:**
```bash
free        # Linux
vm_stat     # macOS
```

**Human-readable format:**
```bash
free -h
```

**Example output:**
```
               total        used        free      shared  buff/cache   available
Mem:            7.7Gi       2.1Gi       3.2Gi       234Mi       2.4Gi       5.1Gi
Swap:           2.0Gi          0B       2.0Gi
```

**Show memory in different units:**
```bash
free -m    # Megabytes
free -g    # Gigabytes
free -k    # Kilobytes
```

### `/proc/meminfo` - Detailed Memory Information

**View detailed memory information:**
```bash
cat /proc/meminfo
```

**Show specific memory information:**
```bash
grep MemTotal /proc/meminfo
grep MemAvailable /proc/meminfo
grep SwapTotal /proc/meminfo
```

## Storage Information

### `df` - Disk Space Usage

**Show disk usage:**
```bash
df
```

**Human-readable format:**
```bash
df -h
```

**Example output:**
```
Filesystem      Size  Used Avail Use% Mounted on
/dev/sda3       489G  156G  308G  34% /
/dev/sda1       511M  6.1M  505M   2% /boot/efi
/dev/sda2       974M  203M  704M  23% /boot
```

**Show filesystem type:**
```bash
df -T
```

**Show inodes usage:**
```bash
df -i
```

### `du` - Directory Space Usage

**Show directory sizes:**
```bash
du -h
```

**Show summary of current directory:**
```bash
du -sh
```

**Show sizes of subdirectories:**
```bash
du -h --max-depth=1
```

**Show largest directories:**
```bash
du -h | sort -hr | head -10
```

### `fdisk` - Disk Partitions

**List disk partitions (requires root):**
```bash
sudo fdisk -l
```

**Show partition table for specific disk:**
```bash
sudo fdisk -l /dev/sda
```

## Process Information

### `ps` - Process Status

**Show all processes:**
```bash
ps aux
```

**Show processes in tree format:**
```bash
ps auxf
```

**Show specific user processes:**
```bash
ps -u username
```

**Show processes with specific command:**
```bash
ps aux | grep firefox
```

### `top` - Real-time Process Monitor

**Interactive process monitor:**
```bash
top
```

**Key shortcuts in top:**  
- `q` - quit  
- `k` - kill process  
- `r` - renice process  
- `M` - sort by memory usage  
- `P` - sort by CPU usage  
- `h` - help  

### `htop` - Enhanced Process Monitor

**Enhanced interactive process monitor:**
```bash
htop
```

*Note: May need to install with `sudo apt install htop` (Linux) or `brew install htop` (macOS; need to [install homebrew](https://brew.sh) first)*

### `pstree` - Process Tree

**Show process tree:**
```bash
pstree
```

**Show with PIDs:**
```bash
pstree -p
```

**Show specific user's processes:**
```bash
pstree username
```

## System Load and Performance

### `uptime` - System Uptime and Load

**Show system uptime and load:**
```bash
uptime
```

**Example output:**
```
 14:23:45 up 2 days,  3:45,  2 users,  load average: 0.15, 0.20, 0.18
```

### `vmstat` (Linux) / `vm_stat` (macOS) - Virtual Memory Statistics

**Show system statistics:**
```bash
vmstat
```

**Continuous monitoring (every 2 seconds):**
```bash
vmstat 2
```

**Show statistics 5 times with 2-second intervals:**
```bash
vmstat 2 5
```

### `iostat` - I/O Statistics

**Show I/O statistics:**
```bash
iostat
```

**Continuous monitoring:**
```bash
iostat 2
```

*Note: May need to install with `sudo apt install sysstat`*

### `sar` (Linux) - System Activity Reporter

**Show CPU usage:**
```bash
sar -u
```

**Show memory usage:**
```bash
sar -r
```

**Show I/O statistics:**
```bash
sar -b
```

## System Services and Processes

### `systemctl` (Linux) / `launchctl` (macOS) - System Service Manager
On macOS, there is no direct equivalent to systemctl because :  
	- macOS uses a different init system, based on launchd.  
	- Services (daemons) are managed via launchctl, which is the command-line tool for interacting with launchd.  

**List all services:**
```bash
systemctl list-units --type=service
launchctl list
```

**Check service status:**
```bash
systemctl status nginx
systemctl status ssh
```

**Show active services:**
```bash
systemctl list-units --type=service --state=active
```

**Show failed services:**
```bash
systemctl list-units --type=service --state=failed
```

### `jobs` - Background Jobs

**Show current jobs:**
```bash
jobs
```

**Show jobs with PIDs:**
```bash
jobs -l
```

## System Files and Configuration

### `/proc` (Linux) Filesystem - System Information

The `/proc` filesystem provides real-time system information:

**CPU information:**
```bash
cat /proc/cpuinfo
```

**Memory information:**
```bash
cat /proc/meminfo
```

**System version:**
```bash
cat /proc/version
```

**Current processes:**
```bash
cat /proc/loadavg
```

**Uptime:**
```bash
cat /proc/uptime
```

**Mounted filesystems:**
```bash
cat /proc/mounts
```

### System Files

**Operating system information:**
```bash
cat /etc/os-release     # Linux
sw_vers                 # macOS
```

**Distribution information:**
```bash
cat /etc/issue          # Linux
sw_vers -productName    # macOS
```

**System timezone:**
```bash
cat /etc/timezone           # Linux
systemsetup -gettimezone    # macOS
```

**DNS configuration:**
```bash
cat /etc/resolv.conf        # Linux
scutil --dns                # macOS
```

**Network interfaces:**
```bash
cat /etc/network/interfaces             # Linux
networksetup -listallnetworkservices    # macOS    
```

## Environment and Shell Information

### `env` - Environment Variables

**Show all environment variables:**
```bash
env
```

**Show specific variable:**
```bash
echo $PATH
echo $HOME
echo $USER
echo $SHELL
```

**Show shell information:**
```bash
echo $0
echo $BASH_VERSION
```

### `which` and `whereis` - Command Locations

**Find command location:**
```bash
which python3
which gcc
which ls
```

**Find command and manual locations:**
```bash
whereis python3
whereis gcc
```

### `history` - Command History

**Show command history:**
```bash
history
```

**Show last 10 commands:**
```bash
history 10
```

## System Monitoring Scripts

### Creating a System Info Script

Here's a simple script to gather system information:

```bash
#!/bin/bash
# system_info.sh - System Information Script

echo "=== SYSTEM INFORMATION ==="
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo "User: $(whoami)"
echo "Uptime: $(uptime)"
echo

echo "=== HARDWARE INFORMATION ==="
echo "CPU: $(lscpu | grep 'Model name' | cut -d':' -f2 | xargs)"
echo "Memory: $(free -h | grep 'Mem:' | awk '{print $2}')"
echo "Disk: $(df -h / | tail -1 | awk '{print $2}')"
echo

echo "=== SYSTEM LOAD ==="
echo "CPU Usage: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)"
echo "Memory Usage: $(free | grep Mem | awk '{printf("%.2f%%"), $3/$2 * 100.0}')"
echo "Disk Usage: $(df -h / | tail -1 | awk '{print $5}')"
echo

echo "=== NETWORK INFORMATION ==="
echo "IP Address: $(hostname -I | awk '{print $1}')"
echo "Active Connections: $(netstat -an | grep ESTABLISHED | wc -l)"
```

Save this as `system_info.sh`, make it executable, and run it:

```bash
chmod +x system_info.sh
./system_info.sh
```

## Practical Examples

### Example 1: Check System Resources Before Installing Software

```bash
# Check available memory
free -h

# Check available disk space
df -h

# Check CPU information
lscpu | grep -E 'Model name|CPU\(s\)|Architecture'

# Check if system is 64-bit
uname -m
```

### Example 2: Monitor System Performance

```bash
# Check current system load
uptime

# Monitor processes using most CPU
ps aux --sort=-%cpu | head -10

# Monitor processes using most memory
ps aux --sort=-%mem | head -10

# Check I/O usage
iostat 1 3
```

### Example 3: Troubleshoot High Resource Usage

```bash
# Find processes using high CPU
top -o %CPU

# Find processes using high memory
top -o %MEM

# Check disk I/O
iotop  # if available

# Check network connections
netstat -tuln
```

### Example 4: System Health Check

```bash
# Check disk usage
df -h | grep -E '(8[0-9]|9[0-9])%'

# Check memory usage
free -h | grep -E 'Mem:|Swap:'

# Check system load
uptime

# Check failed services
systemctl --failed

# Check system logs for errors
journalctl -p err --since today
```

## Common System Information One-liners

**Quick system overview:**
```bash
echo "$(uname -n) running $(uname -s) $(uname -r) on $(uname -m)"
```

**Memory usage percentage:**
```bash
free | grep Mem | awk '{printf("Memory Usage: %.2f%%\n", $3/$2 * 100.0)}'
```

**CPU usage:**
```bash
grep 'cpu ' /proc/stat | awk '{usage=($2+$4)*100/($2+$4+$5)} END {print usage "%"}'
```

**Disk usage of root partition:**
```bash
df -h / | tail -1 | awk '{print "Root partition: " $5 " used of " $2}'
```

**Number of processes:**
```bash
ps aux | wc -l
```

**System architecture:**
```bash
uname -m
```

## Tips for Beginners

1. **Start with basic commands:** Begin with `uname`, `whoami`, `free`, and `df`
2. **Use help options:** Most commands support `-h` or `--help` flags
3. **Read man pages:** Use `man command` for detailed documentation
4. **Practice regularly:** Set up a cron job to run system info commands daily
5. **Create aliases:** Make shortcuts for frequently used command combinations
6. **Monitor trends:** Run the same commands at different times to understand patterns
7. **Combine commands:** Use pipes (`|`) to filter and format output

## Useful Aliases

Add these to your `.bashrc` file:

```bash
alias sysinfo='uname -a && free -h && df -h'
alias meminfo='free -h && cat /proc/meminfo | head -20'
alias cpuinfo='lscpu | grep -E "Model name|CPU\(s\)|Architecture"'
alias diskinfo='df -h && lsblk'
alias netinfo='hostname -I && netstat -tuln'
```

## Security Considerations

- Some commands require root privileges (use `sudo`)
- Be careful when running commands that can modify system settings
- Monitor system logs regularly for security issues
- Don't share detailed system information publicly
- Use system monitoring to detect unusual activity

## Conclusion

Understanding system information is fundamental for effective system administration and troubleshooting. These commands provide comprehensive insights into your system's hardware, software, performance, and configuration. Practice these commands regularly to become proficient in system monitoring and maintenance.

Remember that system information gathering is an ongoing process - the more you understand about your system's normal behavior, the better you'll be at identifying and resolving issues when they arise.