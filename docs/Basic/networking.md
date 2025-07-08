# Networking in Bash

## Introduction

Networking is a fundamental aspect of modern computing, and the bash shell provides powerful tools for network troubleshooting, monitoring, and administration. This guide introduces essential networking commands that every system administrator and developer should know.

## Table of Content

[TOC]

## Basic Network Information Commands

### `ifconfig` - Network Interface Configuration

The `ifconfig` command displays and configures network interfaces on your system.

**Basic usage:**
```bash
ifconfig
```

This shows all active network interfaces with their IP addresses, subnet masks, and other configuration details.

**Example output:**
```
eth0: flags=4163<UP,BROADCAST,RUNNING,MULTICAST>  mtu 1500
        inet 192.168.1.100  netmask 255.255.255.0  broadcast 192.168.1.255
        ether 08:00:27:12:34:56  txqueuelen 1000  (Ethernet)
```

**Specific interface:**
```bash
ifconfig eth0
```

**Configure an interface (requires root privileges):**
```bash
sudo ifconfig eth0 192.168.1.50 netmask 255.255.255.0
```


### `ip` - Modern Network Configuration Tool

The `ip` command is the modern replacement for several networking tools including `ifconfig`.

**Show all interfaces:**
```bash
ip addr show
# or shorter:
ip a
```

**Show specific interface:**
```bash
ip addr show eth0
```

**Show routing table:**
```bash
ip route show
# or shorter:
ip r
```

**Add a static route:**
```bash
sudo ip route add 192.168.2.0/24 via 192.168.1.1
```


## Network Connectivity Testing

### `ping` - Test Network Connectivity

The `ping` command sends ICMP echo requests to test connectivity to remote hosts.

**Basic ping:**
```bash
ping google.com
```

**Ping with count limit:**
```bash
ping -c 4 google.com
```

**Ping with specific interval:**
```bash
ping -i 2 google.com  # ping every 2 seconds
```

**Ping IPv6:**
```bash
ping6 google.com
```

**Example output:**
```
PING google.com (142.250.191.14) 56(84) bytes of data.
64 bytes from lhr35s10-in-f14.1e100.net (142.250.191.14): icmp_seq=1 ttl=119 time=12.3 ms
64 bytes from lhr35s10-in-f14.1e100.net (142.250.191.14): icmp_seq=2 ttl=119 time=11.8 ms
```


### `traceroute` - Trace Network Path

Shows the route packets take to reach a destination.

**Basic traceroute:**
```bash
traceroute google.com
```

**Traceroute with specific interface:**
```bash
traceroute -i eth0 google.com
```

**Example output:**
```
traceroute to google.com (142.250.191.14), 30 hops max, 60 byte packets
 1  192.168.1.1 (192.168.1.1)  1.234 ms  1.198 ms  1.167 ms
 2  10.0.0.1 (10.0.0.1)  5.432 ms  5.401 ms  5.378 ms
 3  203.0.113.1 (203.0.113.1)  12.567 ms  12.534 ms  12.501 ms
```



## DNS and Name Resolution

### `nslookup` - DNS Lookup Tool

Query DNS servers to resolve domain names.

**Basic lookup:**
```bash
nslookup google.com
```

**Reverse lookup (IP to domain):**
```bash
nslookup 8.8.8.8
```

**Query specific DNS server:**
```bash
nslookup google.com 8.8.8.8
```

**Query specific record type:**
```bash
nslookup -type=MX google.com  # Mail exchange records
nslookup -type=NS google.com  # Name server records
```


### `dig` - Advanced DNS Lookup

More powerful and flexible than `nslookup`.

**Basic dig:**
```bash
dig google.com
```

**Short answer:**
```bash
dig +short google.com
```

**Query specific record types:**
```bash
dig MX google.com      # Mail exchange records
dig NS google.com      # Name server records
dig AAAA google.com    # IPv6 addresses
```

**Reverse lookup:**
```bash
dig -x 8.8.8.8
```


### `host` - Simple DNS Lookup

**Basic lookup:**
```bash
host google.com
```

**All record types:**
```bash
host -a google.com
```



## Port and Service Information

### `netstat` - Network Statistics

Shows network connections, routing tables, and network interface statistics.

**Show all connections:**
```bash
netstat -a
```

**Show listening ports:**
```bash
netstat -l
```

**Show TCP connections:**
```bash
netstat -t
```

**Show UDP connections:**
```bash
netstat -u
```

**Show processes using ports:**
```bash
netstat -p
```

**Combined useful options:**
```bash
netstat -tulpn
# t = TCP, u = UDP, l = listening, p = process, n = numerical addresses
```

**Example output:**
```
Proto Recv-Q Send-Q Local Address           Foreign Address         State       PID/Program name
tcp        0      0 0.0.0.0:22              0.0.0.0:*               LISTEN      1234/sshd
tcp        0      0 127.0.0.1:3306          0.0.0.0:*               LISTEN      5678/mysqld
```


### `ss` - Socket Statistics (Modern Alternative)

**Show all sockets:**
```bash
ss -a
```

**Show listening sockets:**
```bash
ss -l
```

**Show TCP sockets:**
```bash
ss -t
```

**Show processes:**
```bash
ss -p
```

**Combined options:**
```bash
ss -tulpn
```


### `lsof` - List Open Files

Can show network connections by process.

**Show network connections:**
```bash
lsof -i
```

**Show connections on specific port:**
```bash
lsof -i :80
lsof -i :22
```

**Show connections by specific process:**
```bash
lsof -i -p 1234
```

## Network Testing and Troubleshooting

### `telnet` - Test Port Connectivity

**Test if a port is open:**
```bash
telnet google.com 80
telnet localhost 22
```

If successful, you'll see:
```
Trying 142.250.191.14...
Connected to google.com.
Escape character is '^]'.
```


### `nc` (netcat) - Network Swiss Army Knife

**Test port connectivity:**
```bash
nc -zv google.com 80
```

**Listen on a port:**
```bash
nc -l 8080
```

**Send data to a port:**
```bash
echo "Hello" | nc localhost 8080
```

**Port scanning:**
```bash
nc -zv google.com 80-90
```


### `curl` - Transfer Data from Servers

**Basic HTTP request:**
```bash
curl http://google.com
```

**Show headers:**
```bash
curl -I http://google.com
```

**Follow redirects:**
```bash
curl -L http://google.com
```

**Download a file:**
```bash
curl -O http://example.com/file.txt
```

**POST data:**
```bash
curl -X POST -d "key=value" http://example.com/api
```


### `wget` - Download Files

**Download a file:**
```bash
wget http://example.com/file.txt
```

**Download recursively:**
```bash
wget -r http://example.com/
```

**Download in background:**
```bash
wget -b http://example.com/largefile.zip
```


## Network Monitoring

### `iftop` - Real-time Network Usage

Shows bandwidth usage by connection.

```bash
sudo iftop
```


### `netstat` for Monitoring

**Show network statistics:**
```bash
netstat -s
```

**Monitor connections continuously:**
```bash
watch netstat -tulpn
```


### `arp` - ARP Table Management

**Show ARP table:**
```bash
arp -a
```

**Add static ARP entry:**
```bash
sudo arp -s 192.168.1.100 00:11:22:33:44:55
```

**Delete ARP entry:**
```bash
sudo arp -d 192.168.1.100
```

## Common Network Troubleshooting Workflow

When troubleshooting network issues, follow this systematic approach:

1. **Check local network configuration:**
   ```bash
   ip addr show
   ip route show
   ```

2. **Test local connectivity:**
   ```bash
   ping 127.0.0.1        # Test loopback
   ping 192.168.1.1      # Test gateway
   ```

3. **Test DNS resolution:**
   ```bash
   nslookup google.com
   dig google.com
   ```

4. **Test external connectivity:**
   ```bash
   ping 8.8.8.8          # Test external IP
   ping google.com       # Test with DNS
   ```

5. **Check service ports:**
   ```bash
   netstat -tulpn | grep :80
   telnet localhost 80
   ```

6. **Trace network path:**
   ```bash
   traceroute google.com
   ```

## Practical Examples

### Example 1: Check if a Web Server is Running

```bash
# Check if port 80 is listening
netstat -tulpn | grep :80

# Test connectivity to the port
telnet localhost 80

# Check with curl
curl -I http://localhost
```


### Example 2: Troubleshoot SSH Connection

```bash
# Check if SSH service is running
netstat -tulpn | grep :22

# Test SSH port connectivity
telnet server.example.com 22

# Check SSH with verbose output
ssh -v user@server.example.com
```


### Example 3: Monitor Network Usage

```bash
# Show current connections
netstat -i

# Monitor bandwidth usage
sudo iftop

# Show network statistics
netstat -s
```


## Security Considerations

- Always be cautious when using networking commands on production systems
- Some commands require root privileges (use `sudo`)
- Network scanning tools like `nmap` should only be used on systems you own or have permission to test
- Be aware that some network tools generate significant traffic


## Tips for Beginners

1. **Start with basic commands:** Begin with `ping`, `ifconfig`, and `netstat`
2. **Use help options:** Most commands support `-h` or `--help` flags
3. **Read man pages:** Use `man command` for detailed documentation
4. **Practice in safe environments:** Use virtual machines or isolated networks for learning
5. **Combine commands:** Use pipes (`|`) and redirection (`>`) to combine tools effectively


## Conclusion

These networking commands form the foundation of network troubleshooting and administration in Unix-like systems. Practice these commands regularly, and gradually incorporate more advanced options as you become comfortable with the basics. Remember that networking is both an art and a science - systematic troubleshooting combined with experience will make you proficient in diagnosing and resolving network issues.