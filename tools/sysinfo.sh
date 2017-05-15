#!/bin/sh

# This script displays a series of system information
# that can be used to troubleshoot the system settings 
# prior to MFIX installation.
# To run the script , type the following at the prompt: sh ./sysinfo.sh
# To redirect the output to a file, type: sh ./sysinfo.sh > sysinfo.txt
# When requesting help related to MFIX installation, please send the
# file sysinfo.txt to the mfix-help mailing list at mfix-help@mfix.netl.doe.gov
# along with a description of the issue.


Print_Header()
{
# Displays a title centered between two lines
    echo ''
    echo ''
    echo ''   | awk '{d="#";for (i=1;i<80;i++) d=d "#" ;print d;}'
    echo "$1" | awk '{s="";for (i=1;i<(80-length($0))/2;i++) s=s " " ;print s $0;}'
    echo ''   | awk '{d="#";for (i=1;i<80;i++) d=d "#" ;print d;}'
}


Print_Header "System Information" 
cat /etc/issue
echo "Kernel name       : " `uname -s`
echo "Node hostname     : " `uname -n`
echo "Kernel release    : " `uname -r`
echo "Kernel version    : " `uname -v`
echo "Machine name      : " `uname -m`
echo "Processor type    : " `uname -p`
echo "Hardware platform : " `uname -i`
echo "Operating system  : " `uname -o`

Print_Header "CPU Info"
cat /proc/cpuinfo

Print_Header "Memory Info"
cat /proc/meminfo

Print_Header "Environment"
env

Print_Header "Linker"
ld -v

# Print_Header "GPU Info"
# cat /proc/driver/nvidia/version
# pgaccelinfo


Print_Header "End of script"
