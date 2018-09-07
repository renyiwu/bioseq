#!/bin/bash 
#  v 0.3
# Renyi Wu, 2018
#  
# Usage: in terminal, cd to the path of fastq files, then run "bash /path/to/fastqc-0.x.sh"
# This file needs not to be in the same directory of your fastq files.

# cd "$(dirname "$0")" #or cd "${0%/*}" or cd "${0%/*}"
#Run FastQC,optional

[ -d FastQC_reports_auto ] || mkdir FastQC_reports_auto  # make directory if not exists.
if [ -z "$@" ]
then
 fastqc --extract --outdir=./FastQC_reports_auto/ *.gz
else
 fastqc --extract --outdir=./FastQC_reports_auto/ $@
fi
 
#echo "done"
#echo "counting reads"
date >fastqc_reads_for_all.txt
for i in ./FastQC_reports_auto/*/fastqc_data.txt
do
grep "Filename" "$i" >>fastqc_reads_for_all.txt
grep "Total Sequences" "$i" >>fastqc_reads_for_all.txt
grep ">>Sequence Duplication Levels" "$i" >>fastqc_reads_for_all.txt
grep "#Total Deduplicated Percentage" "$i" >>fastqc_reads_for_all.txt
grep "#Duplication Level" "$i" >>fastqc_reads_for_all.txt
#grep "^1\t" "$i" >>fastqc_reads_for_all.txt
awk '/^1\t/{n++}n==5{print; exit}' "$i" >>fastqc_reads_for_all.txt
echo "\n" >>fastqc_reads_for_all.txt
done
#read -n1 -r -p "Press any key to quit" key
