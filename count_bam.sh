#!/usr/bin/env bash
# R Wu Sep 2018
date
SECONDS=0
date > bam_file_reads.txt
for i in $@
do
samtools view $i | echo -n $(wc -l) >> bam_file_reads.txt
echo -e "\t$i" >> bam_file_reads.txt
done
echo -e "\n\n" >> bam_file_reads.txt
date
echo
echo "Finished in $SECONDS seconds."
