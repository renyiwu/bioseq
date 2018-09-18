#!/usr/bin/env bash


date >> bam_file_reads.txt
for i in $@
do
samtools view $i | echo -n $(($(wc -l)/4)) >> bam_file_reads.txt
echo -e "\t$i" >> bam_file_reads.txt
done
echo -e "\n\n" >> bam_file_reads.txt
