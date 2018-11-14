#RNAs-seq workflow. 
#R Wu Sept 2017
# Revised on April 24, 2018. Added [-k 1] and [--mm] flag for hisat2. and added [-p] for featureCounts

#Hisat2 -> picard -> feartureCounts



#1, sequencing adapter trimming
trim_galore --paired --trim1  sample1_R1.fastq.gz sample1_R2.fastq.gz # for paired-end reads
# or below for single end reads
trim_galore --trim1 sample1_R1.fastq.gz

#2 hisat2
# Single end reads
for i in *.fastq.gz; do hisat2 -p 6 -x /path/to/hisat2/index/ -U $i | samtools view -bh -o ${i%.fastq.gz}.bam -; done

# For paired end reads,
for i in 1 2 3 4 5; do
hisat2 -p 6 -x ~/genomes/Mus_musculus/UCSC/mm10/Hisat2_Genome/genome -k 1 --mm -t --un-gz ./ht2reports/C"$i".k1_un.sam.gz --al-gz ./ht2reports/C"$i".k1_al.sam.gz --summary-file ./ht2reports/C"$i".k1_sum.txt --met-file ./ht2reports/C"$i".k1_met.txt -1 ${i}_R1.fastq.gz -2 ${i}_R2.fastq.gz | samtools view -bh -o ./bam/${i}.bam; done
# [-k 1] limits the number of maximum alignments to 1 for each read.
# [--mm] sets the memory mapping so many hisat2 instances can share the genome reference index. 




#3 sorting
for i in *.bam; do samtools sort -@ 6 -o ${i%.bam}.sorted.bam $i; done

#4 deduplicate by picard
for i in *.sorted.bam; do java -jar /path/to/picard.jar MarkDuplicates I=$i O=${i%.sorted.bam}.dedup.bam M=${i%.sorted.bam}.dedup.txt REMOVE_DUPLICATES=true; done

#or samtools
for i in *.sorted.bam; do samtools rmdup -s $i ${i%.sorted.bam}.rmdup.bam; done #Do NOT use -s for paired end data
for i in *.sorted.bam; do samtools rmdup $i ${i%.sorted.bam}.rmdup.bam; done # For paired end sequences.


#5 FeatureCounts --Single end
featureCounts --primary -T 8 -a /path/to/genes.gtf -o featurecounts.results.csv *dup.bam

#for paired end sequences, 
featureCounts --primary -p -T 8 -a /path/to/genes.gtf -o featurecounts.results.csv *dup.bam
# [-p] counts fragments instead of reads. i.e., a pair (of reads) will be counted as 1 (fragment) instead of 2 (reads).
# [-T 8] uses 8 threats. change the number accordingly.




