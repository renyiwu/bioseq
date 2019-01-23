#RNAs-seq workflow. 
#R Wu Sept 2017
# Revised on April 24, 2018. Added [-k 1] and [--mm] flag for hisat2. and added [-p] for featureCounts

#Hisat2 -> picard -> feartureCounts



#1, sequencing adapter trimming
trim_galore --paired --trim1  sample1_R1.fastq.gz sample1_R2.fastq.gz # for paired-end reads
# parallel
time parallel -j 8 trim_galore --paired --trim1 {}R1_001.fastq.gz {}R2_001.fastq.gz ::: $(ls *R1*fastq.gz | sed 's/R1_001.fastq.gz//') 
# output file name examples:
# sample1-R[12]_001.fastq.gz --->
# sample1_R1_001_val_1.fq.gz
# sample1_R2_001_val_2.fq.gz



# or below for single end reads
trim_galore --trim1 sample1_R1.fastq.gz
# sample1.fastq.gz --->
# sample1_trimmed.fq.gz


#2 hisat2
# Single end reads
for i in *.fastq.gz; do hisat2 -p 6 -x /path/to/hisat2/index/ -k 1 -U $i | samtools view -bh -o ${i%.fastq.gz}.bam; done

# For paired end reads,
for i in 1 2 3 4 5; do
hisat2 -p 6 -x ~/genomes/Mus_musculus/UCSC/mm10/Hisat2_Genome/genome -k 1 --mm -t --un-gz ./ht2reports/C"$i".k1_un.sam.gz --al-gz ./ht2reports/C"$i".k1_al.sam.gz --summary-file ./ht2reports/C"$i".k1_sum.txt --met-file ./ht2reports/C"$i".k1_met.txt -1 ${i}_R1.fastq.gz -2 ${i}_R2.fastq.gz | samtools view -bh -o ./bam/${i}.bam; done
# [-k 1] limits the number of maximum alignments to 1 for each read.
# [--mm] sets the memory mapping so many hisat2 instances can share the genome reference index. 

# or for trimmed fastq files,
# do
for i in *val_1.fq.gz; do ( hisat2 -p 6 -x ~/genomes/Mus_musculus/UCSC/mm10/Hisat2_Genome/genome -k 1 -1 ${i} -2 ${i%R1_001_val_1.fq.gz}R2_001_val_2.fq.gz | samtools view -bh -o ../bam/${i%_R1_001_val_1.fq.gz}.bam ); done


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

############################
# one line for PE reads on SOP-1482
time ( parallel -j 8 trim_galore --paired --trim1 {}R1_001.fastq.gz {}R2_001.fastq.gz ::: $(ls *R1*fastq.gz | sed 's/R1_001.fastq.gz//') && for i in *val_1.fq.gz; do ( hisat2 -p 6 -x ~/genomes/Mus_musculus/UCSC/mm10/Hisat2_Genome/genome -k 1 -1 ${i} -2 ${i%R1_001_val_1.fq.gz}R2_001_val_2.fq.gz | samtools view -bh -o ${i%_R1_001_val_1.fq.gz}.bam ); done && for i in *.bam; do samtools sort -@ 6 -o ${i%.bam}.sorted.bam $i; done && for i in *.sorted.bam; do java -jar ~/tools/picard/picard.jar MarkDuplicates I=$i O=${i%.sorted.bam}.dedup.bam M=${i%.sorted.bam}.dedup.txt REMOVE_DUPLICATES=true; done && featureCounts --primary -p -T 8 -a ~/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf  -o featurecounts.results.csv *dup.bam )


# one line for SE on SOP-1482

time ( parallel -j 8 trim_galore --trim1 {} ::: *fastq.gz && for i in *_trimmed.fq.gz; do ( hisat2 -p 6 -x ~/genomes/Mus_musculus/UCSC/mm10/Hisat2_Genome/genome -U $i | samtools view -bh -o ${i%_trimmed.fq.gz}.bam ); done && for i in *.bam; do samtools sort -@ 6 -o ${i%.bam}.sorted.bam $i; done && for i in *.sorted.bam; do java -jar ~/tools/picard/picard.jar MarkDuplicates I=$i O=${i%.sorted.bam}.dedup.bam M=${i%.sorted.bam}.dedup.txt REMOVE_DUPLICATES=true; done && featureCounts --primary -p -T 8 -a ~/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf  -o featurecounts.results.csv *dup.bam )


