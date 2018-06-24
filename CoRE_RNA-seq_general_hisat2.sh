# RNAs-seq workflow for Konglab, revised for nextflow implementation 
# R Wu and others
# June 2018

#Hisat2 -> picard -> feartureCounts

#1, FastQC
fastqc --extract --outdir=./FastQC_reports/ *.gz

#2, cutadapt

#3, Hisat2
#3.1 single end reads
for i in *.fastq.gz; do hisat2 -p 6 -x /path/to/hisat2/index/ -U $i | samtools view -bh -o ${i%.fastq.gz}.bam -; done

#3.2 paired-end reads
for i in {1..10};
do hisat2 -p 16 -x ~/genomes/Mus_musculus/UCSC/mm10/Hisat2_Genome/genome -k 1 --mm -t --un-gz ./ht2reports/"$i"_un.sam.gz --al-gz ./ht2reports/"$i"_al.sam.gz --summary-file ./ht2reports/"$i"_sum.txt --met-file ./ht2reports/"$i"_met.txt -1 "$i"_R1.fastq.gz -2 "$i"_R2.fastq.gz | \
samtools view -bh -o ./bam/"$i".bam; \
done
# [-k 1] limits the number of maximum alignments to 1 for each read.
# [--mm] sets the memory mapping so many hisat2 instances can share the genome reference index. 

#3.3 Tophat2 


#4 sort
for i in *.bam; do samtools sort -@ 6 -o ${i%.bam}.sorted.bam $i; done

#5 deduplicate with picard
for i in *.sorted.bam; do java -jar ~/tools/picard-2.10.5/picard.jar MarkDuplicates I=$i O=${i%.sorted.bam}.dedup.bam M=${i%.sorted.bam}.dedup.txt REMOVE_DUPLICATES=true; done

#6 FeatureCounts
#6.1 single end reads
featureCounts --primary -T 8 -a /path/to/genes.gtf -o featurecounts.results.csv *dup.bam


#6.2 for paired end reads, 
featureCounts --primary -p -T -a /path/to/genes.gtf -o featurecounts.results.csv *dup.bam
# [-p] counts fragments instead of reads. i.e., a pair (of reads) will be counted as 1 (fragment) instead of 2 (reads).






### references

# R1 concatnate bam files
samtools cat -o RW04.bam RW04?.bam

# R2 concatnate fastq files. eg, RW041 and RW042.
zcat RW04?.fastq.gz | gzip >RW04.fastq.gz

# R3 Rename (shorten) fastq file names if needed.
rename 's/_S.*gz/.fastq.gz/' *.fastq.gz
#eg, C75_S1_R1_001_.fastq.gz to C75.fastq.gz


# R4 Count reads on fastq
for i in RW??.fastq.gz; do zcat $i | echo $((`wc -l`/4)); done



