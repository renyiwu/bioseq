#!/usr/bin/Rscript

# R Wu (renyiwu@gmail.com)
# April 2018

# A wrapper for adding gene features to DMR data with Chipseeker package.

version <- '0.2'
copyright <- 'Copyright (C) 2018 R Wu'

changelog <- function() {
  cat('---Change log---
      v0.2:
      set level default to "gene"
      ')
  q()
}



printVersion <- function() {
  cat('ChipAnno.r Version', version, '\n')
  cat(copyright, '\n')
  q()
}

usage <- function() {
  cat('Usage: Rscript ChipAnno.r  -x <genome>  -i <input>  -o <output>
	-x <genome>  Gene annotation file with file name ending in .gtf
	-i <input>    input file to read.
 	-o <output>   Output file to save. Keep the same as -i to override.
')
  q()
}



# default args/parameters
infile <- outfile <- genome <- NULL

# get CL args
args <- commandArgs(trailingOnly=T)
i <- 1
while (i <= length(args)) {
  if (substr(args[i], 1, 1) == '-') {
    if (args[i] == '-h' || args[i] == '--help') {
      usage()
    } else if (args[i] == '--version') {
      printVersion()
    } else if (args[i] == '-v') {
      printVersion()
    } else if (i < length(args)) {
      if (args[i] == '-i') {
        infile <- args[i + 1]
      } else if (args[i] == '-o') {
        outfile <- args[i + 1]
      } else if (args[i] == '-x') {
        genome <- args[i + 1]
      } else {
        cat('Error! Unknown parameter:', args[i], '\n')
        usage()
      }
      i <- i + 1
    } else {
      cat('Error! Unknown parameter with no arg:', args[i], '\n')
      usage()
    }
  } else {
    temp <- strsplit(args[i], '[ ,]')[[1]]
    names <- c(names, list(temp[nchar(temp) > 0]))
  }
  i <- i + 1
}

# check for parameter errors
if (is.null(infile) || is.null(outfile)) {
  cat('Error! Must specify input and output files\n')
  usage()
}
if (is.null(genome)) {
  cat('Error! Must specify genome annotation file\n')
  usage()
}

#####
# genome <- "~/genomes/Homo_Sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
# infile <- "~/Methyl-seq_rucdr/combined_fr5c5_Yen.csv"
# outfile <- "~/Methyl-seq_rucdr/combined_fr5c5_yen_test1.csv"
# #####



if (!suppressMessages(suppressWarnings(require("GenomicFeatures")))) {
  stop('Required package "GenomicFeatures" not installed.\n',
       '  For installation information, please see:\n',
       '  https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html\n')
}
if (!suppressMessages(suppressWarnings(require("ChIPseeker")))) {
  stop('Required package "ChIPseeker" not installed.\n',
       '  For installation information, please see:\n',
       '  https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html\n')
}

txdb <- makeTxDbFromGFF(genome, format='gtf')
peak <- readPeakFile(infile)

peakAnno <- annotatePeak(peak,
                         TxDb=txdb,
                         level = "gene")
out <- as.data.frame(peakAnno)
out$strand <- NULL
out$geneChr <- NULL

# head(peakAnno@anno, 1)
# # add annotations
# peak2 <- read.table(infile, header=T)
# peak2 <- peak2[peak2$chr!="chrM",] #This is necessary for some files...
# peak2$feature <- peakAnno@anno$annotation
# peak2$distance <- peakAnno@anno$distanceToTSS
# peak2$gene <- peakAnno@anno$geneId

write.table(out, outfile, sep='\t', quote=F, row.names = F)
