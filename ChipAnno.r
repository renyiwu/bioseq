#!/usr/bin/Rscript

# R Wu (renyiwu@gmail.com)
# April 2018

# A wrapper for adding gene features to DMR data with Chipseeker package.

version <- '0.1'
copyright <- 'Copyright (C) 2018 R Wu'

printVersion <- function() {
  cat('ChipAnno.r Version', version, '\n')
  cat(copyright, '\n')
  q()
}

usage <- function() {
  cat('Usage: Rscript ChipAnno.r  -x <genome>  -i <input>  -o <output> 
    -i <input>    File to be read
    -o <output>   Output file to save. Keep the same as -i to override.
    -x <genome>  Gene annotation file with file name ending in .gtf
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

##
library(GenomicFeatures)
library(ChIPseeker)
txdb <- makeTxDbFromGFF(genome, format='gtf')
peak <- readPeakFile(infile)
peakAnno <- annotatePeak(peak, TxDb=txdb)
# add annotations
peak2 <- read.table(infile, header=T)
peak2 <- peak2[peak2$chr!="chrM",] #This is necessary for some files...
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$gene <- peakAnno@anno$geneId


write.table(peak2, outfile, sep='\t', quote=F, row.names=F)
