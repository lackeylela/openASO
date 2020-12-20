# Title     : TODO
# Objective : TODO
# Created by: james
# Created on: 12/20/2020


library("stringr")
library(insect)

# defining here for now but don't know that this is actually the right place for these to be defined (seems like they should actually be environment variables)
# create an EnsDB object containing genomic position information.
dbfile <- system.file("extdata/EnsDb.Hsapiens.v86.sqlite", package = "EnsDb.Hsapiens.v86")
db <- EnsDb(dbfile)

# this is actually something used down below but felt I should move it to a more "global" section
dna <- ensembldb::getGenomeTwoBitFile(db)


# returns the cds of a transcript given the dna, db, and transcript name
findCds <- function(dna, db, transcript) {

  # return null for transcripts that do not return a valid cds
  tryCatch({

    # set up filter for specified transcript
    txflt <- AnnotationFilter(~ tx_id == transcript)

    # extract the cds ( coding region only ) for the specified transcript above
    cds <- cdsBy(db, filter=txflt)
    cds <- extractTranscriptSeqs(dna, cds)
    cds <- toString(cds)

    return(cds)

  }, error=function(e){})
}


# return the cdna of a transcript given the dna, db, and transcript name
findCdna <- function(dna, db, transcript) {

  # return null for transcripts that do not return a valid cdna
  tryCatch({

    # set up filter for specified transcript
    txflt <- AnnotationFilter(~ tx_id == transcript)

    # extract the whole transcript sequence for the specified transcript above
    cdna <- extractTranscriptSeqs(dna, db, filter=txflt)
    cdna <- toString(cdna)

    return(cdna)

  }, error=function(e){})
}


# identify the location of the start codon within the cdna
findStartCodon <- function(cds, cdna){

  start_codon_location <- as.data.frame(str_locate(cdna, cds))[1, 1]

  return(start_codon_location)
}


# identify the location of the last codon of the coding sequence of a cdna sequence
findLastNucleotide <- function(cds, cdna){

  last_nucleotide_location <- as.data.frame(str_locate(cdna, cds))[1, 2]

  return(last_nucleotide_location)
}


# identify whether the aso/transcript combo spans multiple exons
exonExonAsoCheck <- function(dna, db, transcript, aso){

  # find cdna of transcript
  txflt <- AnnotationFilter(~ tx_id == transcript)
  cdna <- extractTranscriptSeqs(dna, db, filter=txflt)

  # get reverse complement of aso
  rc <- rc(aso)

  # identify location of aso within transcript
  start <- as.data.frame(str_locate(cdna, rc))[1, 1]
  end <- as.data.frame(str_locate(cdna, rc))[1, 2]

  # define the irange
  ir <- IRanges(start = start, end = end, names = transcript)

  # identify the GRange
  gr <- transcriptToGenome(ir, db)

  # determine number of exons the aso spans
  numberOfExons <- length(unlist(gr))

  # if aso spans multiple exons, return a 1, else return a 0
  if (numberOfExons > 1) {
    return(1)
  } else {
    return(0)
  }
}