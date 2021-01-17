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
dna <- ensembldb::getGenomeTwoBitFile(db)


# must include if statement to check whether the transcript is non-protein coding
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


# must include if statement to check whether the transcript is non-protein coding
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

  cds <- toString(cds)
  cdna <- toString(cdna)

  start_codon_location <- as.data.frame(str_locate(cdna, cds))[1, 1]

  return(start_codon_location)
}



# identify the location of the last codon of the coding sequence of a cdna sequence
findLastNucleotide <- function(cds, cdna){

  cds <- toString(cds)
  cdna <- toString(cdna)

  last_nucleotide_location <- as.data.frame(str_locate(cdna, cds))[1, 2]

  return(last_nucleotide_location)
}



# identify how much of the aso is binding to the 5' utr, if it binds at all
# returns two values (binds, length_of_bind)
fivePrimeUtrBind <- function(start_codon_location, aso_start_location, aso) {

  # determine the length of the aso
  aso_length <- nchar(aso)

  # check if the aso binds upstream of the 5' utr
  five_prime_utr_nucleotides <- start_codon_location - aso_start_location

  # return 1 for bind or 0 for no bind, and return number of nucletoide crossover of utr bind
  if (five_prime_utr_nucleotides > 0){

    # return 1 for bind
    binds <- as.integer(1)

    # define length of crossover
    if (five_prime_utr_nucleotides > aso_length){

      # whole aso binds to the 5' utr so return length of aso
      length_of_bind <- as.integer(aso_length)

    } else {

      # only portion of aso binds to 5' utr so return the bind length
      length_of_bind <- as.integer(five_prime_utr_nucleotides)
    }
  } else {

    # return 0 for no bind
    binds <- as.integer(0)
    length_of_bind <- as.integer(0)
  }

  my_list <- list("binds" = binds, "length_of_bind" = length_of_bind)
  return(my_list)
}



# identify how much of the aso is binding to the 3' utr, if it binds at all
# returns two values (binds, length_of_bind)
threePrimeUtrBind <- function(last_nucleotide_location, aso_end_location, aso) {

  # determine the length of the aso
  aso_length <- nchar(aso)

  # check if the aso binds upstream of the 5' utr
  three_prime_utr_nucleotides <- aso_end_location - last_nucleotide_location

  # return 1 for bind or 0 for no bind, and return number of nucletoide crossover of utr bind
  if (three_prime_utr_nucleotides > 0){

    # return 1 for bind
    binds <- as.integer(1)

    # define length of crossover
    if (three_prime_utr_nucleotides > aso_length){

      # whole aso binds to the 5' utr so return length of aso
      length_of_bind <- as.integer(aso_length)

    } else {

      # only portion of aso binds to 5' utr so return the bind length
      length_of_bind <- as.integer(three_prime_utr_nucleotides)
    }
  } else {

    # return 0 for no bind
    binds <- as.integer(0)
    length_of_bind <- as.integer(0)
  }

  my_list <- list("binds" = binds, "length_of_bind" = length_of_bind)
  return(my_list)
}



# identify whether the aso/transcript combo spans multiple exons
exonExonAsoCheck <- function(dna, db, transcript, aso, cdna, reverse_complement, start_location, end_location){

  # find cdna of transcript
  txflt <- AnnotationFilter(~ tx_id == transcript)

  if (missing(cdna)){
    cdna <- extractTranscriptSeqs(dna, db, filter=txflt)
  }

  # get reverse complement of aso
  if (missing(reverse_complement)){
    reverse_complement <- rc(aso)
  }

  # identify location of aso within transcript
  if (missing(start_location)){
    start_location <- as.data.frame(str_locate(cdna, reverse_complement))[1, 1]
  }

  if (missing(end_location)){
    end_location <- as.data.frame(str_locate(cdna, reverse_complement))[1, 2]
  }

  # define the irange
  ir <- IRanges(start = start_location, end = end_location, names = transcript)

  # identify the GRange
  gr <- transcriptToGenome(ir, db)

  # determine number of exons the aso spans
  numberOfExons <- length(unlist(gr))

  # if aso spans multiple exons, return a 1, else return a 0
  if (numberOfExons > 1) {
    return(as.integer(1))
  } else {
    return(as.integer(0))
  }
}