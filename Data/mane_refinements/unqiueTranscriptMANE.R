# identify singlue unique tx_id for each gene of interest
# we lose MAP4K4 and MAPK8 using the MANE info. may want to go manually check.

# read in the 2 csvs
mane_tx_id = read.csv('geneToTranscript.csv') # MANE table
data = read.csv('../Complete_ASOtoTranscriptSeq.tsv',
                sep='',stringsAsFactors = FALSE) # primary project dataset

# generate new column for eventual unique tx id
data$tx_id = NA

# identify unique genes only
unique_genes = unique(data$GeneID)

# populate tx_id column with unique transcript IDs by indexing MANE table
for (gene in unique_genes){
  tryCatch({
    data$tx_id[which(data$GeneID == gene)] = mane_tx_id$transcript_id[which(mane_tx_id$gene_name == gene)]
  }, error = function(e){})
}

# substring transcript name to get rid of .[0-9]
data$hg38.knownGene.name <- substr(data$hg38.knownGene.name, 0, 15)
data$tx_id <- substr(data$tx_id, 0, 15)

# create unique ids
data$unique_id <- paste(data$GeneID, " _ ", data$ASOseq)

# drop transcripts that don't match the mane transcripts for each gene
mane_refine = data[which(data$hg38.knownGene.name == data$tx_id),]

# drop all duplicates (ie PIK3CD was dropped completely)
mane_refine = mane_refine[which(table(mane_refine$unique_id) == 1),]

# write table to tsv
write.table(mane_refine, file='C:/Users/james/R/openASO/openASOcorrect/Data/uniqueIds.tsv',
            sep='\t',
            row.names=FALSE)
