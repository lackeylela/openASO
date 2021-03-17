# identify singlue unique tx_id for each gene of interest
# we lose MAP4K4 and MAPK8 using the MANE info. may want to go manually check.

# read in the 2 csvs
mane_tx_id = read.csv('geneToTranscript.csv') # MANE table
data = read.csv('Data/Complete_ASOtoTranscriptSeq.tsv',
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

# write table to tsv
write.table(data, file='completeASO_w_maneRefine.tsv', sep='\t', col.names=NA)
