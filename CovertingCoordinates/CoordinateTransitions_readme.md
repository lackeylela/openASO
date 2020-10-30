We want to identify the location of ASOs based on sequence within a transcript and covert these numbers to genomic coordinates. We also want to identify elements from external datasources in genomic coordinates and convert these to transcript coordinates.

Our initial script (spear-headed by Kim) is 20200523_openASO_position.Rmd
It takes "Complete_ASOtoTranscriptSeq.tsv" as input, which is a combination of Processed_IDT_ASO_Data, hgnc-symbol-check-manual-additions-3.tsv and result.txt. The R commands used to generate this file are below. 

```
library(seqinr)
library(dplyr)

Processed_IDT_ASO_Data <- read.delim("Processed_IDT_ASO_Data", header=FALSE)
hgnc.symbol.check.manual.additions.3 <- read.delim("hgnc-symbol-check-manual-additions-3.tsv") %>% select(Input, HGNC.ID)
results <- read.delim("results.txt") %>% select(Approved.symbol, HGNC.ID)
Key <- merge(hgnc.symbol.check.manual.additions.3, results)
names(Processed_IDT_ASO_Data) <- c("GeneID", "ASOseq", "ASOeffective")
Combo <- merge(Processed_IDT_ASO_Data, Key, by.x = "GeneID", by.y = "Input")

output <- read.delim("output.csv", header=FALSE, comment.char="#")

names(output) <- c("hg38.knownGene.name","hg38.knownGene.chrom","hg38.knownGene.exonCount","hg38.knownGene.exonStarts","hg38.knownGene.exonEnds","hg38.knownGene.alignID","hg38.kgXref.kgID","hg38.kgXref.mRNA","hg38.kgXref.geneSymbol","Sequence","End")

SequenceWithASO <- merge(Combo, output, by.x = "Approved.symbol", by.y = "hg38.kgXref.geneSymbol") %>% 
  select(GeneID = Approved.symbol, OriginalName = GeneID, ASOseq, ASOeffective, HGNC.ID, hg38.knownGene.name, Sequence)

#write.table(SequenceWithASO, file = "Complete_ASOtoTranscriptSeq.tsv", quote = FALSE, row.names = FALSE)
```
