import gffpandas.gffpandas as gffpd


file = 'utr_features/MANE.GRCh38.v0.93.select_ensembl_genomic.gff'

annotation = gffpd.read_gff3(file)
attr_to_columns = annotation.attributes_to_columns()

gene_to_transcript_df = attr_to_columns[['gene_name', 'transcript_id']]
gene_to_transcript_df.drop_duplicates(inplace=True)
gene_to_transcript_df.dropna(inplace=True)
gene_to_transcript_df.to_csv('geneToTranscript.csv', index=False)