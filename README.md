**openASO**

**Design effective antisense oligonucleotides based on target RNA structural ensembles**

RNAs are dynamic molecules that exist as populations of structures. Can population-based structural models inform ASO design? Help develop a tool to advance ASO design for experimental and clinical applications!

Antisense oligonucleotides (ASOs) have a high potential therapeutic value as they can be used to target any transcript and regulate its expression. One factor in choosing an effective ASO is identifying an open, unstructured region of the target RNA. However, RNA is a flexible molecule best represented by an ensemble of structures.  Does using ensemble modeling to identify open regions within the target RNA promote selection of more suitable ASO sequences? To answer this question we will use machine learning to ask whether ensemble guided structure models perform better than traditional minimum free energy structure using publicly available datasets of ASOs and their targeting efficiency. Other known attributes that affect ASO targeting will be incorporated in both algorithms including RNA binding protein data (eCLIP) data.  Finally, we will build a visualization tool where users can input their gene of interest and map regions suitable for further ASO design and testing.

Example files are generated for PTEN, one of the genes in the ASO test dataset. We will be generating data for the approximately 100 additional genes to complete the full ASO project. We will need to integrate the disparate data sources (sequence, transcript and genomic coordinates), develop a structure metric, create a machine learning algorithm and evaluate its performance.

Files
1. sequence_hg38_ncbiRefSeq_NM_000314.8.fa - A fasta sequence file for PTEN from NCBI (NM_000314).

2. gene_architecture_NM_000314.8.gtf - A gtf file describing the architecture of the PTEN gene for transcript NM_000314.8.
 
3. structure_ensembles1-10_NM_000314.8.ct - A ct format including 1000 structures representative of the ensemble. Ensemble structures were generated with the RNAstructure commands "partition PTEN.fa PTEN.pfs" followed by "stochastic PTEN.pfs PTEN_ensemble.ct". Only the first ten structures are included as the pfs and full ct files are very large.

4. structure1_ensemble_NM_000314.8.dot - A single structure (one of 5) in dot-bracket format generated with ct2dot from the RNAstructure software.

5. structure_mfe_NM_000314.8.ct - A ct minimum free energy structure generated with the RNAstructure command "Fold PTEN.fa PTEN.ct -mfe". A similar minimum free energy structure with a max pairing distance of 200 nucleotides is also included as PTEN_200md_mfe.ct.

6. structure_mfe_NM_000314.8.dot - The minimum free energy structure in dot-bracket format generated with ct2dot from the RNAstructure software.

These structures will have to be converted into a metric, possibly through forgi element representation.

7. aso_efficacy_PTEN.tab - A tab delimited file of tested PTEN ASOs and their efficiencies (gene, sequence, efficiency).

Further references and resources:

[RNA structure modeling](https://rna.urmc.rochester.edu/Text/index.html)

[ENCODE – eCLIP data](https://www.encodeproject.org/eclip/)

[PolyASite – alternative polyadenylation](https://polyasite.unibas.ch/)

[Breakdown structures by base-pairing](https://viennarna.github.io/forgi/graph_tutorial.html)

[Gene sequence and architecture resource](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)

[REDIportal – A-to-I editing sites](http://srv00.recas.ba.infn.it/atlas/)


McQuisten KA, Peek AS. Identification of sequence motifs significantly associated with antisense activity. BMC Bioinformatics. 2007;8:184. doi:10.1186/1471-2105-8-184 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/17555590)

Johnson E, Srivastava R. Volatility in mRNA secondary structure as a design principle for antisense. Nucleic Acids Res. 2013; 41(3):e43. doi: 10.1093/nar/gks902, [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/23161691)

Further questions:
Are there additional features we should include? Some ideas are splice junctions, single nucleotide polymorphisms and distance from the start or stop codon.

