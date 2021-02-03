# openASO

**Design effective antisense oligonucleotides**

Can machine learning be used to identify aspects of RNA targets that make antisense oligonucleotide (ASOs) more effective? This project seeks to develop datasets and tools to understand important features of effective ASO targets to improve ASO design for experimental and clinical applications!

Antisense oligonucleotides (ASOs) have a high potential therapeutic value as they can be used to target any transcript and regulate its expression. One factor in choosing an effective ASO is identifying an open, unstructured region of the target RNA. However, RNA is a flexible molecule best represented by an ensemble of structures.  Does using ensemble modeling to identify open regions within the target RNA promote selection of more suitable ASO sequences? To answer this question we will use machine learning to ask whether ensemble guided structure models perform better than traditional minimum free energy structure using publicly available datasets of ASOs and their targeting efficiency. Other known attributes that affect ASO targeting will be incorporated in both algorithms including RNA binding protein data (eCLIP) data.  Finally, we will build a visualization tool where users can input their gene of interest and map regions suitable for further ASO design and testing.

## Team members

Lela Lackey was the team lead and coordinator. Team members worked in three major groups including:

**Data processing**

Axel Hauduc, Ivan Jimenez, Jayashree Kumar, James Adler and Kimberly Wellman are members of the data processing group of team openASO. This group took ASO sequences and incomplete gene target information, and built a data frame from which new features could be analyzed. This entailed data clean up and standardization of the gene names to HGNC nomenclature. Transcript sequences of ASO targets were also included so that ASOs could be mapped to transcript coordinates as well as genome coordinates. Some of the features included in this data frame were secondary structure, RNA-binding protein occupancy, and frequency of SNP variation.

**Machine Learning**

Shawn Hsueh, Alex Sweeten and Tony Shen are members of the machine learning team of openASO. THey worked on the machine learning algorithm. They used sci-kit learn to develop a SVM Regression algorithm. This will be used on the final dataset including target RNA secondary structure, RNA-binding protein occupancy, and SNP variation.

**Visualization**

Sophia Shen, Owen Tsai, Helena Catena Sánchez and Kelly Wei are members of the visualization team of openASO. They made an interactive R Shiny table, where you can load the ASO data in tsv format and visualize it. There is an interactive table where you can select and filter by table attributes. We plan to expand this and show a representation of each transcript and where the different ASO bind and their features.


# Data sources:

[ENCODE – eCLIP data](https://www.encodeproject.org/eclip/)

[Gene sequence and architecture resource](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)

[REDIportal – A-to-I editing sites](http://srv00.recas.ba.infn.it/atlas/)

[Genomic variation](https://www.ncbi.nlm.nih.gov/snp/)

McQuisten KA, Peek AS. Identification of sequence motifs significantly associated with antisense activity. BMC Bioinformatics. 2007;8:184. doi:10.1186/1471-2105-8-184 [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/17555590)
