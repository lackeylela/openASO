import pandas as pd
import numpy as np
import sys

# Intersection data between bed file with coords and common variants
intersectFileInput=str(sys.argv[1])
# Yes or no string to only use single variants for common variants
onlyIncludeSNPs = str(sys.argv[2])

# Read in the intersecting sites
intersectingSites = pd.read_csv(intersectFileInput,sep="\t",header=None)

# Group the intersecting Sites by the ASO id
intersectingSitesGroupedByASO = intersectingSites.groupby([3])

# Collect all the values for common variant feature
feat_commonVar = []

# Go through each group by value
for aso, frame in intersectingSitesGroupedByASO:
    #print(aso)
    if onlyIncludeSNPs=="Yes":
        frame=frame[(frame[11].isin(['A','G','T','C','N']))&(frame[12].isin(['A','G','T','C','N']))]
    # Get unique SNPs for ASO
    frame=frame.drop_duplicates(subset=[9])
    # Get the length of the aso
    aso_length = frame.iloc[0,2]-frame.iloc[0,1]
    # Get the MAFs for the SNPs from column 14
    x=frame[14].str.split(";")
    # Go through every SNP's column 14 which contains the MAF information
    all_MAFs_ASO= []
    for i in x.values:
        # Get the allele frequency values for that SNP
        k=[j for j in i if j is not None and "CAF=" in j ][0]
        # Split up the allele frequencies in readable format 
        z= k.split("=")[1].split(",")
        # Make the frequencies floats since they are strings
        mafs = [float(n) for n in z if n!="."]
        all_MAFs_ASO.append(np.min(mafs))
    feat_commonVar.append([aso,str(np.median(all_MAFs_ASO)),str((sum([1 for i in all_MAFs_ASO if i > 0.01])/aso_length))])
    

with open("CommonVariation_features_ASOs_bindingSites.txt","w") as fw:
    for i in feat_commonVar:
        fw.write("\t".join(i))
        fw.write("\n")    


