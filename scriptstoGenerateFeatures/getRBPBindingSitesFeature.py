#!/bin/python

import sys
import pandas as pd


asoCoordinateFile = sys.argv[1]
intersectingFile = sys.argv[2]

# Get file with intersecting Sites
intersectingSites = pd.read_csv(intersectingFile,sep="\t",header=None)

# Group the intersecting Sites file by the ASO unique identifier
# Get a pandas series that has ASO unique identifer name and the number of RBP binding sites overlapping  the ASO 
intersectingSitesGroupedByASO = intersectingSites.groupby([3]).count()[0]

# Get all the ASO sites
with open(asoCoordinateFile) as f:
    aso_sites = [line.strip().split('\t')[3] for line in f]

# Go through all the ASO unique identifiers, if they are found in the intersecting Sites file, then append 1 otherwise append 0
# Just asking if there is any RBP binding site within the ASO binding region
with open("../tmp2/RBPbindingSites_Features_DoesItContainRBPSite.txt","w") as fw:
    for i in aso_sites:
        if i in list(intersectingSitesGroupedByASO.index.values):
            fw.write(i+"\t"+str(1)+"\n")
        else:
            fw.write(i+"\t"+str(0)+"\n")

# Go through all ASO unique identifiers, and put in the number of RBP binding sites, if they are not found, then append zero
with open("../tmp2/RBPbindingSites_Features_NumberOfRBPsites.txt","w") as fw:
    for i in aso_sites:
        if i in list(intersectingSitesGroupedByASO.index.values):
            numOfSites = intersectingSitesGroupedByASO[intersectingSitesGroupedByASO.index==i].values[0]
            fw.write(i+"\t"+str(numOfSites)+"\n")
        else:
            fw.write(i+"\t"+str(0)+"\n")