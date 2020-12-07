#!/bin/python

import sys
import pandas as pd
import numpy as np

#print "Number of arguments: ", len(sys.argv)
#print "The arguments are: " , str(sys.argv)

# Get the base pairing file
inputfile=str(sys.argv[1])
length_rna = int(sys.argv[2])
writefile=str(sys.argv[3])

# Read in base pairing prob values
bp_values = pd.read_csv(inputfile,sep="\t",skiprows=1)

# Get the actual probability, value given is the -log10 
bp_values = bp_values.assign(prob=np.power(10,-bp_values["-log10(Probability)"]))

#Define matrix that will contain the base pairing probabilities
global overall_matrix
overall_matrix = np.zeros((length_rna,length_rna))
# Definie function which will allocate the probabilities from dataframe to matrix using index values from dataframe
def fm(z):
    #print(z)
    overall_matrix[int(z["i"])-1,int(z["j"])-1]=z["prob"]
    overall_matrix[int(z["j"])-1,int(z["i"])-1]=z["prob"]
# Apply function to dataframe of basepairing probabilities
y = bp_values.apply(fm,axis=1)

overall_bpProbs = np.sum(overall_matrix,axis=0)

with open(writefile,"w") as fw:
    for p in overall_bpProbs:
        fw.write(str(p)+"\n")
    