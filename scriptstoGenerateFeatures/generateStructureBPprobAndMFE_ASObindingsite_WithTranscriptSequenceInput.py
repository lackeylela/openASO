import pandas as pd
import subprocess
import os
import sys

# tab separated val file with ASO id, sequence around ASO and start pos of ASO binding site
seqFileInput=str(sys.argv[1])

# Simple function that writes a fasta file
def writefasta(seqname,seq,file_name):#function that writes a seq, seqname and file_name into a fasta file
	''' writes seqname and seq into a file of name file_name '''
	file = open(file_name,"w")
	file.write('>'+seqname+"\n")
	file.write(seq+'\n')
	file.close()
	return

# open up the sequence file that is tab separate
with open(seqFileInput) as f:
    data = [line.strip().split("\t") for line in f]
print(len(data))
# this list will store the data for all the mfe features
mfe_feats_overall=[]

# this list will store the data for the bp probs features
bpProbs_feats_overall=[]
    
# Go through each line, with each sequence and create a fasta file for that seq, run the partition function on fasta file, and then calculate base pairing probabilities 
for data_i in data:
    print(data_i)
    # Get name of sequence, the sequence itself and the length of sequence
    name=data_i[0]
    seq=data_i[1].upper()
    startASO=int(data_i[2])
    seqlen=len(seq)
    
    # Generate a fasta file for that sequence
    writefasta(name,seq,"ASOBindingSiteToFold.fa")
    
    # First lets get the base pairing probabilities
    subprocess.check_output(['partition-smp','ASOBindingSiteToFold.fa','ASOBindingSiteToFold.pfs'])
    subprocess.check_output(['ProbabilityPlot','ASOBindingSiteToFold.pfs','ASOBindingSiteToFold_ProbabilityPlot.txt','--text'])
    subprocess.check_output(['python','calculateBasePairingProbability.py','ASOBindingSiteToFold_ProbabilityPlot.txt',str(seqlen),'ASOBindingSiteToFold_BasePairingProbability.txt'])
    
    with open("ASOBindingSiteToFold_BasePairingProbability.txt") as f:
        bp_probs = [float(line.strip()) for line in f]
        
    bp_probs_feat = [str(i) for i in bp_probs[startASO:startASO+18]]
    
    bpProbs_feats_overall.append([name]+bp_probs_feat)
    
    # Let's get the MFE structure
    subprocess.check_output(['Fold','ASOBindingSiteToFold.fa','ASOBindingSiteToFold.ct','--MFE'])
    subprocess.check_output(['ct2dot','ASOBindingSiteToFold.ct','ALL','ASOBindingSiteToFold.db','--format','simple'])
    
    # Read in the folded data
    with open("../tmp/ASOBindingSiteToFold.db") as f:
        folded_data = [line.strip() for line in f]
        
    # Get the dotbracket structure
    db_data = folded_data[2]
    # Get the dotbracket for just the binding site
    db_bindingSite=db_data[startASO:startASO+18]
    db_bindingSite=db_bindingSite.replace("(","1")
    db_bindingSite=db_bindingSite.replace(")","1")
    db_bindingSite=db_bindingSite.replace(".","0")
    
    mfe_feat = '\t'.join(db_bindingSite)
    
    mfe_feats_overall.append(name+"\t"+mfe_feat)
    
# Write the mfe feats overall 
with open("MFE_features_ASOs_bindingSites.txt","w") as fw:
    for i in mfe_feats_overall:
        fw.write(i+"\n")

# Write the mfe feats overall 
with open("BpProbs_features_ASOs_bindingSites.txt","w") as fw:
    for i in bpProbs_feats_overall:
        fw.write("\t".join(i))
        fw.write("\n")
    


