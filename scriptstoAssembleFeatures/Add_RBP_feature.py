# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:19:46 2021

@author: Shawn
"""

RBP_data = pd.read_csv('../Data/Features/RBPbindingSites_Features_NumberOfRBPsites.txt', delimiter="\t", header=None)
data = pd.read_csv('../Data/openASO_gr_All.tsv', delimiter="\t")

#create two new columns for gene name and the RBP binding sequence
RBP_data['gene'] = [i.split(' _ ')[0] for i in RBP_data.iloc[:, 0]]
RBP_data['seq'] = [i.split(' _ ')[1] for i in RBP_data.iloc[:, 0]]

#data = ASO_score_data
#create two new columns that are self explainatory
data['DoesItContainRBPSite'] = 0
data['NumberOfRBPsites'] = 0

for i in range(len(data)):
    print('processed', i, '/', len(data))
    for j in range(len(RBP_data)):
        RBP_seq = RBP_data.loc[j, 'seq']
        RBP_number = RBP_data.loc[j, 1]
        if data.loc[i, 'ASOseq'].find(RBP_seq) != -1:
            data.loc[i, 'DoesItContainRBPSite'] = 1
            data.loc[i, 'NumberOfRBPsites'] = RBP_number
            break
        
data.to_csv('../Data/openASO_gr_All+RBP.tsv', sep='\t')