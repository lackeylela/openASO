# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 21:53:13 2021

@author: Chung
"""
# -*- coding: utf-8 -*-

import pandas as pd

with open('../Data/Features/BpProbs_features_ASOs_bindingSites.txt') as f:
    txt = f.read()
txt = txt.split('\n')[:-1]
data_new = pd.DataFrame()
for i in range(len(txt)):
    print(i, '/', len(txt))
    toks = txt[i].split('\t')
    if len(toks) < 22:
        toks = toks + [-1] * (22-len(toks)) #assign -1 if the ASO is shorter than 21 nt
    data_new = data_new.append([toks], ignore_index = True)
#define column names 
column_names = ['gene+seq']
for i in range(21):
    column_names += ['BpProbs' + str(i+1)]
data_new.columns = column_names

#data_new = pd.read_csv('../Data/Features/BpProbs_features_ASOs_bindingSites.txt', delimiter="\t", header=None)
data = pd.read_csv('../Data/openASO_gr_All+RBP+UTR+MFE.tsv', delimiter="\t")

#create two new columns for gene name and the RBP binding sequence
data_new['gene'] = [i.split(' _ ')[0] for i in data_new.iloc[:, 0]]
data_new['seq'] = [i.split(' _ ')[1] for i in data_new.iloc[:, 0]]

#create new columns 
for i in range(1, len(data_new.columns)):
    data[data_new.columns[i]] = 0

for i in range(len(data)):
    print('processed', i, '/', len(data))
    for j in range(len(data_new)):
        seq = data_new.loc[j, 'seq']
        if data.loc[i, 'ASOseq'].find(seq) != -1:
            for k in range(1, 22):
                    data.loc[i, data_new.columns[k]] = data_new.loc[j, data_new.columns[k]]
            break
        
data.to_csv('../Data/openASO_gr_All+RBP+UTR+MFE+BpProbs.tsv', sep='\t')