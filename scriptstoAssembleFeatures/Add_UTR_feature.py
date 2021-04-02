# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 19:59:41 2021

@author: Chung
"""
import pandas as pd

data_new = pd.read_csv('../Data/utr_features/utr_features.tsv', delimiter="\t", header=None)
data = pd.read_csv('../Data/openASO_gr_All+RBP.tsv', delimiter="\t")

#create two new columns for gene name and the RBP binding sequence
data_new['gene'] = [i.split('_')[0] for i in data_new.iloc[:, 0]]
data_new['seq'] = [i.split('_')[1] for i in data_new.iloc[:, 0]]

#data = ASO_score_data
#create two new columns that are self explainatory
data['threePrimeUtrBind'] = 0
data['lengthThreePrimeBind'] = 0
data['fivePrimeUtrBind'] = 0
data['lengthFivePrimeBind'] = 0
for i in range(len(data)):
    print('processed', i, '/', len(data))
    for j in range(len(data_new)):
        seq = data_new.loc[j, 'seq']
        RBP_number = data_new.loc[j, 1]
        if data.loc[i, 'ASOseq'].find(seq) != -1:
            data.loc[i, 'threePrimeUtrBind'] = data_new.loc[j, 'threePrimeUtrBind']
            data.loc[i, 'lengthThreePrimeBind'] = data_new.loc[j, 'lengthThreePrimeBind']
            data.loc[i, 'fivePrimeUtrBind'] = data_new.loc[j, 'fivePrimeUtrBind']
            data.loc[i, 'lengthFivePrimeBind'] = data_new.loc[j, 'lengthFivePrimeBind']
            break
        
data.to_csv('../Data/openASO_gr_All+RBP+UTR.tsv', sep='\t')