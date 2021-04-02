# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 22:55:34 2021

@author: Chung
"""
# -*- coding: utf-8 -*-

import pandas as pd

data_new = pd.read_csv('../Data/Features/conservationSummaryData.tsv', delimiter="\t")
data = pd.read_csv('../Data/openASO_gr_All+RBP+UTR+MFE+BpProbs.tsv', delimiter="\t")

#create two new columns for gene name and the RBP binding sequence
data_new['gene'] = [i.split(' _ ')[0] for i in data_new.iloc[:, 0]]
data_new['seq'] = [i.split(' _ ')[1] for i in data_new.iloc[:, 0]]

#data = ASO_score_data
#create two new columns that are self explainatory
data['averageConservationScore'] = 0
data['countConserveAboveZero'] = 0

for i in range(len(data)):
    print('processed', i, '/', len(data))
    for j in range(len(data_new)):
        seq = data_new.loc[j, 'seq']
        if data.loc[i, 'ASOseq'].find(seq) != -1:
            data.loc[i, 'averageConservationScore'] = data_new.loc[j, 'averageConservationScore']
            data.loc[i, 'countConserveAboveZero'] = data_new.loc[j, 'countConserveAboveZero']
            break
        
data.to_csv('../Data/openASO_gr_All+RBP+UTR+MFE+BpProbs+PhyloP.tsv', sep='\t')