# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 21:36:55 2021

@author: Chung
"""
import numpy as np
import pandas as pd
from sklearn import svm
import sys
import argparse

from sklearn.linear_model import Perceptron
from sklearn.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
#from sklearn.ensemble import VotingRegressor
from sklearn.inspection import permutation_importance

from itertools import combinations_with_replacement
from itertools import permutations
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPRegressor
#import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder

def parse_args():
    parser = argparse.ArgumentParser(description = "Run regression on tab file")
    required = parser.add_argument_group("Required arguments to run")
    required.add_argument("input", type=argparse.FileType('r'))
    args = parser.parse_args()
    return args

def get_input_file():
    try:
        ASO_file = pd.read_csv(args.input.name, delimiter="\t")
        return(ASO_file)
    except Exception as err:
        print("Error opening file.")
        print(err)
        sys.exit(1)

def split(word):
    return list(word)

def one_hot(nuc_list):
    
    #equalize the ASO length
    patch_len = 21-len(nuc_list)
    nuc_list = nuc_list + ['X']*patch_len
    
    one_hot = []
    for nucelotide in nuc_list:
        if nucelotide == 'A':
            one_hot.append([1,0,0,0])
        elif nucelotide == 'C':
            one_hot.append([0,1,0,0])
        elif nucelotide == 'G':
            one_hot.append([0,0,1,0])
        elif nucelotide == 'T':
            one_hot.append([0,0,0,1])
        else:
            one_hot.append([0,0,0,0])
    return one_hot

def one_hot_gene(gene, gene_pool):
    
    one_hot = [0]*len(gene_pool)
    for i in range(len(gene_pool)):
        if gene == gene_pool[i]:
            one_hot[i] = 1
    return one_hot

def get_kmers(kmer_length): 

	all_kmers = []
	
	kmers = list(combinations_with_replacement("ATCG", kmer_length))
	for index in range(len(kmers)):
		kmers[index] = "".join(kmers[index])
		
	for kmer in kmers:
		permut = list(permutations(kmer))
		permut = list(set(["".join(x) for x in permut]))
		all_kmers = all_kmers + permut

	return all_kmers

def get_features(data, kmer_list):

	features = []

	seq_len = len(data[0])
	kmer_length = len(kmer_list[0])

	for entry in data:

		feature_count = [0] * len(kmer_list)
		start = 0
		end = kmer_length

		while end <= seq_len + 1:
			seq_slice = entry[start:end]
			start += 1
			end += 1

			if len(seq_slice) == kmer_length:
				#append to list of counts
				index = kmer_list.index(seq_slice)
				feature_count[index] += 1 

		#get the gc content and add this to the list of features for this sequence  
		#gc_content = calc_gc(entry)
        
        #feature_count.append(gc_content)
        #append to list of all features 
		features.append(feature_count)
		
	return features

def combine_features(*args):
    
    combined_features = []
    
    num_features = len(args[0])
    
    for index in range(num_features):
        
        new_list = []
        
        for a in args:

            if isinstance(a[index], list):
                new_list = new_list + a[index]
            else:
                new_list = new_list + list(a)
        
        combined_features.append(new_list)
        
    return combined_features
 
#sys.argv = ['RF.py', '../example_data/Processed_IDT_ASO_Data']
print('Process: loading data')
sys.argv = ['RF_v2.py', '../Data/openASO_gr_All+RBP+UTR+MFE+BpProbs+PhyloP.tsv']
args = parse_args()
ASO_score_data = get_input_file()
#ASO_score_data = ASO_score_data[:1000] ########Tempararily reduce the size of the input. 

'''
######Clean redundant data
ASO_score_data_clean = pd.DataFrame(ASO_score_data.iloc[0:2,:])
overlap = 0
for i in range(len(ASO_score_data)):
    print('processed', i, '/', len(ASO_score_data))
    overlap = 0
    seq = ASO_score_data.loc[i, 'ASOseq']
    for j in range(len(ASO_score_data_clean)):
        if ASO_score_data_clean.loc[j, 'ASOseq'].find(seq) != -1:
            overlap = 1
            break
    if overlap == 0:
        ASO_score_data_clean = ASO_score_data_clean.append(ASO_score_data.iloc[i,:], ignore_index=True)
ASO_score_data_clean.to_csv('../Data/openASO_gr_All+RBP+UTR_clean.tsv', sep='\t')
'''

#####build location info features
print('Process: construct features')
location_info = ASO_score_data["ASOseq"].apply(lambda x: one_hot(split(x)))
location_info = location_info.apply(lambda l: [item for sublist in l for item in sublist])

'''
#####build gene info features
structure_features = []
for i in range(len(ASO_score_data.columns)):
    if ASO_score_data.columns[i].find('RNAstructScore') !=-1:
        structure_features.append(ASO_score_data.columns[i])
'''

#build gene info features
gene_pool = list(set(ASO_score_data["GeneID"]))[:]
gene_info = ASO_score_data["GeneID"].apply(lambda x: one_hot_gene(x, gene_pool))

#build chromosome info features
chrom_pool = list(set(ASO_score_data["chrom"]))[:]
chrom_info = ASO_score_data["chrom"].apply(lambda x: one_hot_gene(x, chrom_pool))

#build transcript info features
transcript_pool = list(set(ASO_score_data["transcript_id"]))[:]
transcript_info = ASO_score_data["transcript_id"].apply(lambda x: one_hot_gene(x, transcript_pool))

#####build kmers(2-5 mer) features
symbol_2mers = get_kmers(2)
symbol_3mers = get_kmers(3)
symbol_4mers = get_kmers(4)
symbol_5mers = get_kmers(5)
features_2mer = get_features(ASO_score_data["ASOseq"], symbol_2mers)
features_3mer = get_features(ASO_score_data["ASOseq"], symbol_3mers)
features_4mer = get_features(ASO_score_data["ASOseq"], symbol_4mers)
features_5mer = get_features(ASO_score_data["ASOseq"], symbol_5mers)

#####build intron feature
ASO_score_data['is_across_intron'] = 0
ASO_length = [len(i) for i in ASO_score_data['ASOseq']]
ASO_width = ASO_score_data['width']
for i in range(len(ASO_length)):
    if ASO_length[i] != ASO_width[i]:
        ASO_score_data.loc[i, 'is_across_intron'] = 1

#####build strand feature
ASO_score_data['strand_boolean'] = 0
for i in range(len(ASO_score_data)):
    if ASO_score_data.loc[i, 'strand'] == '+':
        ASO_score_data.loc[i, 'strand_boolean'] = 1


#####construct X
X = location_info.to_list()
X = combine_features(X, gene_info.to_list())
X = combine_features(X, chrom_info.to_list())
X = combine_features(X, transcript_info.to_list())
X = combine_features(X, features_2mer)
X = combine_features(X, features_3mer)
#X = combine_features(X, features_4mer)
#X = combine_features(X, features_5mer)
X = combine_features(X, [ [i] for i in ASO_score_data["is_across_intron"].to_list()])
X = combine_features(X, [ [i] for i in ASO_score_data["strand_boolean"].to_list()])

added_features = ASO_score_data.columns[23:]
for i in range(len(added_features)):
    if added_features[i] == 'gene' or added_features[i] == 'seq':
        continue
    X = combine_features(X, [ [i] for i in ASO_score_data[added_features[i]].to_list()])

#####construct Y
Y = ASO_score_data["ASOeffective"]

#####constrcut feature names
max_length = max([len(i) for i in ASO_score_data['ASOseq']])
ATCG_identity = ['A','C','G','T']
features = ["tmp"] * 4 * max_length
for i in range(len(features)):
    base = i%4
    features[i] = 'if position ' + str(i//4) +' is ' + ATCG_identity[base]


features = features + gene_pool + ['chrom ' + str(i) for i in chrom_pool] + transcript_pool\
        + symbol_2mers +symbol_3mers\
        + ['is_across_intron', 'strand_boolean']\
        + list(added_features)
 #  features + gene_pool + ['chrom ' + str(i) for i in chrom_pool] + transcript_pool\
 #   + symbol_2mers + symbol_3mers + symbol_4mers + symbol_5mers\
 #   + ['DoesItContainRBPSite', 'NumberOfRBPsites', 'is_across_intron', 'strand_boolean']
features = np.array(features)

#####split test/train
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=0)
X_train = np.array(X_train)
X_test = np.array(X_test)
Y_train = np.array(Y_train)
Y_test = np.array(Y_test)

#####Find if any test ASOseq is identical as tranining
N_test = len(X_test)
ASO_test = X_test[:, : max_length*4]
ASO_train = X_train[:, : max_length*4]
index_identical = np.zeros(N_test)
for i in range(N_test):
    #find if there are training sequence entry that is exactly the same as i-th test entry
    ASO_diff = np.sum(abs(ASO_train - ASO_test[i]), axis = 1) 
    identical_hit = np.sum(np.where(ASO_diff == 0))
#    print(identical_hit)
    if identical_hit:
        index_identical[i] = 1
X_test_new = X_test[index_identical == 0]
Y_test_new = Y_test[index_identical == 0]

#####training model
print('Process: training')
#tfbs_classifier = RandomForestRegressor(n_estimators=100, max_features = "sqrt")
tfbs_classifier = MLPRegressor(hidden_layer_sizes=(500,500), activation='tanh', solver='adam', alpha=0.0001, early_stopping=True)
tfbs_classifier.fit(X_train, Y_train)

######predict Y
print('Process: predicting')
Y_pred = tfbs_classifier.predict(X_test)
Y_pred_new = tfbs_classifier.predict(X_test_new)

#####evaluate the model
error = mean_squared_error(Y_test, Y_pred)
print("root_mean_squared_error: %.3f" % error**0.5)
error_new = mean_squared_error(Y_test_new, Y_pred_new)
print("root_mean_squared_error_new : %.3f" % error_new**0.5)


#####analyze the importance of each feature
importances = tfbs_classifier.feature_importances_
std = np.std([tree.feature_importances_ for tree in tfbs_classifier.estimators_],
             axis=0)
indices = np.argsort(importances)
indices = indices[-50:]
print("Feature ranking:")
for f in range(len(indices)):
    print("%d. feature %s (%f)" % (f + 1, features[indices[f]], importances[indices[f]]))


#####Plot feature importance: premutation-based
perm_importance = permutation_importance(tfbs_classifier, X_test_new, Y_test_new)
sorted_idx = perm_importance.importances_mean.argsort()
fig, ax = plt.subplots(figsize=(15,15))
plt.barh(features[sorted_idx[-50:]], perm_importance.importances_mean[sorted_idx[-50:]])
plt.xlabel("Permutation Importance")
fig.savefig('../figure/feature_importance_permutation.png', bbox_inches='tight', dpi=200)

##### Plot feature importance: Gini
fig, ax = plt.subplots(figsize=(15,15))
plt.title("Feature importances")
#ax.barh(range(len(indices)), importances[indices], color="r", xerr=std[indices], align="center")
ax.barh(range(len(indices)), importances[indices], color="r", align="center")
ax.set_yticks(range(len(indices)))
ax.set_yticklabels(features[indices])
fig.savefig('../figure/feature_importance_new.png', bbox_inches='tight', dpi=200)
#plt.show()


#visual inspection to the predicted data
fig, ax = plt.subplots()
#plt.plot(Y_pred, Y_test, '.')
plt.plot(Y_pred_new, Y_test_new, '.')
plt.plot([0.2, 1], [0.2, 1], '--')
plt.xlabel('Predicted ASO effect', fontsize=20)
plt.ylabel('Real ASO effect', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()
fig.savefig('../figure/pred_vs_real_ASOeffect.png', dpi=500, bbox_inches='tight')


#pre addition of MFE
'''
1, 0.233, Include every features
2, 0.230, Include every features, clean
3, 0.228, Include every features, clean, 200trees
4, 0.227, Include every features, clean, 1000trees
5, 0.226, 
6, 0.227, Include every features, clean, 1000trees, sqrt
7, 0.228, No UTR features, clean, 100trees, sqrt
8, 0.224, No 3mer features, clean, 100trees, sqrt
9, 0.261, No 3mer features, NN 500x500
10, 0.227, No 3mer features, NN 500x500, early stopping
11, 0.222, No 3mer features, NN 500x500, early stopping, tanh
'''
'''
12, 0.235,Include every features, clean, 100trees, sqrt
'''