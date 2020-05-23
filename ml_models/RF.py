# -*- coding: utf-8 -*-
"""
Created on Sat May 23 12:02:26 2020

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

from itertools import combinations_with_replacement
from itertools import permutations
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
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
    patch_len = 22-len(nuc_list)
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
                new_list = new_list + list(a[index])
        
        combined_features.append(new_list)
        
    return combined_features



#####args.input.name = 'example_data/Processed_IDT_ASO_Data'
sys.argv = ['RF.py', '../example_data/Processed_IDT_ASO_Data']
args = parse_args()
ASO_score_data = get_input_file()

#####build location info features
location_info = ASO_score_data["ASO"].apply(lambda x: one_hot(split(x)))
location_info = location_info.apply(lambda l: [item for sublist in l for item in sublist])

#build gene info features
gene_pool = list(set(ASO_score_data["gene"]))[1:]
gene_info = ASO_score_data["gene"].apply(lambda x: one_hot_gene(x, gene_pool))

#####build kmers(2-5 mer) features
symbol_2mers = get_kmers(2)
symbol_3mers = get_kmers(3)
symbol_4mers = get_kmers(4)
symbol_5mers = get_kmers(5)
features_2mer = get_features(ASO_score_data["ASO"], symbol_2mers)
features_3mer = get_features(ASO_score_data["ASO"], symbol_3mers)
features_4mer = get_features(ASO_score_data["ASO"], symbol_4mers)
features_5mer = get_features(ASO_score_data["ASO"], symbol_5mers)

#####construct X
X = location_info.to_list()
X = combine_features(X, gene_info.to_list())
X = combine_features(X, features_2mer)
X = combine_features(X, features_3mer)
#X = combine_features(X, features_4mer)
#X = combine_features(X, features_5mer)

#####construct Y
Y = ASO_score_data["score"]

#####feature names
ATCG_identity = ['A','C','G','T']
features = ["1"]*4*22
for i in range(len(features)):
    base = i%4
    features[i] = 'if position ' + str(i//4) +' is ' + ATCG_identity[base]

features = features + gene_pool + symbol_2mers + symbol_3mers
features = np.array(features)

#####split test/train
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=0)
X_train = np.array(X_train)
X_test = np.array(X_test)
Y_train = np.array(Y_train)
Y_test = np.array(Y_test)

#####training model
tfbs_classifier = RandomForestRegressor(n_estimators=1000)
tfbs_classifier.fit(X_train, Y_train)
Y_pred = tfbs_classifier.predict(X_test)

#####evaluate the model
error = mean_squared_error(Y_test, Y_pred)
print('mean_squared_error: ', error)


#####analyze the importance of each feature
importances = tfbs_classifier.feature_importances_
std = np.std([tree.feature_importances_ for tree in tfbs_classifier.estimators_],
             axis=0)
indices = np.argsort(importances)
indices = indices[-50:]
print("Feature ranking:")
for f in range(len(indices)):
    print("%d. feature %s (%f)" % (f + 1, features[indices[f]], importances[indices[f]]))

##### Plot feature importance
fig, ax = plt.subplots()
plt.title("Feature importances")
#ax.barh(range(len(indices)), importances[indices], color="r", xerr=std[indices], align="center")
ax.barh(range(len(indices)), importances[indices], color="r", align="center")
ax.set_yticks(range(len(indices)))
ax.set_yticklabels(features[indices])
fig.savefig('figure/feature_importance.png', bbox_inches='tight', dpi=200)
#plt.show()

'''
#visual inspection to the predicted data
plt.plot(Y_test, Y_pred, '.')
plt.show()
'''