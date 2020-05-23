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
from sklearn.ensemble import VotingRegressor
from ml_models.conv_net import OneDConvNet


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


def main():

    ASO_score_data = get_input_file()

    # Split ASO string into a list.
    ASO_score_data["ASO"] = ASO_score_data["ASO"].apply(lambda x: one_hot(split(x)))
    ASO_score_data["ASO"] = ASO_score_data["ASO"].apply(lambda l: [item for sublist in l for item in sublist])

    X = pd.DataFrame(ASO_score_data["ASO"].to_list())
    X = X.replace(np.nan, 0)
    Y = ASO_score_data["score"]

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=0)
    X_train = X_train.to_numpy()
    X_test = X_test.to_numpy()
    Y_train = Y_train.to_numpy()
    Y_test = Y_test.to_numpy()

    # Conv net sample code
    one_d_convnet = OneDConvNet()
    trained_model = one_d_convnet.fit(X_train, X_test, Y_train, Y_test, epochs=1, batch_size=32)

    clf = svm.SVR()
    clf.fit(X_train, Y_train)
    Y_predicted = clf.predict(X_test)

    mse = np.mean(np.square(Y_predicted - Y_test))

    for i in range(len(Y_predicted)):
        print("Predicted Value: " + str(Y_predicted[i]))
        print("Actual Value: " + str(Y_test[i]))
        print("-------------------------------------------------")

    print("Mean Squared Error: " + str(mse))

if __name__ == '__main__':
    args = parse_args()
    if args.input:
        main()
    else:
        print("Enter input file")
