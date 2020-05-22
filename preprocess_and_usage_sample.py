import numpy as np
import pandas as pd
from sklearn import svm
import sys
import argparse

from sklearn.linear_model import Perceptron
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import VotingRegressor

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

def main():

    test = get_input_file()
    print(test["ASO"])
    base2idx = {"A":0, "T":1, "G":2, "C":3, "X":4} #X is padding for varying gene segment length

    entry_count = 0
    max_sequence_len = -1
    with open(aso_datapath) as datafile:
        for line in datafile:
            columns = line.split()

            effectiveness = columns[-1]
            sequence = columns[-2]

            curr_sequence_len = len(sequence)

            entry_count += 1
            max_sequence_len = curr_sequence_len if curr_sequence_len > max_sequence_len else max_sequence_len


    X = np.zeros(shape=(entry_count, max_sequence_len*5))
    Y = np.zeros(shape=(entry_count, ))

    with open(aso_datapath) as datafile:
        for i, line in enumerate(datafile):
            columns = line.split()

            effectiveness = columns[-1]
            sequence = columns[-2]

            curr_sequence_len = len(sequence)

            if curr_sequence_len < max_sequence_len:
                padding_count = max_sequence_len - curr_sequence_len
                sequence += "X" * padding_count # append X to end to achieve uniform length

            one_hot_sequence = np.zeros(shape=(5, max_sequence_len), dtype=np.int32)
            for j, base in enumerate(sequence):
                one_hot_sequence[base2idx[base], j] = 1

            input_features = one_hot_sequence.flatten()
            output = effectiveness

            X[i, :] = input_features
            Y[i] = effectiveness

    regr = svm.SVR()
    regr.fit(X, Y)

    abs_error_sum = 0
    for sample_id in range(entry_count):
        x = X[sample_id, :].reshape(1, -1)
        res = regr.predict(x)[0]

        abs_error_sum += abs(res - Y[sample_id])

    print(f"SVR avg error is {abs_error_sum/entry_count}")

if __name__ == '__main__':
    args = parse_args()
    if args.input:
        main()
    else:
        print("Enter input file with --i")
