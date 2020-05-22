import numpy as np

aso_datapath = "example_data\Processed_IDT_ASO_Data"

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

from sklearn import svm
regr = svm.SVR()
regr.fit(X, Y)

abs_error_sum = 0
for sample_id in range(entry_count):
    x = X[sample_id, :].reshape(1, -1)
    res = regr.predict(x)[0]

    abs_error_sum += abs(res - Y[sample_id])

print(f"SVR avg error is {abs_error_sum/entry_count}")
