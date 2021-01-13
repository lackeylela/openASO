# import necessary modules
import pandas as pd
from random import sample


def subsample(file='Data/Complete_ASOtoTranscriptSeq.tsv',
              sep=' ',
              n=100):
    """
    subsamples a tsv file by randomly selecting 100 rows

    Parameters:
        file = tsv file
        sep = separator of file type. default is ' ' of primary tsv
        n = number of rows to randomly sample

    Returns:
        100-row randomly sampled pandas dataframe
    """

    # read in the tsv
    df = pd.read_csv(file, sep=' ')

    # randomly subsample 100 rows
    subsample_df = df.sample(n=n)

    return subsample_df


def distributed_sample(file='Data/Complete_ASOtoTranscriptSeq.tsv',
                       sep=' ',
                       intervals=None):

    """
    return a subsample dataframe with even distribution of a range of intervals

    Parameters:
        file = file to subsample (default is primary ASO tsv)
        sep = type of separator (default is the ' ' sep that is present in the primary tsv)
        intervals (list)= desired range of intervals. Default will be intervals of .1 from 0 to 1
    """

    # define intervals if not defined in function call
    if intervals is None:
        intervals = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1]

    # initiate empty variables
    i = 0
    lengths = []

    # read tsv
    df = pd.read_csv(file, sep=sep)

    # determine lowest represented interval
    while i < len(intervals) - 1:
        low = intervals[i]
        high = intervals[i + 1]
        index_list = df.index[(df.ASOeffective > low) & (df.ASOeffective < high)].tolist()
        lengths.append(len(index_list))
        i += 1

    lowest_rep = min(lengths)

    # reset empty variable
    i = 0
    indeces = []

    # develop list of indices to subsample
    while i < len(intervals) - 1:
        low = intervals[i]
        high = intervals[i + 1]
        index_list = df.index[(df.ASOeffective >= low) & (df.ASOeffective <= high)].tolist()
        s = sample(index_list, lowest_rep)
        indeces.extend(s)
        i += 1

    # subsample the dataframe
    subsample_df = df.iloc[indeces]

    return subsample_df
