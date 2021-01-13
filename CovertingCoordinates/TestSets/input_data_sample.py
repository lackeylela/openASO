# import necessary modules
import pandas as pd


# # read in the tsv
# df = pd.read_csv('../Data/Complete_ASOtoTranscriptSeq.tsv', sep=' ')
#
# # random 100 row sample
# subsample_df = df.sample(n=100)
#

# define function that can be imported
def subsample(file='../Data/Complete_ASOtoTranscriptSeq.tsv'):
    """
    subsamples a tsv file by randomly selecting 100 rows

    Parameters:
        file = tsv file

    Returns:
        100-row randomly sampled pandas dataframe
    """

    df = pd.read_csv(file, sep=' ')
    subsample_df = df.sample(n=100)

    return subsample_df
