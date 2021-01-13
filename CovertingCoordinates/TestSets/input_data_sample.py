# import necessary modules
import pandas as pd


# define function that can be imported
def subsample(file='../Data/Complete_ASOtoTranscriptSeq.tsv'):
    """
    subsamples a tsv file by randomly selecting 100 rows

    Parameters:
        file = tsv file

    Returns:
        100-row randomly sampled pandas dataframe
    """

    # read in the tsv
    df = pd.read_csv(file, sep=' ')

    # randomly subsample 100 rows
    subsample_df = df.sample(n=100)

    return subsample_df
