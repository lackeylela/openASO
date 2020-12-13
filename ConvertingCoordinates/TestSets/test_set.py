# import necessary modules
from sklearn.model_selection import train_test_split
import pandas as pd


# read in the dataframe
df = pd.read_csv('../../deeper_utr_scope.tsv', sep='\t')

# drop the single value effectives to get an even split of effective
df = df.loc[(df.effective != 0.06) & (df.effective != 0.03)]

# split (stratify on effectiveness)
y = df.effective
train_df, test_df = train_test_split(df, test_size=0.5, random_state=111, stratify=y)

# write to csv
train_df.to_csv('train.csv', index=False)
test_df.to_csv('test.csv', index=False)