# takes a look at the no match aso-to-transcript problem for BIRC2
# evaluates potential that BIRC2 is not the correct gene for
#   'InhibitorofApoptosis-2' and explores option that 'BIRC3'
#   might actually be the correct gene

# import necessary modules
import pandas as pd
import os


# define working directory
os.chdir("C:/Users/james/PycharmProjects/openASOtroubleshoot/unique_checks/birc2_checks")

# read in the tables
aso_table = pd.read_csv('birc2_asos.tsv', sep='\t')
birc2 = pd.read_csv('birc2_transcripts.tsv', sep='\t')
birc3 = pd.read_csv('birc3_transcripts.tsv', sep='\t')

# isolate aso and transcript sequences
rc_asos = aso_table.rc.values
birc2_seqs = birc2.cdna.values
birc3_seqs = birc3.cdna.values

# initiate empty counters
birc2_counter = 0
birc3_counter = 0

# count number of aso matches in each set of sequences
for rc in rc_asos:
    for cdna in birc2_seqs:
        if rc in cdna:
            birc2_counter += 1

for rc in rc_asos:
    for cdna in birc3_seqs:
        if rc in cdna:
            birc3_counter += 1

print(f'BIRC2 match count: {birc2_counter}\n'
      f'BIRC3 match count: {birc3_counter}')
