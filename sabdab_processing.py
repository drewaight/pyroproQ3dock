import numpy as np
import pandas as pd

sabdab = '/home/drewaight/hdd1/sabdab/protein_structues/20190811_0954601_summary.tsv'

df1 = pd.read_csv(sabdab, delimiter='\t', header=0)
df2 = df1[pd.notnull(df1['Lchain'])]
print(df2.light_ctype.unique())
print(len(df2[(df2['light_ctype'] == 'unknown')]))
print(len(df2.pdb.unique()))
print(df2.describe())
print(df2)