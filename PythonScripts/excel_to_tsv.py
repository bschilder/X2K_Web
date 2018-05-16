import pandas as pd
import os

def convert(in_file, out_file):
    df = pd.read_excel(in_file).applymap(lambda s: s.replace('\n', ' ') if type(s) == str else s)
    cols = list(df)
    cols.insert(0, cols.pop(cols.index('Database')))
    cols.insert(1, cols.pop(cols.index('Filename(s)')))
    df = df.loc[:, cols]
    df.to_csv(out_file, sep='\t', index=False)

files=!ls *.xlsx
for file in files:
    base, ext = os.path.splitext(file)
    convert(file, '../X2K/src/main/webapp/resources/' + base.split('_')[0] + '.tsv')
