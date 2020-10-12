import os
import re
import sys
import numpy as np
import pandas as pd
import pathlib
from biopandas.pdb import PandasPdb
import glob
import json

complexes = []
hdir = "/media/hdd1/roproQ3drew"
dirs = os.listdir(hdir)
master_blaster = pd.DataFrame()
for dr in dirs:
    if os.path.isdir(os.path.join(hdir, dr)):
        if re.match(r'\b\w{4}\b', dr):
            name = re.match(r'\b\w{4}\b', dr).group(0)
            print(name)
            if os.path.isdir(os.path.join(hdir, dr, 'megadock')):
                os.chdir(os.path.join(hdir, dr, "megadock"))
                files = os.listdir()
                r_list = glob.glob('*_relaxed.pdb')
                ncmplx = len(r_list)
                if len(r_list) > 0:                   
                    for file in r_list:
                        num = re.search('[0-9]+', file).group(0) 
                        ppdb = PandasPdb()
                        ppdb.read_pdb(''.join(['complex_', str(num), '.pdb']))
                        df = ppdb.df['ATOM']
                        ppdb2 = PandasPdb()
                        ppdb2.read_pdb(file)
                        df2 = ppdb.df['ATOM']
                        df2.add_suffix('_rlx')
                        df3 = df.merge(df2)
                        df3['complex_name'] = name
                        master_blaster = master_blaster.append(df3)
                        print("done: " + file)
        
master_blaster.to_json('master_blaster.json')