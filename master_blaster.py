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
master_blaster = pd.DataFrame(columns=['key', 'complex_name', 'complex', 'complex_relaxed', 'directory'])
for dr in dirs:
    if os.path.isdir(os.path.join(hdir, dr)):
        if re.match(r'\b\w{4}\b', dr):
            name = re.match(r'\b\w{4}\b', dr).group(0)
            print(name)
            if os.path.isdir(os.path.join(hdir, dr, 'megadock')):
                os.chdir(os.path.join(hdir, dr, "megadock"))
                cwd = os.getcwd()
                files = os.listdir()
                r_list = glob.glob('*_relaxed.pdb')
                if len(r_list) > 0:                   
                    for file in r_list:
                        d = {}
                        num = re.search('[0-9]+', file).group(0)
                        base = ''.join(['complex_', str(num), '_aho.pdb'])
                        d['key'] = ''.join([name, "_", base])
                        d['complex_name'] = name
                        d['complex'] = base
                        d['complex_relaxed'] = file
                        d['directory'] = cwd
                        print(d)
                        master_blaster = master_blaster.append(d, ignore_index=True)
                        print("done: " + file)
        
master_blaster.to_json(os.path.join(hdir, 'master_blaster.json'))