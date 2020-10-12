import drewdock
import drewdock_biopandas
import natural_sort
from pathlib import Path
import numpy as np
import pandas as pd
import os
import re

sab_dir = '/home/drewaight/hdd1/sabdab_short/'
exec = '/media/hdd1/proQ3drewdock_big/drewdock_exec/'
ProQ3_dir = '/home/drewaight/proq3/'
decoynum = 200

cwd = os.getcwd()
print(cwd)
score = pd.DataFrame()
master_list = []
for file in os.listdir(sab_dir):
    if file.endswith(".pdb"):
        master_list.append(file)

master_list_sorted = natural_sort.main(master_list)
for structure in master_list_sorted:
    if os.path.isdir(cwd+'/'+ Path(structure).stem):   
        print('Looks like you already did ' + Path(structure).stem + ' bud.')
        continue
    else: 
        print('You didnt do ' + Path(structure).stem + ' yet')
        try:
            print("############################# NOW WORKING ON " + structure + " ##########################################")
            sc = drewdock_biopandas.main(structure, sab_dir, exec, ProQ3_dir, decoynum)
            score = score.append(sc)
            score.to_csv('TOTAL.csv', sep='\t', mode='a')
        except:
            dir = exec +"../"
            os.chdir(dir)
            with open(dir + "ERROR.txt", "a+") as e:
                e.write("failed "+ structure+ '\n')
            continue

    
score.to_csv('MASTER_TOTAL.csv', sep='\t', mode='a')