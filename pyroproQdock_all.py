import pyroproQdrewdock
import natural_sort
from pathlib import Path
import numpy as np
import pandas as pd
import os
import re
import subprocess

sab_dir = '/home/drewaight/hdd1/structure_finder/output/output_folder/'
exec = '/media/hdd1/roproQ3drew/drewdock_exec/'
ProQ3_dir = '/home/drewaight/proq3/'
decoynum = 95


cwd = os.getcwd()
print(cwd)
score = pd.DataFrame()
master_list = []
for file in os.listdir(sab_dir):
    if file.endswith(".pdb"):
        master_list.append(file)

# master_list_sorted = natural_sort.main(master_list)
for structure in master_list:
    if os.path.isdir(cwd+'/'+ Path(structure).stem):   
        print('Looks like you already did ' + Path(structure).stem + ' bud.')
        continue
    else: 
        print('You didnt do ' + Path(structure).stem + ' yet')
        try:
            print("############################# NOW WORKING ON " + structure + " ##########################################")
            pyroproQdrewdock.main(structure, sab_dir, exec, ProQ3_dir, decoynum)
               
        except:
            dir = exec +"../"
            os.chdir(dir)
            with open(dir + "ERROR.txt", "a+") as e:
                e.write("failed "+ structure+ '\n')
            continue

