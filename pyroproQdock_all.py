import pyroproQdrewdock
import natural_sort
from pathlib import Path
import numpy as np
import pandas as pd
import os
import re
import subprocess

sab_dir = '/home/drewaight/hdd1/sabdab_short/'
exec = '/media/hdd1/roproQ3drew/drewdock_exec/'
ProQ3_dir = '/home/drewaight/proq3/'
decoynum = 96
retries = 2
thresh = 8

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
            sc = pyroproQdrewdock.main(structure, sab_dir, exec, ProQ3_dir, decoynum)
            
            attempt = 1
            if sc.iloc[1]['rmsd'] > thresh:
                
                while attempt < retries:
                    os.chdir(cwd)
                    subprocess.run(['mv', Path(structure).stem, Path(structure).stem+"_"+str(attempt)])
                    print("######################### RMSD under threshold " + structure + "Starting attempt: " + str(attempt+1) + " #######################")
                    sc = pyroproQdrewdock.main(structure, sab_dir, exec, ProQ3_dir, decoynum)
                    attempt +=1
                    if sc.iloc[1]['rmsd'] <= thresh:
                        break
                    else:
                        continue

            rmsd = sc.iloc[1]['rmsd']                
            with open("attempts.txt", 'w') as a:
                a.write("Best RMSD: %f"%rmsd)
                a.write("Number of attempts: %d"%attempt)        
            score = score.append(sc)
            score.to_csv('TOTAL.csv', sep='\t', mode='a')

            
        except:
            dir = exec +"../"
            os.chdir(dir)
            with open(dir + "ERROR.txt", "a+") as e:
                e.write("failed "+ structure+ '\n')
            continue

    
score.to_csv('MASTER_TOTAL.csv', sep='\t', mode='a')