import os
import re
import sys
import fileinput
from pathlib import Path
import numpy as np
import pandas as pd
import natural_sort

def main():
    complex_list = []
    for file in os.listdir():
        if file.endswith(".pdb"):
            complex_list.append(file)
    complex_list_sorted = natural_sort.main(complex_list)

    df_pQ3sc = pd.DataFrame(columns=['name','ProQ2D', 'ProQRosCenD', 'ProQRosFAD', 'ProQ3D']) 
    
    for i, cmplx in enumerate(complex_list_sorted):
        name = cmplx
        with open(name+'.proq3.global', 'r') as pq_sc:
            lines = pq_sc.readlines()
            scores = lines[1].split()
            df_pQ3sc.loc[i] = [name, scores[0], scores[1], scores[2], scores[3]]                                                                                                                       
    
    return(df_pQ3sc)

if __name__ == '__main__':
    main()