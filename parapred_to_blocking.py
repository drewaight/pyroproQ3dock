#! 
import os
import sys
import fileinput
from pathlib import Path
import shutil

def main(file,receptor):
    path = '../2vir/'
    outfile = path+Path(file).stem+'_blocking.pdb'
    shutil.copyfile(file, outfile)
    for line in fileinput.input(outfile, inplace=False):
        field = line.strip().split()
        if field[0] == 'ATOM':
            print(field)
            if float(field[9]) == 0.00:
                line = line.replace(field[3], 'BLK', 2).rstrip('\n')
                print(line)
            else:
                line = line.rstrip('\n')
                print(line)
        else:
            line = line.rstrip('\n')
            print(line)

    f = open("SAMPLE.table", "w")
    f.write("TITLE= sample jobs\n")
    f.write("PARAM= -R $1 -L $2\n")
    f.write(outfile+"\t"+receptor)
    f.close()

if __name__ == "__main__":
    

    file = "../2vir/2vir_1_rand_LH.pdb"
    receptor = "../2vir/2vir_1_rand_A.pdb"
    main(file, receptor)