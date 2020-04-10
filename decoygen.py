import re
import os
import sys
import subprocess
from pathlib import Path
import fileinput

decoygen_path = '/home/drewaight/drewdock'

def main(lig_file, infile, num, catfile):
     num_rank = int(num)
     list = []
     outlist = []
     prefix = Path(lig_file).stem
     f = open(infile)
     for i, line in enumerate(f):
          if i < num_rank:
               ligout_file = 'lig_'+prefix+'.'+str(i+1)+'.pdb'
               subprocess.run([decoygen_path +'/decoygen', ligout_file, lig_file, infile, str(i)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
               list.append(ligout_file)
          else:
               continue

     else:
          for i, l in enumerate(list):
               outfilename = 'complex_'+str(i+1)+'.pdb'
               outfile = open(outfilename, "w")
               out = subprocess.run(['cat', catfile, l], stdout=outfile)
               outlist.append(outfilename)


     return(outlist)

if __name__ == "__main__":
     main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
     if len(sys.argv) < 3:
          print(
               "usage:\n\t decoygen.py <ligand file> <megadock output file> <return number>\noptional: <receptor file>\n\t ex. decoygen.py EGFR.pdb dock.out 10 panitumumab.pdb")

     try:
          lig_file = sys.argv[1]
     except IndexError:
          quit()
     try:
          infile = sys.argv[2]
     except IndexError:
          quit()
     try:
          num_rank = int(sys.argv[3])
     except IndexError:
          quit()
     try:
          catfile = sys.argv[4]
     except NameError:
          print("ligand only")