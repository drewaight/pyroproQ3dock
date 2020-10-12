import sys
import fileinput
from pathlib import Path
from shutil import copyfile

def main(file,receptor):
    outfile = copyfile(file, Path(file).stem+'_blocking.pdb')
    #with open(outfile, 'r') as f:
    with fileinput.input(outfile, inplace=True) as f:
        for line in f:
            field = line.strip().split()
            if line.startswith('ATOM'):
                if field[4] == 'H' or field[4] =='L':
                    if int(field[5]) in range(1,6) or int(field[5]) in range(54,89) or int(field[5]) in range(104, 142):
                        line = line.rstrip('\n')
                        print(line)
                    else:
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
    main(sys.argv[1],sys.argv[2])
    if len(sys.argv) < 2:
        print(
            "Creates the megadock blocking file and writes output to SAMPLE.table\nusage: parapred_to_blocking_2.py <antibody file> <receptor file>")
        sys.exit()
        file = sys.argv[1]
        receptor = sys.argv[2]