import os
import sys
import re
import subprocess
from pathlib import Path
from Bio.PDB import *
from Bio.PDB.Polypeptide import *
from Bio.PDB.PDBIO import PDBIO


def main(pdbfile,scheme,outfile):
    pdb_io = PDBIO()
    parser = PDBParser()
    structure = parser.get_structure('self',pdbfile)
    model = structure[0]
    chain = ' '
    anarci_dict = {}
    fastafile = Path(pdbfile).stem+'.fasta'
    with open(fastafile, 'w+') as ff:
        subprocess.run(['pdb_tofasta', '-multi', pdbfile], stdout=ff)
    out =subprocess.run(['anarci', '-i', fastafile, '--scheme', scheme], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = out.stdout.decode("utf-8").splitlines()
    for line in output:
        l = line.strip()
        li = line.split()
        if line.startswith('#'):
            if re.search("PDB", l):
                ch = l.split('|')
                chain = ch[1]
                anarci_dict[chain] = []
        elif line.startswith(chain):
            if len(li) == 3:
                if li[2] == '-':
                    continue
                else:
                    ch = li[0]
                    resseq = int(li[1])
                    inscode = ' '
                    residue = li[2]
                    het = ' '
                    tuple = (het, resseq, inscode)
                    anarci_dict[chain].append(tuple)
            elif len(li) == 4:
                ch = li[0]
                resseq = int(li[1])
                inscode = str(li[2])
                residue = li[3]
                het = ' '
                tuple = (het, resseq, inscode)
                anarci_dict[chain].append(tuple)
        else:
            continue

    for residue in structure.get_residues():
        residue.id = (' ', residue.id[1] + 900, residue.id[2])

    for d in anarci_dict:
        count = 0
        for i, residue in enumerate(model[d]):
            if i < len(anarci_dict[d]):
                residue.id = anarci_dict[d][i]
                count = residue.id[1]
            else:
                count = count + 1
                residue.id = ( ' ', count, ' ')


    pdb_io.set_structure(structure)
    pdb_io.save(outfile)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])
    if len(sys.argv) < 3:
        print("Writes the Anarci scheme into the pdb file\nusage: anarci_to_pdb.py <pdbfile> <scheme> <outfile>")
        sys.exit()