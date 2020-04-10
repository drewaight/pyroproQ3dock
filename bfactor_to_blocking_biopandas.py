import sys
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
from biopandas.pdb import PandasPdb

def main(pdbfile):

    from biopandas.pdb import PandasPdb
    outfile = Path(pdbfile).stem+'_blocking.pdb'
    ppdb = PandasPdb()
    ppdb.read_pdb(pdbfile)
    ppdb.df['ATOM'].loc[ppdb.df['ATOM'].b_factor == 0.0, 'residue_name'] = 'BLK'
    ppdb.to_pdb(outfile, records=None)

if __name__ == "__main__":
    pdbfile = "../2vir/2vir_1_rand_LH.pdb"
    receptor = "../2vir/2vir_1_rand_A.pdb"
    main(pdbfile)

