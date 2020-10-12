import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb


def main(pdbfile):
    ppdb = PandasPdb()
    ppdb.read_pdb('./'+pdbfile)
    pd = ppdb.df['ATOM'][(ppdb.df['ATOM']['chain_id'] != 'H') & (ppdb.df['ATOM']['chain_id'] != 'L') & (ppdb.df['ATOM']['chain_id'] != 'A')]
    if len(pd) != 0:
        for i in pd.index:
            ppdb.df['ATOM'].loc[i,'chain_id'] = ppdb.df['ATOM'].loc[i,'chain_id'] = 'A'

    print(ppdb.df['ATOM'].loc[pd.index]['chain_id'])
    ppdb.to_pdb(path=pdbfile)

if __name__ == "__main__":
    main(sys.argv[1])
