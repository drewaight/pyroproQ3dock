import numpy as np
import pandas as pd
import fileinput
import sys
from pathlib import Path
import subprocess
import re
import zrank_mpi

def main(megadockfile,ligfile,recfile,exec):
    ligh = ligfile+'.h'
    rech = recfile+'.h'
    with open(ligh, 'w') as redlig:
        subprocess.run(['reduce', ligfile], stdout=redlig)
    with open(rech, 'w') as redrec:
        subprocess.run(['reduce', recfile], stdout=redrec)
    for line in fileinput.input(megadockfile, inplace=1):
        if ligfile in line:
            line = line.replace(ligfile, ligh)
        if recfile in line:
            line = line.replace(recfile,rech)
        sys.stdout.write(line)
    subprocess.run([exec+'zrank', megadockfile, '1', '10000'])
    zrankfile = megadockfile+'.zr.out'
    outfile = Path(megadockfile).stem+'_zranked.out'
    df1 = pd.read_csv(megadockfile, header= 3, sep='\t', names=['rot1','rot2','rot3','vox1','vox2','vox3','score'])
    df2 = pd.read_csv(zrankfile, sep='\t', names=['id','zrank'])
    df1 = df1.join(df2)
    df1['id'] = df1['id'].values.astype(np.int64)
    df1['score'] = df1['score'].round(2)
    df1.sort_values(by=['zrank'], inplace=True)
    outfp = open(outfile, 'w')
    with open(megadockfile, 'r') as fp:
        for i, line in enumerate(fp):
            if i < 4:
                outfp.write(line)
    df1['rot1'] = df1['rot1'].map(lambda x: '%.6f' % x)
    df1['rot2'] = df1['rot2'].map(lambda x: '%.6f' % x)
    df1['rot3'] = df1['rot3'].map(lambda x: '%.6f' % x)
    df1['score'] = df1['score'].map(lambda x: '%.2f' % x)
    df1.drop(['id','zrank'],axis=1, inplace=True)
    #df1['zrank'] = df1['zrank'].map(lambda x: '%.5f' % x)
    df1.to_csv(outfp, sep='\t', mode='a', header=False, index=False)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

