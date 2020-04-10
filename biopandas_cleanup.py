import sys
import numpy as np
import pandas as pd
import re
import subprocess
from pathlib import Path
import fileinput
#outfile='5veb_renum.pdb'

# ----------------------------------------------------------Anarci--------------------------------------------------

def anarci(seq,scheme):
    out =subprocess.run(['anarci', '--sequence', seq, '--scheme', scheme, '--assign_germline'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = out.stdout.decode("utf-8").splitlines()
    stats = {}
    hd = []
    anarci_dict = {'chain':[], 'new_residue_number':[], 'new_insertion':[]}
    for line_i, line in enumerate(output):
        l = line.strip()
        li = line.split()
        if line.startswith('#'):
            hspl= line.split('|')
            hd.append(hspl)
        else:
            if len(li) == 3:
                if li[2] == '-':
                    continue
                else:
                    anarci_dict['chain'].append(li[0])
                    anarci_dict['new_residue_number'].append(int(li[1]))
                    anarci_dict['new_insertion'].append('')
                    residue = li[2]
                    
            elif len(li) == 4:
                anarci_dict['chain'].append(li[0])
                anarci_dict['new_residue_number'].append(int(li[1]))
                anarci_dict['new_insertion'].append(str(li[2]))
                residue = li[2]   

    stats = {'chain' : hd[5][2], 'start': hd[5][5], 'stop': hd[5][6], 'species' : hd[8][1], 
            'v_gene' : hd[8][2], 'v_id' : hd[8][3], 'j_gene' : hd[8][4], 'j_id' : hd[8][5]}     
    df_renum = pd.DataFrame.from_dict(anarci_dict)
    return(stats, df_renum)
    
# ------------------------------------------------BioPandas Cleanup--------------------------------------------------

def main(pdbfile, scheme, outfile):

    from biopandas.pdb import PandasPdb
    
    ppdb = PandasPdb()
    ppdb.read_pdb('./'+pdbfile)

    chains = ppdb.df['ATOM']['chain_id'].unique()
    chain_num = len(chains)
    print("Number of chains:", chain_num)
    print('\nRaw PDB file contents:\n\n%s\n...' % ppdb.pdb_text[:1000])

    header  = ppdb.pdb_text[:1000].split('\n')
    complex = []
    num = 1
    for row in header:
        if re.search(r'\bPAIRED_HL\b', row):
            fields = row.split()
            hchain = fields[3].split('=')[1]
            lchain = fields[4].split('=')[1]
            agchain = fields[5].split('=')[1]
            agtype = fields[6].split('=')[1]
            complex.append([num, hchain, lchain, agchain, agtype])
            num += 1

    del ppdb.df['HETATM']
    del ppdb.df['ANISOU']
    del ppdb.df['OTHERS']

    for cmplx in complex:
        compchains = ppdb.df['ATOM'][((ppdb.df['ATOM']['chain_id'] == cmplx[0]) | 
                                    (ppdb.df['ATOM']['chain_id'] == cmplx[1]) |
                                    (ppdb.df['ATOM']['chain_id'] == cmplx[2])) &
                                    ((ppdb.df['ATOM']['atom_name'] == 'C') | 
                                    (ppdb.df['ATOM']['atom_name'] == 'O') | 
                                    (ppdb.df['ATOM']['atom_name'] == 'N') | 
                                    (ppdb.df['ATOM']['atom_name'] == 'CA'))]
        bfact_avg = compchains['b_factor'].mean()
        cmplx.append(bfact_avg)  
    bestcmplx = complex[0]
    for cmplx in complex:
        if cmplx[5] < bestcmplx[5]:
            bestcmplx = cmplx

    print("Best Complex is:", bestcmplx)

    ppdb.df['ATOM'] = ppdb.df['ATOM'][(ppdb.df['ATOM']['chain_id'] == bestcmplx[1]) |
                                      (ppdb.df['ATOM']['chain_id'] == bestcmplx[2]) | 
                                      (ppdb.df['ATOM']['chain_id'] == bestcmplx[3])]

    ppdb.df['ATOM'].loc[ppdb.df['ATOM'].chain_id == bestcmplx[1], ['chain_id']] = 'H'   
    ppdb.df['ATOM'].loc[ppdb.df['ATOM'].chain_id == bestcmplx[2], ['chain_id']] = 'L'   
    ppdb.df['ATOM'].loc[ppdb.df['ATOM'].chain_id == bestcmplx[3], ['chain_id']] = 'A'   
    df =ppdb.df['ATOM'].copy()
    hdf = df[(df['chain_id'] == 'H')].reset_index(drop=True)
    ldf = df[(df['chain_id'] == 'L')].reset_index(drop=True)
    adf = df[(df['chain_id'] == 'A')].reset_index(drop=True)
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 'A']

    def seq_order(df):
        from collections import OrderedDict
        df['residue_insertion'] = df['residue_number'].astype(str)+df['insertion'].astype(str)
        ordered_seq = list(OrderedDict.fromkeys(df['residue_insertion']))
        seq_dict = {ordered_seq[i]: i+1 for i in range(0, len(ordered_seq))}
        df['residue_insertion'] = df['residue_insertion'].map(seq_dict)
        df['residue_number'] = df['residue_insertion']
        df.drop(['residue_insertion'], axis=1, inplace = True)
        df['insertion'] = None
        return(df)
    
    hdf = seq_order(hdf)
    ldf = seq_order(ldf)
    ppdb.df['ATOM'] = ppdb.df['ATOM'].append(hdf)
    ppdb.df['ATOM'] = ppdb.df['ATOM'].append(ldf).reset_index(drop=True)

    def df_to_seq(chain):
        seq = ppdb.amino3to1()
        seq = str(''.join(seq.loc[seq['chain_id'] == chain, 'residue_name']))
        return(seq)

    h_seq = df_to_seq('H')
    l_seq = df_to_seq('L')

    hdf_stats, hdf_renum = anarci(h_seq, scheme)
    hdf_renum.drop(['chain'], inplace = True, axis=1)
    hdf_renum['residue_number'] = np.arange(len(hdf_renum))
    hdf.drop(hdf[hdf['residue_number'] > int(hdf_stats['stop'])].index, inplace=True)
    hdf = hdf.merge(hdf_renum, how='left') \
    .drop(columns=['residue_number', 'insertion']) \
    .rename(columns={'new_residue_number': 'residue_number','new_insertion': 'insertion'})

    ldf_stats, ldf_renum = anarci(l_seq, scheme)
    ldf_renum.drop(['chain'], inplace = True, axis=1)
    ldf_renum['residue_number'] = np.arange(len(ldf_renum))
    ldf.drop(ldf[ldf['residue_number'] > int(ldf_stats['stop'])].index, inplace=True)
    ldf = ldf.merge(ldf_renum, how='left') \
    .drop(columns=['residue_number', 'insertion']) \
    .rename(columns={'new_residue_number': 'residue_number','new_insertion': 'insertion'})

    def atom_renum(df):
        from collections import OrderedDict
        atom_list = list(OrderedDict.fromkeys(df['atom_number']))
        atom_dict = {atom_list[i]: i+1 for i in range(0, len(atom_list))}
        df['atom_number'] = df['atom_number'].map(atom_dict)
        return(df)
    
    ppdb.df['ATOM'].drop(ppdb.df['ATOM'].index, inplace=True)
    hdf = hdf.append(ldf, sort=False)
    hdf = hdf.append(adf, sort=False)
    hdf = atom_renum(hdf).reset_index(drop=True)
    ppdb.df['ATOM'] = ppdb.df['ATOM'].append(hdf, sort=False).reset_index(drop=True)
    ppdb.to_pdb(outfile, records=['ATOM'])
    print(ppdb.df['ATOM'].head())

    with open(Path(outfile).stem + '_HL.pdb', 'w') as f1:
        subprocess.run(['pdb_selchain','-H,L', outfile],stdout=f1)
    
    with open(Path(outfile).stem + '_A.pdb', 'w') as f2:
        subprocess.run(['pdb_selchain','-A', outfile],stdout=f2)

    with open(outfile, 'w') as f3:
        subprocess.run(['pdb_merge', Path(outfile).stem + '_HL.pdb', Path(outfile).stem + '_A.pdb'],stdout=f3)

    h_stats = pd.DataFrame.from_records([hdf_stats])
    l_stats = pd.DataFrame.from_records([ldf_stats])

    return(h_stats, l_stats)

if __name__ == "__main__":
    main(pdbfile=sys.argv[1], scheme='chothia', outfile=sys.argv[2])


