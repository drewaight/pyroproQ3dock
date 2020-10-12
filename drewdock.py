import os
import re
import sys
import fileinput
import subprocess
from pathlib import Path
import docking_structure_randomizer
import anarci_to_pdb
import parapred_to_blocking
import decoygen
import zrank_processing
import rmsd
import numpy as np
import pandas as pd
import chain_check
import collector
import natural_sort
import timeit

def main(structure, exec, ProQ3_dir, decoynum):
    start_time = timeit.default_timer()
    topdir = os.getcwd()
    chain = subprocess.run(['pdbcount', structure], stdout=subprocess.PIPE)
    counts = chain.stdout.split()
    if int(counts[1]) > 3:
        print("too many chains bud")
        return()
    else:
        print(chain.stdout.decode('utf-8'))
    chain_check.main(structure)
    subprocess.run(['mkdir', Path(structure).stem])
    subprocess.run(['mv', structure, Path(structure).stem])
    os.chdir(Path(structure).stem)
    docking_structure_randomizer.main(structure)
    test_str = Path(structure).stem + '_rand_LH.pdb'
    print(test_str)
    recep_str = Path(structure).stem + '_rand_A.pdb'
    out_str = Path(test_str).stem + '_cho.pdb'
    print(out_str)
    scheme = 'chothia'
    anarci_to_pdb.main(test_str, scheme, out_str)
    subprocess.run(['parapred', 'pdb', out_str])
    parapred_to_blocking.main(out_str,recep_str)
    subprocess.run(['/usr/bin/mpiexec', '-n', '16', '--use-hwthread-cpus', 'megadock-gpu-dp',
                   '-v', '1.0', '-D', '-t', '3', '-N', '54000', '-tb', 'SAMPLE.table'])
    #subprocess.run(['/usr/bin/mpiexec', '-n', '16', '--use-hwthread-cpus', 'megadock-gpu-dp',
    #                '-tb', 'SAMPLE.table'])
    subprocess.run(['mkdir', 'megadock'])
    subprocess.run(['mv', Path(structure).stem + '_rand_LH_cho_blocking-' + Path(structure).stem + '_rand_A.out',
                    'megadock'])
    subprocess.run(['cp', recep_str, 'megadock'])
    subprocess.run(['cp', test_str, 'megadock'])
    subprocess.run(['cp', Path(structure).stem+ '.clean.pdb', 'megadock'])
    os.chdir('megadock')
    subprocess.run(['mv', Path(structure).stem+ '.clean.pdb', 'complex_0.pdb'])
    subprocess.run(['mv', Path(structure).stem + '_rand_LH_cho_blocking-' + Path(structure).stem + '_rand_A.out',
                    Path(structure).stem + '_complex.out'])
    infile = Path(structure).stem + '_complex.out'
    df_zsc = zrank_processing.main(infile, recep_str, test_str, exec)
    print(df_zsc.describe())
    zinfile = Path(structure).stem + '_complex_zranked.out'
    complex_list = decoygen.main(recep_str, zinfile, decoynum, test_str)

    subprocess.run([exec+'ProQDock_prepare.sh', 'complex_0.pdb'])

    with open('proQ3_list.txt', 'w') as p:
        p.write('complex_0_merge.pdb'+'\n')
        for cmplx in complex_list:
            subprocess.run([exec+'ProQDock_prepare.sh', cmplx])
            p.write(Path(cmplx).stem + '_merge.pdb'+'\n')

    try:
        print('ProQ3 processing')
        subprocess.run([ProQ3_dir + 'proq3_all.sh', './complex.fasta', 'proQ3_list.txt', topdir+'/'+Path(structure).stem+'_profile',
                                    './run_all', '16'], stdout=subprocess.PIPE)
    except:
        print('ProQ3 parallel error')
    
    os.chdir('run_all')
    ProQ_list = []
    for file in os.listdir():
        if file.endswith(".pdb"):
            ProQ_list.append(file)
    ProQ_list_sorted = natural_sort.main(ProQ_list)

    total_sc = collector.main()
    rmsd_df = pd.DataFrame(columns=['rmsd','all_rmsd', 'sc_score', 'p2sc']) 
    for i, p_cmplx in enumerate(ProQ_list_sorted):
        rmd_sc = rmsd.main('../../' + structure, p_cmplx)
        rmsd_df.loc[i] = [rmd_sc[0][0], rmd_sc[0][1], rmd_sc[0][2], rmd_sc[0][3]]

    total_sc = total_sc.join(rmsd_df)    
    os.chdir('../../')
    outfp = Path(structure).stem + '_TOTAL.csv'
    total_sc.to_csv(outfp, sep='\t', mode='a')
    stop_time = timeit.default_timer()
    time_tot = str(start_time-stop_time)
    with open ('time.txt', 'w') as t:
        t.write(time_tot)
    os.chdir('..')
    return(total_sc)

if __name__ == '__main__':
    main(sys.argv[1], exec = '/media/hdd1/proQ3drewdock_big/drewdock_exec', ProQ3_dir = '/home/drewaight/proq3/', decoynum = 10)
    
