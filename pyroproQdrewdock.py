#########################################################################################################
#   pyroproQdrewdock.py 0.2 - pyRosetta ProQ3D docking program 04/09/2020 Andrew Waight
#   takes chothia labelled input crystal structures from SabDab and returns relaxed and scored decoys
#   using pyRosetta and ProQ3D 
# 
# 
# 
# 
#########################################################################################################

import os
import re
import sys
import fileinput
import subprocess
from pathlib import Path
import antibody_info_and_dock_prep
import biopandas_cleanup
import pandas_anarci_numberer
import bfactor_to_blocking_biopandas
import decoygen
import zrank_processing
import numpy as np
import pandas as pd
import chain_check
import collector
import natural_sort
import timeit

def main(structure, sab_dir, exec, ProQ3_dir, decoynum):
    # The .theano directory gets very large RAM usage over time using ProQ3D
    subprocess.run(['rm', '-rf', '/home/drewaight/.theano'])
    start_time = timeit.default_timer()

    # Working directory where you want all the results stored
    topdir = os.getcwd()
    scheme = 'chothia'

    # Make a folder and copy the target structure into it
    subprocess.run(['mkdir', Path(structure).stem])
    subprocess.run(['cp', sab_dir + structure, Path(structure).stem])
    os.chdir(Path(structure).stem)

    # The biopandas_cleanup module selects the best vL vH Epitope complex from a structure it uses Bfactors as criteria
    # the best complex is relabelled 'L' 'H' 'A' and the antibodies are renumbered by Anarci with the input scheme
    # and returns dataframes with the heavy and light stats. See biopandas_cleanup.py for more information
    struct_clean = Path(structure).stem + '_1.pdb'
    h_stats, l_stats = biopandas_cleanup.main(structure, scheme, struct_clean)

    # RaBdaB takes IgG Kappa definition pulled here from the anarci stats
    if re.search(r"IGK", str(l_stats['v_gene'])):
        light_chain = 'kappa'
    elif re.search(r"IGL", str(l_stats['v_gene'])):
        light_chain = 'lambda'
    else:
        light_chain = 'undefined'

    cwd = os.getcwd()

    # The program antibody_info_and_dock.prep.py calls structure_extract_relax.py which uses pyrosetta to run packing and
    # minimization on the input structure and extracts the releveant antibody CDR statistics and Dunbrack clusters into a returned dataframe.
    # After this the Aho numbered CDRs + 3 stem residues are set to 1.0 in the Bfactor column and all other residues are set to 0.0
    # the antibody_info_and_dock.prep then separates the LH and A molecules randomizes the orientation and slides them back together
    # in preparation for docking. 
    abdf = antibody_info_and_dock_prep.main(cwd, struct_clean, light_chain)
        
    h_stats.to_json('h_stats.json')
    l_stats.to_json('l_stats.json')
    abdf.to_json('abdf.json')

    # At this point we want to fork the two prepared structures into the broad CDR definition above, and a targeted defintion using
    # Parapred. The randomized LH structure is copied with the ppd extension to denote Parapred scoring
    test_str = Path(struct_clean).stem + '_rand_LH.pdb'
    test2_str = Path(struct_clean).stem + '_ppd_rand_LH.pdb'
    subprocess.run(['cp', test_str, test2_str])
    recep_str = Path(struct_clean).stem + '_rand_A.pdb'
    init_str  = Path(struct_clean).stem + '_init.pdb'
    print([test_str, recep_str, init_str])
    
    # Execute parapred on the copied structure 
    subprocess.run(['parapred', 'pdb', test2_str])
    
    # The bfactor_to_blocking_biopandas.py program simply converts all residues with bfactor == 0 to the residue BLK which avoids rigid
    # body docking on these residues (ZDOCK and Megadock) standard
    bfactor_to_blocking_biopandas.main(test_str)
    bfactor_to_blocking_biopandas.main(test2_str)
    block_cdr = Path(test_str).stem + "_blocking.pdb"
    parap_cdr = Path(test2_str).stem + "_blocking.pdb"

    # Writing the jobs table for the megadock program
    with open("SAMPLE.table", "w") as f:
        f.write("TITLE= sample jobs\n")
        f.write("PARAM= -R $1 -L $2\n")
        f.write(block_cdr+"\t"+recep_str+"\n")
        f.write(parap_cdr+"\t"+recep_str)

    # Execute megadock GPU docking on both the broader CDR and targeted Parapred treated structures
    subprocess.run(['/usr/bin/mpiexec', '-n', '16', '--use-hwthread-cpus', 'megadock-gpu-dp',
                    '-v', '1.0', '-D', '-t', '3', '-N', '54000', '-tb', 'SAMPLE.table'])
    # subprocess.run(['/usr/bin/mpiexec', '-n', '16', '--use-hwthread-cpus', 'megadock-gpu-dp',
    #                '-tb', 'SAMPLE.table'])
    
    # Everything gets moved into the megadock folder
    subprocess.run(['mkdir', 'megadock'])
    subprocess.run(['mv', Path(struct_clean).stem + '_rand_LH_blocking-' + Path(struct_clean).stem + '_rand_A.out',
                    'megadock'])
    subprocess.run(['mv', Path(struct_clean).stem + '_ppd_rand_LH_blocking-' + Path(struct_clean).stem + '_rand_A.out',
                    'megadock'])
    subprocess.run(['cp', recep_str, 'megadock'])
    subprocess.run(['cp', test_str, 'megadock'])
    subprocess.run(['cp', test2_str, 'megadock'])
    subprocess.run(['cp', init_str, 'megadock'])
    os.chdir('megadock')
    subprocess.run(['mv', init_str, 'complex_0.pdb'])
    subprocess.run(['mv', Path(struct_clean).stem + '_rand_LH_blocking-' + Path(struct_clean).stem + '_rand_A.out',
                     Path(struct_clean).stem + '_complex.out'])
    subprocess.run(['mv', Path(struct_clean).stem + '_ppd_rand_LH_blocking-' + Path(struct_clean).stem + '_rand_A.out',
                     Path(struct_clean).stem + '_ppd_complex.out'])                     
    
    # The zrank_processing.py module rescores the megadock output files using the zrank program, it calls the mpi enabled module zrank_mpi.py 
    # and returns a dataframe of the results sorted by zscore
    infile = Path(struct_clean).stem + '_complex.out'
    infile2 = Path(struct_clean).stem + '_ppd_complex.out'
    df_zsc = zrank_processing.main(infile, recep_str, test_str, exec)
    df2_zsc = zrank_processing.main(infile2, recep_str, test2_str, exec)

    # We want to keep trach of which method performs better, so add a identifier column before merging and resorting the zranked datasets
    df_zsc['method'] = "CDR"
    df2_zsc['method'] = "Parapred"
    df_zsc = pd.concat([df_zsc, df2_zsc], ignore_index=True)
    df1 = df_zsc.sort_values(by=['zrank'])
    df1 = df1.reset_index(drop=True)

    # Reformat the merged zranked dataframes from back into megadock format
    outfile = Path(struct_clean).stem + '_complex_zranked.out'
    outfp = open(outfile, 'w')
    with open(infile, 'r') as fp:
        for i, line in enumerate(fp):
            if i < 4:
                outfp.write(line)
    df1['rot1'] = df1['rot1'].map(lambda x: '%.6f' % x)
    df1['rot2'] = df1['rot2'].map(lambda x: '%.6f' % x)
    df1['rot3'] = df1['rot3'].map(lambda x: '%.6f' % x)
    df1['score'] = df1['score'].map(lambda x: '%.2f' % x)
    df1['zrank'] = df1['zrank'].map(lambda x: '%.5f' % x)

    # We want to keep track of the zranked stats, this data is split off as dfout, the rest gets re-written to the outfile for decoy generation
    df3 = df1.drop(['id','zrank','method'],axis=1)
    df3.to_csv(outfp, sep='\t', mode='a', header=False, index=False)
    dfout = df1.drop(['rot1','rot2','rot3','vox1','vox2','vox3','id','score'],axis=1)
    dfout['complex'] = dfout.index + 1
    print(dfout.head())
    dfout.to_csv('zrank.csv')
    dfout.to_json('zrank.json')

    # The decoygen.py program converts the zranked list of rotational coordinates to back into pdb complexes. The number of complexes
    # returned is set by the decoynum variable, complex_0.pdb (native complex) is added to the beginning of the list as the reference.
    zinfile = Path(struct_clean).stem + '_complex_zranked.out'
    complexes = decoygen.main(recep_str, zinfile, decoynum, test_str)
    complexes.append('complex_0.pdb')
    complex_list = natural_sort.main(complexes)
    with open('relax_list.txt', 'w') as r:
        for cmplx in complex_list:
            r.write(cmplx+'\n')

    cwd = os.getcwd()
    jsonfile="relaxed_list.json"
    reference='complex_0_relaxed.pdb'

    # The relax_mpi.py program uses mpi to run the drew_relax.py program over 16 cores. The output is a json file of the relaxed and interface
    # parameters of each complex
    print('Running Relax:')
    subprocess.run(['mpiexec', '-n', '16', 'python', exec+'relax_mpi.py', 'relax_list.txt', cwd], 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # The rmsd_align_mpi.py program uses mpi to run rmsd.py which aligns all of the relaxed complexes by their Antigen chain, and calculates 
    # the RMSD value from the reference pose
    print('Align and Rmsd:')
    subprocess.run(['mpiexec', '-n', '16', 'python', exec+'rmsd_align_mpi.py', cwd, jsonfile, reference], 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Read back in the output json file from the above mpi programs
    ros_df = pd.read_json(cwd + '/relaxed_list_rmsd.json')
    ros_df = ros_df.drop(['dG_cross', 'dG_cross/dSASAx100', 'packstat'], axis=1)

    # We want ProQ3D to use our relaxed structures write the list for it to process
    with open('proQ3_list.txt', 'w') as p:
        for cmplx in complex_list:
            pandas_anarci_numberer.main(Path(cmplx).stem + '_relaxed.pdb', 'chothia', Path(cmplx).stem + '_relaxed_cho.pdb' )
            subprocess.run(['bash', exec+'ProQDock_prepare.sh', Path(cmplx).stem + '_relaxed_cho.pdb'])
            p.write(Path(cmplx).stem + '_relaxed_cho_reres.pdb'+'\n')

    # ProQ3D is a big complicated program. It is called using GNU parallel through the proq3_all.sh script, if there is not yet a profile
    # one is created in the top level working directory. All ProQ3D output is put into the run_all directory
    try:
        print('ProQ3 processing')
        subprocess.run([ProQ3_dir + 'proq3_all.sh', './complex.fasta', 'proQ3_list.txt', topdir+'/'+Path(structure).stem+'_profile',
                                    'run_all', '16'], stdout=subprocess.PIPE)
    except:
        print('ProQ3 parallel error')
    
    os.chdir('run_all')
    ProQ_list = []
    for file in os.listdir():
        if file.endswith(".pdb"):
            ProQ_list.append(file)
    ProQ_list_sorted = natural_sort.main(ProQ_list)
 
    # The collector.py program scrapes the scores from the proq3.global files and returns a dataframe
    total_sc = collector.main()
    # Join the rosetta scores
    total_sc = total_sc.join(ros_df)
    total_sc.to_json('total_sc.json')
    # Merge the zrank data into the total scores
    total_sc = total_sc.merge(dfout, how='left', on='complex') 
    total_sc['structure'] = Path(structure).stem  
    abdf['structure'] = Path(structure).stem
    # Merge the antibody information the single statistics for the antibody including CDR classifications are populated to all the complexes
    total_sc = total_sc.merge(abdf, how='right')
    # Shuffle important columns to the front
    col = total_sc.pop("zrank")
    total_sc.insert(0, col.name, col)
    col = total_sc.pop("method")
    total_sc.insert(0, col.name, col)
    col = total_sc.pop("structure")
    total_sc.insert(0, col.name, col)
    print(total_sc)
    os.chdir('../../')
    outfp = Path(structure).stem + '_TOTAL.csv'
    # Sort by rmsd and output
    total_sc = total_sc.sort_values(by=['rmsd'])
    total_sc.to_csv(outfp, sep='\t', mode='a')
    total_sc.to_json(Path(structure).stem + '_TOTAL.json')
    stop_time = timeit.default_timer()
    time_tot = str(stop_time-start_time)
    with open ('time.txt', 'w') as t:
        t.write(time_tot)
    os.chdir('..')
    # Return the dataframe and thats it!
    return(total_sc)  

if __name__ == '__main__':
    main(sys.argv[1], exec = '/media/hdd1/roproQ3drew/drewdock_exec/', sab_dir = '../sabdab_short/',
    ProQ3_dir = '/home/drewaight/proq3/', decoynum = 96)
    
