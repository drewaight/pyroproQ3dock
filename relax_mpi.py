import sys
import fileinput
import concurrent.futures
import timeit
import subprocess
from random import random
from mpi4py import MPI
import numpy as np
import pandas as pd
from pathlib import Path
from communication import *
import drew_relax
import json

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def relax(path, loc_cmplx, rank):
  start_time = timeit.default_timer()
  with concurrent.futures.ThreadPoolExecutor() as executor:
    future = executor.submit(drew_relax.main, path, loc_cmplx)
    sc = future.result()
    executor.shutdown(wait=True)
  stop_time = timeit.default_timer()
  time_tot = int(start_time-stop_time)
  print('Processor %d finished'%rank)
  print('Time: %d'%time_tot)
  return(sc)

def main(listfile, path):
  st_time = timeit.default_timer()
  with open(path+'/'+listfile) as f:
    list = []
    for line in f:
      list.append(line.strip())

  z = None
  nproc = size
  ncmplx = np.int64(len(list))

  if rank == 0:
    x = np.arange(0, ncmplx, dtype=np.int32)
    print('Test started with %d processors.'%size)
    z = np.zeros([ncmplx, 41], dtype=np.float64)
    for x in x:
      zr = z[x]
      zr[0] = x

  z_scatt = scatter_array(z)

  for z in z_scatt:
    zi = z[0]
    zin = int(zi)
    loc_cmplx = list[zin]
    sc = relax(path, loc_cmplx, rank)
    for i in range (0, z.size-1):
      z[i+1] = sc[0][i]

  gath = np.zeros([ncmplx, 41], dtype=np.float64)
  gather_array(gath, z_scatt)
 
  if rank == 0:
    df = pd.DataFrame(gath, columns=['complex','complex_normalized', 'dG_cross', 'dG_cross/dSASAx100', 'dG_separated',
       'dG_separated/dSASAx100', 'dSASA_hphobic', 'dSASA_int', 'dSASA_polar',
       'delta_unsatHbonds', 'hbond_E_fraction', 'hbonds_int', 'nres_all',
       'nres_int', 'packstat', 'per_residue_energy_int', 'sc_value',
       'side1_normalized', 'side1_score', 'side2_normalized', 'side2_score',
       'fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 'fa_intra_sol_xover4',
       'lk_ball_wtd', 'fa_elec', 'pro_close', 'hbond_sr_bb', 'hbond_lr_bb',
       'hbond_bb_sc', 'hbond_sc', 'dslf_fa13', 'omega', 'fa_dun', 'p_aa_pp',
       'yhh_planarity', 'ref', 'rama_prepro', 'total_score'], dtype='object')
    df.sort_values(by=['complex'])
    print(df)
    df.to_csv(path+'/'+"test.csv")
    print("finished!")
    stp_time = timeit.default_timer()
    tim_tot = int(stp_time-st_time)
    print('Total Time: %d'%tim_tot)
    df.to_json(path +'/relaxed_list.json')

  # # MPI.Finalize()

if __name__ == '__main__':
  main(sys.argv[1],sys.argv[2])