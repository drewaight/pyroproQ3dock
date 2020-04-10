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
import rmsd
import json

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def align_rmsd(path, cmplx, reference):
  start_time = timeit.default_timer()
  with concurrent.futures.ThreadPoolExecutor() as executor:
    future = executor.submit(rmsd.main, path, cmplx, reference)
    rms = future.result()
    executor.shutdown(wait=True)
  print('Processor %d finished'%rank)
  return(rms)

def main(path, jsonfile, reference):
  st_time = timeit.default_timer()
  ros_df = pd.read_json(path + '/' + jsonfile)
  print(ros_df)
  z = None
  nproc = size
  ncmplx = np.int64(len(ros_df.index))

  if rank == 0:
    x = np.arange(0, ncmplx, dtype=np.int32)
    print('Test started with %d processors.'%size)
    z = np.zeros([ncmplx, 2], dtype=np.float64)
    for x in x:
      zr = z[x]
      zr[0] = x

  z_scatt = scatter_array(z)

  for z in z_scatt:
    zi = int(z[0])
    zin = str(zi)
    loc_cmplx = ('complex_'+zin+'_relaxed.pdb')
    rms = align_rmsd(path, loc_cmplx, reference)
    z[1] = rms

  gath = np.zeros([ncmplx, 2], dtype=np.float64)
  gather_array(gath, z_scatt)
 
  if rank == 0:
    df = pd.DataFrame(gath, columns=['complex','rmsd'], dtype='object')
    df.sort_values(by=['complex'])
    df = pd.merge(df, ros_df, on='complex', how='left')
    stp_time = timeit.default_timer()
    tim_tot = int(stp_time-st_time)
    print('Total Time: %d'%tim_tot)
    print(df)
    df.to_json(path +'/relaxed_list_rmsd.json')

  # MPI.Finalize()

if __name__ == '__main__':
  # main(path = '/media/hdd1/roproQ3drew/2xwt_nopack/megadock', jsonfile="relaxed_list.json", reference='complex_0_relaxed.pdb')
  main(sys.argv[1],sys.argv[2],sys.argv[3])