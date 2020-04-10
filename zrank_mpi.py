from mpi4py import MPI
import numpy as np
import time
import subprocess
from pathlib import Path
from communication import *
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def main(infile,exec,num):
  infile = infile
  z = None
  nproc = size
  nmega = np.int64(num)

  if rank == 0:
    x = np.arange(0, nmega, dtype=np.int32)
    print('Test started with %d processors.'%size)
    z = np.zeros([nmega, 2], dtype=np.float32)
    for x in x:
      zr = z[x]
      zr[0] = x+1
    #print('Total Array:\n%s\n'%z)

  comm.Barrier()

  z_scatt = scatter_array(z)

  for z in z_scatt:
    zi = z[0]
    zis = str(int(zi))
    loc_infile = Path(infile).stem+'_'+zis+'.out'
    subprocess.run(['cp', infile, loc_infile])
    subprocess.run(['zrank', loc_infile, zis, zis])
    subprocess.run(['rm', loc_infile])
    loc_zrfile = Path(infile).stem + '_' + zis + '.out.zr.out'
    with open(loc_zrfile, 'r') as zin:
      for line in zin:
        l= line.split()
        z[1] = l[1]
    subprocess.run(['rm', loc_zrfile])

  gath = np.zeros([nmega,2], dtype=np.float32)
  x_gath = gather_array(gath, z_scatt)
  # print(x_gath)
  if rank == 0:
    sorted = gath[np.argsort(gath[:, 1])]
    gathout = Path(infile).stem + '_unsorted.zr.out'
    sortedout = Path(infile).stem + '_sorted.zr.out'
    np.savetxt(gathout, gath, fmt= '%5d \t %5f ')
    np.savetxt(sortedout, sorted, fmt = '%5d \t %5f')

  MPI.Finalize()

if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3])