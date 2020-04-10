from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.docking import *
from rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.simple_moves import *
import numpy as np
import pandas as pd


def main(path, cmplx,reference):
    scores = []
    init()
    pose = pose_from_pdb(path + '/' + cmplx)
    reference = pose_from_pdb(path + '/' + reference)
    acm = AlignChainMover()
    acm.source_chain(3)
    acm.target_chain(3)
    acm.pose(reference)
    acm.apply(pose)
    pose.dump_pdb(path + '/' + cmplx)
    rmsd = CA_rmsd(reference, pose)
    return(rmsd)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
