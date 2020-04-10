import fileinput
import tempfile
import sys
import subprocess
import structure_extract_relax
from pathlib import Path
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta.rosetta.protocols.docking import *
import pyrosetta.rosetta.protocols.rigid as rigid_moves

def main(path, pdb, light_chain):
    
    abdf = structure_extract_relax.main(path, pdb, light_chain)

    init()
    outfile = Path(pdb).stem+'_rand.pdb'
    pose = pose_from_pdb(pdb)
    init_pose = pose.clone()
    setup_foldtree(pose, "LH_A", Vector1([1]))
    print(pose.fold_tree())
    jump_num = 1
    randomize1 = rigid_moves.RigidBodyRandomizeMover(pose, jump_num, rigid_moves.partner_downstream)
    randomize1.apply(pose)
    slide = DockingSlideIntoContact(jump_num)
    slide.apply(pose)
    # docking_pre_prot = DockingPrepackProtocol()
    # docking_pre_prot.apply(pose)
    # docking_pre_prot.apply(init_pose)
    pose.dump_pdb(outfile)
    init_pose.dump_pdb(Path(pdb).stem+'_init.H.pdb')

    readout = open(outfile)
    abfile = open(Path(outfile).stem+'_LH.pdb', "w")
    with tempfile.TemporaryFile(mode='w+t') as temp:
        temp.writelines(readout)
        temp.seek(0)
        p1 = subprocess.Popen(['pdb_selchain', '-L,H'], stdin=temp, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['pdb_delelem','-H'], stdin=p1.stdout, stdout=abfile)
    abfile.close()
    temp.close()

    readout = open(outfile)
    ligfile = open(Path(outfile).stem+'_A.pdb', "w")
    with tempfile.TemporaryFile(mode='w+t') as temp:
        temp.writelines(readout)
        temp.seek(0)
        p1 = subprocess.Popen(['pdb_selchain', '-A'], stdin=temp, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['pdb_delelem','-H'], stdin=p1.stdout, stdout=ligfile)
    ligfile.close()
    temp.close()
    readout.close()

    readout = open(Path(pdb).stem+'_init.H.pdb')
    initfile = open(Path(pdb).stem+'_init.pdb', "w")
    with tempfile.TemporaryFile(mode='w+t') as temp:
        temp.writelines(readout)
        temp.seek(0)
        p2 = subprocess.Popen(['pdb_delelem','-H'], stdin=temp, stdout=initfile)
    initfile.close()
    temp.close()

    return(abdf)

if __name__ == '__main__':
    main('/home/drewaight/hdd1/proQ3drewdock_big/4ydl/4ydl_1.pdb')