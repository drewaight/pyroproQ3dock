import os 
import sys
import numpy as np
import pandas as pd
import pandas_anarci_numberer
from pathlib import Path
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.utility import *
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.analysis import *
from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.antibody.design import *
from pyrosetta.rosetta.protocols.loops import *
from pyrosetta.rosetta.core.scoring import *

import timeit

def main(path, structure):

    init('-input_ab_scheme AHo -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false')
    # init('-input_ab_scheme AHo')

    pandas_anarci_numberer.main(path +'/'+ structure, 'aho', Path(structure).stem+'_aho.pdb')
    pose = pose_from_pdb(path +'/'+ Path(structure).stem+'_aho.pdb')
    original_pose = pose.clone()
    scorefxn = get_score_function()

    iam = InterfaceAnalyzerMover("LH_A")
    iam.set_pack_separated(True)
    iam.set_compute_packstat(False)
    iam.set_compute_interface_sc(True)
    iam.set_compute_separated_sasa(True)
    dg_term = "dG_separated"
    iam.apply(original_pose)

    abinfo = AntibodyInfo(pose, AHO_Scheme, North)
    abinfo.set_light_chain_type(kappa)

    ft = FoldTree()
    original_ft = pose.fold_tree()
    cdr_loops = Loops()
    for i in range (1,7):
        start = abinfo.get_CDR_start(CDRNameEnum(i), pose)
        stop =  abinfo.get_CDR_end(CDRNameEnum(i), pose)
        cutpoint = int((stop-start)/2+start)
        cdr_loop = Loop(start, stop, cutpoint)
        cdr_loops.add_loop(cdr_loop)
    fold_tree_from_loops(pose, cdr_loops, ft)
    pose.fold_tree(ft)
    add_cutpoint_variants(pose)

    FR_protocol = """
    <ROSETTASCRIPTS>
        <SCOREFXNS>
            <ScoreFunction name="r15_cart" weights="ref2015_cart">
                <Reweight scoretype="chainbreak" weight="100"/>
            </ScoreFunction>
            <ScoreFunction name="r15" weights="ref2015" />
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
            <AntibodyRegion name="cdrs" region="cdr_region" />
            <AntibodyRegion name="antigen" region="antigen_region" />
            <Neighborhood name="cdr_neighbor" selector="cdrs" distance="10.0" />
        </RESIDUE_SELECTORS>
        <MOVE_MAP_FACTORIES>
            <MoveMapFactory name="movemap_cdrs" bb="1" chi="1">
                <Backbone residue_selector="cdrs" />
                <Backbone residue_selector="cdr_neighbor" />
                <Chi residue_selector="cdrs" />
                <Chi residue_selector="cdr_neighbor" />
            </MoveMapFactory>
        </MOVE_MAP_FACTORIES>
        <TASKOPERATIONS>
            <RestrictToCDRsAndNeighbors name="restrict"/>
        </TASKOPERATIONS>
        <FILTERS>
        </FILTERS>
        <MOVERS>
            <FastRelax name="fast_relax" scorefxn="r15_cart" task_operations="restrict" 
            dualspace="1" cartesian="0" repeats="1" movemap_factory="movemap_cdrs"/>     
        </MOVERS>
        <PROTOCOLS>
            <Add mover_name="fast_relax" />
        </PROTOCOLS>
        <OUTPUT />
    </ROSETTASCRIPTS>
    """
    xml = XmlObjects.create_from_string(FR_protocol)
    protocol = xml.get_mover("ParsedProtocol")
    protocol.apply(pose)

    pose.fold_tree(original_ft)

    iam.apply(original_pose)
    iam.apply(pose)

    print("original", scorefxn(original_pose), original_pose.scores[dg_term])
    print("relaxed", scorefxn(pose), pose.scores[dg_term])
    pose.dump_pdb(path + '/' + Path(structure).stem+"_relaxed.pdb")

    score_dict = {}
    score_dict.update(pose.scores)
    df = pd.DataFrame(score_dict, index=[0])
    sc = np.array(df.to_numpy())
    return(sc)

if __name__ == '__main__':
    main(path = '/media/hdd1/proQ3drewdock_big/2xwt/megadock/' , structure='complex_0.pdb')