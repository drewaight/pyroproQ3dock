from pyrosetta import *
from pyrosetta.rosetta import *
from rosetta.protocols.rosetta_scripts import *
from rosetta.protocols.antibody import *
from rosetta.protocols.antibody.design import *
from pyrosetta.rosetta.core.scoring import *
import pandas_anarci_numberer
from rosetta.utility import *
import os
import sys
import timeit
import numpy as np
import pandas as pd
import re
from pathlib import Path


def main(path, structure, light_chain):
    aho_structure = Path(structure).stem + '_aho.pdb'
    pandas_anarci_numberer.main(path +'/'+ structure, 'aho', path +'/'+ aho_structure)
    init('-input_ab_scheme AHo -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -allow_omega_mismatches_for_north_clusters')
    pose = pose_from_pdb(path +'/'+ aho_structure)
    original_pose = pose.clone()
    scorefxn = get_score_function()
    start = timeit.default_timer()
    min_packer = """
        <ROSETTASCRIPTS>
            <SCOREFXNS>
                <ScoreFunction name="r15" weights="ref2015" />
                <ScoreFunction name="molmech" weights="mm_std_fa_elec_dslf_fa13" />
                <ScoreFunction name="r15_cart" weights="ref2015" >
                <Reweight scoretype="pro_close" weight="0.0" />
                <Reweight scoretype="cart_bonded" weight="0.625" />
                </ScoreFunction>
            </SCOREFXNS>
            <RESIDUE_SELECTORS>
            </RESIDUE_SELECTORS>
            <TASKOPERATIONS>
                <RestrictToRepacking name="no_design" /> 
                <ExtraRotamersGeneric name="extrachi" ex1="1" ex2="1" ex1_sample_level="1" ex2_sample_level="1" /> 
            </TASKOPERATIONS>
            <FILTERS>
            </FILTERS>
            <MOVERS>
                <MinMover name="min_torsion" scorefxn="molmech" chi="true" bb="1" cartesian="F" />
                <MinMover name="min_cart" scorefxn="r15_cart" chi="true" bb="1" cartesian="T" />
                <PackRotamersMover name="pack1" scorefxn="r15" task_operations="no_design,extrachi" />
            </MOVERS>
            <APPLY_TO_POSE>
            </APPLY_TO_POSE>
            <PROTOCOLS>
                <Add mover="pack1" />
                <Add mover="min_torsion" />
                <Add mover="min_cart" />
            </PROTOCOLS>
            <OUTPUT scorefxn="r15" />
        </ROSETTASCRIPTS>
        """
    xml = XmlObjects.create_from_string(min_packer)
    min_pack_protocol = xml.get_mover("ParsedProtocol")
    min_pack_protocol.apply(pose)
    
    abinfo = AntibodyInfo(original_pose, AHO_Scheme, North)
    if light_chain == "kappa":
        abinfo.set_light_chain_type(kappa)
    elif light_chain == "lambda":
        abinfo.set_light_chain_type(LightChainTypeEnum_start)
    print(abinfo)    

    pdbinfo = pose.pdb_info()
    chains = ("H","L")
    for res in range(1, pose.total_residue()+1):
        for atom in range(1, pdbinfo.natoms(res)+1):
            pdbinfo.temperature(res,atom, 0.0)
                      
    for num in range(1,3):        
        for res in range(pose.conformation().chain_begin(num), pose.conformation().chain_begin(num)+4):
            for atom in range(1, pdbinfo.natoms(res)+1):
                pdbinfo.temperature(res,atom, 1.0)
        
    for chain in chains:
        for res in range(pdbinfo.pdb2pose(chain, 73), pdbinfo.pdb2pose(chain, 93)):
            for atom in range(1, pdbinfo.natoms(res)+1):
                pdbinfo.temperature(res,atom, 1.0)
    
    for i in range(1, 7):
        CDR_range = abinfo.get_CDR_start(CDRNameEnum(i), pose), abinfo.get_CDR_end(CDRNameEnum(i), pose)
        for r in range(CDR_range[0]-3, CDR_range[1]+4):
            for atom in range(1, pdbinfo.natoms(r)+1):
                pdbinfo.temperature(r,atom, 1.0)

    pose.dump_pdb(path +'/'+ Path(structure).stem + '_min.pdb')

    pandas_anarci_numberer.main(path +'/'+ Path(structure).stem + '_min.pdb', 'chothia', path +'/'+ structure)
    print(scorefxn(original_pose))
    print(scorefxn(pose),CA_rmsd(pose, original_pose))
    
    list = []
    for i in range(1, 7):
        list.append(abinfo.get_cluster_name(abinfo.get_CDR_cluster(CDRNameEnum(i)).cluster()))
        list.append(abinfo.get_CDR_cluster(CDRNameEnum(i)).normalized_distance_in_degrees())
        list.append(abinfo.get_CDR_sequence_with_stem(CDRNameEnum(i), pose, 3, 3))
    abdf_T = pd.DataFrame(list)
    abdf = pd.DataFrame(abdf_T.T)
    abdf.columns=['H1_cluster', 'H1_distance', 'H1_sequence', 'H2_cluster', 'H2_distance', 'H2_sequence',
                'H3_cluster', 'H3_distance', 'H3_sequence', 'L1_cluster', 'L1_distance', 'L1_sequence',
                'L2_cluster', 'L2_distance', 'L2_sequence', 'L3_cluster', 'L3_distance', 'L3_sequence']
    
    end = timeit.default_timer()
    print(end-start)
    print(abdf)
    return(abdf)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])