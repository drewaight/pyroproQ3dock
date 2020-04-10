#!/bin/bash

pdb=$1
pdb_selchain -A,H,L $1 | pdb_tidy | pdb_reres | pdb_delinsertion | pdb_chain -A  > "${1%.*}_reres.pdb"
if [ $1 == "complex_0_cho.pdb" ]
then
    sed -i '/^\(ATOM\|TER\)/!d' "${1%.*}_reres.pdb"
    pdb_tofasta "${1%.*}_reres.pdb" > "complex.fasta" 
fi
