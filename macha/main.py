#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:03:58 2022

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""
import argparse
from functions import checkInput, Preparation, CharmmManipulation

################################################################################
### VARIABLES/SETTINGS
################################################################################

parent_dir = "."
original_dir = "data/original"
input_ext = "pdb"  # exclusive support for PDB
protein_name = "protein" # -> protein.pdb
cgenff_path = "/site/raid2/johannes/programs/silcsbio/silcsbio.2022.1/cgenff/cgenff"

################################################################################
# MAIN (WORK)
################################################################################

# Check for protein.pdb and if none is found, assume complexes are provided.
# If protein.pdb is found, assume ligands alone.
protein_id, ligand_ids = checkInput(parent_dir = parent_dir, original_dir = original_dir, protein_name = protein_name)

# ARGUMENT HANDLING
# Initialize
parser = argparse.ArgumentParser(description='Options')
# Add expected/possible arguments
parser.add_argument('-nc', '--nocomplex', action='store_true', help="Toggle production of waterboxes only (no complexes)")
# Parse arguments
args = parser.parse_args()

# Determine which environments/run types will be used
if args.nocomplex == True:
    no_complex = True
    envs = ["waterbox"]
    print("Attention: Only waterboxes will be produced!")
else:
    no_complex = False
    envs = ["waterbox", "complex"]

# ITERATE THROUGH LIGANDS
for ligand_id in ligand_ids:

    print("\n")
    print(f"Processing ligand {ligand_id}")

    for env in envs:
        # Preparation of ligands
        # Instantiate class
        preparation = Preparation(
            parent_dir=parent_dir,
            original_dir=original_dir,
            ligand_id=ligand_id,
            protein_id=protein_id,
            env=env,
        )

        # Check input types
        segids, df = preparation.checkInputType()

        # Make a Transformato style folder structure
        preparation.makeTFFolderStructure()

        # Create CHARMM Coordinate files
        segids, used_segids = preparation.createCRDfiles(segids, df)
        print("The following segment IDs were found/assigned:")
        print(*segids)
        print("The following segment IDs were used/not excluded:")
        print(*used_segids)

        # Get the toppar stream from a local CGenFF binary
        preparation.getTopparFromLocalCGenFF(cgenff_path=cgenff_path)

        # Edit the CHARMM-GUI scripts
        # Instantiate class
        charmmManipulation = CharmmManipulation(
            parent_dir=parent_dir,
            ligand_id=ligand_id,
            original_dir=original_dir,
            resname=preparation.resname,
            env=env,
        )
        
        # Copy Files from the template folder
        charmmManipulation.copyFiles()

        # Modify step1_pdbreader.inp to read in correct amount of chains/residues
        charmmManipulation.modifyStep1(used_segids)

        # Run Charmm giving the correct executable path
        charmmManipulation.executeCHARMM(charmm_exe="charmm")
        charmmManipulation.createOpenMMSystem()
        charmmManipulation.applyHMR()
        charmmManipulation.createTFYamlFile(dt=0.002)