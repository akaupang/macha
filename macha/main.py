#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated in February 2023

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""
import argparse
from functions import checkInput, Preparation, CharmmManipulation

################################################################################
### VARIABLES/SETTINGS
################################################################################

parent_dir = "data"
original_dir = "original"
input_ext = "pdb"         # exclusive support for PDB files for now
protein_name = "protein"  # -> protein.pdb
cgenff_path = ""          # MUST BE SET BY USER

################################################################################
# MAIN (WORK)
################################################################################

# Check for protein.pdb and if none is found, assume complexes are provided.
# If protein.pdb is found, assume ligands alone.
protein_id, ligand_ids = checkInput(
    parent_dir = parent_dir, 
    original_dir = original_dir, 
    protein_name = protein_name
    )

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

    print(f"---------------------------------------------------------------\n")
    print(f"Processing ligand {ligand_id}")

    for env in envs:
        # Preparation of ligands
        # Instantiate class
        preparation = Preparation(
            parent_dir=parent_dir,
            original_dir=original_dir,
            ligand_id=ligand_id,
            env=env,
            protein_id=protein_id,
            small_molecule=False,
            rna=False,
        )
        
        # Make a Transformato style folder structure
        preparation.makeTFFolderStructure()
        
        # Check input types
        segids, pm_obj_df = preparation.checkInputType()

        # Create CHARMM Coordinate files
        segids, used_segids = preparation.createCRDfiles(segids, pm_obj_df)
        print(f"The following segment IDs were found/assigned:    " + " ,".join(segids)) 
        print(f"The following segment IDs were used/not excluded: " + " ,".join(used_segids))

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
        
        # Convert the CHARMM system to an OpenMM system
        charmmManipulation.createOpenMMSystem()
        
        # Apply HMR to the OpenMM system
        charmmManipulation.applyHMR()
        
        # Create YAML file for ASFE simulations using TF (requires ./parent/config directory)
        charmmManipulation.createTFYamlFile(dt=0.002, nstep=2500000)
