#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:03:58 2022

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""

from functions import check_ligands, Preparation, CharmmManipulation

################################################################################
### VARIABLES/SETTINGS
################################################################################

parent_dir = "."
original_dir = "original"
ligands_dir = "../ligands"
input_ext = "pdb"  # for testing - should be pdb
cgenff_path = "/site/raid2/johannes/programs/silcsbio/silcsbio.2022.1/cgenff/cgenff"

################################################################################
# MAIN (WORK)
################################################################################

# Check for ligand as argument, and if found, run in single ligand mode
ligand_ids = check_ligands(parent_dir = parent_dir, original_dir = original_dir, ligands_dir = ligands_dir)


# ITERATE THROUGH LIGANDS
for ligand_id in ligand_ids:
    
    print(f"\n")
    print(f"Processing ligand {ligand_id}")

    for env in ["waterbox", "complex"]:
        # Preparation of ligands
        # Instantiate class
        preparation = Preparation(
            parent_dir=parent_dir,
            ligand_id=ligand_id,
            original_dir=original_dir,
            env=env,
        )

        # Make a Transformato style folder structure
        preparation.makeTFFolderStructure()
        # Create CHARMM Coordinate files
        segids, used_segids = preparation.createCRDfiles()
        print("The following segment IDs were found/created:")
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
        charmmManipulation.applyHMR()
