#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated in February 2023

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""
import argparse
from functions import checkInput, Preparation, CharmmManipulation
#import logging

################################################################################
### VARIABLES/SETTINGS
################################################################################

parent_dir = "data"
original_dir = "original"
input_ext = "pdb"         # exclusive support for PDB files for now
protein_name = "protein"  # -> protein.pdb
cgenff_path = ""          # MUST BE SET BY USER

# ################################################################################
# ### LOGGING
# ################################################################################
# # create logger
# logger = logging.getLogger('macha')
# logger.setLevel(logging.DEBUG)

# # create console handler and set level to debug
# ch = logging.StreamHandler()
# ch.setLevel(logging.DEBUG)

# # create formatter
# #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# formatter = logging.Formatter(%(levelname)s:%(message)s')

# # add formatter to ch
# ch.setFormatter(formatter)

# # add ch to logger
# logger.addHandler(ch)

################################################################################
# MAIN (WORK)
################################################################################
# Introductions are in order
print(f"{str(' '):>90}")
print(f"======================================================================   m   a   c   h   a")
print(f"{str(' '):>90}")

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
    # Pretty printing
    print(f"__________________________________________________________________________________________")
    print(f"{str(' '):>90}")
    print(f"Processing ligand {ligand_id}")

    for env in envs:
        # Announce the current environment
        print(f"_________________________________________________________________________________ {env.upper():>8}")
        print(f"{str(' '):>90}")

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
            ligand_input_sanitation=True,   # If ligand input sanitation is selected and
            system_ph=7.4,                  # the ligand has no hydrogens, it will be
                                            # protonated by OpenBabel according to this pH
        )
        
        # THESE LINES COULD BE MOVED INTO THE CLASS AND CALLED IN ___INIT___ METHOD
        # Make a Transformato style folder structure
        #preparation.makeTFFolderStructure() # ALREADY MOVED
        
        segids = set(i.residue.segid for i in preparation.pm_obj)
        pm_obj_df = preparation.pm_obj.to_dataframe()

        # Create CHARMM Coordinate files
        segids, used_segids = preparation.createCRDfiles(segids, pm_obj_df)
        print(f"The following segment IDs were found/assigned:    " + ", ".join(segids)) 
        print(f"The following segment IDs were used/not excluded: " + ", ".join(used_segids))

        # Get the toppar stream from a local CGenFF binary
        preparation.getTopparFromLocalCGenFF(cgenff_path=cgenff_path)
        ############################################################################

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
        
        # Space to pause and consider
        print(f"{str(' '):>90}")
        
