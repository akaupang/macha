#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated in February 2023

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""
import argparse
from macha.functions import checkInput, Preparation
from macha.charmm_factory import CharmmManipulation

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
    print(f"Processing {ligand_id}")

    for env in envs:
        # Announce the current environment
        print(f"_________________________________________________________________________________ {env.upper():>8}")
        print(f"{str(' '):>90}")


        ########################################################################
        # Preparation of ligands
        # Instantiate class
        preparation = Preparation(
            parent_dir=parent_dir,
            original_dir=original_dir,
            ligand_id=ligand_id,
            env=env,
            cgenff_path=cgenff_path,
            protein_id=protein_id,
            rna=False,
            raw_pdb_to_mol2=False,          # There is use for this as OpenBabel fails (silently) to treat certain structures e.g. triazoles, and creates a different molecule
            system_ph=7.4,                  # If a ligand is missing all hydrogens, it will be protonated by OpenBabel according to this pH
            ligand_input_sanitation=True,   # if False, ligands will not be checked for hydrogens, repeated atom names/numbers and double-uppercase element names
            relax_input_segid_filter=False, # if False, only PROA and HETA will be used
            #include_xray_water=False,       # Whether to include x-ray water molecules
        )
        
        # Create CHARMM Coordinate files
        used_segids = preparation.createCRDfiles()
        print(
            f"Input parsing finished. Segids {', '.join(used_segids)} will be used "\
            f"for CRD generation and related tasks"
        )

        ########################################################################
        
        # Edit the CHARMM-GUI scripts
        # Instantiate class
        charmmManipulation = CharmmManipulation(
            parent_dir=parent_dir,
            ligand_id=ligand_id,
            original_dir=original_dir,
            het_resnames=preparation.het_resnames,
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
        
