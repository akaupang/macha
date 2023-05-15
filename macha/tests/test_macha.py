from random import randint
import pytest
from macha.functions import Preparation, checkInput
from macha.charmm_factory import CharmmManipulation
import os
import parmed as pm

# Paths for GitHub tests
parent_dir = "macha/data"
original_dir = "original"

# # Imports/Paths for local tests
# parent_dir = "data"
# original_dir = "original/external_protein"
# cgenff_path = "/home/manny/Documents/Work/UiO/Modeling/wien/ligands/silcsbio/silcsbio.2022.1/cgenff/cgenff"
# protein_name = "protein"

def test_import():
    import parmed as pm
    import openbabel
    import natsort
    import re

def test_checkInput():
    protein_id, ligand_ids = checkInput(
        parent_dir=parent_dir,
        original_dir=original_dir,
    )

    # Assure that all ligands are found by the function
    if not protein_id:
        assert len(os.listdir(f"{parent_dir}/{original_dir}")) == len(ligand_ids)

def test_createFolders():
    ligand_id = "cdk2_l32"

    preparation = Preparation(
        parent_dir=parent_dir,
        ligand_id=ligand_id,
        original_dir=original_dir,
        env="complex",
        cgenff_path="",
        ligand_input_sanitation=True
    )

    pdb_macha = preparation.pm_objs[0]
    pdb_orig = pm.load_file(f"{parent_dir}/{original_dir}/cdk2_l32.pdb")
 
    # Check if PDB file is read in correctly by ParmEd and compare a random atom to the original file
    number = randint(0, 4818)#the length of chain A
    assert str(pdb_macha.atoms[number]) == str(pdb_orig.atoms[number])

@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_macha_external_protein():
    envs = ['waterbox', 'complex']
    
    protein_id, ligand_ids = checkInput(
        parent_dir=parent_dir,
        original_dir=original_dir,
        protein_name = protein_name, # IMPORTANT FOR EXTERNAL PROTEIN
    )
    
    for ligand_id in ligand_ids:
        print(f"Working on {ligand_id}")
        for env in envs:
            preparation = Preparation(
                parent_dir=parent_dir,
                original_dir=original_dir,
                ligand_id=ligand_id,
                env=env,
                cgenff_path=cgenff_path,
                protein_id=protein_id,
                rna=False,
                raw_pdb_to_mol2=False,
                system_ph=7.4,
                ligand_input_sanitation=True,
                relax_input_segid_filter=False,
            )
            
            used_segids = preparation.createCRDfiles()

            charmmManipulation = CharmmManipulation(
                parent_dir=parent_dir,
                ligand_id=ligand_id,
                original_dir=original_dir,
                het_resnames=preparation.het_resnames,
                env=env,
            )
            
            charmmManipulation.copyFiles()
            charmmManipulation.modifyStep1(used_segids)
            charmmManipulation.executeCHARMM(charmm_exe="charmm")
            charmmManipulation.createOpenMMSystem()
            charmmManipulation.applyHMR()
            charmmManipulation.createTFYamlFile(dt=0.002, nstep=2500000)


@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_macha_for_rna():
    original_dir = "original/rna"
    envs = ["waterbox", "complex"]
    cgenff_path = "" # Not needed for RNA

    protein_id, ligand_ids = checkInput(
        parent_dir=parent_dir,
        original_dir=original_dir,
    )
    
    for ligand_id in ligand_ids:
        for env in envs:
            preparation = Preparation(
                parent_dir=parent_dir,
                original_dir=original_dir,
                ligand_id=ligand_id,
                env=env,
                cgenff_path=cgenff_path,
                protein_id=protein_id,
                rna=True,
                raw_pdb_to_mol2=False,
                system_ph=7.4,
                ligand_input_sanitation=True,
                relax_input_segid_filter=False,
            )
            
            used_segids = preparation.createCRDfiles()

            charmmManipulation = CharmmManipulation(
                parent_dir=parent_dir,
                ligand_id=ligand_id,
                original_dir=original_dir,
                het_resnames=preparation.het_resnames,
                env=env,
            )
            
            charmmManipulation.copyFiles()
            charmmManipulation.modifyStep1(used_segids)
            charmmManipulation.executeCHARMM(charmm_exe="charmm")
            charmmManipulation.createOpenMMSystem()
            charmmManipulation.applyHMR()
            charmmManipulation.createTFYamlFile(dt=0.002, nstep=2500000)
        
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_macha_for_multiple_proteins_and_ligands():
    # THIS TEST CAN BE PERFORMED E.G. ON PDB ID: 2XYJ
    # AND WILL PRODUCE THE PROTEIN DIMER AND ITS FOUR(!) LIGANDS
    # THE LARGE STRUCTURE WILL TAKE A VERY LONG TIME TO FINISH
    # note: relax_input_segid_filter=True
    original_dir = "original/complex"
    envs = ["waterbox", "complex"]
   
    protein_id, ligand_ids = checkInput(
        parent_dir = parent_dir, 
        original_dir = original_dir, 
        protein_name = protein_name, # IMPORTANT FOR EXTERNAL PROTEIN
        )

    for ligand_id in ligand_ids:
        for env in envs:
            preparation = Preparation(
                parent_dir=parent_dir,
                original_dir=original_dir,
                ligand_id=ligand_id,
                env=env,
                cgenff_path=cgenff_path,
                protein_id=protein_id,
                rna=False,
                raw_pdb_to_mol2=False,
                system_ph=7.4,
                ligand_input_sanitation=True,
                relax_input_segid_filter=True,
            )
            
            used_segids = preparation.createCRDfiles()

            charmmManipulation = CharmmManipulation(
                parent_dir=parent_dir,
                ligand_id=ligand_id,
                original_dir=original_dir,
                het_resnames=preparation.het_resnames,
                env=env,
            )
            
            charmmManipulation.copyFiles()
            charmmManipulation.modifyStep1(used_segids)
            charmmManipulation.executeCHARMM(charmm_exe="charmm")
            charmmManipulation.createOpenMMSystem()
            charmmManipulation.applyHMR()
            charmmManipulation.createTFYamlFile(dt=0.002, nstep=2500000)
