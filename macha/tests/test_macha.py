from random import randint
import pytest
from macha.functions import Preparation, checkInput, CharmmManipulation
import os
import parmed as pm

ligand_id = "cdk2_l32"
# for github
parent_dir = "."
original_dir = "macha/data/original"
# for local
#parent_dir = "../data"
#original_dir = "original"


def test_import():
    import parmed as pm
    import openbabel
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

    preparation = Preparation(
        parent_dir=parent_dir,
        ligand_id=ligand_id,
        original_dir=original_dir,
        env="complex",
        cgenff_path="",
        ligand_input_sanitation=True
    )

    pdb_macha = preparation.pm_objs[0]
    pdb_macha_df = pdb_macha.to_dataframe()
    #print(f"We are looking at segid {pdb_macha_df.segid.unique()}, chain {pdb_macha_df.chain.unique()}")
    pdb_orig = pm.load_file(f"{parent_dir}/{original_dir}/cdk2_l32.pdb")
 
    # check if pdb file is read in correctly by parmed and compare a random atom to the original file
    number = randint(0, 4818)#the length of chain A
    assert str(pdb_macha.atoms[number]) == str(pdb_orig.atoms[number])



@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_macha_for_rna():

    ligand_id = "ino5"
    parent_dir = "macha/data"
    original_dir = "original"
    input_ext = "pdb"  # for testing - should be pdb
    cgenff_path = "/site/raid2/johannes/programs/silcsbio/silcsbio.2022.1/cgenff/cgenff"

    for env in ["single_strand", "double_strand"]:
        preparation = Preparation(
            parent_dir=parent_dir,
            ligand_id=ligand_id,
            original_dir=original_dir,
            env=env,
            rna=True,
        )
        segids, df = preparation.checkInputType()
        print(segids)
        # # Make a Transformato style folder structure below a folder bearing
        # # the name of the ligand
        preparation.makeTFFolderStructure()
        segids, used_segids = preparation.createCRDfiles(segids, df)

        print(segids, used_segids)
        # # # # Get the toppar stream from a local CGenFF binary
        # preparation.getTopparFromLocalCGenFF(cgenff_path=cgenff_path)

        charmmManipulation = CharmmManipulation(
            parent_dir=parent_dir,
            ligand_id=ligand_id,
            original_dir=original_dir,
            resname=preparation.resname,
            env=env,
            ion_name="SOD",
            ion_conc=1.0,
        )
        # # # Copy Files from the template folder
        charmmManipulation.copyFiles()
        # Modify step1_pdbreader.inp to read in correct amount of chains/residues
        charmmManipulation.modifyStep1(used_segids)
        # # Run Charmm giving the correct executable path
        charmmManipulation.executeCHARMM(charmm_exe="charmm")
        # # charmmManipulation.applyHMR()
        charmmManipulation.createOpenMMSystem()
        # charmmManipulation.createTFYamlFile()
        
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_macha_for_multiple_proteins_and_ligands():
    # THIS TEST CAN BE PERFORMED ON PDB ID: 2XYJ
    # AND WILL PRODUCE THE PROTEIN DIMER AND ITS FOUR(!) LIGANDS
    # THE LARGE STRUCTURE WILL TAKE A VERY LONG TIME TO FINISH
    # note: relax_input_segid_filter=True

    parent_dir = "/home/manny/Documents/Work/UiO/Modeling/wien/proteins/ppar/rbfe/test/data"
    original_dir = "original"
    protein_name = "protein" # -> protein.pdb
    cgenff_path = "/home/manny/Documents/Work/UiO/Modeling/wien/ligands/silcsbio/silcsbio.2022.1/cgenff/cgenff"

    # Introductions are in order
    print(f"{str(' '):>90}")
    print(f"=======================================================================  m   a   c   h   a")
    print(f"{str(' '):>90}")

    # Check for protein.pdb and if none is found, assume complexes are provided.
    # If protein.pdb is found, assume ligands alone.
    protein_id, ligand_ids = checkInput(
        parent_dir = parent_dir, 
        original_dir = original_dir, 
        protein_name = protein_name
        )

    envs = ["waterbox",
            "complex"
            ]

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
                relax_input_segid_filter=True,  # if False, only PROA and HETA will be used
                #include_xray_water=False,        # Whether to include x-ray water molecules
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

            ########################################################################
