from random import randint
import pytest
from macha.functions import Preparation, checkInput, CharmmManipulation
import os
import parmed as pm

ligand_id = "cdk2_l32"
parent_dir = "macha/data"
original_dir = "original"


def test_import():
    import parmed as pm
    import openbabel


def test_checkInput():

    protein_id, ligand_ids = checkInput(
        parent_dir=parent_dir,
        original_dir=original_dir,
    )

    # Assure that all ligands are found by the function
    if not protein_id:
        assert len(os.listdir(f"macha/data/{original_dir}")) == len(ligand_ids)


def test_createFolders():

    preparation = Preparation(
        parent_dir=parent_dir,
        ligand_id=ligand_id,
        original_dir=original_dir,
        env="waterbox"
        ligand_input_sanitation=True
    )

    pdb_macha = preparation.pm_obj
    pdb_orig = pm.load_file("macha/data/original/cdk2_l32.pdb")

|   # I think this test will fail since the pm_obj is a now a dynamical entity which changes size 
    # quite a bit during the preparation run. The new statics are;
    # preparation.ligand_pm_pbj and preparation.protein_pm_obj 
 
    # check if pdb file is read in correctly by parmed and compare a random atom to the original file
    number = randint(0, 9065)
    assert str(pdb_macha.atoms[number]) == str(pdb_orig.atoms[number])


# Test for handling a small molecule


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
