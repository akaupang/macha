#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:50:05 2022

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""
import sys
import os
import glob
import shutil
import subprocess
import parmed as pm
import pandas as pd
from openbabel import openbabel

################################################################################
### FUNCTIONS
################################################################################
# HET BLOCK GENERATION
################################################################################
def genPROblock(protein_id, segment_id="PROA"):

    block = (
        f"!-------------------------------------------------------------------------------\n"
        f"! PROTEIN:\n"
        f"open read card unit 10 name {protein_id.lower()}_{segment_id.lower()}.crd\n"
        f"read sequence coor card unit 10 resid\n"
        f"generate {segment_id.upper()} setup warn first NTER last CTER\n"
        f"\n"
        f"open read unit 10 card name {protein_id.lower()}_{segment_id.lower()}.crd\n"
        f"read coor unit 10 card resid\n"
        f"!-------------------------------------------------------------------------------\n"
    )

    return block


def genHETblock(ligand_id, segment_id="HETA"):

    block = (
        f"!-------------------------------------------------------------------------------\n"
        f"open read card unit 10 name {ligand_id.lower()}.crd\n"
        f"read sequence coor card unit 10 resid\n"
        f"generate {segment_id.upper()} setup warn first none last none\n"
        f"\n"
        f"open read unit 10 card name {ligand_id.lower()}.crd\n"
        f"read coor unit 10 card resid\n"
        f"!-------------------------------------------------------------------------------\n"
    )

    return block


def genWATblock(protein_id, segment_id="WATA"):

    block = (
        f"!-------------------------------------------------------------------------------\n"
        f"open read card unit 10 name {protein_id.lower()}_{segment_id.lower()}.crd\n"
        f"read sequence coor card unit 10 resid\n"
        f"generate {segment_id.upper()} setup warn noangle nodihedral\n"
        f"\n"
        f"open read unit 10 card name {protein_id.lower()}_{segment_id.lower()}.crd\n"
        f"read coor unit 10 card resid\n"
        f"!-------------------------------------------------------------------------------\n"
    )

    return block


################################################################################
# HBUILD CONTROL
################################################################################
def hbuild_preserve_explicit_H(segids):

    segid_head = f"hbuild sele hydr .and. .not. ("
    segid_vars = ""
    for idx, segid in enumerate(segids):
        if idx == len(segids) - 1:
            segid_add = f"segid {segid}"
        else:
            segid_add = f"segid {segid} .or. "
        segid_vars = segid_vars + segid_add
    segid_tail = f") end"

    block = (
        f"!-------------------------------------------------------------------------------\n"
        f"! MOD: HBUILD control - preserve explicit H-coordinates\n"
        f"prnlev 5\n"
        f"echo START_HBUILD\n"
        f"{segid_head + segid_vars + segid_tail}\n"
        f"echo END_HBUILD\n"
        f"!-------------------------------------------------------------------------------\n"
    )
    return block


################################################################################
# PRINT USED PARAMETERS
################################################################################
def insert_parameter_print():
    block = (
        f"!-------------------------------------------------------------------------------\n"
        f"! MOD: Print used parameters\n"
        f"echo START_PAR\n"
        f"print para used\n"
        f"echo END_PAR\n"
        f"!-------------------------------------------------------------------------------\n"
    )
    return block


################################################################################
# OTHER FUNCTIONS
################################################################################


def streamWriter(work_dir, stream_name, blocks_as_lst):
    with open(f"{work_dir}/{stream_name}.str", "w") as ostr:
        for block in blocks_as_lst:
            ostr.write(block)


def inputFileInserter(inpfile, cases, blocks, inversions):
    assert len(cases) == len(blocks)

    # Check for the presence of a backup, and start from this
    # or make a backup if none exists
    if os.path.exists(f"{inpfile}.premod"):
        shutil.copy(f"{inpfile}.premod", f"{inpfile}")
    else:
        shutil.copy(f"{inpfile}", f"{inpfile}.premod")

    inp_backup = f"{inpfile}.premod"
    outfile = inpfile

    with open(outfile, "w") as out:
        with open(inp_backup, "r") as inf:
            for line in inf:
                for case, block, inversion in zip(cases, blocks, inversions):
                    if case in line:
                        if inversion == True:
                            line = f"{block}\n" f"\n" f"{case}\n"
                        elif inversion == None:
                            line = f"{block}\n"
                        else:
                            line = f"{case}\n" f"\n" f"{block}\n"
                out.write(line)


class Preparation:

    def __init__(self, parent_dir, ligand_id, original_dir):
        self.parent_dir = parent_dir
        self.ligand_id = ligand_id
        self.original_dir = original_dir
        self.resname = str

    def makeFolder(self, path):

        try:
            os.makedirs(path)
            print(f"Creating folder in {path}")
        except OSError:
            if not os.path.isdir(path):
                raise

    def makeTFFolderStructure(self):

        self.makeFolder(f"{self.parent_dir}/{self.ligand_id}")
        self.makeFolder(f"{self.parent_dir}/{self.ligand_id}/complex")
        self.makeFolder(f"{self.parent_dir}/{self.ligand_id}/waterbox")

        print(
            f"Creating folders in {self.parent_dir}/{self.ligand_id} for complex and waterbox"
        )


    def _create_mol2_file(self):

        print(f"Converting the residue pdb file to a mol2 file")

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "mol2")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.pdb")
        #mol.AddHydrogens() TODO: Do we need this?

        assert (mol.NumResidues()) == 1
        obConversion.WriteFile(mol,f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.mol2")
        
    def _modify_resname_in_str(self):
        
        fin = open(f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.str", "rt")
        fout = open(f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}_tmp.str", "wt")
        for line in fin:
            if line.startswith("RESI"):
                fout.write(line.replace(line.split()[1], 'UNK'))
            else:
                fout.write(line)

        fin.close()
        fout.close()

        shutil.copy(fout.name,fin.name)
    
    def copyREStocomplex(self):

        for file in glob.glob(f"{self.ligand_id}/waterbox/{self.resname.upper()}/*"):
            shutil.copy(file, f"{self.ligand_id}/complex/{self.resname.upper()}/")


    def getTopparFromLocalCGenFF(
        self,
        cgenff_path=False,
    ):
        cgenff_bin = None
        cgenff_output = None

        # CGenFF needs a mol2 file as input file
        self._create_mol2_file()

        # If no particular path is given, check whether CGenFF is available
        if cgenff_path == False:
            cgenff_path = shutil.which("cgenff")
            if cgenff_path == None:
                print("This function requires cgenff.")
                print(
                    "Please install it in the active environment or point the routine"
                )
                print("to the right path using the key cgenff_path='/path/to/cgenff' .")
            else:
                cgenff_bin = cgenff_path
        else:
            cgenff_bin = cgenff_path

        # CGenFF exists - start program
        if cgenff_bin != None:

            # Run CGenFF
            ligand_path = f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.mol2"

            cgenff_output = subprocess.run(  #
                [cgenff_bin]
                + [ligand_path]
                + ["-v"]
                + ["-f"]
                + [
                    f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.str"
                ]
                + ["-m"]
                + [
                    f"{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.log"
                ],
                # [cgenff_bin] + [ligand_path] + ["-v"] + ["-m"] + [f"{ligand_id}.log"],#
                capture_output=True,
                text=True,  #
            )

            # Evaluate the subprocess return code
            if cgenff_output.returncode == 1:
                print(f"CGenFF returned an error after being called with:\n")
                print(" ".join(cgenff_output.args))
                print(cgenff_output.stdout)
                print(cgenff_output.stderr)
            else:
                print(f"CGenFF executed successfully")
                # print(cgenff_output.stdout)
                # print(cgenff_output.stderr)
        
        self._modify_resname_in_str()

        return cgenff_output

    def createCRDfiles(self):

        pdb_file = pm.load_file(f"{self.original_dir}/{self.ligand_id}.pdb")
        df = pdb_file.to_dataframe()
        segids = set(i.residue.chain for i in pdb_file)
        if (
            len(segids) > 1
        ):  # for pdb file created with MAESTRO containing the chaing segment
            print(f"Processing a Maestro based pdb file")
            aa = [
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLN",
                "GLU",
                "GLY",
                "HIS",
                "ILE",
                "LEU",
                "LYS",
                "MET",
                "PHE",
                "PRO",
                "SER",
                "THR",
                "TPO",
                "TRP",
                "TYR",
                "VAL",
                "HSD",
                "HSE",
            ]  # we need to finde the residue of the ligand which should be the only one beeing not an aa
            for (
                chain
            ) in (
                segids
            ):  # rename the chain names (a,b, ...) to segnames (proa,prob,...)
                for i in pdb_file.view[df.chain == f"{chain}"]:
                    if i.residue.name not in aa:
                        i.residue.chain = f"HETA"
                        self.resname = i.residue.name
                    else:
                        if i.residue.name == "HIS":
                            i.residue.name = "HSD" # ATTENTION!! Here we make all HIS to HSD
                            i.residue.chain = f"PRO{chain}"
                        i.residue.chain = f"PRO{chain}"

            self.makeFolder(
                f"{self.parent_dir}/{self.ligand_id}/complex/{self.resname.upper()}"
            )
            self.makeFolder(
                f"{self.parent_dir}/{self.ligand_id}/waterbox/{self.resname.upper()}"
            )

            df = pdb_file.to_dataframe()
            segids = set(i.residue.chain for i in pdb_file)
            for segid in segids:  # now we can save the crd files
                if segid not in ["SOLV", "IONS"]:
                    if segid == "HETA":
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/waterbox/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.pdb",
                            overwrite=True,
                        )
                    else:
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.pdb",
                            overwrite=True,
                        )

        else:  # CHARMM-GUI generated pdb files
            print(f"Processing a CHARMM-GUI based pdb file")
            segids = set(i.residue.segid for i in pdb_file)
            for segid in segids:
                if segid not in ["SOLV", "IONS"]:
                    if segid == "HETA":
                        self.resname = pdb_file[df.segid == f"{segid}"].residues[0].name
                        self.makeFolder(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{self.resname.upper()}"
                        )
                        self.makeFolder(
                            f"{self.parent_dir}/{self.ligand_id}/waterbox/{self.resname.upper()}"
                        )
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/waterbox/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/waterbox/{self.resname.upper()}/{self.resname.lower()}.pdb",
                            overwrite=True,
                        )

                    else:
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.pdb",
                            overwrite=True,
                        )
        
        # resname should be a 3 or 4 letters code                
        assert len(self.resname) < 5


class CharmmManipulation():
    
    def __init__(self, parent_dir, ligand_id, original_dir):

        self.parent_dir = parent_dir
        self.ligand_id = ligand_id
        self.original_dir = original_dir
        self.default_path = f"{self.parent_dir}/templates/default/"
        self.resname: str = None

    def manipulateToppar(self, resname):
        
        self.resname = resname
        shutil.copy(f"{self.default_path}/toppar.str", f"{self.ligand_id}/waterbox/toppar.str")
        shutil.copy(f"{self.default_path}/toppar.str", f"{self.ligand_id}/complex/toppar.str")
        try:
            shutil.copytree(f"{self.default_path}/toppar", f"{self.ligand_id}/complex/toppar")
            shutil.copytree(f"{self.default_path}/toppar", f"{self.ligand_id}/waterbox/toppar")                  
        except:
            print(f"Toppar directory is already available")

        # manipulate toppar.str file
        file = open(f"{self.ligand_id}/waterbox/toppar.str","a")
        file.write(f"stream {self.resname.upper()}/{self.resname.lower()}.str")
        file.close()
        file = open(f"{self.ligand_id}/complex/toppar.str","a")
        file.write(f"stream {self.resname.upper()}/{self.resname.lower()}.str")
        file.close()

    def copyINPfiles(self):
        for file in glob.glob(f"{self.default_path}/*.inp"):
            shutil.copy(file, f"{self.ligand_id}/waterbox/")
            shutil.copy(file, f"{self.ligand_id}/complex/")        

    def prepareStep1(self):

        correction_on = False
        
        fout = open(f"{self.ligand_id}/waterbox/step1_pdbreader_tmp.inp", "wt")
        with open(f"{self.ligand_id}/waterbox/step1_pdbreader.inp", "r+") as f:
            for line in f:
                if line.startswith("! Read PROA"):
                    correction_on = True
                    f.readline()
                    fout.write('Hallihallo \n')
                    fout.write('und das hier noch \n')
                    fout.write('und das \n')

                if line.startswith("!Print heavy atoms with "):
                    fout.write(f"{line}")
                    correction_on = False
                    
                if correction_on == True:
                    pass
                else:
                    fout.write(line)


