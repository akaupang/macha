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
from openbabel import openbabel
# from charmm_factory import CharmmFactory

def check_ligands(parent_dir= ".", original_dir="original", ligands_dir = "ligands", input_ext="pdb"):

    if not os.path.isdir(original_dir):
        sys.exit(f"Input file {original_dir} not found!")
    
    try:
        ligand_id = sys.argv[1]
        # Check for existence of ligands/"ligandid"
        if os.path.exists(f"{parent_dir}/{ligands_dir}/{ligand_id}.{input_ext}"):
            pass
        else:
            sys.exit(
                f"Input file: {parent_dir}/{ligands_dir}/{ligand_id}.{input_ext} not found!"
            )

    # If no argument is given, run in multiple ligand mode
    except IndexError:
        ligand_id = None  # None in particular

    # Create the master list of the ligand/system names and fill it based on
    # the findings above
    ligand_ids = []
    if ligand_id == None:
        for ifile in glob.glob(f"{parent_dir}/{original_dir}/*.{input_ext}"):
            ligand_ids.append(os.path.splitext(os.path.basename(ifile))[0])  #
    elif ligand_id != "":
        ligand_ids.append(ligand_id)
    else:
        sys.exit(f"Unknown ligand (file) name error with name: {ligand_id}")

    return ligand_ids


class Preparation:
    def __init__(self, parent_dir, ligand_id, original_dir, env):
        """
        This class prepares everything for further use with CHARMM. The pdb files are sliced into pieces 
        and the ligand is converted to a mol2 file.
        A local version of CGenFF creates a stream file for the ligand.
        """
        self.parent_dir = parent_dir
        self.ligand_id = ligand_id
        self.original_dir = original_dir
        self.resname = str
        self.env: str = env

    def makeFolder(self, path):

        try:
            os.makedirs(path)
            print(f"Creating folder in {path} for the {self.env}")
        except OSError:
            print(f"There exists a folder in {path} for {self.env} we will use it")
            if not os.path.isdir(path):
                raise

    def makeTFFolderStructure(self):

        self.makeFolder(f"{self.parent_dir}/{self.ligand_id}")
        self.makeFolder(f"{self.parent_dir}/{self.ligand_id}/{self.env}")
        self.makeFolder(f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/")

    def _create_mol2_file(self):

        print(f"Converting the residue pdb file to a mol2 file")

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "mol2")
        mol = openbabel.OBMol()
        obConversion.ReadFile(
            mol,
            f"{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
        )
        # mol.AddHydrogens() TODO: Do we need this?

        assert (mol.NumResidues()) == 1
        obConversion.WriteFile(
            mol,
            f"{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.mol2",
        )

    def _modify_resname_in_str(self):

        fin = open(
            f"{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str",
            "rt",
        )
        fout = open(
            f"{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_tmp.str",
            "wt",
        )
        for line in fin:
            if line.startswith("RESI"):
                fout.write(line.replace(line.split()[1], "UNK"))
            else:
                fout.write(line)

        fin.close()
        fout.close()

        shutil.copy(fout.name, fin.name)

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
            
            # remove str file otherwise cgenff will append the new one at the end
            stream_file = f"{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str"
            if os.path.isfile(stream_file):
                os.remove(stream_file)

            # Run CGenFF
            ligand_path = f"{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}"

            cgenff_output = subprocess.run(
                [cgenff_bin]
                + [f"{ligand_path}.mol2"]
                + ["-v"]
                + ["-f"]
                + [f"{stream_file}"]
                + ["-m"]
                + [f"{ligand_path}.log"],
                capture_output=True,
                text=True,
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

    def _remove_lp(self, pdb_file_orig):
        # Will check if there are lone pairs and remove them
        # CGenFF will add them later on
        pdb_file = pm.load_file(pdb_file_orig, structure=True)
        lps = []
        for atom in pdb_file:
            if atom.name.startswith("LP"):
                print(f"We will remove {atom}, {atom.idx} ")
                lps.append(atom.idx)
        for i in range(len(lps)):
            pdb_file.strip(f"@{lps[i]+1-i}")

        return pdb_file

    def createCRDfiles(self):

        pdb_file_orig = f"{self.original_dir}/{self.ligand_id}.pdb"
        pdb_file = self._remove_lp(pdb_file_orig)
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
                            i.residue.name = (
                                "HSD"  # ATTENTION!! Here we make all HIS to HSD
                            )
                            i.residue.chain = f"PRO{chain}"
                        i.residue.chain = f"PRO{chain}"

            self.makeFolder(
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}"
            )

            df = pdb_file.to_dataframe()
            segids = set(i.residue.chain for i in pdb_file)
            for segid in segids:  # now we can save the crd files
                if segid not in ["SOLV", "IONS"]:
                    if segid == "HETA":
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.chain == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                            overwrite=True,
                        )
                    else:
                        try:
                            pdb_file[df.chain == f"{segid}"].save(
                                f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.crd",
                                overwrite=True,
                            )
                            pdb_file[df.chain == f"{segid}"].save(
                                f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.pdb",
                                overwrite=True,
                            )
                        except:
                            pass

        else:  # CHARMM-GUI generated pdb files
            print(f"Processing a CHARMM-GUI based pdb file")
            segids = set(i.residue.segid for i in pdb_file)
            for segid in segids:
                if segid not in ["SOLV", "IONS"]:
                    if segid == "HETA":
                        self.resname = pdb_file[df.segid == f"{segid}"].residues[0].name
                        self.makeFolder(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}"
                        )
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                            overwrite=True,
                        )

                    else:
                        try:
                            pdb_file[df.segid == f"{segid}"].save(
                                f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.crd",
                                overwrite=True,
                            )
                            pdb_file[df.segid == f"{segid}"].save(
                                f"{self.parent_dir}/{self.ligand_id}/complex/{segid.lower()}.pdb",
                                overwrite=True,
                            )
                        except:
                            pass

        # resname should be a 3 or 4 letters code
        assert len(self.resname) < 5

        return segids


class CharmmManipulation:
    def __init__(self, parent_dir, ligand_id, original_dir, resname, env):
        """
        CHARMM related files like the toppar file are modified, later the CHARMM executable is
        executed
        """ 
        self.parent_dir: str = parent_dir
        self.ligand_id: str = ligand_id
        self.original_dir: str = original_dir
        self.default_path: str = f"{self.parent_dir}/../data/templates/default/"
        self.resname: str = resname
        self.env: str = env

    def _manipulateToppar(self):

        # manipulate toppar_charmm.str file
        file = open(f"{self.ligand_id}/{self.env}/toppar.str", "a")
        file.write(f"stream {self.resname.lower()}/{self.resname.lower()}.str")
        file.close()
        # manipulate toppar.str file for use with openmm
        file = open(f"{self.ligand_id}/{self.env}/openmm/toppar.str", "a")
        file.write(f"stream ../{self.resname.lower()}/{self.resname.lower()}.str")
        file.close()

    def copyFiles(self):

        # copy CHARMM related files
        for file in glob.glob(f"{self.default_path}/[!toppar]*"):
            shutil.copy(file, f"{self.ligand_id}/{self.env}/")

        shutil.copy(
            f"{self.default_path}/toppar_charmm.str", f"{self.ligand_id}/{self.env}/toppar.str"
        )
        shutil.copy(
            f"{self.default_path}/toppar.str", f"{self.ligand_id}/{self.env}/openmm/toppar.str"
        )

        # copy files for OpenMM
        for file in glob.glob(f"{self.default_path}/*py"):
            shutil.copy(file, f"{self.ligand_id}/{self.env}/openmm/")

        try:
            shutil.copytree(
                f"{self.default_path}/toppar", f"{self.ligand_id}/{self.env}/toppar"
            )
        except:
            print(f"Toppar directory is already available")

    def modifyStep1(self, segids):

        self._manipulateToppar()

        correction_on = False

        fout = open(f"{self.ligand_id}/{self.env}/step1_pdbreader_tmp.inp", "wt")
        with open(f"{self.ligand_id}/{self.env}/step1_pdbreader.inp", "r+") as f:
            for line in f:
                if line.startswith("! Read PROA"):
                    correction_on = True
                    f.readline()
                    string = CharmmFactory.createHeader(segids, self.env)
                    fout.write(f"{string} \n")
                if line.startswith("!Print heavy atoms with "):
                    correction_on = False

                if correction_on == True:
                    pass
                else:
                    fout.write(line)
        fout.close()

        shutil.copy(fout.name, f"{self.ligand_id}/{self.env}/step1_pdbreader.inp")
        os.remove(fout.name)

    def _runCHARMM(self, step, charmm_exe):
        # We need to go to the specific directory to run CHARMM
        os.chdir(f"{self.ligand_id}/{self.env}/")
        output = subprocess.run(
            [f"{charmm_exe}"] + ["-i"] + [f"{step}.inp"] + ["-o"] + [f"{step}.out"],
            text=True,
            capture_output=True,
        )
        os.chdir("../../")
        if output.returncode:
            print(
                f"Something went wrong in step1 please check the outputfile in {self.ligand_id}/{self.env}/{step}.out"
            )
            sys.exit
        else:
            print(f"CHARMM process finished for {step}")

    def executeCHARMM(self, charmm_exe):

        steps = [
            "step1_pdbreader",
            "step2.1_waterbox",
            "step2.2_ions",
            "step2_solvator",
            "step3_pbcsetup_mod",
        ]
        for step in steps:
            self._runCHARMM(step, charmm_exe)
