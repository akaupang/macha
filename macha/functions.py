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
import string
import numpy as np
from .charmm_factory import CharmmFactory

import warnings

# Supress ParmEd warnings
warnings.filterwarnings("ignore", module="parmed")


def checkInput(
    parent_dir=".", original_dir="original", protein_name=None, input_ext="pdb"
):

    # Make sure the input directory exists
    if not os.path.isdir(f"{parent_dir}/{original_dir}"):
        sys.exit(f"Input directory: {parent_dir}/{original_dir} not found!")

    # Look for an input protein for ligand exchange
    if protein_name != None:
        protein_path = f"{parent_dir}/{original_dir}/{protein_name}.{input_ext}"
        if not os.path.isfile(protein_path):
            print(f"No protein for ligand exchange found.")
            protein_id = None
        else:
            print(f"Protein for ligand exchange found at: {protein_path}")
            protein_id = os.path.splitext(os.path.basename(protein_path))[0]
            print(f"Protein ID: {protein_id}")
    else:
        protein_id = None

    if protein_id == None:
        print(
            """
        No protein was specified for ligand exchange with multiple ligands.
        The input is thus assumed to consist either of complexes from which to create
        waterboxes + complexes, OR of single ligands from which to create waterboxes
        (use the -nc switch for efficiency).
        """
        )

    # Create the master list of the ligand/system names and fill it
    ligand_ids = []
    # Add multiple ligands
    for ifile in glob.glob(f"{parent_dir}/{original_dir}/*.{input_ext}"):
        lig_id = os.path.splitext(os.path.basename(ifile))[0]
        # But don't add the protein/complex
        if lig_id == protein_id:
            pass
        else:
            ligand_ids.append(lig_id)  #
    return protein_id, ligand_ids


class Preparation:
    def __init__(
        self,
        parent_dir,
        original_dir,
        ligand_id,
        env,
        protein_id=None,
        small_molecule=False,
        rna=True,
    ):
        """
        This class prepares everything for further use with CHARMM. The PDB
        files are sliced into pieces and the ligand is converted to a mol2 file.
        A local version of CGenFF creates a stream file for the ligand.
        """
        self.parent_dir = parent_dir
        self.original_dir = original_dir
        self.ligand_id = ligand_id
        self.protein_id = protein_id
        self.resname = str
        self.env: str = env
        self.small_molecule: bool = small_molecule
        self.rna: bool = rna
        self.input_type = None

        # Distinguished treatment if the current ligand is to be merged with a protein

        if self.protein_id == None:
            # For complexes->waterbox+complex or single ligands->waterbox
            # Load the PDB file into ParmEd
            self.pdb_file_orig = (
                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
            )
            self.pdb_file = pm.load_file(self.pdb_file_orig, structure=True)
            if self.env == "single_strand":
                self.pdb_file_orig = (
                    f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
                )
                self.pdb_file = pm.load_file(self.pdb_file_orig, structure=True)
                self.pdb_file = self.pdb_file["A", :, :]  # select only CHAIN A
        else:
            if self.env == "waterbox" or self.env == "double_strand":
                # Run like normal single ligand
                self.pdb_file_orig = (
                    f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
                )
                self.pdb_file = pm.load_file(self.pdb_file_orig, structure=True)
            elif self.env == "single_strand":
                self.pdb_file_orig = (
                    f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
                )
                self.pdb_file = pm.load_file(self.pdb_file_orig, structure=True)
                self.pdb_file = self.pdb_file["A", :, :]  # select only CHAIN A
                print(f"Wir sind hier")
            elif self.env == "complex":
                # Start the merger function
                self.mergeToComplex()

            else:
                sys.exit(f"Unrecognized environment: {self.env}")

    def mergeToComplex(self):
        print("* * * Start merging * * *")
        # Here, it is practical to borrow the "self.pdb_file" temporarily
        # for both the input protein and later for the ligand. This solution
        # benefits from the fact that most of the internal functions in this
        # class already work on this object. The merge function will in the end
        # restore "self.pdb_file" as the complex.

        # Load the protein
        print(f"Processing protein (input typing, segid addition, HETx removal)")
        self.pdb_file = pm.load_file(
            f"{self.original_dir}/{self.protein_id}.pdb", structure=True
        )

        # What kind of protein file has been provided?
        # Remember: segids will be added if there are none
        segids, df = self.checkInputType()

        # Get the parts of the parmed object (protein) that does not contain
        # "HET" in the segid column
        # Note the tilde to invert the boolean mask!
        apo_pdb = self.pdb_file[~df.segid.str.contains("HET")]

        # Load the ligand
        print(f"Processing ligand: {self.ligand_id} (input typing, segid addition)")
        self.pdb_file = pm.load_file(
            f"{self.original_dir}/{self.ligand_id}.pdb", structure=True
        )

        # Add segids to ligand PDB object if there are none
        segids, df = self.checkInputType()

        lig_pdb = self.pdb_file

        # Overwrite the self.pdb_file by joining the structures
        # (no reordering to a typical PROA, HETA, IONS, WAT/SOLV)
        self.pdb_file = apo_pdb + lig_pdb

        print("* * * End merging * * *")

    def _make_folder(self, path):

        try:
            os.makedirs(path)
            print(f"Creating folder in {path} for the {self.env}")
        except OSError:
            print(f"There exists a folder in {path} for {self.env} - we will use it")
            if not os.path.isdir(path):
                raise

    def makeTFFolderStructure(self):
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/")
        self._make_folder(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/restraints/"
        )

    def _create_mol2_file(self):

        print(f"Converting ligand {self.ligand_id} to MOL2 and SDF files")

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "mol2")
        mol = openbabel.OBMol()

        if self.small_molecule:
            obConversion.ReadFile(
                mol,
                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb",
            )
        else:
            obConversion.ReadFile(
                mol,
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
            )
        # Useful for H-less ligands directly from x-ray crystallography PDBs
        if self.input_type == "allhydrogens":
            pass
        elif self.input_type == "missinghydrogens":
            mol.AddHydrogens()

        assert (mol.NumResidues()) == 1
        obConversion.WriteFile(
            mol,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.mol2",
        )
        obConversion.SetInAndOutFormats("pdb", "sdf")
        # mol = openbabel.OBMol()
        # obConversion.ReadFile(
        #     mol,
        #     f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
        # )
        obConversion.WriteFile(
            mol,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.sdf",
        )

    def _modify_resname_in_stream(self):

        fin = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str",
            "rt",
        )
        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_tmp.str",
            "wt",
        )
        for line in fin:
            if line.startswith("RESI"):
                fout.write(line.replace(line.split()[1], self.resname))
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

            # Remove stream file otherwise CGenFF will append the new one at
            # the end of the old one
            stream_file = f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str"
            if os.path.isfile(stream_file):
                os.remove(stream_file)

            # Run CGenFF
            ligand_path = f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}"

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

        self._modify_resname_in_stream()

        return cgenff_output

    def _remove_lp(self):
        # Will check if there are lone pairs and remove them
        # CGenFF will add them later on
        lps = []
        for atom in self.pdb_file:
            if (atom.name.startswith("LP")) or (atom.name.startswith("Lp")):
                print(f"We will remove {atom}, {atom.idx} ")
                lps.append(atom.idx)

        for i in range(len(lps)):
            self.pdb_file.strip(f"@{lps[i]+1-i}")

        return self.pdb_file

    def _check_ionizable(self):
        # A preliminary list of residue names for titrable residues that indicate that the user has set them
        # AMBER names will be converted to CHARMM names, no questions asked.
        # TODO This list should be expanded to be comprehensive
        ionized_res = [
            "HID",
            "HIE",
            "HIP",
            "HSD",
            "HSE",
            "HSP",
            "CYX",
            "CYM",
            "ASPP",
            "GLUP",
        ]

        for (
            res
        ) in (
            self.pdb_file.residues
        ):  # produce a select view of this chain using a boolean mask
            if res.name in ionized_res:
                # Harmonize to CHARMM names (from AMBER names)
                # TODO This section should be expanded to be comprehensive
                if res.name == "HID":
                    res.name = "HSD"
                elif res.name == "HIE":
                    res.name = "HSE"
                elif res.name == "HIP":
                    res.name = "HSP"
                elif res.name == "CYX":
                    res.name = "CYM"
                else:
                    pass

            elif res.name == "HIS":
                # TODO A lookup using PROPKA or another solution could be of interest here
                res.name = "HSD"  # ATTENTION!! Here we make all HIS to HSD

        return self.pdb_file

    def _check_for_hydrogens(self):
        # Do hydrogens cover all residues?
        h_cont_res = [
            res.name
            for res in self.pdb_file.residues
            if any([True for atm in res.atoms if atm.element_name == "H"])
        ]
        no_h_res = [
            res.name
            for res in self.pdb_file.residues
            if any([True for atm in res.atoms if atm.element_name != "H"])
        ]

        h_miss_res = list(np.setdiff1d(no_h_res, h_cont_res))

        # if any([True for ele in self.pdb_file.atoms if ele.element_name == "H"]):
        if len(h_cont_res) == len(no_h_res):
            self.input_type = "allhydrogens"
            print("All residues appear to be fully hydrogenated")
        else:
            self.input_type = "missinghydrogens"
            print(f"Some residues appear to lack hydrogens: {*h_miss_res,}")

    def _create_tlc_rna(self):

        for atom in self.pdb_file.atoms:
            if atom.residue.name == "G":
                atom.residue.name = "GUA"
            elif atom.residue.name == "C":
                atom.residue.name = "CYT"
            elif atom.residue.name == "A":
                atom.residue.name = "ADE"
            elif atom.residue.name == "U":
                atom.residue.name = "URA"
            elif atom.residue.name == "T":
                atom.residue.name = "THY"
            elif atom.residue.name == "I":
                atom.residue.name = "INO"

        return self.pdb_file

    def _add_segids_rna(self, df):

        chids = set(i.residue.chain for i in self.pdb_file)
        print(f"Found chain IDs: {chids}")
        for chain in chids:
            for res in self.pdb_file.view[
                df.chain == f"{chain}"
            ].residues:  # produce a select view of this chain's residues using a boolean mask
                res.segid = f"RNA{chain}"

        return self.pdb_file

    def _add_segids(self, df):
        # This function adds segids to the Parmed object from a PDB that does
        # not contain segids, but only chain ids, ensuring that segids can be
        # used in making of crd files later.
        # The function attempts to establish residue names that lie outside the
        # range of "normal" amino acids names and will give these the name HETx
        # in which x is the chain ID

        # TODO This list should be expanded to be comprehensive
        aa_res = [
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
            "HID",
            "HIE",
            "HIP",
            "HSD",
            "HSE",
            "HSP",
            "CYX",
            "CYM",
            "ASPP",
            "GLUP",
        ]  # we need to find the residue of the ligand which should be the only one being not an aa

        exclude_res = [  # TODO Many more common crystallization additives should be added to this list
            "1PE"
        ]

        # ions_res =[] #TODO A list of ions could be added for their particular handling

        # Make a list of ABCD... for use as x in HETx
        het_letters = list(string.ascii_uppercase)

        # TODO Make this loop handle multiple ligands in one chain
        chids = set(i.residue.chain for i in self.pdb_file)
        # DEBUG
        # print(f"Found chain IDs: {chids}")
        for (
            chain
        ) in chids:  # rename the chain names (a,b, ...) to segnames (proa,prob,...)
            # DEBUG
            # print(f"Working on chain {chain}.")

            # Set a temporary previous residue number
            prev_resnum = -1

            # Loop over the residues in the current chain
            for res in self.pdb_file.view[
                df.chain == f"{chain}"
            ].residues:  # produce a select view of this chain's residues using a boolean mask
                resnum = res.number
                if res.name in aa_res:
                    res.segid = f"PRO{chain}"
                elif res.name == "HOH":
                    res.segid = f"WAT{chain}"
                elif res.name in exclude_res:
                    res.segid = f"CRST"
                else:
                    if resnum == prev_resnum:
                        pass
                    else:
                        res.segid = f"HET{het_letters.pop(0)}"
                    # self.resname = i.residue.name
                prev_resnum = resnum

        return self.pdb_file

    def checkInputType(self):

        # Remove lone pairs (if any)
        self.pdb_file = self._remove_lp()

        # Check/rename ionizable residues (if the user has not set their protonation states)
        # ATTENTION: All occurrences of HIS will become HSD
        self.pdb_file = self._check_ionizable()

        # Store the input PDB as a dataframe for later use
        df = self.pdb_file.to_dataframe()

        # We will cater for two types of PDB files
        # 1: Canonical/MAESTRO-derived PDB files
        # 2: CHARMM-derived PDB files
        # This distinction is made by looking for segids.
        # Canonical/MAESTRO-derived PDB files will be void of segids, but will
        # have chain IDs, while CHARMM PDB files contain the segid column, but
        # lack chain IDs.

        # # Count the segids in the PDB file (if any)
        segids = set(i.residue.segid for i in self.pdb_file)

        if self.rna:
            self.pdb_file = self._add_segids_rna(df)
            self.pdb_file = self._create_tlc_rna()
            segids = set(i.residue.segid for i in self.pdb_file)
            df = self.pdb_file.to_dataframe()
        # For canonical/MAESTRO-derived PDB files containing chain IDs
        elif segids == {""}:  # empty set of segids
            print(f"Processing a canonical/Maestro based pdb file (without segids)")
            # Check for hydrogens
            # set self.input_type
            self._check_for_hydrogens()

            # Add segids to the Parmed object
            self.pdb_file = self._add_segids(df)
            # Update the dataframe and segid list
            df = self.pdb_file.to_dataframe()
            segids = set(i.residue.segid for i in self.pdb_file)
        # For CHARMM-GUI generated PDB files
        else:
            print(f"Processing a CHARMM PDB file (with segids)")
            # Check for hydrogens
            # set self.input_type
            self._check_for_hydrogens()

        return segids, df

    def createCRDfiles(self, segids, df):

        exclude_segids = [
            "SOLV",
            "IONS",
            "WATA",
            "WATB",
            "WATC",
            "CRST",
            "HETB",
            "HETC",
        ]
        used_segids = []
        for segid in segids:
            if segid not in exclude_segids:  # multiple ligands are excluded for now
                # Store in the segids to be given to CharmmManipulation
                used_segids.append(segid)

                # WATERBOX ENVIRONMENT
                if self.env != "complex":
                    if segid.startswith("HET"):

                        # Note the residue name for checks
                        self.resname = (
                            self.pdb_file[df.segid == f"{segid}"].residues[0].name
                        )
                        # resname should be a 3 or 4 letters code
                        # assert len(self.resname) < 5

                        self._make_folder(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}"
                        )
                        self.pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        if self.small_molecule:
                            # We copy the pdf file since there shouldn't be any changes
                            shutil.copy(
                                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb",
                                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                            )
                        else:
                            self.pdb_file[df.segid == f"{segid}"].save(
                                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                                overwrite=True,
                            )
                    elif segid.startswith("RNA"):

                        self.pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )

                # COMPLEX ENVIRONMENT
                elif self.env == "complex":
                    if segid.startswith("HET"):

                        # TODO: Should this check be more generally applied?
                        # Note the residue name for checks
                        self.resname = (
                            self.pdb_file[df.segid == f"{segid}"].residues[0].name
                        )
                        # resname should be a 3 or 4 letters code
                        assert len(self.resname) < 5

                        self._make_folder(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}"
                        )
                        self.pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        self.pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                            overwrite=True,
                        )

                    else:

                        self.pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        self.pdb_file[df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.pdb",
                            overwrite=True,
                        )
                else:
                    sys.exit(f"Unrecognized environment: {self.env}")

        return segids, used_segids


class CharmmManipulation:
    def __init__(
        self,
        parent_dir,
        ligand_id,
        original_dir,
        resname,
        env,
        include_ions=True,
        ion_name="POT",
        ion_conc=0.15,
        default_path="",
    ):
        """
        CHARMM related files like the toppar file are modified, later the CHARMM executable is
        executed
        """
        self.parent_dir: str = parent_dir
        self.ligand_id: str = ligand_id
        self.original_dir: str = original_dir
        self.resname: str = resname
        self.env: str = env
        self.include_ions: bool = include_ions
        self.ion_name: str = ion_name
        self.ion_conc: float = ion_conc
        if not default_path:
            self.default_path: str = self.get_default_path()
        else:
            self.default_path = default_path

    def get_default_path(self):

        import macha

        return f"{macha.__path__[0]}/data/templates/default/"

    def _appendToppar(self):

        # Add line to CHARMM toppar stream to read ligand parameters
        file = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str", "a")
        file.write(f"stream {self.resname.lower()}/{self.resname.lower()}.str")
        file.close()

    def _convertCharmmTopparStreamToOpenmm(self, charmmstr_path, openmmstr_path):
        new_toppar = ""
        # Read the CHARMM toppar stream from the provided location
        with open(charmmstr_path, "r") as topparstream:
            for line in topparstream:
                if line.startswith("*"):  # DROP: Header
                    pass
                elif line.startswith("!"):  # DROP: Comments
                    pass
                elif line.strip(" ") == "\n":  # DROP: Empty lines
                    pass
                elif line.startswith("read"):  # DROP: Legacy read statements
                    pass
                else:
                    strpath = line.split(" ")[-1]
                    new_toppar += f"./{strpath}"
                    # Split on whitespace and take the last item in the
                    # resulting list of words (the path to each stream file,
                    # including the newline characters.

        # Write the new OpenMM toppar stream to the requested location
        with open(openmmstr_path, "w") as topparstream:
            topparstream.write(new_toppar)

    def _getExternalTopparFromTopparStream(self, topparstream):
        external_toppars = []
        with open(topparstream, "r") as tp_str:
            for line in tp_str:
                if "toppar" in line:
                    pass
                else:
                    external_toppars.append(line.split("/")[1])

        return external_toppars

    def _modifyIonsCountFile(self):

        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count_tmp.str",
            "wt",
        )
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count.str",
            "r+",
        ) as f:
            for line in f:
                if "set pos" in line:
                    new_line = line.replace("POT", str(self.ion_name))
                    fout.write(new_line)
                elif "set conc" in line:
                    new_line = line.replace("0.15", str(self.ion_conc))
                    fout.write(new_line)
                else:
                    fout.write(line)
        fout.close()

        shutil.copy(
            fout.name,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count.str",
        )
        os.remove(fout.name)

        pass

    def copyFiles(self):
        # Copy files from template path
        # CHARMM related files
        for file in glob.glob(f"{self.default_path}/*[!omm_*][!toppar][!__pycache__]*"):
            shutil.copy(file, f"{self.parent_dir}/{self.ligand_id}/{self.env}/")

        if self.ion_name != "POT" or self.ion_conc != 0.15:
            self._modifyIonsCountFile()

        # CheckFFT script
        shutil.copy(
            f"{self.default_path}/checkfft.py",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/checkfft.py",
        )

        # The CHARMM toppar.str (which later will be converted to OpenMM format)
        shutil.copy(
            f"{self.default_path}/toppar.str",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str",
        )
        # The whole toppar directory
        try:
            shutil.copytree(
                f"{self.default_path}/toppar",
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar",
            )
        except:
            print(f"Toppar directory is already available")

    def modifyStep1(self, segids):

        # Append the ligand toppar to the CHARMM toppar stream
        try:
            self._appendToppar()
        except TypeError:
            print("Won't add a ligand stream file, assuming it's a RNA strand")

        correction_on = False

        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader_tmp.inp",
            "wt",
        )
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader.inp", "r+"
        ) as f:
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

        shutil.copy(
            fout.name,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader.inp",
        )
        os.remove(fout.name)

    def _runCHARMM(self, step, charmm_exe):
        # We need to go to the specific directory to run CHARMM
        main_dir = os.getcwd()
        os.chdir(f"{self.parent_dir}/{self.ligand_id}/{self.env}/")
        output = subprocess.run(
            [f"{charmm_exe}"] + ["-i"] + [f"{step}.inp"] + ["-o"] + [f"{step}.out"],
            text=True,
            capture_output=True,
        )
        os.chdir(main_dir)
        if output.returncode:
            print(
                f"Something went wrong, please check the outputfile in {self.parent_dir}/{self.ligand_id}/{self.env}/{step}.out"
            )
            sys.exit("Processing terminated")
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

        if not self.include_ions:
            # avoid processing the step2.2_ions.inp CHARMM file
            # and remove the lines reading in these results in the
            # next CHARMM input file (step2_solvator)
            steps.remove("step2.2_ions")

            dir = f"{self.parent_dir}/{self.ligand_id}/{self.env}"
            correction_on = False
            fout = open(f"{dir}/step2_solvator.inp", "wt")
            with open(f"{self.default_path}/step2_solvator.inp", "r+") as f:
                for line in f:
                    if line.startswith("! Add ions?"):
                        correction_on = True

                    if line.startswith("! Remove water molecules"):
                        correction_on = False

                    if correction_on == True:
                        pass
                    elif line.startswith("* set"):
                        pass
                    elif line.startswith("* stream step2.2_ions.str"):
                        pass
                    else:
                        fout.write(line)
            fout.close()

        for step in steps:
            returncode = self._runCHARMM(step, charmm_exe)

    def createOpenMMSystem(self):
        print("Creating a stand-alone OpenMM system")

        # Copy files for OpenMM
        # OpenMM Python scripts
        for file in glob.glob(f"{self.default_path}/[!checkfft.py]*py"):
            shutil.copy(file, f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/")

        # Equilibration input file
        shutil.copy(
            f"{self.default_path}/omm_step4_equilibration.inp",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step4_equilibration.inp",
        )

        # Production input file
        shutil.copy(
            f"{self.default_path}/omm_step5_production.inp",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step5_production.inp",
        )

        # Toppar directory from CHARMM base directory to the OpenMM directory
        try:
            shutil.copytree(
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar",
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar",
            )
        except:
            print(f"Toppar directory is already available in the OpenMM directory")

        # Convert CHARMM toppar.str to OpenMM toppar.str
        self._convertCharmmTopparStreamToOpenmm(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str",  #
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar.str",
        )

        # Inspect OpenMM toppar.str to identify ligands whose parameter
        # folders should be copied
        external_toppars = self._getExternalTopparFromTopparStream(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar.str"
        )

        if external_toppars != []:
            for ext_top in external_toppars:
                try:
                    shutil.copytree(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{ext_top}",
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/{ext_top}",
                    )
                except FileNotFoundError:
                    print(
                        f"ERROR: External topology/parameter folder {ext_top} not found!"
                    )
                except FileExistsError:
                    print(
                        f"The topology/parameters folder {ext_top} is already present."
                    )

        # Create sysinfo.dat file, which contains information about the box
        file = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/step3_pbcsetup.str")
        for line in file:
            if line.startswith(f" SET A "):
                a = line.split(" ")[-1].strip()
            elif line.startswith(f" SET B "):
                b = line.split(" ")[-1].strip()
            elif line.startswith(f" SET C "):
                c = line.split(" ")[-1].strip()
            elif line.startswith(f" SET ALPHA "):
                alpha = line.split(" ")[-1].strip()
            elif line.startswith(f" SET BETA "):
                beta = line.split(" ")[-1].strip()
            elif line.startswith(f" SET GAMMA "):
                gamma = line.split(" ")[-1].strip()

        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/sysinfo.dat", "w"
        ) as f:
            f.write(
                '{"dimensions":[' + f"{a}, {b}, {c}, {alpha}, {beta}, {gamma}" + "]}"
            )

    def applyHMR(self):

        # try:
        # If input_orig does not exist (the system has not been updated/overwritten)
        if not os.path.isfile(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input_orig.psf"
        ):

            # Load parameters
            # parms = ()
            # for file in glob.glob(f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar/*[!tip216.crd]*"):
            #     parms += ( file, )

            # parms += (f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str",)

            # Copy the original PSF to a backup (input_orig)
            shutil.copy(
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input.psf",
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input_orig.psf",
            )
            print("Backed up OpenMM PSF file")

            input_psf = "step3_input.psf"

        # A backup of the orginal system exists (input_orig)
        else:
            input_psf = "step3_input_orig.psf"
            print("Reading PSF file from backup")

        ## Load parameter into ParmEd
        # PARAMETERS ARE LOADED DUE TO A BUG IN PARMED*
        #
        # This way of reading the parameters preserves the order of the
        # original/CHARMM-GUI-derived toppar.str, and ensures that
        # duplicate parameters are read in their "intended" order - the
        # last parameters read will overwrite older parameters and be
        # those used.
        #
        parms = ()
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar.str", "r"
        ) as ommtopparstream:
            for line in ommtopparstream:
                parms += (line.strip("\n"),)

        # Change directory to for ParmEd to find the parameter files
        # load the parameters and change back to the original run directory
        cur_dir = os.getcwd()
        os.chdir(f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/")
        params = pm.charmm.CharmmParameterSet(*parms)
        os.chdir(cur_dir)

        # Load the PSF into ParmEd
        psf = pm.charmm.CharmmPsfFile(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/{input_psf}"
        )

        # *Since ParmEd forgot about the masses when loading the PSF above,
        # we need to load the parameters to update the masses in the CharmmPsfFile class.
        psf.load_parameters(params)

        # Apply HMR and save PSF
        pm.tools.actions.HMassRepartition(psf).execute()
        psf.save(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input.psf",
            overwrite=True,
        )

        # Assure that mass is greater than one
        for atom in psf:
            if atom.name.startswith("H") and atom.residue.name != "TIP3":
                assert atom.mass > 1.5
        print("Successfully applied HMR.")

        # except:
        #    print("Masses were left unchanged! HMR not possible, check your output! ")

    def _apply_constraints(self, fout):

        with open(fout, "r") as file:
            filedata = file.read()
        filedata = filedata.replace("cons: None", "cons: HBonds")
        with open(fout, "w") as file:
            file.write(filedata)

    def createTFYamlFile(self, dt=0.001, nstep=2500000):

        try:
            os.makedirs(f"{self.parent_dir}/config")
        except:
            pass

        if self.env == "waterbox":

            # Creates a yaml file for ASFE simulations using TF

            fin = open(f"{self.default_path}/../temp_asfe.yaml", "r")
            fout = open(f"{self.parent_dir}/config/{self.ligand_id}.yaml", "w")
            for line in fin:
                if line.strip().startswith("name"):
                    new_line = line.replace("NAME", f"{self.ligand_id}")
                    fout.write(new_line)
                elif line.strip().startswith("tlc"):
                    new_line = line.replace("UNK", f"{self.resname}")
                    fout.write(new_line)
                elif line.strip().startswith("dt:"):
                    new_line = line.replace("TIMESTEP", str(dt))
                    fout.write(new_line)
                elif line.strip().startswith("nstep:"):
                    new_line = line.replace("NUMBEROFSTEPS", str(nstep))
                    fout.write(new_line)
                else:
                    fout.write(line)

            fin.close()
            fout.close()

            if dt != 0.001:
                self._apply_constraints(fout.name)
