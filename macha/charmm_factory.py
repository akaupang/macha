#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in March 2023

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""

import sys
import os
import glob
import shutil
import subprocess
import parmed as pm
#from openbabel import openbabel
#from openbabel import pybel
#import string
#import numpy as np
#from .charmm_factory import CharmmFactory # SEPARATE PREPARATION AND CHARMMANIPULATION CLASSES AS SEPARATE FILES
#import re
import natsort

import warnings

# Supress ParmEd warnings
warnings.filterwarnings("ignore", module="parmed")


###############################################################################
class CharmmManipulation:
###############################################################################
    def __init__(
        self,
        parent_dir,
        ligand_id,
        original_dir,
        het_resnames,
        env,
        include_ions=True,
        ion_name="POT", # EXPAND FOR CATION AND ANION?
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
        self.het_resnames: list = het_resnames
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
    
    def copyFiles(self):
        # Copy files from template path
        # CHARMM related files
        for file in glob.glob(f"{self.default_path}/*[!omm_*][!toppar][!__pycache__]*"):
            shutil.copy(file, f"{self.parent_dir}/{self.ligand_id}/{self.env}/")

        # CHECK FOR USER INPUT ON ION TYPE AND ION STRENGTH
        # If someething other than the default is requested, edit the step2.2 ions 
        # stream file on the fly.
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
            print(f"A CHARMM toppar directory is already available")
    
    def _modifyIonsCountFile(self):

        # Create a target for a new step2.2 ions stream file
        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count_tmp.str",
            "wt",
        )
        
        # Read the existing template step2.2 ions stream file and replace the occurrences
        # of "POT" and "0.15", with the cation type and concentration requested by the user
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

    def _create_header(self, segids):      
        header = ""
        
        for segi in natsort.natsorted(segids):
            # RNA
            if segi.startswith("RNA"):
                header += f"""
! Read {segi.upper()}
open read card unit 10 name {segi.lower()}.crd
read sequence coor card unit 10 resid
generate {segi.upper()} setup warn first 5TER last 3TER

open read unit 10 card name {segi.lower()}.crd
read coor unit 10 card resid
        """
            # PROTEIN
            elif segi.startswith("PRO"):
                header += f"""
! Read {segi.upper()}
open read card unit 10 name {segi.lower()}.crd
read sequence coor card unit 10 resid
generate {segi.upper()} setup warn first NTER last CTER

open read unit 10 card name {segi.lower()}.crd
read coor unit 10 card resid
        """
            # HET SEGIDS
            elif segi.startswith("HET"):
                header += f"""
bomlev -1  ! not ideal but necessary for taking three-membered rings into account
! Read {segi.upper()}
open read card unit 10 name {segi.lower()}.crd
read sequence coor card unit 10 resid
generate {segi.lower()} setup warn first none last none

open read unit 10 card name {segi.lower()}.crd
read coor unit 10 card resid
bomlev 0
        """
            else:
                pass

        return header
    
    def modifyStep1(self, segids):     
        # Append the ligand toppar to the CHARMM toppar stream
        try:
            self._appendLigandTopparToTopparStream()
        except TypeError:
            print("Won't add a ligand stream file - assuming that the ligand is an RNA strand")

        # Open a a target for the new input file
        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader_tmp.inp",
            "wt",
        )
        # Read the template input file and 

        header_section_active = False # A switch to mark whether we are in the header 
                                      # section of step1_pdbreader as we read the file
                                      # and replace/copy its lines
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader.inp", "r+"
        ) as f:
            for line in f:
                if line.startswith("! Read PROA"): # Marks the start of the header section
                                                   # with the read/generate statements
                    header_section_active = True
                    f.readline() # Go to the next line
                    
                    # Create read/generate statements from segids
                    # NOTE WHY IS THIS NOT A FUNCTION IN THE CHARMMMANIPULATION CLASS, 
                    # BUT RATHER IN ITS OWN CLASS (CHARMMFACTORY)?
                    #header_block = CharmmFactory.createHeader(segids, self.env)
                    header_block = self._create_header(segids)

                    fout.write(f"{header_block}\n")                    
                    
                if line.startswith("!Print heavy atoms with "): # Marks the end of the header section
                    header_section_active = False

                # If we are in the header section, we should not record the lines in the
                # template file since these have been replaced by new ones
                if header_section_active == True:
                    pass
                else: # The rest of the step1_pdbreader.inp however, we copy verbatim
                    fout.write(line)
                    
                    
        fout.close()
        
        # Copy the finished modified file to the environment system directory
        shutil.copy(
            fout.name,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader.inp",
        )
        os.remove(fout.name)
        
    def _appendLigandTopparToTopparStream(self):
        # Add line to CHARMM toppar stream to read ligand parameters
        file = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str", "a")
        file.write(f"\n"\
                   f"! Ligand stream files\n"
                   )
        for resname in self.het_resnames:
            file.write(f"stream {resname.lower()}/{resname.lower()}.str\n")
        file.close()
        
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

            # Open a target for the modified step2_solvator input file
            fout = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2_solvator.inp", "wt")
            
            read_ions_section_active = False

            with open(f"{self.default_path}/step2_solvator.inp", "r+") as f:
                for line in f:
                    if line.startswith("! Add ions?"):              # Marks the start of the 
                        read_ions_section_active = True             # read ions section

                    if line.startswith("! Remove water molecules"): # Marks the end of the
                        read_ions_section_active = False            # read ions section

                    # If we are in the read ions section, we should not record these lines
                    # in the new file
                    if read_ions_section_active == True:
                        pass
                    elif line.startswith("* set"):                     # Also exclude
                        pass
                    elif line.startswith("* stream step2.2_ions.str"): # Also exclude
                        pass
                    else: # The rest of the step2_solvator.inp however, we copy verbatim
                        fout.write(line)
            fout.close()
            
            
        # Now, run CHARMM-GUI input scripts
        for step in steps:
            returncode = self._runCHARMM(step, charmm_exe)

    def _runCHARMM(self, step, charmm_exe):
        # Record the directory where we start
        main_dir = os.getcwd()
        
        # We need to go to the specific directory to run CHARMM
        os.chdir(f"{self.parent_dir}/{self.ligand_id}/{self.env}/")
        
        # Start a subprocess
        output = subprocess.run(
            [f"{charmm_exe}"] + ["-i"] + [f"{step}.inp"] + ["-o"] + [f"{step}.out"],
            text=True,
            capture_output=True,
        )
        
        # Return to the main directory
        os.chdir(main_dir)
        
        # A basic evaluation of the subprocess output codes
        if output.returncode:
            print(
                f"Something went wrong, please check the outputfile in {self.parent_dir}/{self.ligand_id}/{self.env}/{step}.out"
            )
            if step == "step1_pdbreader":
                print(
                    f"\n"\
                    f"INFO on common errors:\n"\
                    f"If GENIC is saying it cannot find residue XYZ, please check that CGenFF actually\n"\
                    f"produced a stream file with meaningful content from the input MOL2 file\n"\
                    f"(which was converted with OpenBabel from the input PDB file).\n"\
                    f"\n"\
                    f"If BILDC reports that it was called with a null IC table, please check the MOL2 file\n"\
                    f"and verify that atom numbers are present"
                )
            
            sys.exit("Processing terminated")
        else:
            print(f"CHARMM process finished for {step}")
            
            
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
                    
                # NOTE Here we are not allowing the overwriting of the ligand toppar folder
                # though that could be a useful behaviour. Consider adding
                # shutil.copytree(from, to, dirs_exist_ok=True) which works since Python 3.8
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
        
    def _convertCharmmTopparStreamToOpenmm(self, charmmstr_path, openmmstr_path):
        new_toppar = ""
        # Read the CHARMM toppar stream from the provided location
        with open(charmmstr_path, "r") as topparstream:
            for line in topparstream:
                if line.startswith("*"):       # DROP: Header
                    pass
                elif line.startswith("!"):     # DROP: Comments
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


    def applyHMR(self):

        # If input_orig does not exist, the system has not been updated/overwritten yet, and
        # this is likely the first pass. Thus, we should make a backup.
        if not os.path.isfile(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input_orig.psf"
        ):

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
        # we need to load the parameters now to update the masses in the CharmmPsfFile class.
        psf.load_parameters(params)

        # Apply HMR and save PSF
        pm.tools.actions.HMassRepartition(psf).execute()
        psf.save(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input.psf",
            overwrite=True,
        )

        # Assure that hydrogen mass in a TIP3P water is greater than one
        for atom in psf:
            if atom.name.startswith("H") and atom.residue.name != "TIP3":
                assert atom.mass > 1.5
        print("Successfully applied HMR")



    def createTFYamlFile(self, dt=0.001, nstep=2500000):

        try:
            os.makedirs(f"{self.parent_dir}/config")
        except:
            pass

        if self.env == "waterbox":
            for resname in self.het_resnames:
                # Creates a yaml file for ASFE simulations using TF
                fin = open(f"{self.default_path}/../temp_asfe.yaml", "r")
                fout = open(f"{self.parent_dir}/config/{self.ligand_id}_{resname.lower()}.yaml", "w")
                for line in fin:
                    if line.strip().startswith("name"):
                        new_line = line.replace("NAME", f"{self.ligand_id}")
                        fout.write(new_line)
                    elif line.strip().startswith("tlc"):
                        new_line = line.replace("UNK", f"{resname}")
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

                print(f"Created YAML file {self.parent_dir}/config/{self.ligand_id}_{resname.lower()}.yaml for ASFE simulations using Transformato")
                
    def _apply_constraints(self, fout):
        with open(fout, "r") as file:
            filedata = file.read()
        filedata = filedata.replace("cons: None", "cons: HBonds")
        with open(fout, "w") as file:
            file.write(filedata)
