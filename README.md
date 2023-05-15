[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/akaupang/macha/workflows/CI/badge.svg)](https://github.com/akaupang/macha/actions)
[![codecov](https://codecov.io/gh/akuapang/macha/branch/main/graph/badge.svg)](https://codecov.io/gh/akaupang/macha/branch/main)


## maNUAL chaRMM (macha)

#### Description
macha is meta-tool designed to (semi-)automate certain tasks related to setting up MD simulation systems. It is divided into a CLI front-end, written in bash, and a backend written in Python. The latter is accessible either directly or through the CLI frontend. 

Typical use cases include rerunning CHARMM-GUI scripts locally (after edits or changes) and setting up water boxes and complexes (e.g. for further setup of free energy calculations with [transformato](https://github.com/wiederm/transformato)), starting from PDB files of separate proteins and/or ligands and/or complexes, e.g. from protein data banks such as the [PDBe](https://www.ebi.ac.uk/pdbe/) or the [RCSB PDB](https://www.rcsb.org/).

#### Installation

macha can be obtained by cloning this repository, issuing: 

`git clone git@github.com:akaupang/macha.git`

macha is also available as a conda package and can be installed in any conda environment with:

`conda install -c johanneskarwou -c conda-forge macha`

<!-- TODO: Update the conda package. Add packages parmed, openbabel, natsort to conda package requirements. -->

##### Post-installation setup
After installation, the user must set the path to `cgenff` in `main.py` to be able to use CGenFF parameterization of ligands (not required for single- or double-stranded RNA, which are natively parameterized in CHARMM).

For the CLI interface to work, a few variables should be set at the top of `manual_charmm_system_setup.sh`. The `macha_py_base` *must* be set to point to the base directory of the python backend for this to be usable. A user can also choose to set the `charmm_bin_man` if a particular CHARMM binary should be used (by default, `charmm` is assumed to be in the $PATH). 

For convenient usage of the CLI frontend, we recommend setting an alias in your `~/.bashrc` to `macha`, like so:

`alias macha="/path/to/manual_charmm_system_setup.sh"`

Make sure that the `charmm` binary is in the $PATH or you have set a location in the script.

`export PATH="/path/to/charmm/bin:$PATH"`

#### Disclaimer:
Please note: This an early incarnation of macha, which may or may not be suitable for general use. We assume no responsibility for your use of the provided code.
The macha CLI interface assumes that you are working in a folder containing input scripts 
generated by CHARMM-GUI (https://charmm-gui.org/) and was designed to work with these scripts
in the state they were provided in the years 2022-2023.

The CLI interface intends to keep up with developments in the Python backend, but may not
always be up to date. Please inspect and/or use/modify the python scripts directly if you
encounter any problems.

#### Usage:

The main task of the CLI frontend is to enable quick access to rerunning a particular step of the CHARMM-GUI-derived CHARMM input scripts, or all of them consecutively. The CLI also gives access to the main run types of the Python backend. The menus function by single-key selections and their operation should be fairly self-explanatory to anyone familiar with the CHARMM-GUI scripts and the function of each step.

##### System creation for transformato
The submenu "System creation for transformato" is accessed by clicking "t". Here, four run types are exposed; normal run (for proteins and ligands), RNA run, a run with exclusive generation of water boxes (no complexes) and a run with exclusive generation of complexes.

macha should be called in a working directory of choice. This directory must contain a subdirectory called "data", which in turn contains a subdirectory called "original", in which PDB files of ligands, complexes or RNA should reside. These directory names can be changed by editing main.py.

Before running system creation, the `main.py` script must be copied to the working directory. This is done by choosing menu option "1" in the system creation for transformato submenu. This will also copy the path to the local cgenff binary to the python script, if it has been set at the top of the bash script `manual_charmm_system_setup.sh`, AND if the `cgenff_path` in `main.py` has not been set in the package version of this file `../macha/macha/main.py`.

One may then select the run type of choice by clicking "2", "3", "4" or "5".

Note that mixed batches of input files are not supported - a default run (for proteins and/or ligands) will fail to process input files with RNA as the ligand/guest, and vice versa; small-molecules will not be handled correctly by an RNA run. 

