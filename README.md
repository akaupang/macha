[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/akaupang/macha/workflows/CI/badge.svg)](https://github.com/akaupang/macha/actions)
[![codecov](https://codecov.io/gh/akuapang/macha/branch/main/graph/badge.svg)](https://codecov.io/gh/akaupang/macha/branch/main)


## maNUAL chaRMM (macha)

#### Description
macha is a meta-tool designed to (semi-)automate certain tasks related to setting up MD simulation systems. It is divided into a CLI front-end, written in bash, and a backend written in Python. The latter is accessible either directly or through the CLI frontend. 

Typical use cases include rerunning CHARMM-GUI scripts[^1] locally (after edits or changes) and preparing water boxes and complexes for further setup of free energy calculations with **transformato**[^2][^3], starting from PDB files of separate proteins and/or ligands and/or complexes, e.g. from protein data banks such as the [PDBe](https://www.ebi.ac.uk/pdbe/) or the [RCSB PDB](https://www.rcsb.org/). 

macha relies on CHARMM-GUI input scripts[^1] for CHARMM[^4], as well as on ParmEd[^5] and OpenBabel[^6]. Please, see the references at the bottom.

#### Installation

macha can be obtained by cloning this repository, issuing: 

`git clone git@github.com:akaupang/macha.git`

macha will soon be available as a conda package and can be installed in a conda environment with:

`conda install -c conda-forge macha`

In your python environment, you will also need the packages parmed, openbabel and [natsort](https://github.com/SethMMorton/natsort).

##### Post-installation setup
After installation, the user must set the path to a local `cgenff` binary in `main.py` to be able to use CGenFF parameterization of ligands (not required for single- or double-stranded RNA, which are natively parameterized in CHARMM).

For the CLI to work, a few variables should be set at the top of `manual_charmm_system_setup.sh`. The `macha_py_base` *must* point to the base directory of the python backend for this to be usable. If not set in the package version of `main.py`, the path to `cgenff` can be set here using the variable `cgenff_bin`. A user can also choose to set the variable `charmm_bin_man`, if a particular CHARMM binary should be used (by default, `charmm` is assumed to be in the $PATH). 

For convenient usage of the CLI frontend, we recommend setting an alias in your `~/.bashrc` to `macha`, like so:

`alias macha="/path/to/manual_charmm_system_setup.sh"`

Make sure that the `charmm` binary is in the $PATH or you have set a location in the script.

`export PATH="/path/to/charmm/bin:$PATH"`

#### Disclaimer:
Please note: This an early incarnation of macha, which may or may not be suitable for general use. We assume no responsibility for your use of the provided code. The macha CLI assumes that you are working in a folder containing input scripts generated by CHARMM-GUI[[1]](https://charmm-gui.org/)  and was designed to work with these scripts in the state they were provided in the years 2022-2023.

The CLI intends to keep up with developments in the Python backend, but may not always be up to date in terms of exposed functionality. Please inspect and/or use/modify the Python scripts directly if you encounter any problems.

#### Usage:

A main purpose for the CLI frontend is to give quick access to rerunning particular steps of the CHARMM-GUI-derived CHARMM input scripts, or all of them consecutively. The CLI also gives access to the main run types of the Python backend. The menus function by single-key selections and their operation should be fairly self-explanatory. 

Having set an alias as suggested above, issue `macha` in the working directory of choice. This directory should contain the CHARMM-GUI input scripts. Python3 needs to be in the $PATH or available through an active conda environment.

##### System creation for transformato
Call `macha` in the working directory of choice. This directory must contain a subdirectory called "data", which in turn contains a subdirectory called "original", in which PDB files of ligands, complexes or RNA should reside. These directory names can be changed by editing `main.py`. More advanced options, such as a segment ID filter to allow production of multimeric complexes and adjustable system pH for addition of hydrogens, are exposed in `main.py` and users with less typical use cases/input structures are encouraged to explore these. Proper documentation may follow at some point - for now, users are referred to the source code (found in `functions.py` and in 'charmm_factory.py`)

The submenu "System creation for transformato" is accessed by clicking "t". Here, five run types are exposed (the corresponding direct Python calls are shown in parentheses); 
- Make water boxes/complexes from ligands/proteins/complexes (`python3 main.py`)
- Make complexes (no water boxes) from ligands/proteins/complexes (`python3 main.py --nowaterbox`)
- Make water boxes (no complexes) from ligands/complexes (`python3 main.py --nocomplex`)
- Make water boxes/complexes from double-stranded RNA (`python3 main.py --rna`)
- Make water boxes (no complexes) from double- or single-stranded RNA (`python3 main.py --rna --nocomplex`)

Before running system creation, the `main.py` script must be copied to the working directory. This is done by choosing menu option "1" in the system creation for transformato submenu (*note that this overwrites any existing local copy*). This will also copy the path to the local cgenff binary to this local copy of `main.py`, if the variable `cgenff_bin` has been set at the top of the bash script `manual_charmm_system_setup.sh`, AND if the `cgenff_path` in `main.py` has not been set in the package version of this file `../macha/macha/main.py`.

One may then select the run type of choice by clicking a number from 2 - 6 (if using the CLI).

**Note** that mixed batches of input files are *not supported*, e.g. a typical run (for ligands/proteins/complexes) will fail to process input files with RNA as the ligand/guest, and vice versa; small-molecules will not be handled correctly by an RNA run. 


#### References
[^1]: Jo, S.; Kim, T.; Iyer, V. G.; Im, W. CHARMM-GUI: A Web-Based Graphical User Interface for CHARMM. *Journal of Computational Chemistry* **2008**, *29*(11), 1859–1865. https://doi.org/10.1002/jcc.20945 and https://charmm-gui.org/

[^2]: Wieder, M.; Fleck, M.; Braunsfeld, B.; Boresch, S. Alchemical Free Energy Simulations without Speed Limits. A Generic Framework to Calculate Free Energy Differences Independent of the Underlying Molecular Dynamics Program. *Journal of Computational Chemistry* **2022**, *43*(17), 1151–1160. https://doi.org/10.1002/jcc.26877 and https://github.com/wiederm/transformato

[^3]: Karwounopoulos, J.; Wieder, M.; Boresch, S. Relative Binding Free Energy Calculations with Transformato: A Molecular Dynamics Engine-Independent Tool. *Frontiers in Molecular Biosciences* **2022**, *9*, 954638. https://doi.org/10.3389/fmolb.2022.954638.

[^4]: Brooks, B. R.; Bruccoleri, R. E.; Olafson, B. D.; States, D. J.; Swaminathan, S.; Karplus, M. CHARMM: A Program for Macromolecular Energy, Minimization, and Dynamics Calculations. *Journal of Computational Chemistry* **1983**, *4*(2), 187–217. https://doi.org/10.1002/jcc.540040211 and https://academiccharmm.org/

[^5]: Shirts, M. R.; Klein, C.; Swails, J. M.; Yin, J.; Gilson, M. K.; Mobley, D. L.; Case, D. A.; Zhong, E. D. Lessons Learned from Comparing Molecular Dynamics Engines on the SAMPL5 Dataset. *Journal of Computer-Aided Molecular Design* **2017**, *31*(1), 147–161. https://doi.org/10.1007/s10822-016-9977-1 and https://github.com/ParmEd/ParmEd

[^6]: O’Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch, T.; Hutchison, G. R. Open Babel: An Open Chemical Toolbox. *Journal of Cheminformatics* **2011**, *3*(1), 33. https://doi.org/10.1186/1758-2946-3-33 and https://github.com/openbabel/openbabel




