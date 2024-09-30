# Changelog

This file documents all notable changes to the Cloud-J repository since the initial commit under git version control.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [8.0.1] - 2024-09-30
### Added
- Added error handling to exit the model if JVN_ is less than number of entries in FJX_j2j.dat

### Changed
- Updated usage of variable JVN_ to be max # of J-values rather than exact # of J-values

## [8.0.0] - 2024-08-29
### Added
- Added M. Prather's new Cloud-J v8 feature of UV absorption by water (source of differences between v7.7 and v8)
- Added Cloud-J input to specify whether to turn on UV absorption by water
- Added MODEL_GEOSCHEM C-preprocessor blocks to skip reading dat-files not used
- Added MODEL_STANDALONE C-preprocessor block to only read T and O3 climatology and H2O and CH4 profiles if using standalone
- Added H2O cross-sections to tables/FJX_spec.dat

### Changed
- Applied no-diff changes made by M. Prather between v7.7 and v8
- Moved most parameters set in CJ77_inp.dat to arguments passed to CLOUD_JX
- Hard-coded NRANDO where it is used rather than pass as configurable input
- Updated all files in tools/AddXs to be a copy of M. Prather's files included in his Cloud-J v8 release
- Renamed cldj_osa_sub_mod.F90 to cldj_fjx_osa_mod.F90
- Hardcoded physical constants in the cldj_cmn_mod module
- Renamed CLOUDJ_STANDALONE C-preprocessor switch to MODEL_STANDALONE
- Changed default standalone model setting for spherical correction from flat earth to spherical
- Changed default standalone model temperature climatology to be read from file
- Renamed Cloud-J standalone driver file from CJ77.F90 to cldj_standalone.F90

### Removed
- Removed subroutines SOLAR_JX and private subroutines it used
- Removed unused parameters and code related to RRTMG
- Removed configuration file CJ77_inp.dat

## [7.7.3] - 2024-08-12
### Added
- Added gcc 14 to github action to build Cloud-J on mac

### Changed
- Minimized usage and checks of error codes to speed up the model

### Removed
- Removed gcc 11 from github action to build Cloud-J on mac
- Removed safety initialization of outputs where assigning all values is obvious to improve performance

## [7.7.2] - 2024-07-12
### Added
- `.zenodo.json` file for auto-creation of DOI upon issuing a new release

### Changed
- Changed hard-coding LWEPAR in cldj_cmn_mod to be passed from parent model unless using standalone

### Added
- Added CLOUDJ_STANDALONE c-proprocessor switch to generalize code to use instead within a parent model
- Added integer status flag RC to most subroutines and pass it back up to parent model with error messages
- Added subroutine CLOUDJ_ERROR to print error message and set integer status flag RC

### Removed
- Removed subroutine EXITC and replaced it with CLOUDJ_ERROR_STOP in standalone and CLOUDJ_ERROR elsewhere

## [7.7.1] - 2024-04-02
### Changed
- Changed arguments of Init_Cldj to include root thread logical and LUT data directory
- Renamed global variable FP to FPAR to avoid name conflict in other models
- Changed FJX_scat-aer.dat format to expect 1 character greater width for R-eff and rho
- Changed FJX_scat-aer.dat format to expect 2 characters greater width for QAA (data column 2)
- Changed FJX_scat-aer.dat format to expect 1 character greater width for PAA (data columns 3-9)
- Limit prints to single core
- Renamed all modules to use .F90 suffix and cldj_ prefix
- Modified the cmake files for Standalone to build on mac successfully

## Added
- Added CMake support
- Added .gitignore, authors list, and contribution and support guidelines
- Added module for error handling
- Added C pre-processor blocks for compatibility with GEOS-Chem offline CTM, GCHP, GEOS, and CESM
- Added new global variables RNAMES and BRANCH for species info
- Added GitHub config and PR/issue templates
- Added GitHub action to test builds on Ubuntu
- Added GitHub action to test builds on Mac

### Fixed
- Fixed bugs for compatibility with gfortran compilers
- Fixed bugs for compatibility with intel compilers
- Fixed bug in standalone model for computing ice water path input to Cloud_JX

### Removed
- Removed unnecessary whitespace printed to log
- Removed CLDFLAG, LNRG, and NRANDO from CLOUD_JX argument list since they are global variables
- Removed unused arguments in ICA_* subroutines
- Removed unused common variables JXL_, JXL1_, and JXL2_
- Removed redundant file close in initialization

## [7.7]  - committed 2021-09-01
### Changed
- Replace files with version 7.7 from M. Prather's ftp site

## [7.4d] - committed 2021-09-01
### Changed
- Replaced files with version 7.4d from M. Prather's ftp site

## [7.3e] - committed 2021-09-01
### Added
- Added files for version 7.3e from M. Prather's ftp site
