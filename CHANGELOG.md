# Changelog

This file documents all notable changes to the Cloud-J repository since the initial commit under git version control.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [7.7.1] - 2024-04-02
### Changed
- Changed arguments of Init_Cldj to include root thread logical and LUT data directory
- Renamed global variable FP to FPAR to avoid name conflict in other models
- Changed FJX_scat-aer.dat format to expect 1 character greater width for R-eff and rho
- Changed FJX_scat-aer.dat format to expect 2 characters greater width for QAA (data column 2)
- Changed FJX_scat-aer.dat format to expect 1 character greater width for PAA (data columns 3-9)
- Limit prints to single core
- Renamed all modules to use .F90 suffix and cldj_ prefix

## Added
- Added CMake support
- Added .gitignore, authors list, and contribution and support guidelines
- Added module for error handling
- Added C pre-processor blocks for compatibility with GEOS-Chem offline CTM, GCHP, GEOS, and CESM
- Added new global variables RNAMES and BRANCH for species info
- Added GitHub config and PR/issue templates
- Added GitHub action for build tests

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
