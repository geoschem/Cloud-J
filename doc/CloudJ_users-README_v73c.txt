Fast-J users please note:

Fast-J code, including developments such as Cloud-J, is distributed under
GNU General Public License version 3.0 (GPLv3).
A copy of the GNU public license is at
http://sourceforge.net/p/webrtcbench/code/ci/default/tree/COPYING

Fast-J ftp site:   ftp://128.200.14.8/public/prather/Fast-J     (a.k.a.  halo.ess.uci.edu)

Recommendation (as of 18 Jul 2015):

        Implement ver 7.3c Cloud-J, skip ver 7.2, use the new data formats.
              v7.3c can be run with only 1 call to Fast-J (CLDFLAG = 1,2,3) if you must.
              v7.3c reduces the correlation  of G6 layers with gaps (Jul 18 fix).
              v7.3b corrects segmentation fault and compiler issues (Jul 3 fix) but does
                  not change any numbers in the 2015 GMDD publication.

        Implement ver 7.1c if you do not want to use Cloud-J in the near future.
              v7.1c has the same core Fast-J modules as Cloud-J.

        Implement ver 6.8d ONLY if you heritage code (e.g., back to ver 5.3)
              and cannot update to F90.  Be sure to update cross sections no matter what.


History of versions:

UCI_cloudJ73c-f90.zip
        Reduces the correlation of G6 layers separated by gaps.
        Updates fjx_spec.dat to version 7.3c: small, <% changes with new wavelength data.
        Results change at the <1% level; no change in the 2015 GMDD tables/figures.

UCI_cloudJ73b-f90.zip
        Corrects indexing error that could cause segment fault.
        No change in the 2015 GMDD numbers.

UCI_cloudJ73-f90.zip
        Is a standalone NEW cloud overlap code based on Prather 2015 GMDD paper.
        Has new v7.3+ formatting for the spectral data that is more flexible (as requested by users)

UCI_fastJX72-f90.zip
        This version contains the Cloud-J code based on the Neu et al 2007 JGR paper.
        It has been implemented and is running in the UCI and Oslo CTMs, but not CAM5.

UCI_fastJX71c-f90.zip
        Is the current 'just-Fast-J' fortran-90 version that is recommended.
        Differences bewteen v71c and v70d are minimal, and made to enable computation of
              J-values using only wavelengths >200nm for WACCM.
        BEWARE, there are two versions of fjx_spec.dat, the one for WACCM with zero solar
              fluxes in the first 4 bins.

UCI_fastJX70d-f90.zip
        Is the first fortran-90 Fast-J version with the latest bug fixes

UCI_fastJX68d-f77.zip
        Is the last fortran-77 version and includes corrected FJX_spec.dat (v6.8d 6/2014)
              where previous FJX_spec.dat had some VOC errors.
        Note that FJX_spec_68d.dat is also on the FastJ site separately

UCI_fastJX_newX_v68-v73.zip
        Has examples of code to generate new cross sections for both v6.8 and v7.3 formats
              Both v6.8 and v7.3 formats contain the same cross section data, but starting
              in v7.3, new formatting and extra data are added (solar flux in Watts, PAR spectrum).
