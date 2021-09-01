Cloud-J v7.3c correction (to 7.3dx)

replace line 740 of cld_sub_mod.f90 version 7.3c (the GMD published version)
         if (ZZZ(L) .lt. Zbin(N)) then
with the following
         if (ZZZ(L)-ZZZ(1) .lt. Zbin(N)) then

This fix is recommended because then the code is insensitive to the
definition of ZZZ(L) as long as it is 'cm'.  ZZZ can then altitude
above surface or altitude above mean sea level.  NOTE that ZZZ adds
the Earth's radius to it in order to calculate the curved geometry of
the sun's rays.  In some cases with the UCI CTM runnning OpenMP,
GMAX(1) was undefined and gave segmentation faults.

The updated file on the web site (128.200.14.8) /public/prather/Fast-J
          UCI_CloudJ73dx.zip
contains an explanation and just the cld_sub_mod.f90 module.

The outputs from the standalone test code are not changed with this fix.

Note that Cloud-J v7.4 (designed for the extra super-bins (19:27) that
will be used in heating rates) has been fixed likewise and this version
is recommended as it handles the new formats for aerosol and cloud
spectral data.  V7.4 can still be run in simple Fast-J modes (including
tropospheric with only 8 out of 18 bins calculated) with lowered computation.

Michael Prather
13 Jan 2016
