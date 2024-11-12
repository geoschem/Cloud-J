!------------------------------------------------------------------------------
!    cldj_cmn_mod.F90 is based on 'fjx_cmn_mod.f90' in Cloud/Fast-J by M. Prather
!    Initial adaptation by E. Lundgren 10/2021
!
!    Notes:
!      Cloud/Fast-J code v 8.0 (prather 04/2023):
!         (1) The last Solar-J version 7.6c ran with RRTMG.
!         (2) Current Cloud-J-only version 8.0 with S_ = 18 drops wavelengths 19:27,
!               leaves them in the datsets
!         (3) If the CLIRAD or LLNL heating codes are invoked then Cloud-J then go
!               back to v7.6c and adapt the Spec, Cloud and SSA datasets need to be
!               replaced with the CLIRAD or LLNL bins and you will need to reinsert
!               and check the CLIRAD and LLNL subroutines.
!         (4) Includes cleanup and removal of dead-end variables
!      Solar/Cloud/Fast-J code v 7.7 (prather 02/2020)
!         (1) Enabled LLNL heating code, fix SJSUB dim's
!
!------------------------------------------------------------------------------

      MODULE CLDJ_CMN_MOD

      implicit none
      public

!------------------------------------------------------------------------------
! Physical constants
!------------------------------------------------------------------------------

      real*8, parameter:: RAD      = 6375.0d5  ! Radius of Earth (cm)
      real*8, parameter:: ZZHT     =    5.0d5  ! Scale height (cm) used above top of CTM ZHL(LPAR+1)
      
!------------------------------------------------------------------------------
! Basic vertical grid
!------------------------------------------------------------------------------

      ! Can be changed as needed. These must be set to exact atmospheric dimensions.
      ! 57 levels is only for the UCI CTM atmospheres in the standalone example.
      ! For other uses, pass the # of levels and # layers with clouds from
      ! the parent model.
#ifdef MODEL_STANDALONE
      integer, parameter :: L_  = 57    !  # of CTM layers, set at build-time
      integer, parameter :: L1_ = L_+1  !  L_+1 = # of CTM layer edges (radii)
      integer, parameter :: L2_ = L_+2  !  L_+2 = total # of layer edges counting top (TAU=0)#
      integer, parameter :: LWEPAR = 34 !  # layers that have clouds (LWEPAR < L_)
#else
      integer :: L_     !  # of CTM layers, set at run-time
      integer :: L1_    !  L_+1 = # of CTM layer edges (radii)
      integer :: L2_    !  L_+2 = total # of layer edges counting top (TAU=0)
      integer :: LWEPAR !  # layers that have clouds (LWEPAR < L_)
#endif

!-----------------------------------------------------------------------
! General configuration settings
!-----------------------------------------------------------------------

      ! ATM0: Option for spherical corrections: 0=flat 1=sphr 2=refr 3=geom
      ! standard spherical atmosphere is OK, but refractive and geometric can be used
      ! see  Prather & Hsu, 2019, A round Earth for climate models,
      ! PNAS, 116(39): 1933019335
      integer :: ATM0

      ! ATAU: Factor increase in cloud optical depth (OD) from layer to next below
      real*8  :: ATAU

      ! ATAU0: Minimum cloud OD in uppermost inserted layer
      real*8  :: ATAU0

      ! USEH2OUV: Whether to use H2O UV absorption
      logical :: USEH2OUV

      ! Cloud flag options:
      !       CLDFLAG = 1  :  Clear sky J's
      !       CLDFLAG = 2  :  Averaged cloud cover
      !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
      !       CLDFLAG = 4  :  ****not used (old direct beam avg)
      !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
      !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of QCA bin)
      !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds in layer within each Q-bin)
      !       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)      
      integer :: CLDFLAG

      ! String description of comments to print to log
      character*25, dimension(8), parameter :: TITCLD =  &
         ['clear sky - no clouds       ', &
          'avg cloud cover             ', &
          'avg cloud cover^3/2         ', &
          'ICAs - avg direct beam*VOID*', &
          'ICAs - random N ICAs        ', &
          'QCAs - midpt of bins        ', &
          'QCAs - avg clouds in bins   ', &
          'ICAs - use all ICAs***      ']

      ! CLDCOR: Cloud decorellation between max-overlap blocks (0.00 = random)
      ! Only used for cloud flags 5 and above
      real*8  :: CLDCOR

      ! LNRG: Number of max-overlap blocks, can be 0 (max-ran @ gaps) or 3 (alt blocks)
      ! Only used for cloud flags ___, etc
      integer :: LNRG

      ! Dimensions of readin spec - NWBIN can zero out strat-wavels
      ! NWBIN = 18 = std full Fast-J, for TROP-ONLY =12 (0% err in trop, 33% savings)
      ! = 08 big savings, but 1-2% error in J-O2 and J-OCS in upper trop
      integer  NWBIN

      ! what is this?
      integer, parameter :: NSBIN = 27

!------------------------------------------------------------------------------
! Additional parameters
!------------------------------------------------------------------------------

      ! JVN_ :  max # of J-values
      integer, parameter :: JVN_ = 200

      ! AN_ :  max # of FJX aerosols in layer (needs NDX for each)
#ifdef MODEL_GEOSCHEM
      integer, parameter :: AN_=37
#elif MODEL_STANDALONE
      integer, parameter :: AN_=25
#endif

      !  Dimensions for wavelength data tables, fjx_spec.dat & others
      
      ! used for table dimensions for cross-sections
      integer, parameter ::  WX_= 18
      
      ! used for table dimensions of broad-bands thru IR
      integer, parameter ::  SX_= 27
      
      ! W_   = dim = no. of Fast-J Wavelength bins:
      ! Currenly 18, should be same as WX_ for now
      ! A TROP-ONLY calc. is done by setting the read-in NWBIN to 8 or 12, zeros the FL fluxes
      integer, parameter ::  W_=18

      !------------------------------------------------------------------------------
      ! this section is hard wired to skip bins 19:27, used for just chemistry & J-values
      ! S_   = dim = number of wavelength sub-bins INCLUDING the Solar-J extensions
      ! (RRTMG value = 27)
      integer, parameter ::  S_=W_  ! does NOT do calculations for w>18, simple J's only
      !     else: integer, parameter ::  S_=27  != # of broad bands in entire calculation,
      !     27 includes the 9 RRTMG bins but Cloud-J only calculates heating for clouds and
      !     aerosols in these bins need full Solar-J to get water vapor lines & total heating
      !     rates.  then S_ = sum(NGC(1:27) = 100

      integer, parameter :: NW1=1

      integer, parameter :: NW2=W_

      integer, parameter :: NS1=1

      integer, parameter :: NS2=S_
      
      integer, parameter ::  W_r = S_-W_  ! # of bins that is added for solar IR on top of W_

      integer, parameter, dimension(27) :: NGC = &
          (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &    ! these are Cloud-J, no sub-bins
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &    ! Note that S_= W_, means no near-IR    
            1, 1, 1, 1, 1, 1, 1/)
      ! SJ (/1,  1,  1, 1, 1,  1, 1, 1, 1,  1, &  !these are used in Solar-J for RRTMG
      ! SJ   1,  1,  1, 1, 1,  1, 1, 5, 10, 2, &
      ! SJ   10, 10, 8, 8, 12, 6, 12 /)
      !

      ! X_   = dim = max no. of X-section data sets (input data)
#ifdef MODEL_GEOSCHEM
      integer, parameter ::  X_=123
#elif MODEL_STANDALONE
      integer, parameter ::  X_=72
#endif

      ! A_   = dim = max no. of Aerosol Mie sets (input data) not including
      !        clouds and SSA
#ifdef MODEL_GEOSCHEM
      integer, parameter ::  A_=56
#elif MODEL_STANDALONE
      integer, parameter ::  A_=40
#endif

      ! SSA_ & GGA_ = dim = no. of strat sulfate aerosol types (input data)
      integer, parameter ::  SSA_=18, GGA_=15

      ! C_   = dim = no. of cld-data sets (input data):
      ! liquid-water, irregular-ice, hexagonal ice
      integer, parameter ::  C_=3

      ! CR_   = dim = no. of effective radii in each cld-data sets
      integer, parameter ::  CR_=6

      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*(L_+1) + 1`+ 2*added-cld-layers
      integer, parameter ::  N_=601

      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      integer, parameter ::  M_=4

      ! M2_ = 2*M_ = 8, replaces MFIT
      integer, parameter ::  M2_=2*M_

!-----------------------------------------------------------------------
! 4 Gauss pts = 8-stream
!------------------------------------------------------------------------------

      real*8, DIMENSION(M_), parameter  ::  &
                         EMU = [.06943184420297d0, .33000947820757d0, &
                                .66999052179243d0, .93056815579703d0]
      real*8, DIMENSION(M_), parameter  :: &
                         WT  = [.17392742256873d0, .32607257743127d0, &
                                .32607257743127d0, .17392742256873d0]

!-----------------------------------------------------------------------
! Radiation Field, Cloud Cover & Other fixed parameters
!------------------------------------------------------------------------------

      ! MASFAC: Conversion factor for pressure to column density
      real*8, parameter:: MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

      ! HeatFac_: convert watt/m2 to K/day
      real*8, parameter:: HeatFac_ = 86400.d0*9.80616d0/1.00464d5


!------------------------------------------------------------------------------
! Variables in file 'FJX_spec.dat' (RD_XXX)
!------------------------------------------------------------------------------

      ! # fast-JX J-values based on number of cross-sections read in
      integer NJX

      ! WL: Centres of wavelength bins - 'effective wavelength'  (nm)
      real*8  WL(SX_)

      ! WBIN: Boundaries of wavelength bins                  (microns)
      real*8  WBIN(SX_+1)

      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      real*8  FL(SX_)

      ! FW: Solar flux in W/m2
      real*8  FW(SX_)

      ! FPAR: PAR quantum action spectrum
      ! NOTE: renamed from FP to avoid conflict with flexible precision naming in GEOS-Chem
      real*8  FPAR(SX_)

      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      ! SJ had this set ant S_+1, it should not needed
      real*8  QRAYL(SX_)

      ! SJSUB:  intended for breakdown of the super-bins (1:27) into
      ! smaller sub-bins.
      real*8  SJSUB(SX_,16)

      ! QH2O: H2O UV-blue cross-sections (290-350 nm)
      real*8  QH2O(WX_)
      
      ! KDOKR: Set for RRTMG =18+82 max set at 100
      integer KDOKR(100)

      ! LDOKR: Set for RRTMG =18+82 max set at 100
      integer LDOKR(100)

      integer NSJSUB(SX_)

      ! QO2: O2 cross-sections
      real*8  QO2(WX_,3)

      ! QO3: O3 cross-sections
      real*8  QO3(WX_,3)

      ! Q1D: O3 => O(1D) quantum yield
      real*8  Q1D(WX_,3)

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      real*8  QQQ(WX_,3,X_)

      ! TQQ: Temperature for supplied cross sections
      real*8  TQQ(3,X_)

      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      integer LQQ(X_)

      ! TITLEJX: Title (short & long) for supplied cross sections,
      ! from 'FJX_spec.dat'
      CHARACTER*6  TITLEJX(X_)

      CHARACTER*16 TITLEJL(X_)

      ! SQQ: Flag (pressure or temperature tables) for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1  SQQ(X_)

!------------------------------------------------------------------------------
! Variables in file 'FJX_scat-aer.dat' (RD_MIE)
!------------------------------------------------------------------------------

      ! TITLAA: Aerosol Mie Titles
      character*12  TITLAA(A_)

      ! QAA: Aerosol scattering phase functions
      real*8  QAA(5,A_)

      ! WAA: 5 Wavelengths for the supplied phase functions
      real*8  WAA(5,A_)

      ! PAA: Phase function: first 8 terms of expansion
      real*8  PAA(8,5,A_)

      ! RAA: Effective radius associated with aerosol type
      real*8  RAA(A_)

      ! SAA: Single scattering albedo
      real*8  SAA(5,A_)

      ! DAA: density (g/cm^3)
      real*8  DAA(A_)

      ! NAA: Number of categories for scattering phase functions
      integer NAA

!------------------------------------------------------------------------------
! Variables in file 'FJX_scat-cld.dat' (RD_CLD)  ***major update for v75
!------------------------------------------------------------------------------

      ! NCC: Number of categories for cloud scattering phase functions,
      ! MCC: no. of R_eff
      integer NCC,MCC

      ! TITLCC: Cloud type titles
      character*12  TITLCC(C_)

      ! RCC: Effective radius associated with cloud type:
      ! should be indep of S-bin
      real*8  RCC(CR_,C_)

      ! GCC: Effective geometric cross section:  should be indep of S-bin
      real*8  GCC(CR_,C_)

      ! DCC: density (g/cm^3):  should be indep of S-bin, and eff radius
      real*8  DCC(C_)

      ! QCC: Cloud Q-ext
      real*8  QCC(SX_,CR_,C_)

      ! WCC: Wavelengths for supplied phase functions, should be std S-bins
      real*8  WCC(SX_,C_)

      ! SCC: Single scattering albedo
      real*8  SCC(SX_,CR_,C_)

      ! PCC: Phase function: first 8 terms of expansion
      real*8  PCC(8,SX_,CR_,C_)

!------------------------------------------------------------------------------
! Variables in file 'FJX_scat-ssa.dat' (RD_SSA)
!------------------------------------------------------------------------------

      ! NSS: Number of categories for Stratospheric Sulfate Aerosol scattering
      ! phase functions
      integer NSS

      ! TITLSS: Cloud type titles
      character*12  TITLSS(SSA_)

      ! RSS: Effective radius associated with cloud type
      real*8  RSS(SSA_)

      ! GSS: Effective geometric cross section
      real*8  GSS(SSA_)

      ! DSS: density (g/cm^3)
      real*8  DSS(SSA_)

      ! TSS: temperature (K)
      real*8  TSS(SSA_)

      ! WSS: weight percent sulfuric acid (%)
      real*8  WSS(SSA_)

      ! QSS: Q-ext      ----begin wavelength dependent quantities
      real*8  QSS(SX_,SSA_)

      ! SSS: Single scattering albedo
      real*8  SSS(SX_,SSA_)

      ! PSS: Phase function: first 8 terms of expansion
      real*8  PSS(8,SX_,SSA_)

!------------------------------------------------------------------------------
! Variables in file 'FJX_scat-geo.dat' (RD_GEO)
!------------------------------------------------------------------------------

      ! NGG: # of categories for Strat Sulf Aerosol scattering phase functions
      integer NGG

      ! RGG: Effective radius associated with cloud type
      real*8  RGG(GGA_)

      ! DGG: density (g/cm^3)
      real*8  DGG(GGA_)

      ! QGG: Q-ext      ----begin wavelength dependent quantities
      real*8  QGG(SX_,GGA_)

      ! SGG: Single scattering albedo
      real*8  SGG(SX_,GGA_)

      ! PGG: Phase function: first 8 terms of expansion
      real*8  PGG(8,SX_,GGA_)

!------------------------------------------------------------------------------
! Variables in file 'FJX_scat-UMa.dat' (RD_UM) (UMich aerosol scattering data)
!------------------------------------------------------------------------------

      ! WMM: U Michigan aerosol wavelengths
      real*8  WMM(6)

      ! UMAER: U Michigan aerosol data sets
      real*8  UMAER(3,6,21,33)

!------------------------------------------------------------------------------
! Variables in file 'atmos_std.dat' (RD_PROF) and 'atmos_h2och4.dat' (RD_TRPROF)
! NOTE: only used in Cloud-J standalone
!------------------------------------------------------------------------------

      ! layer dim. in reference profiles
      integer, parameter ::  LREF=51

      ! latitude dim. in reference profiles
      integer, parameter ::  JREF=18

      ! T, O3, H2O, CH4
      ! NOTE: ref profiles added underscore _ because TREF used in RRTMG_SW

      real*8, DIMENSION(LREF,JREF,12) :: T_REF

      real*8, DIMENSION(LREF,JREF,12) :: O_REF

      real*8, DIMENSION(LREF,JREF,12) :: H2O_REF

      real*8, DIMENSION(LREF,JREF,12) :: CH4_REF

!------------------------------------------------------------------------------
! Reference monthly zonal mean profiles for GEOMIP SSA 'atmos_geomip.dat'
!------------------------------------------------------------------------------

      ! layer dim. in reference profiles
      integer, parameter ::  LGREF=19

      ! R = Reff (microns), X = micro-g-H2SO4/kg-air
      real*8, dimension(64,19,12) :: R_GREF

      real*8, dimension(64,19,12) :: X_GREF

      !  A = 4 pi R^2 = microns^2/cm^3
      real*8, dimension(64,19,12) :: A_GREF

      real*8, dimension(64)       :: Y_GREF

      real*8, dimension(19)       :: P_GREF

!------------------------------------------------------------------------------
! Variables in file 'FJX_j2j.dat' (RD_JS_JX)
!------------------------------------------------------------------------------

      ! multiplication factor for fast-JX calculated J
      real*8  JFACTA(JVN_)

      ! index arrays that map Jvalue(j) onto rates
      integer JIND(JVN_)

      ! number of Photolysis reactions in CTM chemistry, NRATJ <= JVN_
      integer NRATJ

      !label of J-value used to match w/FJX J's
      character*6 JVMAP(JVN_)

      ! label of J-value used in the chem model
      character*50 JLABEL(JVN_)

      ! Branches for photolysis species
      integer BRANCH(JVN_)

      ! Names of photolysis species
      character*10 RNAMES(JVN_) 

!------------------------------------------------------------------------------
! Cloud overlap parameters
!-----------------------------------------------------------------------

      ! # of quantized cloud fration bins
      integer, parameter :: CBIN_ = 10

      ! Max # of indep colm atmospheres
      integer, parameter :: ICA_ = 20000

      ! # of cloud quadrature bins (4)
      integer, parameter :: NQD_ = 4

      real*8,  parameter ::  CPI    = 3.141592653589793d0

      real*8,  parameter ::  C2PI   = 2.d0*CPI

      real*8,  parameter ::  CPI180 = CPI/180.d0

      real*8,  parameter ::  G0     = 9.80665d0

      real*8,  parameter ::  G100 = 100.d0/G0

!------------------------------------------------------------------------------
! Data to set up the random number sequence for use in cloud-JX
!------------------------------------------------------------------------------

      ! dimension for random number
      integer, parameter :: NRAN_ = 10007

      ! Random number set
      real*4   RAN4(NRAN_)

      END MODULE CLDJ_CMN_MOD
