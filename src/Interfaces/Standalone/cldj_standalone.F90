!>>>>>>>>  Cloud-J version 8.0c - H2O uv abs allowed!!! (5/2023)
!      Includes overall cleanup of dead ends and removal of Solar-J parallels
!      Now does average clouds in example, IWP error corrected.
!

      program standalone

      USE CLDJ_CMN_MOD
      USE CLDJ_ERROR_MOD
      USE CLDJ_INIT_MOD
      USE CLDJ_FJX_SUB_MOD
      USE CLDJ_SUB_MOD, ONLY : CLOUD_JX
      USE CLDJ_FJX_OSA_MOD

      implicit none
!---------------key params in/out of CLOUD_J-------------------------
      logical                    :: LPRTJ, LDARK
      integer                    :: IRAN, RC
      integer                    :: JVNU,ANU,L1U
      integer                    :: NICA,JCOUNT
      real*8                     :: U0,SZA,SOLF
      real*8,  dimension(L2_  )  :: PPP,ZZZ
      real*8,  dimension(L1_  )  :: TTT,HHH,DDD,RRR,OOO,CCC, DELZ
      real*8,  dimension(L1_  )  :: O3,CH4,H2O, OD18
      real*8,  dimension(S_+2,L1_):: SKPERD
      real*8,  dimension(6)       :: SWMSQ
      real*8,  dimension(L1_)    :: CLF,LWP,IWP,REFFL,REFFI
      integer, dimension(L1_)    :: CLDIW
      real*8,  dimension(L1_,AN_):: AERSP
      integer, dimension(L1_,AN_):: NDXAER
      real*8,  dimension(L_,JVN_) :: VALJXX
      real*8,  dimension(5,S_)    :: RFL
      real*8,  dimension(NQD_)    :: WTQCA

!-------------local use-----------------------
      integer :: NSZA,MONTH,ILAT
! beware the OSA code uses single R*4 variables
! old
      !real*4  :: OWAVEL,OWIND,OCHLR,OSA_dir(5)
      !real*8  :: YLAT,PSURF,ALBEDO(5),WIND,CHLR
      !real*8, dimension(L2_) :: ETAA,ETAB, ZOFL,RI,TI,AER1,AER2,PPPX
      !integer,dimension(L2_) :: NAA1,NAA2
! new
       real*8  :: ALBEDO(5),WIND,CHLR, WAVEL,ANGLES(5),OSA_dir(5)
       real*8  :: YLAT,PSURF, HHH0,ZKM
 ! used for diagnostics here, regarding O(1D) and OH production
       real*8  :: dAIR,dO3,dH2O,J1D,K1O2,K1N2,K1H2O,K1air,KCH4,dO1d,POH,LCH4,HRuv,HRir
 
       real*8, dimension(L2_) :: ETAA,ETAB
       real*8, dimension(L1_) :: RHINP,TINP,ZOFL, AER1,AER2
       integer,dimension(L1_) :: NAA1,NAA2
! end new
      real*8, dimension(L_)  :: WLC,WIC
      real*8, dimension(LWEPAR) :: CLDFRW,CLDIWCW,CLDLWCW
      real*8  SCALEH,CF,PMID,PDEL,ZDEL,ICWC,F1
      integer I,J,K,L,N
      integer LTOP, NJXX, NLEVELS
      character*6,  dimension(JVN_)  ::  TITLJXX
      character*64                   ::  thisloc
      logical :: amIRoot

      integer, dimension(10) ::  &
        SZAscan = [ 0, 30, 60, 80, 86, 88, 90, 92, 94, 96]

      real*8  :: ATAU_in
      real*8  :: ATAU0_in
      real*8  :: CLDCOR_in
      integer :: NWBIN_in
      integer :: LNRG_in
      integer :: ATM0_in
      integer :: CLDFLAG_in
      logical :: Use_H2O_UV_Abs

!Notes:  L_ (param), sometimes passed as LU (integer), = numer of atmosphere layers
!        L_+1 = top layer added from P(top-of-atmos) to P=0, needs air at least.

!!!! IN variables !!!!
! U0    cosine(SZA)
! SZA   solar zenith angle (0:180 degrees)
! RFL(1:5,1:S_) surface albedo for the 4 quadrature angles, EMU(1:4), plus U0
!      for each photolysis wavelength (1:W_) plus IR RRTMG wavelengths W_+1:S_ IF USED
! SOLF  solar flux factor for seasonal change in sun-earth distance: 1 +- 3.3%
! LPRTJ (logical) .T.= triggers detailed print to stdout of J-values and heating rates
! PPP(1:L_+2)  pressure (hPa) of layer edges, PPP(L_+2) = 0 (should be)
! ZZZ(1:L_+2)  altitude (cm) of layer edges, calculated from hydrostatic equilibrium
! TTT(1:L_+1)  temperature (K) of layer
! HHH(1:L_+1)  H2O molecules per cm2 in each layer (column)
! DDD(1:L_+1)  Air molecules per cm2 in each layer (column)
! RRR(1:L_+1)  Relative Humidity (0 - 1),  used in sub OPTICA & OPTICM
! OOO(1:L_+1)  O3 molecules per cm2 in each layer (column)
! CCC(1:L_+1)  CH4 molecules per cm2 in each layer (column)
! LWP(1:L_+1)  liquid water path (g/m2) in each layer
! IWP(1:L_+1)  ice water path (g/m2) in each layer
! REFFL(1:L_+1) effective radius (microns) of liquid water particles in each layer
! REFFI(1:L_+1) effective radius (microns) of ice particles in each layer
! CLDF(1:L_+1)  cloud fraction (0-1) in each layer

!!!!  IN & OUT !!!!
! CLDCOR   read in (fjx_init) = 0.33 = cloud decorellation between max-overlap blocks
!           reset to 0 in ICA_NR if LNRG =0 or 3

!!!! IN !!!!
! CLDIW  index for each layer: 0=no cloud, 1=liquid loud only, 2=ice cloud only, 3=liquid+ice cloud mix
! AERSP(L_+1,AN_)  aerosol path (g/m2) in each layer of aerosol type, type depends on CTM/CCM aerosol package
! NDXAER(L_+1,AN_)  aerosol index for looking up opacity tables, locally in FJX_SUB this becomes NAER
!          NAER is used to invoke a range of aerosol options and should be customized by the CTM/CCM
!          =1:2  standard stratospheric sulfate aerosols in OPTICS, 2 = larger, Pinatubo-like aerosols
!          =1001:1015  selects GEOMIP like aerosols in OPTICG
!          =2:999  looks up aerosols in OPTICA, there are not 999 values in current tables
!          <0 selects from U. MIchigan aerosol tables (includes relative humidity)
! L1U = L1_ = L_+1
! ANU = AN_, max dimension size for number of aerosl types.
! NJXU = NJX_, max # of J-values


!!!! OUT variables !!!!
! VALJXX(L_,JVN_)  weighted average of J-values (/s) in each layer (not layer edge)
! SKPERD(S_+2,L_+1) heating (K/day) for each wavelength bin and each layer.  The 2 extra wavelength bins are 1=UV-VIS & 2=IR
! SWMSQ(6) net total heating (e.g. W/m2 = WMSQ):
!            1=inc TOTAL   2=rfl outtop  3=abs in atm  4=abs at srf  5=PAR direct  6=PAR diffus
! OD18(L_+1) average (over QCAs) of W=18 optical depth in each layer (~600 nm, includes clouds).

!* IRAN  integer seed for random numbers (e.g., =1) used only with CLDFLAG=5
! * = should not be used in regular CTM/CCM implementation of Cloud-J

!!!! IN !!!!
!!!! OUT !!!!
!# NICA   number of independent column atmospheres generated by cloud overlap, must be .le. ICA_=20,000 currently
!# JCOUNT number of J-value calculations for each Cloud-J call
!# LDARK  = .T. = no J-values, atmosphere considered dark
!# WTQCA(4)  weight of each of up-to-4 quadrature column atmospheres used to calculate J's.  sum(WtQC(1:4)) = 1
! # these are for diagnostics of the Cloud-J calculation, not needed unless there is trouble
      
      write(6,'(a)') '>>>begin Cloud-J v8.0 Standalone'

      RC = CLDJ_SUCCESS
      NLEVELS = 57
      ANU = AN_
      JVNU = JVN_
      L1U = L1_
      amIRoot = .true.
      thisloc = 'standalone program in cldj_standalone.F90'

      ! Global vars previously read from config file but now passed from driver
      !!!!  Key variables/parameters used in Cloud-J that are passed to initialization routine
      ! ATAU  factor increase in cloud optical depth (OD) from layer to next below
      ! ATAU0 minimum cloud OD in uppermost inserted layer
      ! CLDCOR correlation of cloud overlap between blocks (0.00 = random)
      ! NWBIN  # wavelength bins in uv-vis for J's <<< can be 8 or 12 also for trop only J's
      ! LNRG  correlated cloud overlap (bumber of max-overlap blocks) =0=Max-Ran at gaps,
      !       =3=Max-Ran in 3 alt blocks, 6= six overlap blocks
      ! ATM0  spherical correction: 0=flat, 1=sphr, 2=refr, 3=geom
      ! CLDFLAG see comments below
      ! CLDFLAG - type of cloud overlap parameterization:
      !       CLDFLAG = 1  :  Clear sky J's
      !       CLDFLAG = 2  :  Averaged cloud cover
      !       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
      !       CLDFLAG = 4  :  ****not used
      !       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
      !       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
      !       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
      !       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)
      ! Use_H2O_UV_Abs  whether to use UV absorption by water vapor
      ATAU_in    = 1.050d0
      ATAU0_in   = 0.005d0
      CLDCOR_in  = 0.33d0
      NWBIN_in   = 18 
      LNRG_in    = 06
      ATM0_in    = 1
      CLDFLAG_in = 7
      Use_H2O_UV_Abs = .true.
      
!---read in & setup fast-JX data and parameters:   single call at set up
!-----------------------------------------------------------------------      
      call INIT_CLDJ (amIRoot,'./tables/',NLEVELS,LWEPAR,TITLJXX,JVNU,  &
           ATAU_in, ATAU0_in, NWBIN_in, CLDFLAG_in, CLDCOR_in, LNRG_in, &
           ATM0_in, Use_H2O_UV_Abs, NJXX,RC )
      if ( RC /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR_STOP( 'Failure in INIT_CLDJ', thisloc )
      endif
!-----------------------------------------------------------------------

!-- this is column-specific data, usually input from the running CTM/CCM model                                                  !
!--P, T, Cld & Aersl profiles, simple test input case
      open (77,file='tables/atmos_PTClds.dat',status='old',err=91)
        read (77,*)
        read (77,'(2i5)') MONTH, ILAT               ! for monthly 2D climatologies
          YLAT = ILAT
        read (77,'(f5.0)') PSURF
        read (77,'(f5.2)') ALBEDO(5)
        read (77,'(4f5.2)') ALBEDO(1),ALBEDO(2),ALBEDO(3),ALBEDO(4)
        read (77,'(4f5.2)') WIND,CHLR
          write(6,'(a,2i5,5x,a,i5)') 'Atmosphere:', L_, LWEPAR, 'L_ / LWEPAR', L1_
          write(6,'(a,f10.4)') 'P surface', PSURF
          write(6,'(a,3i4)') 'MONTH/ LAT',MONTH,ILAT
          write(6,'(a,5f8.4)') 'Albedos 1:4 & 5=SZA', ALBEDO
          write(6,'(a,2f8.3)') 'OSA: wind & chlor-a',WIND,CHLR
        read (77,*)
       do L = 1,L1_
        read (77,'(i3,1x,2f11.7,2x,f5.1,f5.2,f11.2,2(f7.3,i4))') &
                      J,ETAA(L),ETAB(L),TINP(L),RHINP(L),ZOFL(L) &
                     ,AER1(L),NAA1(L),AER2(L),NAA2(L)
       enddo
        read (77,*)
       do L = LWEPAR,1,-1
        read (77,'(i3,1p,e14.5,28x,2e14.5)') &
                      J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
       enddo
      close(77)

      ETAA(L2_) = 0.d0
      ETAB(L2_) = 0.d0
      do L = 1,L2_
       PPP(L) = ETAA(L) + ETAB(L)*PSURF
      enddo

!---sets climatologies for O3, T, D & Z
!-----------------------------------------------------------------------
      call ACLIM_FJX (YLAT,MONTH,PPP, TTT,O3,CH4, L1_)
      if ( RC /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR_STOP( 'Failure in ACLIM_FJX', thisloc )
      endif
!-----------------------------------------------------------------------
      do L = 1,L1_
       TTT(L) = TINP(L) !! keep climatology read-in T & RHum, use O3 from ACLIM_FJX
       RRR(L) = RHINP(L)
      enddo
       ZZZ(1)  = 16.d5*log10(1013.25d0/PPP(1))        ! zzz in cm
      do L = 1,L_
       DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
!!! geopotential since assumes g = constant
       SCALEH      = 1.3806d-19*MASFAC*TTT(L)
       ZZZ(L+1) = ZZZ(L) -( log(PPP(L+1)/PPP(L)) * SCALEH )
       OOO(L) = DDD(L)*O3(L)*1.d-6
       CCC(L) = DDD(L)*CH4(L)*1.d-9
      enddo
      L = L_+1
       ZZZ(L+1)= ZZZ(L) + ZZHT
       DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
       OOO(L) = DDD(L)*O3(L)*1.d-6
       CCC(L) = DDD(L)*CH4(L)*1.d-9


! set up H2O profile, should be using q (kg-H2O/kg-air) from wind fields,
!        HHH(L) = QL(L)*DDD(L)*(28.987/18.)
! for NOW, just estimate to get 1.5e23 column with scale ht of about 2.2 km, use lower Z edge
        HHH0 = 0.030d0
      do L = 1,L_+1
        HHH(L) = DDD(L)*max(HHH0*exp(-ZZZ(L)/2.2d5), 2.0d-6)  ! produces colm H2O of 1.4e23
      enddo
!      write(6,'(i5,1p, 3e12.4)') (L, ZZZ(L),HHH(L),DDD(L), L=1,L1_)

      do L = 1,L_+1
        DELZ(L) = ZZZ(L+1) - ZZZ(L) ! to calc densities (/cm3)
      enddo


!!! set up clouds and aerosols
       AERSP(:,:)  = 0.d0
       NDXAER(:,:) = 0
      do L = 1,L_
       NDXAER(L,1) = NAA1(L)
       AERSP(L,1)  = AER1(L)
       NDXAER(L,2) = NAA2(L)
       AERSP(L,2)  = AER2(L)
      enddo
       LTOP  = LWEPAR
      if (maxval(CLDFRW) .le. 0.005d0) then
        IWP(:) = 0.d0
        REFFI(:) = 0.d0
        LWP(:) = 0.d0
        REFFL(:) = 0.d0
      endif
      do L = 1,LTOP
          CLDIW(L) = 0
          CF  = CLDFRW(L)
        if (CF .gt. 0.005d0) then
          CLF(L) = CF
          WLC(L) = CLDLWCW(L) / CF
          WIC(L) = CLDIWCW(L) / CF
!  CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
         if (WLC(L) .gt. 1.d-11) CLDIW(L) = 1
         if (WIC(L) .gt. 1.d-11) CLDIW(L) = CLDIW(L) + 2
        else
          CLF(L) = 0.d0
          WLC(L) = 0.d0
          WIC(L) = 0.d0
        endif
      enddo
!---derive R-effective for clouds:  the current UCI algorithm - use your own
      do L = 1,LTOP
!---ice clouds
        if (WIC(L) .gt. 1.d-12) then
            PDEL = PPP(L) - PPP(L+1)
            ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! m
          IWP(L) = 1000.d0*WIC(L)*PDEL*G100    ! g/m2
          ICWC =        IWP(L) / ZDEL          ! g/m3
          REFFI(L) = 164.d0 * (ICWC**0.23d0)
        else
          IWP(L) = 0.d0
          REFFI(L) = 0.d0
        endif
!---water clouds
        if (WLC(L) .gt. 1.d-12) then
            PMID = 0.5d0*(PPP(L)+PPP(L+1))
            PDEL = PPP(L) - PPP(L+1)
          F1   = 0.005d0 * (PMID - 610.d0)
          F1   = min(1.d0, max(0.d0, F1))
          LWP(L) = 1000.d0*WLC(L)*PDEL*G100     ! g/m2
          REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
        else
          LWP(L) = 0.d0
          REFFL(L) = 0.d0
        endif
      enddo
      do L = 1,LTOP
        CLDFRW(L) = CLF(L)
      enddo
!!!  end of atmosphere setup


!!! begin call to Cloud_J

      do I=1,3    !!!! begin of SZA scan
       NSZA = SZAscan(I)
       SZA = NSZA

      do L = 1,LTOP
        CLF(L) = CLDFRW(L)
      enddo
      IRAN = 1
      SOLF = 1.d0
      U0 = cos(SZA*CPI180)

      ANGLES(1) = EMU(1)
      ANGLES(2) = EMU(2)
      ANGLES(3) = EMU(3)
      ANGLES(4) = EMU(4)
      ANGLES(5) = U0

!!!! subroutine OSA = Ocean Surface Albedo, in OSA_SUB_MOD.f90
!!!! IN !!!!
! WAVEL         center of wavelength bin (nm)
! WIND          surface wind speed (m/s)
! CHLR          chlorophyl-A (mg/m3)
! ANGLES(1:5)   (1:4) = the 4 RT quad angles (~cos(86,71,48 & 22 deg))  (5) = U0 = cos(SZA)
!!!! OUT !!!!
! OSA_dir       (1:5) Ocean Surface Albedo (Lambertian) for each of the 5 ANGLES
!
!!!! NOTE that Cloud-J supplies OSA, but NOT an equivalent Land Surface Albedo code

      do K = 1,NS2
        WAVEL = WL(K)   !! in nm
        call FJX_OSA(WAVEL,WIND,CHLR,ANGLES, OSA_dir)
         do J = 1,5
          RFL(J,K) = OSA_dir(J)
! this temporary overwrite the OSA with the readin values above just for these tests
            RFL(J,K) = ALBEDO(J)
        enddo

        write(6,'(a,f6.3,f6.1,5f6.3)') ' U0/wvl/OSA: ',U0,WAVEL,OSA_dir
     enddo

      LPRTJ = .true.
      if (LPRTJ) then
          write(6,'(a,f8.3,3f8.5)')'SZA SOLF U0 albedo' &
                ,SZA,SOLF,U0,RFL(5,18)
        call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
        if ( RC /= CLDJ_SUCCESS ) then
           call CLOUDJ_ERROR_STOP( 'Failure in JP_ATM0', thisloc )
        endif
          write(6,*) ' wvl  albedo u1:u4 & u0'
        do K=1,NS2
          write(6,'(i5,f8.1,5f8.4)') K,WL(K), (RFL(J,K), J=1,5)
        enddo
      endif

      SKPERD(:,:)=0.d0
      SWMSQ(:)= 0.d0
      OD18(:) =0.d0
      WTQCA(:)= 0.d0

!=======================================================================
       call CLOUD_JX (U0,SZA,RFL,SOLF,LPRTJ,PPP,ZZZ,TTT,HHH,DDD,       &
               RRR,OOO,CCC, LWP,IWP,REFFL,REFFI, CLF,CLDCOR,CLDIW,    &
               AERSP,NDXAER,L1U,ANU,JVNU, VALJXX,SKPERD,SWMSQ,OD18,    &
               IRAN,NICA, JCOUNT,LDARK,WTQCA,RC)
       if ( RC /= CLDJ_SUCCESS ) then
          call CLOUDJ_ERROR_STOP( 'Failure in CLOUD_JX', thisloc )
       endif
!=======================================================================



! Example of summary diagnostics from cldj_standalone.F90:  calc. O(1D) and OH primary production
      N=7
      write(N,'(a,2i5)') ' v8.0 CLDFLAG/NSZA=', CLDFLAG,NSZA
      write(N,*) ' LDARK WTQCA',LDARK,WTQCA
      write(N,'(1x,a,72(a6,3x))') 'L=  ',TITLEJX(3)
      write(N,'(a)') 'L Z p T [M] [O3] [H2O] [O1D] P-OH L-CH4 HR(K/day):vis&nir'
      do L = L_,1,-1
         ZKM = (ZZZ(L)+ZZZ(L+1))*0.5d-5
         dAIR = DDD(L)/DELZ(L)
         dO3 =  OOO(L)/DELZ(L)
         dH2O = HHH(L)/DELZ(L)
         J1d = VALJxx(L,3)
          K1O2 = 3.30d-11*exp(+55/TTT(L))
          K1N2 = 2.15d-11*exp(+110/TTT(L))
          K1H2O= 1.63d-10*exp(+60/TTT(L))
          K1air = 0.21*k1O2 + 0.78*K1N2
          KCH4 = 2.45d-12*exp(-1775/TTT(L))
         dO1D = J1D*dO3/(K1air*dAIR + k1H2O*dH2O)
         POH = 2*K1H2O*dO1d*dH2O*86400.d9/dAIR    !  ppb/day
         LCH4 = POH*KCH4     ! sort of, POH wtd by k_CH4+OH
         HRuv = SKPERD(S_+1,L)
         HRir = SKPERD(S_+2,L)
         write(N,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,7e10.3,0p,2f8.3)') &
            L,ZKM,PPP(L),TTT(L),dAIR,dO3,dH2O,J1D,dO1D,POH,LCH4,HRuv,HRir
         
      enddo

      write(N,'(a)') 'heating rate profiles in K/day v7.6  180-778nm '
      write(N, '(a4, 32f7.1)')'wvl ',(WL(K),K=NW1,NW2)
      do L = L_,1,-1
         write(N,'(i4,32f7.2)') L,(SKPERD(K,L), K=NW1,NW2)
      enddo
      write(N,'(a)') 'Fast-J  v8.0: Solar fluxes (W/m2)--'
      write(N,'(a11,f12.4)')    ' inc TOTAL ',SWMSQ(1)
      write(N,'(a11,f12.4)')    ' rfl outtop',SWMSQ(2)
      write(N,'(a11,f12.4)')    ' abs in atm',SWMSQ(3)
      write(N,'(a11,f12.4)')    ' abs at srf',SWMSQ(4)
      write(N,'(a11,1p,e12.4)') ' PAR direct',SWMSQ(5)
      write(N,'(a11,1p,e12.4)') ' PAR diffus',SWMSQ(6)

!    other example of remapping j-values onto the CTM table
!---map the J-values from fast-JX onto CTM (ZPJQUAD) using JIND & JFACTA
!--- from the 'FJX_j2j.dat' tables
!      do J = 1,NRATJ
!         JP = JIND(J)
!         if (JP .gt. 0) then
!            do L = 1,L_
!               ZPJQUAD(L,J) = ZPJQUAD(L,J) + VALJXX(L,JP)*JFACTA(J)
!          enddo
!        endif
!      enddo


      enddo !!!! end of SZA scan

      goto 92
   91 stop 'error in opening .dat file'
   92 stop
      end
 
 ! not used but important connection for running in global model
 
 !-----------------------------------------------------------------------
       subroutine FJX_TOD(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
 !-----------------------------------------------------------------------
 !  time of day, user example:  NOT USED HERE
 !     GMTIME = UT for when J-values are wanted
 !           (for implicit solver this is at the end of the time step)
 !     NDAY   = integer day of the year (used for solar lat and declin)
 !     YGRDJ  = laitude (radians) for grid (I,J)
 !     XGDRI  = longitude (radians) for grid (I,J)
 !
 !     SZA = solar zenith angle in degrees
 !     COSSZA = U0 = cos(SZA)
 !     SOLFLX = sun-earth distance factor for solar flux
 !-----------------------------------------------------------------------
       USE CLDJ_CMN_MOD
 
       implicit none
 
       real*8,  intent(in)  ::  GMTIME,YGRDJ,XGRDI
       integer, intent(in)  ::  NDAY
       real*8,  intent(out) ::  SZA,COSSZA,SOLFX
 !
       real*8  LOCT
       real*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ
 !
       SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*CPI180)
       SOLDEK = asin(SINDEC)
       COSDEC = cos(SOLDEK)
       SINLAT = sin(YGRDJ)
       SOLLAT = asin(SINLAT)
       COSLAT = cos(SOLLAT)
 !
       LOCT   = (((GMTIME)*15.d0)-180.d0)*CPI180 + XGRDI
       COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
       SZA    = acos(COSSZA)/CPI180
 !
       SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*C2PI/365.d0))
 !
       END SUBROUTINE FJX_TOD
