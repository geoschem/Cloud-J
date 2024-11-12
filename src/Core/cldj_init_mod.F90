!<<<<<<<<<<<<<<<<<<fastJX initialization codes:  need to be called only once
! Solar/Cloud/Fast-J
! v 7.7, minor edits, added !SJ! lines
! v 8.0 - added QH2O for UV abs

! INTERFACE:

      MODULE CLDJ_INIT_MOD

      USE CLDJ_CMN_MOD
      USE CLDJ_ERROR_MOD
!SJ!  USE CLDJ_FJX_MOD

      IMPLICIT NONE

      PUBLIC   :: INIT_CLDJ

      PRIVATE  :: RD_XXX
      PRIVATE  :: RD_CLD
      PRIVATE  :: RD_SSA
      PRIVATE  :: RD_MIE
      PRIVATE  :: RD_UM
      PRIVATE  :: RD_PROF
      PRIVATE  :: RD_TRPROF
      PRIVATE  :: RD_JS_JX
      PRIVATE  :: RD_GEO
      PRIVATE  :: RD_SSAPROF
      PRIVATE  :: RANSET

      CONTAINS

!-----------------------------------------------------------------------
      subroutine INIT_CLDJ (AMIROOT,DATADIR,NLEVELS,NLEVELS_WITH_CLOUD, &
           TITLEJXX,NJXU,ATAU_in,ATAU0_in, NWBIN_in, CLDFLAG_in,        &
           CLDCOR_in,LNRG_in,ATM0_in,use_H2O_UV_abs,NJXX,RC )
!-----------------------------------------------------------------------

      character(len=255)           :: thisloc
      logical, intent(in)          :: AMIROOT
      character(LEN=*), intent(in) :: DATADIR
      integer, intent(in)          :: NLEVELS
      integer, intent(in)          :: NLEVELS_WITH_CLOUD
      integer, intent(in)          :: NJXU
      real*8,  intent(in)          :: ATAU_in
      real*8,  intent(in)          :: ATAU0_in
      integer, intent(in)          :: NWBIN_in
      integer, intent(in)          :: CLDFLAG_in
      real*8,  intent(in)          :: CLDCOR_in
      integer, intent(in)          :: LNRG_in
      integer, intent(in)          :: ATM0_in
      logical, intent(in)          :: use_H2O_UV_abs
      integer, intent(out)         :: NJXX
      integer, intent(out)         :: RC
      character*6, intent(out), dimension(NJXU) :: TITLEJXX

      character*120  TIT_SPEC
      integer  JXUNIT,I, J, K, KR, RANSEED, NUN

      thisloc = ' -> at INIT_CLDJ in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS
      if (AMIROOT) write(6,*) ' Solar/Cloud-J  ver-8.0 initialization parameters'

      ! Use channel 8 to read fastJ data files:
      JXUNIT = 8

#ifndef MODEL_STANDALONE
      ! Set cldj_cmn_mod variables based on input numbers of levels
      L_    = NLEVELS
      L1_   = L_ + 1
      L2_   = L_ + 2
      LWEPAR= NLEVELS_WITH_CLOUD
#endif

      ! Set global variables from inputs
      ATAU     =  ATAU_in
      ATAU0    =  ATAU0_in
      CLDCOR   =  CLDCOR_in
      NWBIN    =  NWBIN_in
      LNRG     =  LNRG_in
      ATM0     =  ATM0_in
      CLDFLAG  =  CLDFLAG_in
      USEH2OUV =  use_H2O_UV_abs

! a v7.7 fix was done within sub ICA_NR, this is now moved up front and must be correct in setup
      if (LNRG.ne.6) then
        CLDCOR = 0.d0         ! v7.7 safety fix for LNRG=0 or 3
      endif

      if (AMIROOT) then
         write(6,'(a,1p,2e12.5)') ' params RAD ZZHT',RAD, ZZHT
         write(6,'(a,3f8.4)') ' params CLDCOR',CLDCOR
         write(6,'(a,4i5)') ' params NWBIN LNRG CLDFLAG',NWBIN,LNRG,CLDFLAG
         write(6,'(a,f8.4,a,f8.4,a,i2)') 'params ATAU0=',ATAU0,'  ATAU=',ATAU,'   option(ATM0)= ', ATM0
         write(6,'(a,i3,a,i3,a,i3)') 'params W_=',W_,'  S_=',S_,'  W_r=',W_r
         write(6,'(a,l1)') 'params use_H2O_UV_abs=',use_H2O_UV_abs
      endif

      if (W_ .ne. 18) then
        call CLOUDJ_ERROR('Invalid no. wavelengths', thisloc, rc)
        return
      endif

! Cloud-J default:  can have added near IR bins (if S_ > W_) but no sub bins
      NSJSUB(1:S_)= NGC(1:S_)
      SJSUB(:,:) = 0.d0   ! default set up for wavelengths when no sub-bins
      SJSUB(:,1) = 1.d0
      do I = W_,  S_
         SJSUB(I,1)   = 1.d0
         SJSUB(I,2:16)= 0.d0
      enddo
      if (AMIROOT) write(6,'(a/(30i3))')'NJSUB=', NSJSUB

! Read in Fast/Solar-J X-sections (spectral data)
      call RD_XXX(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_spec.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_XXX', thisloc, rc)
         return
      endif

! KR = 1:sum(NSJSUB), LDOKR = 0 = FLux=0 & skip this bin (sub-bin): primarily for Solar-J
      KDOKR(:)=0
      KR = 0
      do K = 1,S_
         do J = 1,NSJSUB(K)
            KR = KR+1
            KDOKR(KR) = K
            if (FL(K) .gt. 0.d0) then
               LDOKR(KR) = 1
            else
               LDOKR(KR) = 0
            endif
            if (AMIROOT) write(6,'(A,3I5)')'KR/KDOKR(KR)/LDOKR(KR)',KR, KDOKR(KR), LDOKR(KR)
         enddo
      enddo
      if (KR .ne. W_+W_r) then
         CALL CLOUDJ_ERROR('>>>error w/ RRTM sub bins: KDOKR', thisloc, rc)
         return
      endif
      
      ! Read in cloud scattering data
      call RD_CLD(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_scat-cld.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_CLD', thisloc, rc)
         return
      endif

      ! Read in UCI strat sulf aerosols scattering data
      call RD_SSA(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_scat-ssa.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_SSA', thisloc, rc)
         return
      endif

      ! Read in aerosols scattering data
      call RD_MIE(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_scat-aer.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_MIE', thisloc, rc)
         return
      endif

      ! Read in UMich aerosol scattering data
#ifdef MODEL_GEOSCHEM
      ! Not used in GEOS-Chem. Set to zero for safety.
      WMM = 0.0d0
      UMAER = 0.0d0
#else
      call RD_UM (AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_scat-UMa.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_UM', thisloc, rc)
         return
      endif
#endif

      ! Read in GEOMIP aerosol scattering data
#ifdef MODEL_GEOSCHEM
      ! Not used in GEOS-Chem. Set to zero for safety.
      NGG = 0
      RGG = 0.0d0
      DGG = 0.0d0
      QGG = 0.0d0
      SGG = 0.0d0
      PGG = 0.0d0
#else
      call RD_GEO (AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_scat-geo.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_GEO', thisloc, rc)
         return
      endif
#endif

#ifdef MODEL_STANDALONE
      ! Read in T & O3 climatology used to fill e.g. upper layers or if O3 not calc.
      call RD_PROF(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'atmos_std.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_PROF', thisloc, rc)
         return
      endif

      ! Read in H2O and CH4 profiles for Solar-J
      call RD_TRPROF(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'atmos_h2och4.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_TRPROF', thisloc, rc)
         return
      endif
#else
      ! These variables are only used in ACLIM_FJX which is only called in
      ! Cloud-J standalone. Zero them out for safety if using external model.
      T_REF = 0.0d0
      O_REF = 0.0d0
      H2O_REF = 0.0d0
      CH4_REF = 0.0d0
#endif

      ! Read in zonal mean Strat-Sulf-Aerosol monthly data
#ifdef MODEL_GEOSCHEM
      ! Not used in GEOS-Chem. Set to zero for safety.
      R_GREF = 0.0d0
      X_GREF = 0.0d0
      A_GREF = 0.0d0
      Y_GREF = 0.0d0
      P_GREF = 0.0d0
#else
      call RD_SSAPROF(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'atmos_geomip.dat',rc)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_SSAPROF', thisloc, rc)
         return
      endif
#endif

      NJXX = NJX
      do J = 1,NJXX
        TITLEJXX(J) = TITLEJX(J)
      enddo

      ! Read in photolysis rates used in chemistry code and mapping onto FJX J's
      !---CTM call:  read in J-values names and link to fast-JX names
      call RD_JS_JX(AMIROOT,JXUNIT,TRIM(DATADIR)//'/'//'FJX_j2j.dat', TITLEJXX,NJXX,RC)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RD_JS_JX', thisloc, rc)
         return
      endif

!---setup the random number sequence RAN4
      RANSEED = 66
      call RANSET (NRAN_,RAN4,RANSEED,RC)
      if ( rc /= CLDJ_SUCCESS ) then
         call CLOUDJ_ERROR('Error in RANSET', thisloc, rc)
         return
      endif

      goto 1
    4 continue
      call CLOUDJ_ERROR('Error in read', thisloc, rc)
      return

    1 continue

      if (AMIROOT) write(6,*) ' end of Solar/Cloud-J initialization'

      END SUBROUTINE INIT_CLDJ


!-----------------------------------------------------------------------
      subroutine RD_XXX(AMIROOT,NUN,NAMFIL,RC)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-8.0  added xH2O cross sections for near UV absorption        
!>>>>NEW v-7.6+ added Solar-J bins for some to expand to S_
!     NOTE:  W_=18, use NWBIN 8,12,18 to zero flux for wavelengths only in strat
!>>>>NEW v-7.3  expanded input, full names & notes
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QH2O     H2O cross-sections,  UV-blue absorption, 290-350 nm        
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC

      integer  I, J, JJ, K, IW, NQRD, LQ, NWWW, NSSS
      character*120  TIT_SPEC, TIT_J1N
      character*16 TIT_J1L
      character*6  TIT_J1S,TIT_J2S
      real*8  FWSUM

      thisloc = ' -> at RD_XXX in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!   note that X_ = max # Xsects read in
!           NJX = # fast-JX J-values derived from this (.le. X_)
      if (W_ .ne. 18) then
        call CLOUDJ_ERROR(' no. wavelengths wrong: W_ .ne. 18', thisloc, rc)
        return 
      endif

      open (NUN,FILE=trim(NAMFIL),status='old',form='formatted')

      read (NUN,'(a120)',err=4) TIT_SPEC
      read (NUN,*,err=4)
      read (NUN,'(i5,5x,i5)',err=4) NWWW, NSSS ! dimensions for data set

      if (AMIROOT) write(6,'(a)') adjustl(trim(TIT_SPEC))
      if (AMIROOT) write(6,'(i5,A20,i5,A20)')  NWWW, ' photo-chem wl bins ', &
            NSSS, ' solar heating bins '

      !---these are the dimensions for readin of the solar data sets, it does NOT mean that bin 19:27 need
      ! be calculated
      if (NWWW.ne.WX_ .or. NSSS.ne.SX_) then
         call CLOUDJ_ERROR(' WX_ or SX_ -- incompatible data sets w.r.t parameters', thisloc, rc)
         return
      endif

!----w-params:  1=w-eff  2=w-bins, 3=solar(photons), 4=solar(W/m2), 5=Y-PAR,
! 6=Rayleigh, 7=SJ sub-bins
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
      read (NUN,'(5x,6e10.3)',err=4)    (WL(IW),IW=1,NSSS)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
      read (NUN,'(5x,6e10.3)',err=4)    (WBIN(IW),IW=1,NSSS)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
      read (NUN,'(5x,6e10.3)',err=4)    (FL(IW),IW=1,NSSS)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
      read (NUN,'(5x,6e10.3)',err=4)    (FW(IW),IW=1,NSSS)
         FWSUM=0.d0
         do IW=1, NSSS
            FWSUM= FWSUM + FW(IW)
         enddo
         if (AMIROOT) write(6,*) 'total Solar flux=', FWSUM

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
      read (NUN,'(5x,6e10.3)',err=4)    (FPAR(IW),IW=1,NSSS)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                        ' notes: ', adjustl(trim(TIT_J1N))
      read (NUN,'(5x,6e10.3)',err=4)    (QRAYL(IW),IW=1,NSSS)

!SJ! SJ-sub-bins for RRMTG, CLIRAD, GGLLNL
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
      do I= NWWW, NSSS  ! fraction of solar radiation for each sub-bin
         read  (NUN,'(5x,5f10.5)',err=4) (SJSUB(I,IW),IW=1,15)
         if (AMIROOT) write(6,'(5x,5f10.6)') (SJSUB(I,IW),IW=1,15)
      enddo
!SJ! SJSUB needs to be reset to 1's, since in Cloud-J, only the first sub-bin is used
     do I = 1,NSSS
        SJSUB(I,1)    = 1.0d0
        SJSUB(I,2:15) = 0.0d0
     enddo
!---Read H2O X-sects  *** v8.0
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(5x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          (QH2O(IW),IW=1,NWWW)
      if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a120)') TIT_J1S,TIT_J1L,' notes:',TIT_J1N   !print
      ! If H2O UV absorption turned off then zero out QH2O and write note:
      if ( .not. USEH2OUV ) then
         QH2O = 0.0d0
         if (AMIROOT) write(6,'(a60)') 'WARNING: H2O UV absorption is turned off in Cloud-J'
      endif
     
!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,1), (QO2(IW,3),IW=1,NWWW)
         TITLEJX(1) = TIT_J1S
         TITLEJL(1) = TIT_J1L
         LQQ(1) = 3
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,2), (QO3(IW,3),IW=1,NWWW)
        TITLEJX(2) = TIT_J1S
        TITLEJL(2) = TIT_J1L
        LQQ(2) = 3
        if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                          ' notes: ',adjustl(trim(TIT_J1N))

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)
        TITLEJX(3) = TIT_J1S
        TITLEJL(3) = TIT_J1L
        LQQ(3) = 3
        if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                          ' notes: ',adjustl(trim(TIT_J1N))

!---Read remaining species:  X-sections at 1-2-3 T_s
!---read in 1 to 3 X-sects per J-value (JJ)
        JJ = 3
!-- read new Xsection block
    3 continue
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
        if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                          ' notes: ',adjustl(trim(TIT_J1N))
        if (TIT_J1S .eq. 'endofJ') goto 1
!---try to add a new Xsect
    2 continue
       JJ = JJ+1
       LQ = 1
       if (JJ .gt. X_) then
          call CLOUDJ_ERROR('X_ not large enough', thisloc, rc)
          return
       endif
       TITLEJX(JJ) = TIT_J1S
       TITLEJL(JJ) = TIT_J1L
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(JJ),TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
!try to read a 2nd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
        if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                          ' notes: ',adjustl(trim(TIT_J1N))
        if (TIT_J1S .eq. 'endofJ') goto 1
      if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 2
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
!try to read a 3rd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (AMIROOT) write(6,'(1x,a6,1x,a16,a8,a)') trim(TIT_J1S),trim(TIT_J1L), &
                                           ' notes: ',adjustl(trim(TIT_J1N))
         if (TIT_J1S .eq. 'endofJ') goto 1
       if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 3
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
       else
        goto 2
       endif
      else
        goto 2
      endif
      goto 3
    4 continue
      call CLOUDJ_ERROR('Error in read', thisloc, rc)
      return

    1 continue
      NJX = JJ

!---read in complete, process Xsects for reduced wavelengths (Trop-Only)
!---    possibly also for WACCM >200nm-only version.
!---TROP-ONLY (W_ = 12 or 8) then drop the strat Xsects (labeled 'x')
      if (NWBIN .eq. 12 .or. NWBIN .eq. 8) then
         if (AMIROOT) write(6,'(a)')  &
              ' >>>TROP-ONLY reduced wavelengths, drop strat X-sects'
         JJ = 3
         do J = 4,NJX
            if (SQQ(J) .ne. 'x') then
!---collapse Xsects
               JJ = JJ+1
               if (JJ .lt. J) then
                  TITLEJX(JJ) = TITLEJX(J)
                  LQQ(JJ) = LQQ(J)
                  SQQ(JJ) = SQQ(J)
                  do LQ = 1,LQQ(J)
                     TQQ(LQ,JJ) = TQQ(LQ,J)
                     do IW = 1,NWWW
                        QQQ(IW,LQ,JJ) = QQQ(IW,LQ,J)
                     enddo
                  enddo
               endif
            endif
         enddo
         NJX = JJ
      endif

      do J = 1,NJX
         if (AMIROOT) write(6,'(a8,i5,2x,a6,2x,a16,2x,a1,i3,2x,3f6.1)') &
           'X-sects', J,trim(TITLEJX(J)),trim(TITLEJL(J)), &
           SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
      enddo
!---need to check that TQQ (= T(K) or p(hPa)) is monotonically increasing:
      do J = 1,NJX
         if ((LQQ(J) .eq. 3) .and. (TQQ(2,J) .ge. TQQ(3,J))) then
            call CLOUDJ_ERROR('TQQ out of order', thisloc, rc)
            return
         endif
         if ((LQQ(J) .eq. 2) .and. (TQQ(1,J) .ge. TQQ(2,J))) then
            call CLOUDJ_ERROR('TQQ out of order', thisloc, rc)
            return
         endif
      enddo

!---if FL(K) =0, then R.T calc skipped, method for dropping to 8 or 12 trop-only bins
!           Note that the 'X' cross sections are also skipped
      if (NWBIN .eq. 12) then
         do IW = 1,4
            FL(IW) = 0.d0
         enddo
         do IW = 9,10
            FL(IW) = 0.d0
         enddo
      endif
      if (NWBIN .eq. 8) then
         do IW = 1,4
            FL(IW) = 0.d0
         enddo
         FL(5) = FL(5) * 2.d0
         do IW = 6,11
            FL(IW) = 0.d0
         enddo
      endif

      close(NUN)

      END SUBROUTINE RD_XXX


!-----------------------------------------------------------------------
      subroutine RD_CLD(AMIROOT,NUN,NAMFIL,RC)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.5
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-cld.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC

      integer  I,J,K,L, JCC
      character*120 TITLE0
      real*8     GCCJ, XNDR,XNDI

      thisloc = ' -> at RD_CLD in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

        read (NUN,'(a80)',err=4) TITLE0
          if (AMIROOT) write(6,'(a)') trim(TITLE0)                                
        read (NUN,'(i4)')  NCC
        read (NUN,'(i4)')  MCC
          if (AMIROOT) write(6,'(3i6,a)') NCC,MCC,SX_,' types of clouds & #Reff'
        read (NUN,*)
        read (NUN,*)
        read (NUN,*)
        read (NUN,*)
        read (NUN,*)

        do K = 1, NCC
           read (NUN,'(a12,f8.5)',err=4) TITLCC(K),DCC(K)
           if (AMIROOT) write(6,'(a,i4,1x,a12,f8.5)') 'Cloud#',K,trim(TITLCC(K)),DCC(K)
           do J = 12, SX_
              do I = 1,MCC
                 read (NUN, &
                      '(i2, 1x, f5.2, f5.1, f7.1, f5.3, e8.1,f6.3,f8.5,7f6.3)',&
                      err=4) &
                      JCC,WCC(J,K),RCC(I,K),GCC(I,K),XNDR,XNDI,                &
                      QCC(J,I,K),SCC(J,I,K),(PCC(L,J,I,K),L=2,8)
                 if (JCC .ne. J) goto 4
                 PCC(1,J,I,K) = 1.d0
              enddo
              read (NUN,*)
           enddo
        enddo

! replicate all cloud data for w < 295 nm from J=12 (= 295 nm),
! OK since mostly trop clouds.
        do K = 1,NCC
           do J = 1,11
              WCC(J,K) = WCC(12,K)
              do I = 1,MCC
                 QCC(J,I,K) = QCC(12,I,K)
                 SCC(J,I,K) = SCC(12,I,K)
                 do L = 1,8
                    PCC(L,J,I,K) = PCC(L,12,I,K)
                 enddo
              enddo
           enddo
        enddo

        goto 2

    4 continue
        call CLOUDJ_ERROR('Error in read', thisloc, rc)
        return
        
    2 continue
        close(NUN)

         if (AMIROOT) write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0',ATAU,ATAU0

      END SUBROUTINE RD_CLD


!-----------------------------------------------------------------------
      subroutine RD_SSA(AMIROOT,NUN,NAMFIL,RC)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.4
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-ssa.dat)
!     NUN      Channel number for reading data file
!     NSS      Number of categories for cloud scattering phase functions
!     QSS      Cloud scattering phase functions
!     WSS      5 Wavelengths for supplied phase functions
!     PSS      Phase function: first 8 terms of expansion
!     RSS      Effective radius associated with cloud type:
!                                                  Integ(r^3 dr)/Integ(r^2 dr)
!     GSS      Effective geometric cross section:  Integ(pi r^2 dr)
!     SSS      Single scattering albedo
!     DSS      density (g/cm^3)
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC

      integer  I, J, JSS, K, JCC, NSX_
      character*120 TITLE0
      real*8     WJSS,XNDR,XNDI

      thisloc = ' -> at RD_SSA in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)
      read (NUN,'(a120)',err=4) TITLE0
      if (AMIROOT) write(6,'(a)') adjustl(trim(TITLE0))
      read (NUN,*)
      read (NUN,'(i4,i4)')  NSS, NSX_
      read (NUN,*)
      if (AMIROOT) write(6,'(i6,a)') NSS, ' types of strat sulf aerosols'
      do K = 1,NSS
         read (NUN,'(a10, 3f8.4, 2f8.1)')    &
              TITLSS(K),RSS(K),GSS(K),DSS(K),TSS(K),WSS(K)
         if (AMIROOT) write(6,'(i4,1x,a12,2f10.4,2f8.1)') K,TITLSS(K),RSS(K),DSS(K),&
               TSS(K),WSS(K)
         do J = 5, NSX_
            read(NUN,'(i2,2f8.4,e8.1,2f8.5,7f6.3)')      &
                 JSS,WJSS,XNDR,XNDI,QSS(J,K),SSS(J,K),(PSS(I,J,K), I=2,8)
            PSS(1,J,K) = 1.d0
         enddo
      enddo
! reproduce all SSA data for J=1:4 with J=5
      do K = 1,NSS
         do J = 1,4
            QSS(J,K) = QSS(5,K)
            SSS(J,K) = SSS(5,K)
            do I = 1,8
               PSS(I,J,K) = PSS(I,5,K)
            enddo
         enddo
      enddo
      goto 2

    4 continue
      call CLOUDJ_ERROR('Error in read', thisloc, rc)
      return

    2 continue
        close(NUN)

      END SUBROUTINE RD_SSA


!-----------------------------------------------------------------------
      subroutine RD_MIE(AMIROOT,NUN,NAMFIL,RC)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX ver 7.3+
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC

      integer  I, J, K , JAA
      character*120 TITLE0
! TITLAA: Title for scat data NEEDS to be in COMMON
!      character*12 TITLAA(A_) 
      Character*12 TITLAAJ
      real*8   RAAJ, DAAJ

      thisloc = ' -> at RD_MIE in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS

      if (AMIROOT) write(6,'(i5,a)') NUN,trim(NAMFIL)

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

      read (NUN,'(a120)',err=4) TITLE0

      if (AMIROOT) write(6,'(a)') adjustl(trim(TITLE0))
      read (NUN,*)
      read (NUN,*)
      do J = 1, A_
         ! Change width of RAAJ and DAAJ from 6 to 7 to accomodate larger vals
         ! (ewl, 2/27/23)
         !read (NUN,'(i4,1x,a12,1x,2f6.3,1x,a120)',err=4) &
         read (NUN,'(i4,1x,a12,1x,2f7.3,1x,a120)',err=4) &
              JAA,TITLAAJ,RAAJ,DAAJ,TITLE0
         if (JAA.gt.0) then
            TITLAA(J) = TITLAAJ
            RAA(J) = RAAJ
            DAA(J) = DAAJ
            do K = 1, 5
               ! Change width of QAA from 7 to 9, and each PAA from 6 to 7
               ! to accomodate larger values (ewl, 2/27/23)
               !read (NUN,'(f4.0,f7.4,f7.4,7f6.3)',err=4) &
               read (NUN,'(f4.0,f9.4,f7.4,7f7.3)',err=4) &
                    WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
               PAA(1,K,J) = 1.d0
            enddo
            NAA = J
            if (AMIROOT) write(6,'(i5,1x,a12,1x,7f9.3,1x,a)')   &
                  J,TITLAAJ,RAAJ,DAAJ,(QAA(K,J),K=1,5),trim(TITLE0)
         else
            goto 2
         endif
      enddo
      goto 2

    4 continue
      call CLOUDJ_ERROR('Error in read', thisloc, rc)
      return

    2 continue

      close(NUN)

      END SUBROUTINE RD_MIE


!-----------------------------------------------------------------------
      subroutine RD_UM(AMIROOT,NUN,NAMFIL,RC)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC

      integer  I, J, K, L
      character*78 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      thisloc = ' -> at RD_UM in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a78)') TITLE0
      if (AMIROOT) write(6,*) 'UMichigan Aerosols ', adjustl(trim(TITLE0))
      read(NUN,'(5x,10f5.0)') WMM
      if (AMIROOT) write(6,'(a,10f7.1)') ' UMIchigan aerosol wavelengths:',WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
         read(NUN,'(a4)') TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
         do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
            read(NUN,'(18f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,6)
         enddo
      enddo
      close(NUN)
      if (AMIROOT) write(6,'(a)') 'collapse UM wavelengths, drop 550 nm'
      WMM(4) = WMM(5)
      WMM(5) = WMM(6)
      do L=1,33
         do K=1,21
            do I=1,3
               UMAER(I,4,K,L) = UMAER(I,5,K,L)
               UMAER(I,5,K,L) = UMAER(I,6,K,L)
            enddo
         enddo
      enddo

      if (AMIROOT) write(6,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)

      END SUBROUTINE RD_UM


!-----------------------------------------------------------------------
      subroutine RD_PROF(AMIROOT,NJ2,NAMFIL,RC)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK

      character*78 TITLE0
!
      thisloc = ' -> at RD_PROF in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
      if (AMIROOT) write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (T_REF(I,L,M), I=1,41)
! volume mixing ratio from 0km to 60 km for every 2km resolution in 
! pressure altitude z*
        read (NJ2,'(3X,11F7.4)') (O_REF(I,L,M), I=1,31)
      enddo
      close (NJ2)

!  Extend climatology to 100 km
!  LREF =51 in cmn_fjx_mod.f90
      OFAC = exp(-2.d5/5.d5)
      do I = 32,LREF
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          O_REF(I,L,M) = O_REF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,LREF
        T_REF(I,L,M) = T_REF(41,L,M)
      enddo
      enddo
      enddo

      close(NJ2)

 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_PROF


!-----------------------------------------------------------------------
      subroutine RD_TRPROF(AMIROOT,NJ2,NAMFIL,RC)
!-----------------------------------------------------------------------
!  Routine to input H2O and CH4 reference profiles 'atmos_h2och4.dat'
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216

      character*78 TITLE0
!
      thisloc = ' -> at RD_TRPROF in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
      if (AMIROOT) write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
!        write(6,'(1X,I3,3X,I2)')   LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11E9.2)') (H2O_REF(I,L,M), I=1,31)
        read (NJ2,'(3X,11F9.2)') (CH4_REF(I,L,M), I=1,31)
!        write (6,'(3X,11E9.2)') (H2O_REF(I,L,M), I=1,31)
!        write(6,'(3X,11F9.2)') (CH4REF(I,L,M), I=1,31)

      enddo
      close (NJ2)

!  Extend climatology to 100 km
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 32,LREF
        H2O_REF(I,L,M) = H2O_REF(31,L,M)
        CH4_REF(I,L,M) = CH4_REF(31,L,M)
      enddo
      enddo
      enddo

      close(NJ2)
 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_TRPROF


!-----------------------------------------------------------------------
      subroutine RD_JS_JX(AMIROOT,NUNIT,NAMFIL,TITLEJX,NJX,RC)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's JVMAP
!    including scaling factor JFACTA
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JVMAP
!           real*8  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JVMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JVMAP     label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
!
      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in)                    ::  NUNIT, NJX
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) :: TITLEJX
      integer, intent(out)                   :: RC

      integer   J,JJ,K
      character*120 CLINE
      character*50 T_REACT
      character*6  T_FJX
      real*8 F_FJX

      thisloc = ' -> at RD_JS_JX in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS
! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JVMAP(:) = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
         if (AMIROOT) write(6,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ .eq. 9999 ) exit
       if (JJ .gt. JVN_) then
          call CLOUDJ_ERROR('JVN_ must be >= number of entries in FJX_j2j.dat!', thisloc, rc)
          return
       endif
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JVMAP(JJ) = T_FJX
        NRATJ = JJ

       ! ewl: new for GEOS-Chem
       ! SDE 03/31/13: Check number of branches
       ! Note that the order of the branches in
       ! globchem.dat must match the order in
       ! FJX_j2j.dat
       READ (T_REACT(1:10),"(a10)") RNAMES(JJ)
       RNAMES(JJ) = TRIM(RNAMES(JJ))
       BRANCH(JJ) = 1
       DO K=1,(JJ-1)
          IF (RNAMES(JJ) == RNAMES(K)) THEN
             BRANCH(JJ) = BRANCH(K) + 1
          ENDIF
       ENDDO

      enddo

      close(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
        if (JVMAP(K) .eq. TITLEJX(J)) then
         JIND(K) = J
        endif
       enddo
      enddo

      if (AMIROOT) write(6,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
      do K=1,NRATJ
       if (JVMAP(K) .ne. '------' ) then
        J = JIND(K)
        if (J.eq.0) then
         if (AMIROOT) write(6,'(i5,1x,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' no mapping onto fast-JX',JVMAP(K)
        else
         if (AMIROOT) write(6,'(i5,1x,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' mapped to FJX:',J,TITLEJX(J)
        endif
       endif
      enddo

      END SUBROUTINE RD_JS_JX


!-----------------------------------------------------------------------
      subroutine RD_GEO(AMIROOT,NUN,NAMFIL,RC)
!-----------------------------------------------------------------------
!-------GEOMIP SSA scattering data set for fast-JX ver 7.5 ONLY FOR Cloud-J & RRTMG 27 bins
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-geo.dat)
!     NUN      Channel number for reading data file
!     NGG      Number of sequentially increasing R-eff's for GEOMIP aerosols
!     RGG      Effective radius associated with cloud type
!     DGG      density (g/cm^3)
!     QGG      ratio optical to geometric X-section
!     SGG      Single scattering albedo
!     PGG      Phase function: first 8 terms of expansion
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC

      integer  I, J, K
      character*120 TITLE0
      real*8     WGGJ,XNDR,XNDI,G1,G2,G3

      thisloc = ' -> at RD_GEO in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

      read (NUN,'(a120)',err=4) TITLE0
      if (AMIROOT) write(6,'(a)') trim(TITLE0)
      read (NUN,'(i4)')  NGG
      if (AMIROOT) write(6,'(i6,a)')  NGG, ' Reff-s for GEO SSA'
      read (NUN,*)
      read (NUN,*)
      do K = 1,NGG
         read(NUN,'(10x,5f8.4)') RGG(K),G1,DGG(K),G2,G3
         if (AMIROOT) write(6,'(i4,1x,3f8.4,2f8.1)') K,RGG(K),DGG(K), G1,G2,G3
         do J = 5, 27
            read (NUN,'(2x,2f8.4,e8.1,2f8.5,7f6.3)',err=4) &
                 WGGJ,XNDR,XNDI,QGG(J,K),SGG(J,K),(PGG(I,J,K),I=2,8)
            PGG(1,J,K) = 1.d0
         enddo
      enddo
! reproduce all GEO SSA data for w < 202 nm (J=1:4)
      do K = 1,NGG
         do J = 1,4
            QGG(J,K) = QGG(5,K)
            SGG(J,K) = SGG(5,K)
            do I = 1,8
               PGG(I,J,K) = PGG(I,5,K)
            enddo
         enddo
      enddo
      goto 2

    4 continue
      call CLOUDJ_ERROR('Error in read', thisloc, rc)
      return

    2 continue
      close(NUN)

      END SUBROUTINE RD_GEO


!-----------------------------------------------------------------------
      subroutine RD_SSAPROF(AMIROOT,NJ2,NAMFIL,RC)
!-----------------------------------------------------------------------
!  Routine to input SSA-GEO reference profiles for 'atmos_geomip.dat'
!      R_GEO = effective radius (microns)
!      X_GEO = mass fraction (1e-9 kg-H2SO4/kg-air)
!-----------------------------------------------------------------------

      character(len=255)  :: thisloc
      logical, intent(in) :: AMIROOT
      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
      integer, intent(out)     :: RC
!
      integer J,L,M
      character*78 TITLE0
!
      thisloc = ' -> at RD_SSAPROF in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(a)') TITLE0
         if (AMIROOT) write(6,'(1x,a)') TITLE0
      read (NJ2,*)
      read (NJ2,*)
! only specify 19 pressure levels from 2.7 hPa to 340 hPa
      read (NJ2,'(19f7.2)') (P_GREF(L),L=1,19)
      if (AMIROOT) write(6,'(19f7.2)')(P_GREF(L),L=1,19)
      read (NJ2,*)
      read (NJ2,*)
      read (NJ2,*)
! latitude bins 1:36 are Gauss, but approx 1.3953 + (J-33)*2.7906 deg
      read (NJ2,'(32f5.1)') (Y_GREF(L),L=1,32)
      if (AMIROOT) write (6,'(32f5.1)') (Y_GREF(L),L=1,32)

      read (NJ2,*)
      read (NJ2,'(32f5.1)') (Y_GREF(L),L=64,33,-1)
      if (AMIROOT) write (6,'(32f5.1)') (Y_GREF(L),L=33,64)
      read (NJ2,'(a)') TITLE0
         if (AMIROOT) write(6,'(1x,a)') TITLE0
      do M = 1,12
          read (NJ2,*)
        do J = 1,64
          read (NJ2,'(11x,18f6.3)') (R_GREF(J,L,M), L=1,18)
                   R_GREF(J,19,M) = 0.d0
        enddo
      enddo
      read (NJ2,'(a)') TITLE0
         if (AMIROOT) write(6,'(1x,a)') TITLE0
      do M = 1,12
          read (NJ2,*)
        do J = 1,64
          read (NJ2,'(11x,18f6.3)') (X_GREF(J,L,M), L=1,18)
                   X_GREF(J,19,M) = 0.d0
        enddo
      enddo
      read (NJ2,'(a)') TITLE0
         if (AMIROOT) write(6,'(1x,a)') TITLE0
      do M = 1,12
          read (NJ2,*)
        do J = 1,64
          read (NJ2,'(11x,18f6.3)') (A_GREF(J,L,M), L=1,18)
                   A_GREF(J,19,M) = 0.d0
        enddo
      enddo

      close(NJ2)

      END SUBROUTINE RD_SSAPROF


!-----------------------------------------------------------------------
      SUBROUTINE RANSET (ND,RAN4L,ISTART,RC)
!-----------------------------------------------------------------------
!  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND)
!     program RAN3 from Press, based on Knuth

      character(len=255) ::  thisloc
      integer, parameter ::  MBIG=1000000000
      integer, parameter ::  MSEED=161803398
      integer, parameter ::  MZ=0
      real*4 , parameter ::  FAC=1.e-9
      integer,intent(in)    :: ND
      real*4, intent(out)   :: RAN4L(ND)
      integer,intent(inout) :: ISTART
      integer,intent(out)   :: RC

      integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP

      thisloc = ' -> at RANSET in module cldj_init_mod.F90'
      rc = CLDJ_SUCCESS
!---initialization and/or fix of ISEED < 0
        MJ = MSEED - abs(ISTART)
        MJ = mod(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        do I=1,54
          II = mod(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          if (MK.lt.MZ) then
            MK=MK+MBIG
          endif
          MJ = MA(II)
        enddo
        do K=1,4
         do I=1,55
           MA(I)=MA(I)-MA(1+MOD(I+30,55))
           if (MA(I) .lt. MZ) then
             MA(I) = MA(I)+MBIG
           endif
         enddo
        enddo
        INEXT = 0
        INEXTP = 31
        ISTART = 1
!---generate next ND pseudo-random numbers
      do J=1,ND
         INEXT = mod(INEXT,55) +1
         INEXTP = mod(INEXTP,55) +1
         MJ = MA(INEXT) - MA(INEXTP)
        if (MJ .lt. MZ) then
          MJ=MJ+MBIG
        endif
         MA(INEXT) = MJ
         RAN4L(J) = MJ*FAC
      enddo

      END SUBROUTINE RANSET


      END MODULE CLDJ_INIT_MOD
