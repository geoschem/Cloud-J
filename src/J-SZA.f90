!>>>>>>>>cloud-JX code (includes fractional cloud treatments) ver 7.4 (8/2015)<<<<<<<<<<<<
!           uses full Solar-J data sets, i.e., the FJX_---.dat data, not the RRTMG data
!           calculates heating for FastJ bins (1:18) & RRTMG superbins (19:27) w/o gas abs.
!           Solar heating >778nm for clouds and sulfuric acid aerosols is included
!           option to turn off bins 19:27 since not used for J-values
!   With Cloud-J v7.4, the Fast-J data sets & modules are consistent with the Solar-J code
!   Version 7.4 can still be run in a trop-only mode (W_=8 old style) and WACCM mode (no J's <200nm)

      program standalone
! USES:
      USE FJX_CMN_MOD
      USE FJX_SUB_MOD
      USE FJX_INIT_MOD
      USE CLD_SUB_MOD, ONLY : CLOUD_JX

      implicit none
      real*8, dimension(L1_) :: ETAA,ETAB, ZOFL,RI,TI,CLDP,AER1,AER2
      real*8, dimension(L_) :: WLC,WIC
      real*8  GMTAU,PHOTAU,ALBEDO, XLNG,YLAT,XGRD,YGRD,PSURF, SCALEH
      real*8  CF,PMID,PDEL,ZDEL,ICWC,F1,ZKM
      integer, dimension(L1_):: NAA1,NAA2
      integer MYEAR, MONTH, IDAY, JLAT, ILON, NHRMET
      integer I,J,K,L,N,JP, NN
      integer NRAN, RANSEED, LTOP, NJXX
      real*8, dimension(LWEPAR) :: CLDFRW,CLDIWCW,CLDLWCW   ! WCW=Cloud Water Content (g/g)
      character*6, dimension(JVN_)   ::  TITLJXX
      real*8, dimension(L_,8,4)  :: ZPJCLD,ZPJAVG
      real*8, dimension(21,4)    :: ERRJ, ERRJJ,ERRJ2
      integer JP04,JP09,JP11,JP15, ICLD
      integer, dimension(8)      :: JCNT
      character*11, dimension(4) ::  TITJX

      character*8, dimension(8), parameter :: TITERR =  &
      ['clearsky','avgcloud','cldf^3/2','ICA-beam', &
       'ICAs-ran','QCAs-mid','QCAs-avg','all ICAs']

      integer               :: CLDFLDS, CLDFLDX

!!!!!!
      integer :: NSZA


!--------------------key params sent to CLOUD_J-------------------------
      real*8                     :: U0,SZA,REFLB,SOLF,  CLDCOR
      real*8                     :: FG0
      logical                    :: LPRTJ
      real*8,  dimension(L1_+1)  :: PPP,ZZZ
      real*8,  dimension(L1_  )  :: TTT,DDD,RRR,OOO
      real*8,  dimension(L1_)    :: LWP,IWP,REFFL,REFFI
      real*8,  dimension(L1_)    :: CLF
      integer, dimension(L1_)    :: CLDIW
      real*8,  dimension(L1_,AN_):: AERSP
      integer, dimension(L1_,AN_):: NDXAER
      real*8, dimension(L_,JVN_) :: VALJXX
      integer                    :: CLDFLAG,NRANDO,IRAN,LNRG,ICNT
      integer                    :: NICA,JCOUNT
!---U0 = cos (SZA), SZA = solar zenith angle
!---REFLB = Lambertian reflectivity at the Lower Boundary
!---SOLF = solar flux factor for sun-earth distance
!---FG0 = scale for asymmetry factor to get equivalent isotropic (CLDFLAG=3 only)
!---LPRTJ = .true. = turn on internal print in both CLOUD_JX & PHOTO_JX
!--- P = edge press (hPa), Z = edge alt (m), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!--- R = layer rel.hum.(fraction)
!---LWP/IWP = Liquid/Ice water path (g/m2)
!---REFFL/REFFI = R-effective(microns) in liquid/ice cloud
!---CLF = cloud fraction (0.0 to 1.0)
!---CLDIW = integer denoting cloud in layer: 1=water, 2=ice, 3=both
!---AERSP = aerosol path (g/m2) & NDXAER = aerosol index type
!---  aerosols are dimensioned with up to AN_ different types in an ICA layer
!---L1_ = parameter, dim of profile variables, L_+1 for top (non CTM) layer
!---AN_ = parameter, dim of number of aerosols being passed
!---VALJXX = J-values from CLOUD_JX & PHOTO_JX
!---JVN_ = dim of max number of J-s reported out (in the order of fast-JX, not CTM)
!---CLDFLAG = integer index for type of cloud overlap
!---CLOUD_JX:   different cloud schemes
!---CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
!       CLDFLAG = 1  :  Clear sky J's
!       CLDFLAG = 2  :  Averaged cloud cover
!       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
!       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
!       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
!       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
!       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)
!---NRANDO = number of random ICAs to do J's for (CLDFLAG=4)
!---IRAN = starting index for random number selection
!---LNRG = flag for setting max-ran overlap groups:
!---     =0   break max overlap groups at cloud fraction = 0
!---     =3   else fixed 3 layers (1:9, 9:last LWcloud, LWclud+1:LTOP)
!---     else(=6) fixed correlated length max-overlap layers
!---NICA = total number of ICAs

      write(6,'(a)') '>>>begin Cloud-J Driver/Standalone version 7.4'

!--Set up atmosphere for a single column and time for J-values calculation
!--Nominally taken from CTM, but for standalone here is read in
      open (77,file='CTM_GrdCld.dat',status='old',err=91)
      read (77,*)
      read (77,'(i6)') MYEAR
      read (77,'(i6)') IDAY
      read (77,'(i6)') MONTH
      read (77,'(i6)') NHRMET
      read (77,'(2f6.1)') GMTAU,PHOTAU
      read (77,'(f6.1,i4)') YLAT, JLAT
      read (77,'(f6.1,i4)') XLNG, ILON
          YGRD = YLAT*CPI180
          XGRD = XLNG*CPI180
      read (77,'(f6.1)') PSURF
      read (77,'(f6.1)') ALBEDO
      read(77,'(2i5,f5.1)') LNRG,NRANDO,CLDCOR
      read(77,'(f5.2)') FG0
      read(77,'(2i5)')  NWBIN,NSBIN
      read(77,'(2i5)') CLDFLDS,CLDFLDX
        write(6,'(a,2i4,a,f7.4)') 'LNRG = ',LNRG,NRANDO,'    Cloud-Corel',CLDCOR
        write(6,'(a,f10.4)') 'FG0 (assym for direct beam)', FG0
        write(6,'(a,2i5)') '#CLDFLDS(out of) = ',CLDFLDS, CLDFLDX
        write(6,'(a,2i5,5x,a,i5)') 'Atmosphere:',LPAR,LWEPAR, 'LPAR / LWEPAR', L1_
      read (77,*)
      do L = 1,LPAR+1
        read (77,'(i3,1x,2f11.7,2x,f5.1,f5.2,f11.2,2(f7.3,i4))') &
                      J,ETAA(L),ETAB(L),TI(L),RI(L),ZOFL(L) &
                     ,AER1(L),NAA1(L),AER2(L),NAA2(L)
      enddo
      do L = 1,L1_
       PPP(L) = ETAA(L) + ETAB(L)*PSURF
      enddo

!---fast-JX:  INIT_JX is called only once to read in & store all fast-JX data:
!              also sets up random sequence for cloud-JX
!-----------------------------------------------------------------------
      call INIT_FJX (TITLJXX,JVN_,NJXX)
!-----------------------------------------------------------------------

!---ACLIM_FJX sets climatologies for O3, T, D & Z - overwrite with CTM data
!-----------------------------------------------------------------------
      call ACLIM_FJX (YLAT,MONTH,PPP, TTT,OOO, L1_)
!-----------------------------------------------------------------------

      do L = 1,L_
       TTT(L) = TI(L)
       RRR(L) = RI(L)
      enddo
       ZZZ(1)  = 16.d5*log10(1013.25d0/PPP(1))        ! zzz in cm
      do L = 1,L_
       DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
       SCALEH      = 1.3806d-19*MASFAC*TTT(L)
       ZZZ(L+1) = ZZZ(L) -( log(PPP(L+1)/PPP(L)) * SCALEH )
      enddo
       DDD(L_+1) = PPP(L_+1)*MASFAC
       ZZZ(L1_+1) = ZZZ(L1_) + ZZHT
      REFLB = ALBEDO
      LPRTJ = .true.

!---following is readin for cloud data, currently has 160 atmospheres
!   from tropical T319 ECMWF atmosphere used in UCI CTM.
        ZPJAVG(:,:,:) = 0.d0

        ERRJJ(:,:) = 0.d0
        ERRJ2(:,:) = 0.d0
        JCNT(:) = 0



!!!!
      CLDFLDS = 1

      do ICLD=1,CLDFLDS

        read (77,*)
       do L = LWEPAR,1,-1
        read (77,'(i3,1p,e14.5,28x,2e14.5)') &
                      J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
       enddo
!---load CTM-based data on ozone and aerosols on top of climatology
!---to convert kg (STT) in grid cell w/AREA (m2) to # molecules/cm^2
!      D_O3(I,J,L) = 6.023d26*STT(I,J,L,1)/48.d0  *1.d-4 /AREAXY(I,J)
!---to convert kg (STT) in grid cell w/AREA (m2) to PATH (g/m^2)
!      P_AERSL(I,J,L) = STT(I,J,L,1)*1.d-3/AREAXY(I,J)
!---this data should be available from the CTM somewhere.
!---   ZZZ(1:L_+1) is geometric altitude (cm), approx is OK.
!---AERSP = aerosol path (g/m2) and NDXAER() must come from CTM
!    the index must match one in the std or UMich data sets (or add your own).
        AERSP(:,:)  = 0.d0
        NDXAER(:,:) = 0
        do L = 1,L_
          NDXAER(L,1) = NAA1(L)
          AERSP(L,1)  = AER1(L)
          NDXAER(L,2) = NAA2(L)
          AERSP(L,2)  = AER2(L)
        enddo
!---convert cloud data from our EC met fields into cloud data for cloud-JX
!---     init data = cloud fraction, and water content (g/m3) averaged over cell
!---     needs: cloud fraction, ice- and liq-water path (in cloud)
!---     and R-effective of ice and liquid clouds
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
!---cloud input as interpreted by fast_JX
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
!           ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L)
!           ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0)

!---fast-JX:  SOLAR_JX is called  once per grid-cell to set U0, SZA, SOLF
!--- your CTM code may have its own way of calculating and passing these quantities
!-----------------------------------------------------------------------
      call SOLAR_JX(PHOTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!-----------------------------------------------------------------------

!!!!!
      do NSZA = 1,18
      SZA = 79.0 + 1.0*float(NSZA)
      U0 = cos(SZA*CPI180)

      write(6,'(a,2f8.4)') ' >>>sza = ',SZA,U0


      if (LPRTJ) then
        write(6,'(a,f8.3,3f8.5)')'solar zenith angle, solar-f' &
               ,SZA,SOLF,U0,REFLB
        write(6,'(a,f8.3,f8.3)') 'lat/lng',YLAT,XLNG
        call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
      endif

!  locate the position of random number sequence based on year/day/hour
      IRAN = 13+ILON+3*JLAT+5*(MYEAR-1900)+7*IDAY + 11*nint(GMTAU)

          ZPJCLD(:,:,:) = 0.d0

!  >>>>>test mode -- loop across all CLDFLAG's to accumulate errors

!!!!!!!!
      do CLDFLAG = 1,1

       LPRTJ = .true.

!=======================================================================

      call CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT,       &
             DDD,RRR,OOO,   LWP,IWP,REFFL,REFFI,  CLF,CLDCOR,CLDIW, &
             AERSP,NDXAER,L1_,AN_,VALJXX,JVN_,                      &
             CLDFLAG,NRANDO,IRAN,LNRG,NICA,JCOUNT)

          JCNT(CLDFLAG) = JCNT(CLDFLAG) + JCOUNT
!=======================================================================



          if (CLDFLAG .eq. 2) LPRTJ = .false.

        if (CLDFLAG .eq. 8) then
          write(6,'(a,3i8,f8.3)') ' Cloud-J v7.4  ALL ICAs: ICLD/LNRG/NICAs',ICLD,LNRG,NICA,CLDCOR
          write(8,'(a,3i8,f8.3)') ' Cloud-J v7.4  ALL ICAs: ICLD/LNRG/NICAs',ICLD,LNRG,NICA,CLDCOR
          do L=1,L_
          write(8,'(i4,1p,4e10.3)') L,VALJXX(L,JP04),VALJXX(L,JP09),VALJXX(L,JP11),VALJXX(L,JP15)
          enddo
        endif
!   4O3        PHOTON    O2        O(1D)     1.000 mapped to FJX:   3 O3(1D)
!   9NO2       PHOTON    N2        O         1.000 mapped to FJX:   9 NO2
!  11NO3       PHOTON    NO        O2        0.114 mapped to FJX:  10 NO3
!  15HNO3      PHOTON    NO2       OH        1.000 mapped to FJX:  13 HNO3
        JP04 = JIND(4)
        JP09 = JIND(9)
        JP11 = JIND(11)
        JP15 = JIND(15)
        TITJX(1) = 'J(O3>O1D)  '
        TITJX(2) = 'J(NO2)     '
        TITJX(3) = 'J(NO3>all) '
        TITJX(4) = 'J(HNO3)    '
       do L = 1,L_
        ZPJCLD(L,CLDFLAG,1) = VALJXX(L,JP04)
        ZPJCLD(L,CLDFLAG,2) = VALJXX(L,JP09)
        ZPJCLD(L,CLDFLAG,3) = VALJXX(L,JP11)
        ZPJCLD(L,CLDFLAG,4) = VALJXX(L,JP15)

        ZPJAVG(L,CLDFLAG,1) = ZPJAVG(L,CLDFLAG,1) + VALJXX(L,JP04)
        ZPJAVG(L,CLDFLAG,2) = ZPJAVG(L,CLDFLAG,2) + VALJXX(L,JP09)
        ZPJAVG(L,CLDFLAG,3) = ZPJAVG(L,CLDFLAG,3) + VALJXX(L,JP11)
        ZPJAVG(L,CLDFLAG,4) = ZPJAVG(L,CLDFLAG,4) + VALJXX(L,JP15)

       enddo

      enddo   !CLDFLAG


!!!!!!!
      enddo    !NSZA loop, clear sky with LPRTJ!



      goto 92


! look at errors at L=1, L=34 (above clouds), and p-wtd mean (1:33)
       do J=1,4
        do I=1,8
         ZPJCLD(L_,I,J) = 0.d0
         do L = 1,33
         ZPJCLD(L_,I,J) = ZPJCLD(L_,I,J) + ZPJCLD(L,I,J)*(PPP(L)-PPP(L+1))
         enddo
         ZPJCLD(L_,I,J) = ZPJCLD(L_,I,J)/(PPP(1)-PPP(34))
        enddo
       enddo
       L=1
       do J=1,4
        ERRJ(1,J) = ZPJCLD(L,1,J)/ZPJCLD(L,8,J)
        ERRJ(2,J) = ZPJCLD(L,2,J)/ZPJCLD(L,8,J)
        ERRJ(3,J) = ZPJCLD(L,3,J)/ZPJCLD(L,8,J)
        ERRJ(4,J) = ZPJCLD(L,4,J)/ZPJCLD(L,8,J)
        ERRJ(5,J) = ZPJCLD(L,5,J)/ZPJCLD(L,8,J)
        ERRJ(6,J) = ZPJCLD(L,6,J)/ZPJCLD(L,8,J)
        ERRJ(7,J) = ZPJCLD(L,7,J)/ZPJCLD(L,8,J)
       enddo
       L=34
       do J=1,4
        ERRJ( 8,J) = ZPJCLD(L,1,J)/ZPJCLD(L,8,J)
        ERRJ( 9,J) = ZPJCLD(L,2,J)/ZPJCLD(L,8,J)
        ERRJ(10,J) = ZPJCLD(L,3,J)/ZPJCLD(L,8,J)
        ERRJ(11,J) = ZPJCLD(L,4,J)/ZPJCLD(L,8,J)
        ERRJ(12,J) = ZPJCLD(L,5,J)/ZPJCLD(L,8,J)
        ERRJ(13,J) = ZPJCLD(L,6,J)/ZPJCLD(L,8,J)
        ERRJ(14,J) = ZPJCLD(L,7,J)/ZPJCLD(L,8,J)
       enddo
       L=L_
       do J=1,4
        ERRJ(15,J) = ZPJCLD(L,1,J)/ZPJCLD(L,8,J)
        ERRJ(16,J) = ZPJCLD(L,2,J)/ZPJCLD(L,8,J)
        ERRJ(17,J) = ZPJCLD(L,3,J)/ZPJCLD(L,8,J)
        ERRJ(18,J) = ZPJCLD(L,4,J)/ZPJCLD(L,8,J)
        ERRJ(19,J) = ZPJCLD(L,5,J)/ZPJCLD(L,8,J)
        ERRJ(20,J) = ZPJCLD(L,6,J)/ZPJCLD(L,8,J)
        ERRJ(21,J) = ZPJCLD(L,7,J)/ZPJCLD(L,8,J)
       enddo
! accumulate errors for each cloud field
       do J=1,4
       do I=1,21
        ERRJ(I,J) = log(ERRJ(I,J))
        ERRJJ(I,J) = ERRJJ(I,J) + ERRJ(I,J)
        ERRJ2(I,J) = ERRJ2(I,J) + ERRJ(I,J)**2
       enddo
       enddo

! accumulate mean J's

      enddo ! ICLD


! print summary over all CLDFLAGS and all atmospheres (ICLDs)

       do J=1,4
        do L=1,34
         do K=1,8
           ZPJAVG(L,K,J) = ZPJAVG(L,K,J)/float(CLDFLDS)
         enddo
        enddo
       enddo

       do K=1,8
        write(6,'(i5,3x,a)') K,TITCLD(K)
       enddo
       do J=1,4
          write(6,'(a11,i7,8i9)') TITJX(J),(K, K=1,8)
        do L=34,1,-1
          ZKM = 1.d-5*ZZZ(L)
          write(6,'(i3,f6.1,1p,72e9.2)') L,ZKM,(ZPJAVG(L,K,J),K=1,8)
        enddo
       enddo


       do J=1,4
        do I=1,21
         ERRJJ(I,J) = ERRJJ(I,J)/float(CLDFLDS)
         ERRJ2(I,J) = sqrt(ERRJ2(I,J)/float(CLDFLDS))
        enddo
       enddo
       write(6,'(a,2i5)') ' no. atmos / LNRG',ICLD-1,LNRG
       write(6,'(5a)') ' mean error (%)  L=1  &  L=34  & p-avg ',(TITJX(J),J=1,4)
       write(6,'(a11,i8,2p,4f7.0,5x,4f7.0,5x,4f7.0)')  (TITERR(I),JCNT(I),  &
         (ERRJJ(I,J),J=1,4),(ERRJJ(I+7,J),J=1,4),(ERRJJ(I+14,J),J=1,4), I=1,7)
       write(6,'(a11,i8)')  TITERR(8),JCNT(8)
       write(6,'(5a)') ' RMS error (%)   L=1  &  L=34  & p-avg ',(TITJX(J),J=1,4)
       write(6,'(a11,i8,2p,4f7.0,5x,4f7.0,5x,4f7.0)') ,(TITERR(I),JCNT(I),  &
         (ERRJ2(I,J),J=1,4),(ERRJ2(I+7,J),J=1,4),(ERRJ2(I+14,J),J=1,4), I=1,7)
       write(6,'(a11,i8)')  TITERR(8),JCNT(8)

      stop

   91 stop 'error in opening CTM_GrdCld.dat file'
   92 stop
      end
