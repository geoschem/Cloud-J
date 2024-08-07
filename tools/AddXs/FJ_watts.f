c-------FJ_index11d.f    modern code for integrating over REFR INDEX, Rayleigh, YPar & fluxes
!           uses both photons and Watts to do weighting  Example = water liqu & ice
!           starts with pratmo full 77+SR bins and then generates the Solar-J bins (18+9)

      implicit none
      integer, parameter :: NC_ = 199
      integer, parameter :: NB_ = 99
      integer, parameter :: NX_ = 6785
      integer, parameter :: NY_ = 40000

      integer, parameter :: IU_ = 0
      real*8,  parameter :: TL_ = 250.d0

      real*8   SRB(15,8)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
      integer  IBINX(NX_),IBINY(NY_)

      real*8, dimension(NC_) ::
     &          FFBIN,WWBIN,RRBIN,YYBIN,AABIN,BBBIN,CCBIN,DDBIN

      real*8  W1(NB_),W2(NB_),WNM, RAYLAY,YPAR
      real    WL,RNW,CNW,RNI,CNI,A1,A2

      real*8  WX(NX_),FX(NX_)
      real*8  WY(NY_),FY(NY_)

      real*8, dimension(NB_) :: FBIN,WBIN,RBIN,YBIN,ABIN,BBIN,CBIN,DBIN
      real*8, dimension(NB_) ::
     &           FBINw,WBINw,RBINw,YBINw,ABINw,BBINw,CBINw,DBINw
      character*78 TITLE
      character*6 TITLNEW
      integer  I, J,J1,J2

      integer NC1,NC2,NC3,NC4, NB4
      real*8 W11,W22

!!!!! this reads in the full set of wavelength bins needs to map the S-R bands
!     and not-adjacent bins into the Fast-JX 18 bins plus the Solar-J bins
!

      open (1, file='SolarJ_bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(a)') TITLE
!          write(6,'(a)') TITLE
        read(1,'(a)') TITLE
!          write(6,'(a)') TITLE
        read(1,'(4i5)') NC1,NC2,NC3,NC4
!   NC1=1, NC2=38 last strat bin (JX#11), NC3=76 (JX#18) last trop bin
          if (NC4 .gt. NC_) stop
        read(1,'(5x,i5,f8.3)') (IJX(I),WCBIN(I), I=1,NC4+1)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
!  I tracks the 77 or 85 or NC4 (w/Solar-J) bins from the high-res pratmo wavelengths
!  J or IJX(I) tracks the 1:18 (Cloud-J) ro 1:27 (Solar-J) bins
! convert all to microns
        do I = 1,NC4+1
          WCBIN(I) = WCBIN(I)*1.d-3
        enddo
        NB4 = IJX(NC4)
          if (NB4 .gt. NB_) stop
      close (1)


      open (2, file='SolarF_watts.dat', status='OLD')
        read(2,'(a)') TITLE
!            write(6,'(a)') TITLE
       do I = 1,NX_
        read(2,'(5x,f10.4,3e14.3)') WX(I),FX(I)
       enddo
       do I = 2,NX_-1
        FX(I) = FX(I) * (WX(I+1)-WX(I-1))*0.5d0
       enddo
        FX(1) = 0.d0
        FX(NX_) = 0.d0
      close (2)
!---now assign bin #(I=1:NC4) to each ASTM1x microbin J (1:1697)
        IBINX(:) = 0
      do I=1,NC4
         W11 = WCBIN(I)
         W22 = WCBIN(I+1)
        do J=1,NX_
          if (WX(J) .gt. W11) goto 11
        enddo
          J = NX_ + 1
   11     J1 = J
        do J=J1,NX_
          if (WX(J) .gt. W22) goto 12
        enddo
          J = NX_ + 1
   12     J2 = J-1
        do J=J1,J2
         IBINX(J) = I
        enddo
!      write(6,'(i6,2f8.4,2i6,2f8.4)') I,W11,W22,J1,J2,WX(J1),WX(J2)
      enddo
!!!!this binning does not interpolate and is OK for large WX bins


      open (3, file='SolarF_photons.dat', status='OLD')
        read(3,'(a)') TITLE
!            write(6,'(a)') TITLE
        read(3,*)
       do I = 1,NY_
        read(3,'(f10.4,e10.3)') WNM,FY(I)
        WY(I) = 1.d-3*WNM
       enddo
      close (3)
!---now assign bin #(I=1:NC4) to each p05nm microbin J (1:40000)
        IBINY(:) = 0
      do I=1,NC4
         W11 = WCBIN(I)
         W22 = WCBIN(I+1)
        do J=1,NY_
          if (WY(J) .gt. W11) goto 16
        enddo
          J = NY_ + 1
   16     J1 = J
        do J=J1,NY_
          if (WY(J) .gt. W22) goto 17
        enddo
          J = NY_ + 1
   17     J2 = J-1
        do J=J1,J2
         IBINY(J) = I
        enddo
!      write(6,'(i6,2f8.4,2i6,2f8.4)') I,W11,W22,J1,J2,WY(J1),WY(J2)
      enddo
!!!! this binning does not interpolate and is OK for large bins
!     it has 7% error in the very short wavel S-R bins of pratmo.


!!!! integration of refractive index at high-res for photons and then Watts
        FFBIN(:) = 0.d0
        WWBIN(:) = 0.d0
        RRBIN(:) = 0.d0
        YYBIN(:) = 0.d0
        AABIN(:) = 0.d0
        BBBIN(:) = 0.d0
        CCBIN(:) = 0.d0
        DDBIN(:) = 0.d0
       J = 1
       I = 1
      do while (WY(I) .lt. WCBIN(J))
       I = I+1
      enddo
      do J = 1,NC4
       do while (WY(I) .lt. WCBIN(J+1) .and. I .lt. NY_)
        WL = WY(I)
        call NDXWATER(IU_, WL, TL_, RNW, CNW, A1,A2)
        call NDXICE  (IU_, WL, TL_, RNI, CNI, A1,A2)
        FFBIN(J) = FFBIN(J) + FY(I)
        WWBIN(J) = WWBIN(J) + FY(I)*WY(I)
        RRBIN(J) = RRBIN(J) + FY(I)*RAYLAY(WY(I))
        YYBIN(J) = YYBIN(J) + FY(I)*YPAR(WY(I))
        AABIN(J) = AABIN(J) + FY(I)*RNW
        BBBIN(J) = BBBIN(J) + FY(I)*CNW
        CCBIN(J) = CCBIN(J) + FY(I)*RNI
        DDBIN(J) = DDBIN(J) + FY(I)*CNI
        I = I+1
       enddo
        WWBIN(J) = WWBIN(J)/FFBIN(J)
        RRBIN(J) = RRBIN(J)/FFBIN(J)
        YYBIN(J) = YYBIN(J)/FFBIN(J)
        AABIN(J) = AABIN(J)/FFBIN(J)
        BBBIN(J) = BBBIN(J)/FFBIN(J)
        CCBIN(J) = CCBIN(J)/FFBIN(J)
        DDBIN(J) = DDBIN(J)/FFBIN(J)
       if (I .ge. NY_) goto 2
      enddo
    2 continue

!!!! NC4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
      write(6,'(2a)') 'pratm bin#   w-eff  phot/Watt X-ray     Y-PAR',
     &    '     Liq:  nr + ni       Ice:  nr + ni'
      do J = 1,NC4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Phot ',J,
     &   WWBIN(J),FFBIN(J),RRBIN(J),YYBIN(J),AABIN(J),
     &   BBBIN(J),CCBIN(J),DDBIN(J)
      enddo

!!!! Second integration from NC4 bins to NB4 bins (Fast-JX + Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        WBIN(:) = 0.d0
        RBIN(:) = 0.d0
        YBIN(:) = 0.d0
        ABIN(:) = 0.d0
        BBIN(:) = 0.d0
        CBIN(:) = 0.d0
        DBIN(:) = 0.d0
       do I=16,NC4
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          WBIN(J) = WBIN(J) + FFBIN(I)*WWBIN(I)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)
          YBIN(J) = YBIN(J) + FFBIN(I)*YYBIN(I)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)
          DBIN(J) = DBIN(J) + FFBIN(I)*DDBIN(I)
!           write(6,'(a,2i5,f8.3)') ' NC/NB# ', I,J
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)         *SRB(I,J)
          WBIN(J) = WBIN(J) + FFBIN(I)*WWBIN(I)*SRB(I,J)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)*SRB(I,J)
          YBIN(J) = YBIN(J) + FFBIN(I)*YYBIN(I)*SRB(I,J)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)*SRB(I,J)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)*SRB(I,J)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)*SRB(I,J)
          DBIN(J) = DBIN(J) + FFBIN(I)*DDBIN(I)*SRB(I,J)
!           write(6,'(a,2i5,f8.3)') ' NC/NB# ', I,J, SRB(I,J)
        enddo
       enddo
       do J=1,NB4
        if (FBIN(J) .gt. 0.d0) then
          WBIN(J) = WBIN(J)/FBIN(J)
          RBIN(J) = RBIN(J)/FBIN(J)
          YBIN(J) = YBIN(J)/FBIN(J)
          ABIN(J) = ABIN(J)/FBIN(J)
          BBIN(J) = BBIN(J)/FBIN(J)
          CBIN(J) = CBIN(J)/FBIN(J)
          DBIN(J) = DBIN(J)/FBIN(J)
        endif
       enddo

!!!! NB4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
      write(6,'(2a)') ' S-J bin#    w-eff  phot/Watt X-ray     Y-PAR',
     &    '     Liq:  nr + ni       Ice:  nr + ni'
      do J = 1,NB4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Phot ',J,
     &   WBIN(J),FBIN(J),RBIN(J),YBIN(J),ABIN(J),
     &   BBIN(J),CBIN(J),DBIN(J)
      enddo

        TITLNEW = 'solflx'
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',FBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',FBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',FBIN(13:18),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    d',FBIN(19:24)
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    e',FBIN(25:30)

        TITLNEW = 'X-Rayl'
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',RBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',RBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',RBIN(13:18),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    d',RBIN(19:24)
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    e',RBIN(25:30)

        TITLNEW = 'w-eff '
        write(6,'(a)') TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    a',WBIN(1:6),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    b',WBIN(7:12),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    c',WBIN(13:18),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    d',WBIN(19:24)
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    e',WBIN(25:30)

!! squeeze Y-PAR all into bin 18:
        YBIN(18) = (YBIN(18)*FBIN(18) + YBIN(19)*FBIN(19))/FBIN(18)
        YBIN(19) = 0.0

        TITLNEW = 'Y-PAR '
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',YBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',YBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',YBIN(13:18),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    d',YBIN(19:24)
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    e',YBIN(25:30)




!!!! Repeat integration over wavelengths for Watts/m2:  FX() and WX()
        FFBIN(:) = 0.d0
        WWBIN(:) = 0.d0
        RRBIN(:) = 0.d0
        YYBIN(:) = 0.d0
        AABIN(:) = 0.d0
        BBBIN(:) = 0.d0
        CCBIN(:) = 0.d0
        DDBIN(:) = 0.d0
       J = 1
       I = 1
      do while (WX(I) .lt. WCBIN(J))
       I = I+1
      enddo
      do J = 1,NC4
       do while (WX(I) .lt. WCBIN(J+1) .and. I .lt. NX_)
        WL = WX(I)
        call NDXWATER(IU_, WL, TL_, RNW, CNW, A1,A2)
        call NDXICE  (IU_, WL, TL_, RNI, CNI, A1,A2)
        FFBIN(J) = FFBIN(J) + FX(I)
        WWBIN(J) = WWBIN(J) + FX(I)*WX(I)
        RRBIN(J) = RRBIN(J) + FX(I)*RAYLAY(WX(I))
        YYBIN(J) = YYBIN(J) + FX(I)*YPAR(WX(I))
        AABIN(J) = AABIN(J) + FX(I)*RNW
        BBBIN(J) = BBBIN(J) + FX(I)*CNW
        CCBIN(J) = CCBIN(J) + FX(I)*RNI
        DDBIN(J) = DDBIN(J) + FX(I)*CNI
        I = I+1
       enddo
        WWBIN(J) = WWBIN(J)/FFBIN(J)
        RRBIN(J) = RRBIN(J)/FFBIN(J)
        YYBIN(J) = YYBIN(J)/FFBIN(J)
        AABIN(J) = AABIN(J)/FFBIN(J)
        BBBIN(J) = BBBIN(J)/FFBIN(J)
        CCBIN(J) = CCBIN(J)/FFBIN(J)
        DDBIN(J) = DDBIN(J)/FFBIN(J)
       if (I .ge. NX_) goto 4
      enddo
    4 continue

!!!! NC4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
      write(6,'(2a)') ' prat bin#   w-eff  phot/Watt X-ray     Y-PAR',
     &    '     Liq:  nr + ni       Ice:  nr + ni'
      do J = 1,NC4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Watts',J,
     &   WWBIN(J),FFBIN(J),RRBIN(J),YYBIN(J),AABIN(J),
     &   BBBIN(J),CCBIN(J),DDBIN(J)
      enddo

!!!! Second integration from NC4 bins to NB4 bins (Fast-JX + Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        WBIN(:) = 0.d0
        RBIN(:) = 0.d0
        YBIN(:) = 0.d0
        ABIN(:) = 0.d0
        BBIN(:) = 0.d0
        CBIN(:) = 0.d0
        DBIN(:) = 0.d0
       do I=16,NC4
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          WBIN(J) = WBIN(J) + FFBIN(I)*WWBIN(I)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)
          YBIN(J) = YBIN(J) + FFBIN(I)*YYBIN(I)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)
          DBIN(J) = DBIN(J) + FFBIN(I)*DDBIN(I)
!          write(6,'(a,2i5,f8.3)') ' NC/NB# ', I,J
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)         *SRB(I,J)
          WBIN(J) = WBIN(J) + FFBIN(I)*WWBIN(I)*SRB(I,J)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)*SRB(I,J)
          YBIN(J) = YBIN(J) + FFBIN(I)*YYBIN(I)*SRB(I,J)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)*SRB(I,J)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)*SRB(I,J)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)*SRB(I,J)
          DBIN(J) = DBIN(J) + FFBIN(I)*DDBIN(I)*SRB(I,J)
!          write(6,'(a,2i5,f8.3)') ' NC/NB# ', I,J, SRB(I,J)
        enddo
       enddo
       do J=1,NB4
        if (FBIN(J) .gt. 0.d0) then
          WBIN(J) = WBIN(J)/FBIN(J)
          RBIN(J) = RBIN(J)/FBIN(J)
          YBIN(J) = YBIN(J)/FBIN(J)
          ABIN(J) = ABIN(J)/FBIN(J)
          BBIN(J) = BBIN(J)/FBIN(J)
          CBIN(J) = CBIN(J)/FBIN(J)
          DBIN(J) = DBIN(J)/FBIN(J)
        endif
       enddo

!!!! NB4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
      write(6,'(2a)') ' S-J bin#    w-eff  phot/Watt X-ray     Y-PAR',
     &    '     Liq:  nr + ni       Ice:  nr + ni'
      do J = 1,NB4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Watts',J,
     &   WBIN(J),FBIN(J),RBIN(J),YBIN(J),ABIN(J),
     &   BBIN(J),CBIN(J),DBIN(J)
      enddo

        TITLNEW = 'solarW'
        write(6,'(a)') TITLNEW
        write(6,'(a5,   6f10.5,1x,a6)')
     &  '    a',FBIN(1:6),TITLNEW
        write(6,'(a5,   6f10.5,1x,a6)')
     &  '    b',FBIN(7:12),TITLNEW
        write(6,'(a5,   6f10.3,1x,a6)')
     &  '    c',FBIN(13:18),TITLNEW
        write(6,'(a5,   6f10.3,1x,a6)')
     &  '    d',FBIN(19:24)
        write(6,'(a5,   6f10.3,1x,a6)')
     &  '    e',FBIN(25:30)

        TITLNEW = 'X-rayl'
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',RBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',RBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',RBIN(13:18),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    d',RBIN(19:24)
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    e',RBIN(25:30)

        TITLNEW = 'w-eff '
        write(6,'(a)') TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    a',WBIN(1:6),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    b',WBIN(7:12),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    c',WBIN(13:18),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    d',WBIN(19:24)
        write(6,'(a5,6f10.5,1x,a6)')
     &  '    e',WBIN(25:30)




      stop
      end


      function RAYLAY(WAVE)
      REAL*8 WAVE, RAYLAY
C-----CALCULATE RAYLEIGH CROSS-SECTION AT WAVE (microns)
C---RAYLEIGH+RAMAN CROSS-SECTION (INCLUDE FOR ALL WAVELENGTHS)
      if (WAVE .lt. 0.170d0) then
          RAYLAY = 1.d-24
      else
       WSQI = 1.d0/(WAVE*WAVE)
       REFRM1 = 1.0d-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
       RAYLAY = 5.40d-21*(REFRM1*WSQI)**2
      endif
      return
      end


c>>>>>>>>>>>>>>>>>>>>>>>added Xsection<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      function YPAR(WAVE)
c---Photosynthetically Active Radiation: action spectrum (quantum): Y-PAR
c---traced from:
c      McCree, Keith J. (1972a). "The action spectrum, absorptance and
c        quantum yield of photosynthesis in crop plants"
c        Agricultural and Forest Meteorology 9:191-216.
c      McCree, Keith J. (1972b). Agric. & Forest Meteorology 10:443-453.
c---PAR in PAR is normally quantified as µmol photons/m2/s =? µE/m2/s
c        photosynthetic photon flux (area) density, or PPFD.

      real*8 YPAR
      real*8, intent(in) :: WAVE
      integer IWW
      real*8 FWW,WWI
!      real*8, dimension(18), parameter :: W = [325.d0,350.d0,375.d0,
!     &     400.d0,425.d0,450.d0,475.d0,500.d0,525.d0,550.d0,575.d0,
!     &     600.d0,625.d0,650.d0,675.d0,700.d0,725.d0,750.d0]
      real*8, dimension(18), parameter :: Y = [0.d0,15.d-2,45.d-2,
     &     64.d-2,78.d-2,75.d-2,68.d-2,70.d-2,74.d-2,88.d-2,95.d-2,
     &     100.d-2,100.d-2,94.d-2,92.d-2,43.d-2,4.d-2,0.d0]

      WWI = 0.04d0*(WAVE*1000.d0 - 300.d0)
        IWW = WWI
        IWW = max( 1, min( 17, IWW))
        FWW = WWI - float(IWW)
        FWW = max( 0.d0, min( 1.d0, FWW))
      YPAR = Y(IWW) + (Y(IWW+1)-Y(IWW))*FWW
      return
      end



      subroutine NDXWATER(IUNIT, XLAM, T, RN, CN, ABSIND, ABSCOF)

C     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR WATER
C     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM .2 MICRONS TO 1 M  [BTJ]
C     10 CM - 1 M IS UNKNOWN, AND WAS ADDED DUE TO LIEBE MODEL ADDITION [BTJ]
C     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 0.1 CM
C
C     ERIC A. SMITH
C     DEPT OF ATMOSPHERIC SCIENCE
C     COLORADO STATE UNIVERSITY
C     FORT COLLINS,CO  80523
C     TEL   303-491-8533
C
C     <> Modifications by Benjamin T. Johnson : August 1998 : Purdue
C     <> University Department of Earth and Atmospheric Sciences.
C     <> See comments with [BTJ] for changes/notes.
C     <> new email: jbenjam@aos.wisc.edu
C
C     <> Further modifications by Michael A. Walters, October 2002
C     <> Microwave Group, Space Science and Engineering Center
C     <> University of Wisconsin - Madison
C     <> walters@aos.wisc.edu
C
C     REFERENCES
C
C     0.2 UM - 0.69 UM
C
C     HALE,G., AND M. QUERRY,1972.
C     OPTICAL CONSTANTS OF WATER IN THE 200 NM TO 200 UM WAVELENGTH REGI
C     APPLIED OPTICS,12,3,555-563.
C
C     0.69 UM - 2.0 UM
C
C     PALMER,K.F., AND D. WILLIAMS,1974.
C     OPTICAL PROPERTIES OF WATER IN THE NEAR INFRARED.
C     JOURNAL OF THE OPTICAL SOCIETY OF AMERICA,64,8,1107-1110.
C
C     2.0 UM - 1000.0 UM
C
C     DOWNING,H.D., AND D. WILLIAMS,1975.
C     OPTICAL CONSTANTS OF WATER IN THE INFRARED.
C     JOURNAL OF GEOPHYSICAL REVIEW,80,12,1656-1661.
C
C     -------------------------------------------------------------------
C     See comments below for the 1.0 MM - 10.0 CM range. [~line 356] [BTJ]
C     -------------------------------------------------------------------
C     1.0 MM - 10.0 CM
C
C     RAY,P.S.,1972.
C     BROADBAND COMPLEX REFRACTIVE INDICES OF ICE AND WATER.
C     APPLIED OPTICS,11,8,1836-1844.
C
C     INPUT PARAMETERS
C
C     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
C           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
C           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
C           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N)
C     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
C     T = TEMPERATURE ( DEGREES KELVIN )
C
C     OUTPUT PARAMETERS
C
C     RN = REAL PORTION ( SCATTERING )
C     CN = COMPLEX PORTION ( ABSORPTION )
C     ABSIND = ABSORPTIVE INDEX ( CN/RN )
C     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
C
      COMPLEX E, M, eps_swd_l91dd
      DIMENSION WLTABW(518), RNTABW(518), CNTABW(518)
      REAL NU, CC
      DATA NUMWAT /518/
      DATA WLMIN, WLMAX /0.2, 1000000.0/
      DATA CUTWAT /1000.0/
      DATA (WLTABW(I),I=  1, 66)/
     *    .20000,    .22500,    .25000,    .27500,    .30000,    .32500,
     *    .35001,    .37500,    .40000,    .42501,    .45000,    .47499,
     *    .50000,    .52499,    .54999,    .57501,    .59999,    .62500,
     *    .64998,    .67499,    .68966,    .70175,    .71429,    .72464,
     *    .73529,    .74627,    .75188,    .75758,    .76923,    .78125,
     *    .79365,    .80645,    .81301,    .81967,    .83333,    .84746,
     *    .86207,    .87719,    .89286,    .90909,    .92593,    .93458,
     *    .94340,    .95238,    .96154,    .97276,    .98039,    .99010,
     *   1.00000,   1.01010,   1.02041,   1.03093,   1.04167,   1.05263,
     *   1.06952,   1.08696,   1.09890,   1.11111,   1.12360,   1.13636,
     *   1.14943,   1.16279,   1.17647,   1.19048,   1.20482,   1.21951/
      DATA (WLTABW(I),I= 67,132)/
     *   1.23457,   1.25000,   1.26582,   1.28205,   1.29870,   1.31579,
     *   1.33333,   1.35135,   1.36986,   1.38889,   1.40845,   1.42857,
     *   1.44300,   1.47059,   1.49254,   1.51515,   1.53846,   1.56250,
     *   1.58730,   1.61290,   1.63934,   1.66667,   1.69492,   1.72414,
     *   1.75439,   1.78571,   1.80180,   1.81818,   1.85185,   1.88679,
     *   1.92678,   1.96078,   2.00000,   2.02020,   2.04082,   2.06186,
     *   2.08333,   2.10526,   2.12766,   2.15054,   2.17391,   2.19780,
     *   2.22222,   2.24719,   2.27273,   2.29885,   2.32558,   2.35294,
     *   2.38095,   2.40964,   2.43902,   2.46914,   2.50000,   2.50627,
     *   2.51256,   2.51889,   2.52525,   2.53165,   2.53807,   2.54453,
     *   2.55102,   2.55754,   2.56410,   2.57069,   2.57732,   2.58398/
      DATA (WLTABW(I),I=133,198)/
     *   2.59067,   2.59740,   2.60417,   2.61097,   2.61780,   2.62467,
     *   2.63158,   2.63852,   2.64550,   2.65252,   2.65957,   2.66667,
     *   2.67380,   2.68097,   2.68817,   2.69542,   2.70270,   2.71003,
     *   2.71739,   2.72480,   2.73224,   2.73973,   2.74725,   2.75482,
     *   2.76243,   2.77008,   2.77778,   2.78552,   2.79330,   2.80112,
     *   2.80899,   2.81690,   2.82486,   2.83286,   2.84091,   2.84900,
     *   2.85714,   2.86533,   2.87356,   2.88184,   2.89017,   2.89855,
     *   2.90698,   2.91545,   2.92398,   2.93255,   2.94118,   2.94985,
     *   2.95858,   2.96736,   2.97619,   2.98507,   2.99401,   3.00300,
     *   3.01205,   3.02115,   3.03030,   3.03951,   3.04878,   3.05810,
     *   3.06748,   3.07692,   3.08642,   3.09598,   3.10559,   3.11526/
      DATA (WLTABW(I),I=199,264)/
     *   3.12500,   3.13480,   3.14465,   3.15457,   3.16456,   3.17460,
     *   3.18471,   3.19489,   3.20513,   3.21543,   3.22581,   3.23625,
     *   3.24675,   3.25733,   3.26797,   3.27869,   3.28947,   3.30033,
     *   3.31126,   3.32226,   3.33333,   3.34448,   3.35570,   3.36700,
     *   3.37838,   3.38983,   3.40136,   3.41297,   3.42466,   3.43643,
     *   3.44828,   3.46021,   3.47222,   3.48432,   3.49650,   3.50877,
     *   3.52113,   3.53357,   3.54610,   3.55872,   3.57143,   3.58423,
     *   3.59712,   3.61011,   3.62319,   3.63636,   3.64964,   3.66300,
     *   3.67647,   3.69004,   3.70370,   3.71747,   3.73134,   3.74532,
     *   3.75940,   3.77358,   3.78788,   3.80228,   3.81679,   3.83142,
     *   3.84615,   3.86100,   3.87597,   3.89105,   3.90625,   3.92157/
      DATA (WLTABW(I),I=265,330)/
     *   3.93701,   3.95257,   3.96825,   3.98406,   4.00000,   4.01606,
     *   4.03226,   4.04858,   4.06504,   4.08163,   4.09836,   4.11523,
     *   4.13223,   4.14938,   4.16667,   4.18410,   4.20168,   4.21941,
     *   4.23729,   4.25532,   4.27350,   4.29185,   4.31034,   4.32900,
     *   4.34783,   4.36681,   4.38596,   4.40529,   4.42478,   4.44444,
     *   4.46429,   4.48430,   4.50450,   4.52489,   4.54545,   4.56621,
     *   4.58716,   4.60829,   4.62963,   4.65116,   4.67290,   4.69484,
     *   4.71698,   4.73934,   4.76190,   4.78469,   4.80769,   4.83092,
     *   4.85437,   4.87805,   4.90196,   4.92611,   4.95050,   4.97512,
     *   5.00000,   5.02513,   5.05051,   5.07614,   5.10204,   5.12821,
     *   5.15464,   5.18135,   5.20833,   5.23560,   5.26316,   5.29101/
      DATA (WLTABW(I),I=331,396)/
     *   5.31915,   5.34759,   5.37634,   5.40541,   5.43478,   5.46448,
     *   5.49451,   5.52486,   5.55556,   5.58659,   5.61798,   5.64972,
     *   5.68182,   5.71429,   5.74713,   5.78035,   5.81395,   5.84795,
     *   5.88235,   5.91716,   5.95238,   5.98802,   6.02410,   6.06061,
     *   6.09756,   6.13497,   6.17284,   6.21118,   6.25000,   6.28931,
     *   6.32911,   6.36943,   6.41026,   6.45161,   6.49351,   6.53595,
     *   6.57895,   6.62252,   6.66667,   6.71141,   6.75676,   6.80272,
     *   6.84932,   6.89655,   6.94444,   6.99301,   7.04225,   7.09220,
     *   7.14286,   7.19424,   7.24638,   7.29927,   7.35294,   7.40741,
     *   7.46269,   7.51880,   7.57576,   7.63359,   7.69231,   7.75194,
     *   7.81250,   7.87402,   7.93651,   8.00000,   8.06452,   8.13008/
      DATA (WLTABW(I),I=397,462)/
     *   8.19672,   8.26446,   8.33333,   8.40336,   8.47458,   8.54701,
     *   8.62069,   8.69565,   8.77193,   8.84956,   8.92857,   9.00901,
     *   9.09091,   9.17431,   9.25926,   9.34579,   9.43396,   9.52381,
     *   9.61538,   9.70874,   9.80392,   9.90099,  10.00000,  10.10101,
     *  10.20408,  10.30928,  10.41667,  10.52632,  10.63830,  10.75269,
     *  10.86957,  10.98901,  11.11111,  11.23596,  11.36364,  11.49425,
     *  11.62791,  11.76471,  11.90476,  12.04819,  12.19512,  12.34568,
     *  12.50000,  12.65823,  12.82051,  12.98701,  13.15789,  13.33333,
     *  13.51351,  13.69863,  13.88889,  14.08451,  14.28571,  14.49275,
     *  14.70588,  14.92537,  15.15152,  15.38462,  15.62500,  15.87302,
     *  16.12903,  16.39344,  16.66667,  16.94915,  17.24138,  17.54386/
      DATA (WLTABW(I),I=463,518)/
     *  17.85714,  18.18182,  18.51852,  18.86792,  19.23077,  19.60784,
     *  20.00000,  20.40816,  20.83333,  21.27660,  21.73913,  22.22222,
     *  22.72727,  23.25581,  23.80952,  24.39024,  25.00000,  25.64103,
     *  26.31579,  27.02703,  27.77778,  28.57143,  29.41176,  30.30303,
     *  31.25000,  32.25806,  33.33333,  34.48276,  35.71429,  37.03704,
     *  38.46154,  40.00000,  41.66667,  43.47826,  45.45455,  47.61905,
     *  50.00000,  52.63158,  55.55556,  58.82353,  62.50000,  66.66667,
     *  71.42857,  76.92308,  83.33333,  90.90909, 100.00000, 111.11111,
     * 125.00000, 142.85714, 166.66667, 200.00000, 250.00000, 333.33333,
     * 500.00000,1000.00000/
      DATA (RNTABW(I),I=  1, 66)/
     *1.396,1.373,1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,
     *1.336,1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.332,1.332,
     *1.332,1.332,1.332,1.332,1.332,1.332,1.331,1.331,1.331,1.331,1.331,
     *1.330,1.330,1.330,1.330,1.330,1.329,1.329,1.329,1.329,1.329,1.328,
     *1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,
     *1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.326,1.325,1.325,1.325/
      DATA (RNTABW(I),I= 67,132)/
     *1.325,1.325,1.324,1.324,1.324,1.324,1.323,1.323,1.323,1.322,1.322,
     *1.321,1.321,1.321,1.320,1.320,1.319,1.319,1.318,1.318,1.317,1.316,
     *1.315,1.314,1.314,1.313,1.312,1.312,1.311,1.310,1.309,1.307,1.306,
     *1.301,1.301,1.300,1.298,1.298,1.296,1.295,1.294,1.293,1.291,1.289,
     *1.287,1.285,1.282,1.280,1.277,1.274,1.270,1.265,1.261,1.260,1.259,
     *1.257,1.256,1.255,1.254,1.252,1.250,1.249,1.247,1.246,1.243,1.241/
      DATA (RNTABW(I),I=133,198)/
     *1.240,1.238,1.235,1.232,1.230,1.227,1.224,1.221,1.218,1.214,1.210,
     *1.205,1.200,1.195,1.191,1.185,1.179,1.172,1.166,1.157,1.149,1.144,
     *1.139,1.138,1.138,1.139,1.141,1.144,1.149,1.154,1.158,1.161,1.165,
     *1.171,1.177,1.183,1.191,1.199,1.212,1.220,1.233,1.246,1.258,1.271,
     *1.282,1.293,1.305,1.317,1.329,1.342,1.353,1.364,1.376,1.386,1.398,
     *1.407,1.417,1.426,1.434,1.442,1.450,1.457,1.465,1.471,1.476,1.480/
      DATA (RNTABW(I),I=199,264)/
     *1.483,1.486,1.487,1.487,1.487,1.486,1.485,1.482,1.479,1.477,1.474,
     *1.472,1.467,1.464,1.461,1.457,1.454,1.451,1.448,1.444,1.441,1.437,
     *1.434,1.431,1.427,1.425,1.421,1.418,1.415,1.413,1.410,1.407,1.405,
     *1.403,1.400,1.398,1.396,1.394,1.392,1.390,1.388,1.387,1.385,1.383,
     *1.382,1.379,1.378,1.377,1.375,1.374,1.372,1.371,1.370,1.369,1.367,
     *1.366,1.365,1.363,1.361,1.361,1.360,1.358,1.358,1.357,1.355,1.354/
      DATA (RNTABW(I),I=265,330)/
     *1.353,1.352,1.351,1.350,1.349,1.348,1.348,1.347,1.346,1.345,1.344,
     *1.344,1.343,1.342,1.341,1.340,1.340,1.338,1.337,1.337,1.335,1.334,
     *1.334,1.333,1.332,1.332,1.331,1.330,1.330,1.330,1.329,1.329,1.329,
     *1.328,1.328,1.327,1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.325,
     *1.325,1.325,1.325,1.325,1.325,1.324,1.324,1.323,1.322,1.322,1.321,
     *1.320,1.319,1.318,1.318,1.317,1.316,1.314,1.313,1.311,1.310,1.308/
      DATA (RNTABW(I),I=331,396)/
     *1.306,1.304,1.302,1.299,1.297,1.294,1.291,1.288,1.285,1.282,1.278,
     *1.275,1.271,1.267,1.262,1.256,1.251,1.247,1.242,1.241,1.241,1.247,
     *1.265,1.289,1.311,1.332,1.349,1.354,1.356,1.354,1.350,1.345,1.341,
     *1.337,1.333,1.330,1.326,1.324,1.322,1.320,1.319,1.318,1.317,1.316,
     *1.315,1.314,1.313,1.311,1.310,1.309,1.308,1.307,1.306,1.305,1.303,
     *1.302,1.301,1.300,1.298,1.296,1.295,1.294,1.293,1.291,1.288,1.286/
      DATA (RNTABW(I),I=397,462)/
     *1.285,1.283,1.281,1.279,1.276,1.274,1.271,1.269,1.267,1.264,1.261,
     *1.259,1.256,1.253,1.249,1.246,1.242,1.238,1.234,1.230,1.224,1.220,
     *1.214,1.208,1.202,1.194,1.189,1.181,1.174,1.168,1.162,1.156,1.149,
     *1.143,1.139,1.135,1.132,1.132,1.131,1.132,1.130,1.130,1.134,1.138,
     *1.142,1.157,1.171,1.182,1.189,1.201,1.213,1.223,1.236,1.249,1.264,
     *1.277,1.289,1.303,1.313,1.324,1.335,1.348,1.361,1.372,1.385,1.396/
      DATA (RNTABW(I),I=463,518)/
     *1.407,1.419,1.431,1.441,1.451,1.462,1.470,1.480,1.488,1.496,1.504,
     *1.510,1.515,1.521,1.527,1.532,1.537,1.541,1.545,1.549,1.552,1.552,
     *1.552,1.550,1.546,1.543,1.541,1.539,1.537,1.534,1.532,1.529,1.525,
     *1.528,1.542,1.567,1.600,1.640,1.689,1.746,1.801,1.848,1.890,1.929,
     *1.960,1.982,1.997,2.000,2.010,2.020,2.040,2.070,2.110,2.150,2.225,
     *2.481/
      DATA (CNTABW(I),I=  1, 66)/
     *1.1000E-07,4.9000E-08,3.4000E-08,2.4000E-08,1.6000E-08,1.1000E-08,
     *6.5000E-09,3.5000E-09,1.9000E-09,1.3000E-09,1.0000E-09,9.4000E-10,
     *1.0000E-09,1.3000E-09,2.0000E-09,3.6000E-09,1.1000E-08,1.4000E-08,
     *1.6000E-08,2.2000E-08,2.7000E-08,3.8000E-08,5.6000E-08,7.7300E-08,
     *1.3900E-07,1.6300E-07,1.6800E-07,1.6400E-07,1.5400E-07,1.4300E-07,
     *1.3300E-07,1.2500E-07,1.2400E-07,1.3000E-07,2.0400E-07,2.6100E-07,
     *2.9400E-07,3.5300E-07,4.3300E-07,5.4300E-07,8.7700E-07,1.1800E-06,
     *1.6100E-06,2.4400E-06,3.6000E-06,3.9800E-06,3.9200E-06,3.7000E-06,
     *3.3100E-06,2.8200E-06,2.3100E-06,1.9000E-06,1.5700E-06,1.3700E-06,
     *1.2600E-06,1.4400E-06,1.6800E-06,2.0500E-06,2.8900E-06,4.9600E-06,
     *8.8700E-06,1.0900E-05,1.1500E-05,1.1800E-05,1.2000E-05,1.1800E-05/
      DATA (CNTABW(I),I= 67,132)/
     *1.1500E-05,1.1000E-05,1.0800E-05,1.1500E-05,1.3800E-05,1.7500E-05,
     *2.3900E-05,4.1600E-05,5.9400E-05,1.0100E-04,2.4100E-04,3.5200E-04,
     *3.6400E-04,3.3400E-04,2.5800E-04,1.8800E-04,1.4800E-04,1.2000E-04,
     *1.0200E-04,8.7300E-05,7.9200E-05,7.4900E-05,7.6200E-05,8.5500E-05,
     *1.0600E-04,1.3000E-04,1.3600E-04,1.3700E-04,1.5900E-04,8.6300E-04,
     *1.9000E-03,1.7000E-03,1.1000E-03,9.0000E-04,7.3100E-04,6.1700E-04,
     *5.1400E-04,4.5200E-04,4.0000E-04,3.5900E-04,3.4100E-04,3.3800E-04,
     *3.4500E-04,3.7600E-04,4.1600E-04,4.6500E-04,5.4200E-04,6.5200E-04,
     *7.9200E-04,9.6800E-04,1.2300E-03,1.5600E-03,1.9000E-03,1.9500E-03,
     *2.0000E-03,2.0500E-03,2.0700E-03,2.1000E-03,2.1200E-03,2.1500E-03,
     *2.1900E-03,2.2400E-03,2.2700E-03,2.3100E-03,2.3400E-03,2.3900E-03/
      DATA (CNTABW(I),I=133,198)/
     *2.4300E-03,2.4800E-03,2.5700E-03,2.7000E-03,2.9800E-03,3.3000E-03,
     *4.0200E-03,4.3700E-03,4.8200E-03,5.3600E-03,6.2700E-03,7.3200E-03,
     *8.5500E-03,1.0500E-02,1.2700E-02,1.4500E-02,1.6400E-02,1.8600E-02,
     *2.0500E-02,2.8200E-02,3.8000E-02,4.6200E-02,5.4800E-02,6.4900E-02,
     *7.4400E-02,8.3600E-02,9.2700E-02,1.0200E-01,1.1200E-01,1.2100E-01,
     *1.3100E-01,1.4200E-01,1.5400E-01,1.6700E-01,1.8000E-01,1.9400E-01,
     *2.0600E-01,2.1800E-01,2.2900E-01,2.3900E-01,2.4900E-01,2.5800E-01,
     *2.6500E-01,2.7100E-01,2.7600E-01,2.8000E-01,2.8100E-01,2.8200E-01,
     *2.8200E-01,2.7900E-01,2.7600E-01,2.7200E-01,2.6700E-01,2.6200E-01,
     *2.5500E-01,2.5000E-01,2.4300E-01,2.3600E-01,2.2800E-01,2.2000E-01,
     *2.1200E-01,2.0400E-01,1.9500E-01,1.8300E-01,1.7300E-01,1.6300E-01/
      DATA (CNTABW(I),I=199,264)/
     *1.5300E-01,1.4400E-01,1.3400E-01,1.2500E-01,1.1700E-01,1.1000E-01,
     *9.9400E-02,9.2000E-02,8.5500E-02,7.8500E-02,7.1600E-02,6.5300E-02,
     *6.0000E-02,5.5000E-02,5.0400E-02,4.6200E-02,4.2200E-02,3.8500E-02,
     *3.4800E-02,3.1500E-02,2.9700E-02,2.7900E-02,2.6200E-02,2.5000E-02,
     *2.2900E-02,2.1000E-02,1.9300E-02,1.7700E-02,1.6300E-02,1.5100E-02,
     *1.3800E-02,1.2800E-02,1.1800E-02,1.1000E-02,1.0100E-02,9.4100E-03,
     *8.6600E-03,8.0700E-03,7.3700E-03,6.8300E-03,6.2500E-03,5.7900E-03,
     *5.3800E-03,5.0600E-03,4.7300E-03,4.4900E-03,4.2400E-03,4.0500E-03,
     *3.8900E-03,3.7600E-03,3.6300E-03,3.5500E-03,3.4700E-03,3.4000E-03,
     *3.3500E-03,3.3600E-03,3.3500E-03,3.3900E-03,3.4000E-03,3.4800E-03,
     *3.5200E-03,3.6300E-03,3.7000E-03,3.7800E-03,3.8900E-03,3.9900E-03/
      DATA (CNTABW(I),I=265,330)/
     *4.1000E-03,4.2200E-03,4.3300E-03,4.5000E-03,4.6500E-03,4.7900E-03,
     *4.9400E-03,5.1200E-03,5.3100E-03,5.4900E-03,5.6800E-03,5.8600E-03,
     *6.0800E-03,6.3100E-03,6.5300E-03,6.7300E-03,6.9600E-03,7.2200E-03,
     *7.4900E-03,7.7900E-03,8.0600E-03,8.3300E-03,8.6400E-03,8.9600E-03,
     *9.2700E-03,9.6600E-03,1.0000E-02,1.0400E-02,1.0800E-02,1.1200E-02,
     *1.1700E-02,1.2200E-02,1.2600E-02,1.3100E-02,1.3600E-02,1.4000E-02,
     *1.4500E-02,1.4900E-02,1.5200E-02,1.5400E-02,1.5600E-02,1.5700E-02,
     *1.5700E-02,1.5700E-02,1.5500E-02,1.5300E-02,1.5100E-02,1.4800E-02,
     *1.4600E-02,1.4300E-02,1.4000E-02,1.3700E-02,1.3300E-02,1.2900E-02,
     *1.2600E-02,1.2200E-02,1.1800E-02,1.1500E-02,1.1000E-02,1.0800E-02,
     *1.0500E-02,1.0300E-02,1.0100E-02,1.0000E-02,9.9300E-03,9.9000E-03/
      DATA (CNTABW(I),I=331,396)/
     *9.9500E-03,1.0000E-02,1.0200E-02,1.0400E-02,1.0700E-02,1.1000E-02,
     *1.1500E-02,1.2000E-02,1.2800E-02,1.3800E-02,1.5000E-02,1.6600E-02,
     *1.8500E-02,2.0500E-02,2.4200E-02,2.9300E-02,3.3200E-02,4.2900E-02,
     *5.4400E-02,6.8800E-02,8.4000E-02,1.0210E-01,1.1700E-01,1.3000E-01,
     *1.3200E-01,1.2400E-01,1.0600E-01,8.8000E-02,7.4000E-02,6.1800E-02,
     *5.3500E-02,4.8400E-02,4.4700E-02,4.2000E-02,3.9800E-02,3.8300E-02,
     *3.7300E-02,3.7000E-02,3.6600E-02,3.6300E-02,3.6000E-02,3.5700E-02,
     *3.5500E-02,3.5200E-02,3.5000E-02,3.4700E-02,3.4600E-02,3.4300E-02,
     *3.4200E-02,3.4200E-02,3.4200E-02,3.4300E-02,3.4200E-02,3.4200E-02,
     *3.4200E-02,3.4200E-02,3.4200E-02,3.4400E-02,3.4500E-02,3.4600E-02,
     *3.4900E-02,3.5100E-02,3.5100E-02,3.5100E-02,3.5200E-02,3.5600E-02/
      DATA (CNTABW(I),I=397,462)/
     *3.5900E-02,3.6100E-02,3.6200E-02,3.6600E-02,3.7000E-02,3.7400E-02,
     *3.7800E-02,3.8300E-02,3.8700E-02,3.9200E-02,3.9800E-02,4.0500E-02,
     *4.1100E-02,4.1700E-02,4.2400E-02,4.3400E-02,4.4300E-02,4.5300E-02,
     *4.6700E-02,4.8100E-02,4.9700E-02,5.1500E-02,5.3400E-02,5.5700E-02,
     *5.8900E-02,6.2200E-02,6.6100E-02,7.0700E-02,7.6400E-02,8.2800E-02,
     *8.9800E-02,9.7300E-02,1.0700E-01,1.1800E-01,1.3000E-01,1.4400E-01,
     *1.5900E-01,1.7600E-01,1.9200E-01,2.0800E-01,2.2600E-01,2.4300E-01,
     *2.6000E-01,2.7700E-01,2.9200E-01,3.0500E-01,3.1700E-01,3.2800E-01,
     *3.3800E-01,3.4700E-01,3.5600E-01,3.6500E-01,3.7300E-01,3.7900E-01,
     *3.8600E-01,3.9200E-01,3.9700E-01,4.0300E-01,4.0800E-01,4.1200E-01,
     *4.1700E-01,4.2000E-01,4.2300E-01,4.2500E-01,4.2700E-01,4.2800E-01/
      DATA (CNTABW(I),I=463,518)/
     *4.2700E-01,4.2700E-01,4.2600E-01,4.2500E-01,4.2300E-01,4.2100E-01,
     *4.1800E-01,4.1500E-01,4.1100E-01,4.0800E-01,4.0400E-01,4.0100E-01,
     *3.9700E-01,3.9400E-01,3.9000E-01,3.8600E-01,3.8200E-01,3.7700E-01,
     *3.7200E-01,3.6800E-01,3.6300E-01,3.5900E-01,3.5600E-01,3.5200E-01,
     *3.5300E-01,3.5700E-01,3.6100E-01,3.6800E-01,3.7500E-01,3.8500E-01,
     *3.9800E-01,4.1400E-01,4.3600E-01,4.6900E-01,5.0500E-01,5.3900E-01,
     *5.7100E-01,5.9700E-01,6.1800E-01,6.2900E-01,6.2200E-01,6.0800E-01,
     *5.9300E-01,5.7700E-01,5.5700E-01,5.3200E-01,5.0700E-01,4.8700E-01,
     *4.6600E-01,4.5000E-01,4.4400E-01,4.3800E-01,4.6000E-01,5.2700E-01,
     *7.1800E-01,8.4657E-01/
      DATA PI/3.14159265/
C
C     FUNCTION FOR TREATING ABSORPTION BANDS NOT CONSIDERED IN THE
C     DEBYE THEOREY
C
C$$$      ABSUM(WL, WLCEN, BET, DEL, GAM)=
C$$$     &     BET * EXP(- ABS(ALOG10(WL/WLCEN)/DEL)**GAM)
C
C     ZERO PARAMETERS
C
      RN = 0.0
      CN = 0.0
      ABSIND = 0.0
      ABSCOF = 0.0
C
C     CONVERT WAVELENGTH TO MICRONS
C
      WL = XLAM
      IF (IUNIT.EQ.1) WL = 1000*WL
      IF (IUNIT.EQ.2) WL = 10000*WL
      IF (IUNIT.EQ.3) WL =  10000*(1.0/WL)
C      IF ((WL .LT. WLMIN) .OR. (WL .GT. WLMAX)) RETURN
      IF (WL .GT. WLMAX) RETURN
C
C     REGION FROM 0.2 MICRON TO 1000.0 MICRON  -  TABLE LOOKUP
C
      IF (WL.GT.CUTWAT) GO TO 3
      DO 1 I = 2, NUMWAT
         IF (WL.GT.WLTABW(I)) GO TO 1
         I1 = I - 1
         I2 = I
         GO TO 2
 1    CONTINUE
      I1 = NUMWAT - 1
      I2 = NUMWAT
 2    FAC = (WL - WLTABW(I1)) / (WLTABW(I2) - WLTABW(I1))
         FAC = max(0.0, min(1.0, FAC))

      RN = RNTABW(I1) + FAC*(RNTABW(I2) - RNTABW(I1))
      CN = CNTABW(I1) + FAC*(CNTABW(I2) - CNTABW(I1))
      GO TO 5
C
C     REGION FROM 0.1 CM TO 10 CM
C
C     EXTENSION OF DEBYE THEOREY BASED ON THE WORK OF
C
C        COLE,K.S.,AND R.H.COLE,1941.JOUR.CHEM.PHYS.,9,P 341.
C
C     DEFINE TEMPERATURE TERMS AND WAVELENGTH IN CM
C
C     ---------------------------------------------------------------------
C     Note, the subroutine EPSW is replacing this section of code.  [BTJ]
C     It is an updated version for the microwave region, using the model
C     of Liebe et al. 1991.
C     Old statements have been commented, yet retained.
C     Some of the old statements were used, and contain lowercase
C     letters.  Compare to refwat.f by Eric A. Smith.
C     July 1998: Benjamin T. Johnson - Purdue University            [BTJ]
C     ---------------------------------------------------------------------

 3    TC = T - 273.15
C      T1=TC+273.0
C      T2=TC-25.0
C     Converts wavelength(WL) from microns to centimeters.          [BTJ]
C     However, to replace this section with EPSW, we need to convert
C     to Gigahertz.                                                 [BTJ]
C      XL=WL/10000.0
C     Conversion to Gigahertz
      XL = WL * 1.0e-6
      CC = 2.99792e+8
      NU = CC / XL
      NU = NU * 1.0e-9
C     Call the eps_swd_l91dd subroutine, passing TC and NU as temperature
C     in degrees celcius, and frequency in gigahertz (i.e. 85.5)    [BTJ]
C     ---------------------------------------------------------------------
      E = eps_swd_l91dd(NU, TC)    ! Replaced old function          [MAW]
C     ---------------------------------------------------------------------
C     ER and EI are carried into some equations below.
C     Most of the old stuff from here on is commented out.          [BTJ]
C
C     DEFINE FREQUENCY INDEPENDENT CONDUCTIVITY(SIGMA) AND
C     SPREAD PARAMETER(ALPHA)
C
C     IN CLASSICAL DEBYE THEOREY THESE TERMS ARE ZERO
C
C     SIGMA GIVEN BY SAXTON,J.A.,1949.WIRELESS ENGINEER,26,P 288.
C     ALPHA GIVEN BY RAY ( EQUATION 7B )
C
C      SIGMA=12.5664E8
C      ALPHA=-16.8129/T1+0.0609265
C
C     DEFINE STATIC DIELECTRIC CONSTANT(ES) - RAY EQN 4
C            HIGH FREQUENCY DIELECTRIC CONSTANT(E00) - RAY EQN 7A
C            RELAXTION WAVELENGTH IN CM(XLAMS) - RAY EQN 7C
C
C     TEMPERATURE DEPENDENCE OF ES GIVEN BY
C
C        WYMAN,J.,AND E.N.INGALLS,1938.JOUR.AM.CHEM.SOC.,60,P 1182.
C
C      ES=78.54*(1.0-4.579E-3*T2+1.19E-5*T2*T2-2.8E-8*T2*T2*T2)
C      E00=5.27137+0.0216474*TC-0.00131198*TC*TC
C      XLAMS=0.00033836*EXP(2513.98/T1)
C
C     CALCULATE EXPRESSIONS USED FOR DIELECTRIC CONSTANT
C
C      TERM=PI*ALPHA/2
C      SINT=SIN(TERM)
C      COST=COS(TERM)
C      XLRAT=XLAMS/XL
C      POWTRM=XLRAT**(1-ALPHA)
C      DENOM=1.0+2*POWTRM*SINT+XLRAT**(2.0*(1-ALPHA))
C
C     CALCULATION OF DIELECTRIC CONSTANT
C
C     REAL PART - RAY EQN 5
C
C      ER=E00+(ES-E00)*(1.0+POWTRM*SINT)/DENOM
C
C     IMAGINARY PART OR LOSS TERM - RAY EQN 6
C
C      EI=(SIGMA*XL/18.8496E10)+(ES-E00)*POWTRM*COST/DENOM
C
C     COMPLEX PERMITTIVITY
C
c$$$
c$$$      E = CMPLX(ER, -EI)
C
C     COMPLEX INDEX OF REFRACTION - RAY EQN 1
C
      M = CSQRT(E)
      RN = REAL(M)
      CN = ABS(AIMAG(M))
C
C     CORRECTION TO IMAGINARY INDEX TO ACCOUNT FOR THE
C     REMAINING ABSORPTION BANDS - RAY EQN 8(TABLE 2)
C
c$$$       IF (WL .GT. 3000.0) GO TO 5
c$$$       CN = CN + ABSUM(WL, 17.0,0.39,0.45,1.3)
c$$$     *         + ABSUM(WL, 62.0,0.41,0.35,1.7)
c$$$     *         + ABSUM(WL,300.0,0.25,0.47,3.0)
C
C     This part no longer needed because of eps_swd_l91dd  [MAW]

C
C     ABSORPTIVE QUANITIES
C
 5    ABSIND = CN / RN
      ABSCOF = 4.0*PI * CN / WL

      RETURN
      END







      SUBROUTINE NDXICE(IUNIT, XLAM, T, RN, CN, ABSIND, ABSCOF)

C Arguments:
      INTEGER IUNIT
      REAL ABSCOF, ABSIND, CN, XLAM, RN, T
C Parameters:
      INTEGER I, LT1, LT2, NWL, NWLT
      PARAMETER (NWL = 468, NWLT = 62)
C Local variables:
      REAL ALAM, CUTICE, PI, T1, T2, TK, WLMAX, WLMIN, X, X1, X2,
     &   Y, Y1, Y2, YLO, YHI

      REAL
     &   TABIM(NWL),TABIMT(NWLT,4),TABRE(NWL),TABRET(NWLT,4),TEMREF(4),
     &   WL(NWL),WLT(NWLT)
C
C     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR ICE.
C     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM 0.045 MICRONS TO 8.6 METER
C     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 167 MICRONS.
C
C     INTERPOLATION IS DONE     RN  VS. LOG(XLAM)
C                               RN  VS.        T
C                           LOG(CN) VS. LOG(XLAM)
C                           LOG(CN) VS.        T
C
C     STEPHEN G. WARREN - 1983
C     DEPT. OF ATMOSPHERIC SCIENCES
C     UNIVERSITY OF WASHINGTON
C     SEATTLE, WA  98195
C
C     BASED ON
C
C        WARREN,S.G.,1984.
C        OPTICAL CONSTANTS OF ICE FROM THE ULTRAVIOLET TO THE MICROWAVE.
C        APPLIED OPTICS,23,1206-1225
C
C     INPUT PARAMETERS
C
C     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
C           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
C           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
C           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS  WAVE N
C     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
C     T = TEMPERATURE ( DEGREES KELVIN )
C
C     OUTPUT PARAMETERS
C
C     RN = REAL PORTION ( SCATTERING )
C     CN = COMPLEX PORTION ( ABSORPTION )
C     ABSIND = ABSORPTIVE INDEX ( CN/RN )
C     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
C
C      DIMENSION WL(NWL),WLT(NWLT)
C      DIMENSION TABRE(NWL),TABRET(NWLT,4),TABIM(NWL),TABIMT(NWLT,4)
C      DIMENSION TEMREF(4)
C
C     REFERENCE TEMPERATURES ARE -1.0,-5.0,-20.0, AND -60.0 DEG CENTIGRA
C
      DATA TEMREF /272.16, 268.16, 253.16, 213.16/

      DATA WLMIN, WLMAX /0.045, 8.6E6/
      DATA CUTICE /167.0/

      DATA (WL(I),I=1,114)/
     +0.4430E-01,0.4510E-01,0.4590E-01,0.4680E-01,0.4770E-01,0.4860E-01,
     +0.4960E-01,0.5060E-01,0.5170E-01,0.5280E-01,0.5390E-01,0.5510E-01,
     +0.5640E-01,0.5770E-01,0.5900E-01,0.6050E-01,0.6200E-01,0.6360E-01,
     +0.6530E-01,0.6700E-01,0.6890E-01,0.7080E-01,0.7290E-01,0.7380E-01,
     +0.7510E-01,0.7750E-01,0.8000E-01,0.8270E-01,0.8550E-01,0.8860E-01,
     +0.9180E-01,0.9300E-01,0.9540E-01,0.9920E-01,0.1033E+00,0.1078E+00,
     +0.1100E+00,0.1127E+00,0.1140E+00,0.1181E+00,0.1210E+00,0.1240E+00,
     +0.1272E+00,0.1295E+00,0.1305E+00,0.1319E+00,0.1333E+00,0.1348E+00,
     +0.1362E+00,0.1370E+00,0.1378E+00,0.1387E+00,0.1393E+00,0.1409E+00,
     +0.1425E+00,0.1435E+00,0.1442E+00,0.1450E+00,0.1459E+00,0.1468E+00,
     +0.1476E+00,0.1480E+00,0.1485E+00,0.1494E+00,0.1512E+00,0.1531E+00,
     +0.1540E+00,0.1550E+00,0.1569E+00,0.1580E+00,0.1589E+00,0.1610E+00,
     +0.1625E+00,0.1648E+00,0.1669E+00,0.1692E+00,0.1713E+00,0.1737E+00,
     +0.1757E+00,0.1779E+00,0.1802E+00,0.1809E+00,0.1821E+00,0.1833E+00,
     +0.1843E+00,0.1850E+00,0.1860E+00,0.1870E+00,0.1880E+00,0.1890E+00,
     +0.1900E+00,0.1910E+00,0.1930E+00,0.1950E+00,0.2100E+00,0.2500E+00,
     +0.3000E+00,0.3500E+00,0.4000E+00,0.4100E+00,0.4200E+00,0.4300E+00,
     +0.4400E+00,0.4500E+00,0.4600E+00,0.4700E+00,0.4800E+00,0.4900E+00,
     +0.5000E+00,0.5100E+00,0.5200E+00,0.5300E+00,0.5400E+00,0.5500E+00/
      DATA (WL(I),I=115,228)/
     +0.5600E+00,0.5700E+00,0.5800E+00,0.5900E+00,0.6000E+00,0.6100E+00,
     +0.6200E+00,0.6300E+00,0.6400E+00,0.6500E+00,0.6600E+00,0.6700E+00,
     +0.6800E+00,0.6900E+00,0.7000E+00,0.7100E+00,0.7200E+00,0.7300E+00,
     +0.7400E+00,0.7500E+00,0.7600E+00,0.7700E+00,0.7800E+00,0.7900E+00,
     +0.8000E+00,0.8100E+00,0.8200E+00,0.8300E+00,0.8400E+00,0.8500E+00,
     +0.8600E+00,0.8700E+00,0.8800E+00,0.8900E+00,0.9000E+00,0.9100E+00,
     +0.9200E+00,0.9300E+00,0.9400E+00,0.9500E+00,0.9600E+00,0.9700E+00,
     +0.9800E+00,0.9900E+00,0.1000E+01,0.1010E+01,0.1020E+01,0.1030E+01,
     +0.1040E+01,0.1050E+01,0.1060E+01,0.1070E+01,0.1080E+01,0.1090E+01,
     +0.1100E+01,0.1110E+01,0.1120E+01,0.1130E+01,0.1140E+01,0.1150E+01,
     +0.1160E+01,0.1170E+01,0.1180E+01,0.1190E+01,0.1200E+01,0.1210E+01,
     +0.1220E+01,0.1230E+01,0.1240E+01,0.1250E+01,0.1260E+01,0.1270E+01,
     +0.1280E+01,0.1290E+01,0.1300E+01,0.1310E+01,0.1320E+01,0.1330E+01,
     +0.1340E+01,0.1350E+01,0.1360E+01,0.1370E+01,0.1380E+01,0.1390E+01,
     +0.1400E+01,0.1410E+01,0.1420E+01,0.1430E+01,0.1440E+01,0.1449E+01,
     +0.1460E+01,0.1471E+01,0.1481E+01,0.1493E+01,0.1504E+01,0.1515E+01,
     +0.1527E+01,0.1538E+01,0.1563E+01,0.1587E+01,0.1613E+01,0.1650E+01,
     +0.1680E+01,0.1700E+01,0.1730E+01,0.1760E+01,0.1800E+01,0.1830E+01,
     +0.1840E+01,0.1850E+01,0.1855E+01,0.1860E+01,0.1870E+01,0.1890E+01/
      DATA (WL(I),I=229,342)/
     +0.1905E+01,0.1923E+01,0.1942E+01,0.1961E+01,0.1980E+01,0.2000E+01,
     +0.2020E+01,0.2041E+01,0.2062E+01,0.2083E+01,0.2105E+01,0.2130E+01,
     +0.2150E+01,0.2170E+01,0.2190E+01,0.2220E+01,0.2240E+01,0.2245E+01,
     +0.2250E+01,0.2260E+01,0.2270E+01,0.2290E+01,0.2310E+01,0.2330E+01,
     +0.2350E+01,0.2370E+01,0.2390E+01,0.2410E+01,0.2430E+01,0.2460E+01,
     +0.2500E+01,0.2520E+01,0.2550E+01,0.2565E+01,0.2580E+01,0.2590E+01,
     +0.2600E+01,0.2620E+01,0.2675E+01,0.2725E+01,0.2778E+01,0.2817E+01,
     +0.2833E+01,0.2849E+01,0.2865E+01,0.2882E+01,0.2899E+01,0.2915E+01,
     +0.2933E+01,0.2950E+01,0.2967E+01,0.2985E+01,0.3003E+01,0.3021E+01,
     +0.3040E+01,0.3058E+01,0.3077E+01,0.3096E+01,0.3115E+01,0.3135E+01,
     +0.3155E+01,0.3175E+01,0.3195E+01,0.3215E+01,0.3236E+01,0.3257E+01,
     +0.3279E+01,0.3300E+01,0.3322E+01,0.3345E+01,0.3367E+01,0.3390E+01,
     +0.3413E+01,0.3436E+01,0.3460E+01,0.3484E+01,0.3509E+01,0.3534E+01,
     +0.3559E+01,0.3624E+01,0.3732E+01,0.3775E+01,0.3847E+01,0.3969E+01,
     +0.4099E+01,0.4239E+01,0.4348E+01,0.4387E+01,0.4444E+01,0.4505E+01,
     +0.4547E+01,0.4560E+01,0.4580E+01,0.4719E+01,0.4904E+01,0.5000E+01,
     +0.5100E+01,0.5200E+01,0.5263E+01,0.5400E+01,0.5556E+01,0.5714E+01,
     +0.5747E+01,0.5780E+01,0.5814E+01,0.5848E+01,0.5882E+01,0.6061E+01,
     +0.6135E+01,0.6250E+01,0.6289E+01,0.6329E+01,0.6369E+01,0.6410E+01/
      DATA (WL(I),I=343,456)/
     +0.6452E+01,0.6494E+01,0.6579E+01,0.6667E+01,0.6757E+01,0.6897E+01,
     +0.7042E+01,0.7143E+01,0.7246E+01,0.7353E+01,0.7463E+01,0.7576E+01,
     +0.7692E+01,0.7812E+01,0.7937E+01,0.8065E+01,0.8197E+01,0.8333E+01,
     +0.8475E+01,0.8696E+01,0.8929E+01,0.9091E+01,0.9259E+01,0.9524E+01,
     +0.9804E+01,0.1000E+02,0.1020E+02,0.1031E+02,0.1042E+02,0.1053E+02,
     +0.1064E+02,0.1075E+02,0.1087E+02,0.1100E+02,0.1111E+02,0.1136E+02,
     +0.1163E+02,0.1190E+02,0.1220E+02,0.1250E+02,0.1282E+02,0.1299E+02,
     +0.1316E+02,0.1333E+02,0.1351E+02,0.1370E+02,0.1389E+02,0.1408E+02,
     +0.1429E+02,0.1471E+02,0.1515E+02,0.1538E+02,0.1563E+02,0.1613E+02,
     +0.1639E+02,0.1667E+02,0.1695E+02,0.1724E+02,0.1818E+02,0.1887E+02,
     +0.1923E+02,0.1961E+02,0.2000E+02,0.2041E+02,0.2083E+02,0.2222E+02,
     +0.2260E+02,0.2305E+02,0.2360E+02,0.2460E+02,0.2500E+02,0.2600E+02,
     +0.2857E+02,0.3100E+02,0.3333E+02,0.3448E+02,0.3564E+02,0.3700E+02,
     +0.3824E+02,0.3960E+02,0.4114E+02,0.4276E+02,0.4358E+02,0.4458E+02,
     +0.4550E+02,0.4615E+02,0.4671E+02,0.4736E+02,0.4800E+02,0.4878E+02,
     +0.5003E+02,0.5128E+02,0.5275E+02,0.5350E+02,0.5424E+02,0.5500E+02,
     +0.5574E+02,0.5640E+02,0.5700E+02,0.5746E+02,0.5840E+02,0.5929E+02,
     +0.6000E+02,0.6100E+02,0.6125E+02,0.6250E+02,0.6378E+02,0.6467E+02,
     +0.6558E+02,0.6655E+02,0.6760E+02,0.6900E+02,0.7053E+02,0.7300E+02/
      DATA (WL(I),I=457,468)/
     +0.7500E+02,0.7629E+02,0.8000E+02,0.8297E+02,0.8500E+02,0.8680E+02,
     +0.9080E+02,0.9517E+02,0.1000E+03,0.1200E+03,0.1500E+03,0.1670E+03/
      DATA  WLT/
     +                                 0.1670E+03,0.1778E+03,0.1884E+03,
     +0.1995E+03,0.2113E+03,0.2239E+03,0.2371E+03,0.2512E+03,0.2661E+03,
     +0.2818E+03,0.2985E+03,0.3162E+03,0.3548E+03,0.3981E+03,0.4467E+03,
     +0.5012E+03,0.5623E+03,0.6310E+03,0.7943E+03,0.1000E+04,0.1259E+04,
     +0.2500E+04,0.5000E+04,0.1000E+05,0.2000E+05,0.3200E+05,0.3500E+05,
     +0.4000E+05,0.4500E+05,0.5000E+05,0.6000E+05,0.7000E+05,0.9000E+05,
     +0.1110E+06,0.1200E+06,0.1300E+06,0.1400E+06,0.1500E+06,0.1600E+06,
     +0.1700E+06,0.1800E+06,0.2000E+06,0.2500E+06,0.2900E+06,0.3200E+06,
     +0.3500E+06,0.3800E+06,0.4000E+06,0.4500E+06,0.5000E+06,0.6000E+06,
     +0.6400E+06,0.6800E+06,0.7200E+06,0.7600E+06,0.8000E+06,0.8400E+06,
     +0.9000E+06,0.1000E+07,0.2000E+07,0.5000E+07,0.8600E+07/
      DATA (TABRE(I),I=1,114)/
     +   0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038,
     +   0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093,
     +   0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089,
     +   0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134,
     +   1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958,
     +   1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221,
     +   1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918,
     +   1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605,
     +   1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302,
     +   1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782,
     +   1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404,
     +   1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269,
     +   1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991,
     +   1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602,
     +   1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216,
     +   1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085,
     +   1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702,
     +   1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349,
     +   1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/
      DATA (TABRE(I),I=115,228)/
     +   1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909,
     +   1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756,
     +   1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625,
     +   1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509,
     +   1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402,
     +   1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299,
     +   1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199,
     +   1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099,
     +   1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997,
     +   1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891,
     +   1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781,
     +   1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666,
     +   1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542,
     +   1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406,
     +   1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252,
     +   1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130,
     +   1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784,
     +   1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137,
     +   1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/
      DATA (TABRE(I),I=229,342)/
     +   1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404,
     +   1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736,
     +   1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815,
     +   1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967,
     +   1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270,
     +   1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507,
     +   1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739,
     +   1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962,
     +   0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597,
     +   1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095,
     +   1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284,
     +   1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139,
     +   1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215,
     +   1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448,
     +   1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053,
     +   1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897,
     +   1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898,
     +   1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156,
     +   1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/
      DATA (TABRE(I),I=343,456)/
     +   1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953,
     +   1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542,
     +   1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934,
     +   1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510,
     +   1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808,
     +   1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779,
     +   1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170,
     +   1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317,
     +   1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807,
     +   1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097,
     +   1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070,
     +   1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955,
     +   1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233,
     +   1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175,
     +   1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877,
     +   1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736,
     +   1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596,
     +   1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713,
     +   1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/
      DATA (TABRE(I),I=457,468)/
     +   1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702,
     +   1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/
      DATA (TABRET(I,1),I=1,NWLT)/
     +                                    1.82961,   1.83258,   1.83149,
     +   1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,
     +   1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,
     +   1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,
     +   1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680,
     +   1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78800/
      DATA (TABRET(I,2),I=1,NWLT)/
     +                         1.82961,   1.83258,   1.83149,   1.82748,
     +   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,
     +   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,
     +   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300,
     +   1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78720/
      DATA(TABRET(I,3),I=1,NWLT)/
     +              1.82961,   1.83258,   1.83149,   1.82748,   1.82224,
     +   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,
     +   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,
     +   1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310,
     +   1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400,
     +   1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78400,   1.78450/
      DATA (TABRET(I,4),I=1,NWLT)/
     +   1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718,
     +   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214,
     +   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201,
     +   1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77800/
      DATA(TABIM(I),I=1,114)/
     +0.1640E+00,0.1730E+00,0.1830E+00,0.1950E+00,0.2080E+00,0.2230E+00,
     +0.2400E+00,0.2500E+00,0.2590E+00,0.2680E+00,0.2790E+00,0.2970E+00,
     +0.3190E+00,0.3400E+00,0.3660E+00,0.3920E+00,0.4160E+00,0.4400E+00,
     +0.4640E+00,0.4920E+00,0.5170E+00,0.5280E+00,0.5330E+00,0.5340E+00,
     +0.5310E+00,0.5240E+00,0.5100E+00,0.5000E+00,0.4990E+00,0.4680E+00,
     +0.3800E+00,0.3600E+00,0.3390E+00,0.3180E+00,0.2910E+00,0.2510E+00,
     +0.2440E+00,0.2390E+00,0.2390E+00,0.2440E+00,0.2470E+00,0.2240E+00,
     +0.1950E+00,0.1740E+00,0.1720E+00,0.1800E+00,0.1940E+00,0.2130E+00,
     +0.2430E+00,0.2710E+00,0.2890E+00,0.3340E+00,0.3440E+00,0.3820E+00,
     +0.4010E+00,0.4065E+00,0.4050E+00,0.3890E+00,0.3770E+00,0.3450E+00,
     +0.3320E+00,0.3150E+00,0.2980E+00,0.2740E+00,0.2280E+00,0.1980E+00,
     +0.1720E+00,0.1560E+00,0.1100E+00,0.8300E-01,0.5800E-01,0.2200E-01,
     +0.1000E-01,0.3000E-02,0.1000E-02,0.3000E-03,0.1000E-03,0.3000E-04,
     +0.1000E-04,0.3000E-05,0.1000E-05,0.7000E-06,0.4000E-06,0.2000E-06,
     +0.1000E-06,0.6377E-07,0.3750E-07,0.2800E-07,0.2400E-07,0.2200E-07,
     +0.1900E-07,0.1750E-07,0.1640E-07,0.1590E-07,0.1325E-07,0.8623E-08,
     +0.5504E-08,0.3765E-08,0.2710E-08,0.2510E-08,0.2260E-08,0.2080E-08,
     +0.1910E-08,0.1540E-08,0.1530E-08,0.1550E-08,0.1640E-08,0.1780E-08,
     +0.1910E-08,0.2140E-08,0.2260E-08,0.2540E-08,0.2930E-08,0.3110E-08/
      DATA(TABIM(I),I=115,228)/
     +0.3290E-08,0.3520E-08,0.4040E-08,0.4880E-08,0.5730E-08,0.6890E-08,
     +0.8580E-08,0.1040E-07,0.1220E-07,0.1430E-07,0.1660E-07,0.1890E-07,
     +0.2090E-07,0.2400E-07,0.2900E-07,0.3440E-07,0.4030E-07,0.4300E-07,
     +0.4920E-07,0.5870E-07,0.7080E-07,0.8580E-07,0.1020E-06,0.1180E-06,
     +0.1340E-06,0.1400E-06,0.1430E-06,0.1450E-06,0.1510E-06,0.1830E-06,
     +0.2150E-06,0.2650E-06,0.3350E-06,0.3920E-06,0.4200E-06,0.4440E-06,
     +0.4740E-06,0.5110E-06,0.5530E-06,0.6020E-06,0.7550E-06,0.9260E-06,
     +0.1120E-05,0.1330E-05,0.1620E-05,0.2000E-05,0.2250E-05,0.2330E-05,
     +0.2330E-05,0.2170E-05,0.1960E-05,0.1810E-05,0.1740E-05,0.1730E-05,
     +0.1700E-05,0.1760E-05,0.1820E-05,0.2040E-05,0.2250E-05,0.2290E-05,
     +0.3040E-05,0.3840E-05,0.4770E-05,0.5760E-05,0.6710E-05,0.8660E-05,
     +0.1020E-04,0.1130E-04,0.1220E-04,0.1290E-04,0.1320E-04,0.1350E-04,
     +0.1330E-04,0.1320E-04,0.1320E-04,0.1310E-04,0.1320E-04,0.1320E-04,
     +0.1340E-04,0.1390E-04,0.1420E-04,0.1480E-04,0.1580E-04,0.1740E-04,
     +0.1980E-04,0.2500E-04,0.5400E-04,0.1040E-03,0.2030E-03,0.2708E-03,
     +0.3511E-03,0.4299E-03,0.5181E-03,0.5855E-03,0.5899E-03,0.5635E-03,
     +0.5480E-03,0.5266E-03,0.4394E-03,0.3701E-03,0.3372E-03,0.2410E-03,
     +0.1890E-03,0.1660E-03,0.1450E-03,0.1280E-03,0.1030E-03,0.8600E-04,
     +0.8220E-04,0.8030E-04,0.8500E-04,0.9900E-04,0.1500E-03,0.2950E-03/
      DATA(TABIM(I),I=229,342)/
     +0.4687E-03,0.7615E-03,0.1010E-02,0.1313E-02,0.1539E-02,0.1588E-02,
     +0.1540E-02,0.1412E-02,0.1244E-02,0.1068E-02,0.8414E-03,0.5650E-03,
     +0.4320E-03,0.3500E-03,0.2870E-03,0.2210E-03,0.2030E-03,0.2010E-03,
     +0.2030E-03,0.2140E-03,0.2320E-03,0.2890E-03,0.3810E-03,0.4620E-03,
     +0.5480E-03,0.6180E-03,0.6800E-03,0.7300E-03,0.7820E-03,0.8480E-03,
     +0.9250E-03,0.9200E-03,0.8920E-03,0.8700E-03,0.8900E-03,0.9300E-03,
     +0.1010E-02,0.1350E-02,0.3420E-02,0.7920E-02,0.2000E-01,0.3800E-01,
     +0.5200E-01,0.6800E-01,0.9230E-01,0.1270E+00,0.1690E+00,0.2210E+00,
     +0.2760E+00,0.3120E+00,0.3470E+00,0.3880E+00,0.4380E+00,0.4930E+00,
     +0.5540E+00,0.6120E+00,0.6250E+00,0.5930E+00,0.5390E+00,0.4910E+00,
     +0.4380E+00,0.3720E+00,0.3000E+00,0.2380E+00,0.1930E+00,0.1580E+00,
     +0.1210E+00,0.1030E+00,0.8360E-01,0.6680E-01,0.5400E-01,0.4220E-01,
     +0.3420E-01,0.2740E-01,0.2200E-01,0.1860E-01,0.1520E-01,0.1260E-01,
     +0.1060E-01,0.8020E-02,0.6850E-02,0.6600E-02,0.6960E-02,0.9160E-02,
     +0.1110E-01,0.1450E-01,0.2000E-01,0.2300E-01,0.2600E-01,0.2900E-01,
     +0.2930E-01,0.3000E-01,0.2850E-01,0.1730E-01,0.1290E-01,0.1200E-01,
     +0.1250E-01,0.1340E-01,0.1400E-01,0.1750E-01,0.2400E-01,0.3500E-01,
     +0.3800E-01,0.4200E-01,0.4600E-01,0.5200E-01,0.5700E-01,0.6900E-01,
     +0.7000E-01,0.6700E-01,0.6500E-01,0.6400E-01,0.6200E-01,0.5900E-01/
      DATA(TABIM(I),I=343,456)/
     +0.5700E-01,0.5600E-01,0.5500E-01,0.5700E-01,0.5800E-01,0.5700E-01,
     +0.5500E-01,0.5500E-01,0.5400E-01,0.5200E-01,0.5200E-01,0.5200E-01,
     +0.5200E-01,0.5000E-01,0.4700E-01,0.4300E-01,0.3900E-01,0.3700E-01,
     +0.3900E-01,0.4000E-01,0.4200E-01,0.4400E-01,0.4500E-01,0.4600E-01,
     +0.4700E-01,0.5100E-01,0.6500E-01,0.7500E-01,0.8800E-01,0.1080E+00,
     +0.1340E+00,0.1680E+00,0.2040E+00,0.2480E+00,0.2800E+00,0.3410E+00,
     +0.3790E+00,0.4090E+00,0.4220E+00,0.4220E+00,0.4030E+00,0.3890E+00,
     +0.3740E+00,0.3540E+00,0.3350E+00,0.3150E+00,0.2940E+00,0.2710E+00,
     +0.2460E+00,0.1980E+00,0.1640E+00,0.1520E+00,0.1420E+00,0.1280E+00,
     +0.1250E+00,0.1230E+00,0.1160E+00,0.1070E+00,0.7900E-01,0.7200E-01,
     +0.7600E-01,0.7500E-01,0.6700E-01,0.5500E-01,0.4500E-01,0.2900E-01,
     +0.2750E-01,0.2700E-01,0.2730E-01,0.2890E-01,0.3000E-01,0.3400E-01,
     +0.5300E-01,0.7550E-01,0.1060E+00,0.1350E+00,0.1761E+00,0.2229E+00,
     +0.2746E+00,0.3280E+00,0.3906E+00,0.4642E+00,0.5247E+00,0.5731E+00,
     +0.6362E+00,0.6839E+00,0.7091E+00,0.6790E+00,0.6250E+00,0.5654E+00,
     +0.5433E+00,0.5292E+00,0.5070E+00,0.4883E+00,0.4707E+00,0.4203E+00,
     +0.3771E+00,0.3376E+00,0.3056E+00,0.2835E+00,0.3170E+00,0.3517E+00,
     +0.3902E+00,0.4509E+00,0.4671E+00,0.4779E+00,0.4890E+00,0.4899E+00,
     +0.4873E+00,0.4766E+00,0.4508E+00,0.4193E+00,0.3880E+00,0.3433E+00/
      DATA(TABIM(I),I=457,468)/
     +0.3118E+00,0.2935E+00,0.2350E+00,0.1981E+00,0.1865E+00,0.1771E+00,
     +0.1620E+00,0.1490E+00,0.1390E+00,0.1200E+00,0.9620E-01,0.8300E-01/
      DATA(TABIMT(I,1),I=1,NWLT)/
     +                                 0.8300E-01,0.6900E-01,0.5700E-01,
     +0.4560E-01,0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,
     +0.1760E-01,0.1665E-01,0.1620E-01,0.1550E-01,0.1470E-01,0.1390E-01,
     +0.1320E-01,0.1250E-01,0.1180E-01,0.1060E-01,0.9540E-02,0.8560E-02,
     +0.6210E-02,0.4490E-02,0.3240E-02,0.2340E-02,0.1880E-02,0.1740E-02,
     +0.1500E-02,0.1320E-02,0.1160E-02,0.8800E-03,0.6950E-03,0.4640E-03,
     +0.3400E-03,0.3110E-03,0.2940E-03,0.2790E-03,0.2700E-03,0.2640E-03,
     +0.2580E-03,0.2520E-03,0.2490E-03,0.2540E-03,0.2640E-03,0.2740E-03,
     +0.2890E-03,0.3050E-03,0.3150E-03,0.3460E-03,0.3820E-03,0.4620E-03,
     +0.5000E-03,0.5500E-03,0.5950E-03,0.6470E-03,0.6920E-03,0.7420E-03,
     +0.8200E-03,0.9700E-03,0.1950E-02,0.5780E-02,0.9700E-02/
      DATA(TABIMT(I,2),I=1,NWLT)/
     +                      0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,
     +0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,0.1760E-01,
     +0.1665E-01,0.1600E-01,0.1500E-01,0.1400E-01,0.1310E-01,0.1230E-01,
     +0.1150E-01,0.1080E-01,0.9460E-02,0.8290E-02,0.7270E-02,0.4910E-02,
     +0.3300E-02,0.2220E-02,0.1490E-02,0.1140E-02,0.1060E-02,0.9480E-03,
     +0.8500E-03,0.7660E-03,0.6300E-03,0.5200E-03,0.3840E-03,0.2960E-03,
     +0.2700E-03,0.2520E-03,0.2440E-03,0.2360E-03,0.2300E-03,0.2280E-03,
     +0.2250E-03,0.2200E-03,0.2160E-03,0.2170E-03,0.2200E-03,0.2250E-03,
     +0.2320E-03,0.2390E-03,0.2600E-03,0.2860E-03,0.3560E-03,0.3830E-03,
     +0.4150E-03,0.4450E-03,0.4760E-03,0.5080E-03,0.5400E-03,0.5860E-03,
     +0.6780E-03,0.1280E-02,0.3550E-02,0.5600E-02/
      DATA(TABIMT(I,3),I=1,NWLT)/
     +           0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,0.3790E-01,
     +0.3140E-01,0.2620E-01,0.2190E-01,0.1880E-01,0.1660E-01,0.1540E-01,
     +0.1470E-01,0.1350E-01,0.1250E-01,0.1150E-01,0.1060E-01,0.9770E-02,
     +0.9010E-02,0.7660E-02,0.6520E-02,0.5540E-02,0.3420E-02,0.2100E-02,
     +0.1290E-02,0.7930E-03,0.5700E-03,0.5350E-03,0.4820E-03,0.4380E-03,
     +0.4080E-03,0.3500E-03,0.3200E-03,0.2550E-03,0.2120E-03,0.2000E-03,
     +0.1860E-03,0.1750E-03,0.1660E-03,0.1560E-03,0.1490E-03,0.1440E-03,
     +0.1350E-03,0.1210E-03,0.1160E-03,0.1160E-03,0.1170E-03,0.1200E-03,
     +0.1230E-03,0.1320E-03,0.1440E-03,0.1680E-03,0.1800E-03,0.1900E-03,
     +0.2090E-03,0.2160E-03,0.2290E-03,0.2400E-03,0.2600E-03,0.2920E-03,
     +0.6100E-03,0.1020E-02,0.1810E-02/
      DATA(TABIMT(I,4),I=1,NWLT)/
     +0.8300E-01,0.6900E-01,0.5700E-01,0.4450E-01,0.3550E-01,0.2910E-01,
     +0.2440E-01,0.1970E-01,0.1670E-01,0.1400E-01,0.1235E-01,0.1080E-01,
     +0.8900E-02,0.7340E-02,0.6400E-02,0.5600E-02,0.5000E-02,0.4520E-02,
     +0.3680E-02,0.2990E-02,0.2490E-02,0.1550E-02,0.9610E-03,0.5950E-03,
     +0.3690E-03,0.2670E-03,0.2510E-03,0.2290E-03,0.2110E-03,0.1960E-03,
     +0.1730E-03,0.1550E-03,0.1310E-03,0.1130E-03,0.1060E-03,0.9900E-04,
     +0.9300E-04,0.8730E-04,0.8300E-04,0.7870E-04,0.7500E-04,0.6830E-04,
     +0.5600E-04,0.4960E-04,0.4550E-04,0.4210E-04,0.3910E-04,0.3760E-04,
     +0.3400E-04,0.3100E-04,0.2640E-04,0.2510E-04,0.2430E-04,0.2390E-04,
     +0.2370E-04,0.2380E-04,0.2400E-04,0.2460E-04,0.2660E-04,0.4450E-04,
     +0.8700E-04,0.1320E-03/

      DATA PI /3.14159265/
C
C     ZERO PARAMETERS
C
      RN = 0.0
      CN = 0.0
      ABSIND = 0.0
      ABSCOF = 0.0
C
C     CONVERT WAVELENGTH TO MICRONS
C
      ALAM = XLAM
      IF (IUNIT .EQ. 1) ALAM = 1000 * ALAM
      IF (IUNIT .EQ. 2) ALAM = 10000 * ALAM
      IF (IUNIT .EQ. 3) ALAM = 10000 * (1.0/ALAM)
      IF ((ALAM .LT. WLMIN) .OR. (ALAM .GT. WLMAX)) RETURN
      IF (ALAM .GT. CUTICE)GO TO 10
C
C     REGION FROM 0.045 MICRONS TO 167.0 MICRONS - NO TEMPERATURE DEPEND
C
      DO 1 I = 2, NWL
         IF (ALAM .LT. WL(I)) GO TO 2
 1    CONTINUE
 2    X1 = ALOG(WL(I-1))
      X2 = ALOG(WL(I))
      Y1 = TABRE(I-1)
      Y2 = TABRE(I)
      X = ALOG(ALAM)
      Y = ((X - X1)*(Y2 - Y1)/(X2 - X1)) + Y1
      RN = Y
      Y1 = ALOG(ABS(TABIM(I-1)))
      Y2 = ALOG(ABS(TABIM(I)))
      Y = ((X - X1)*(Y2 - Y1)/(X2 - X1)) + Y1
      CN = EXP(Y)
      GO TO 20
C
C     REGION FROM 167.0 MICRONS TO 8.6 METERS - TEMPERATURE DEPENDENCE
C
 10   TK = T
      IF (TK .GT. TEMREF(1)) TK = TEMREF(1)
      IF (TK .LT. TEMREF(4)) TK = TEMREF(4)
      DO 11 I = 2, 4
         IF(TK.GE.TEMREF(I)) GO TO 12
 11   CONTINUE
 12   LT1 = I
      LT2 = I - 1
      DO 13 I = 2, NWLT
         IF (ALAM .LE. WLT(I)) GO TO 14
 13   CONTINUE
 14   X1 = ALOG(WLT(I-1))
      X2 = ALOG(WLT(I))
      Y1 = TABRET(I-1,LT1)
      Y2 = TABRET(I,LT1)
      X = ALOG(ALAM)
      YLO = ((X - X1)*(Y2 - Y1)/(X2 - X1)) + Y1
      Y1 = TABRET(I-1,LT2)
      Y2 = TABRET(I,LT2)
      YHI = ((X - X1)*(Y2 - Y1)/(X2 - X1)) + Y1
      T1 = TEMREF(LT1)
      T2 = TEMREF(LT2)
      Y = ((TK - T1)*(YHI - YLO)/(T2 - T1)) + YLO
      RN = Y
      Y1 = ALOG(ABS(TABIMT(I-1,LT1)))
      Y2 = ALOG(ABS(TABIMT(I,LT1)))
      YLO = ((X - X1)*(Y2 - Y1)/(X2 - X1)) + Y1
      Y1 = ALOG(ABS(TABIMT(I-1,LT2)))
      Y2 = ALOG(ABS(TABIMT(I,LT2)))
      YHI = ((X - X1)*(Y2 - Y1)/(X2 - X1)) + Y1
      Y = ((TK - T1)*(YHI - YLO)/(T2 - T1)) + YLO
      CN = EXP(Y)
C
C     ABSORPTIVE QUANITIES
C
 20   ABSIND = CN / RN
      ABSCOF = 4.0*PI * CN / ALAM

      RETURN
      END


      FUNCTION eps_swd_l91sd(f, T)

C     This function returns the complex dielectric constant of suspended
C     (pure) water droplets as a function of frequency f (GHz) and
C     temperature T (deg. C), using the single Debye model of Liebe,
C     Hufford, and Manabe (Int. J. IR & mm Waves V.12 (7) July 1991).
C
C     Note: You may need to remove the underscores from function names if
C     using the g77 compiler.

      REAL f                    ! [GHz]  Frequency (valid from 0 to 1000 GHz)
      REAL T                    ! [�C]   Temperature
      COMPLEX eps_swd_l91sd     !        Dielectric constant
      REAL theta                !        Inverse temperature variable
      REAL eps0                 !        Static dielectric constant
      REAL epsinf               !        High frequency dielectric constant
      REAL gamma                ! [GHz]  Relaxation frequency

      COMPLEX i
      PARAMETER (i = (0.0, 1.0))

      theta = 1.0 - 300.0/(273.15 + T)
      eps0 = 77.66 - 103.3*theta
      epsinf = 0.066 * eps0
      gamma = 20.27 + 146.5 * theta + 314.0 * theta**2

      eps_swd_l91sd = (eps0 - epsinf)/(1.0 - i * f/gamma) + epsinf

      RETURN
      END


C-----------------------------------------------------------------------
      FUNCTION eps_swd_l91dd(f, T)

C     This function returns the complex dielectric constant of suspended
C     (pure) water droplets as a function of frequency f (GHz) and
C     temperature T (deg. C), using the double Debye model of Liebe, Hufford,
C     and Manabe (Int. J. IR & mm Waves V.12 (7) July 1991).
C
C     Note: You may need to remove the underscores from function names if
C     using the g77 compiler.

      REAL f                    ! [GHz]  Frequency (valid from 0 to 1000 GHz)
      REAL T                    ! [�C]   Temperature
      COMPLEX eps_swd_l91dd     !        Dielectric constant
      REAL theta                !        Inverse temperature variable
      REAL eps0                 !        Static dielectric constant
      REAL eps1                 !        First high-frequency constant
      REAL eps2                 !        Second high-frequency constant
      REAL gamma1               ! [GHz]  Primary relaxation frequency
      REAL gamma2               ! [GHz]  Secondary relaxation frequency

      theta = 1.0 - 300.0/(273.15 + T)
      eps0 = 77.66 - 103.3*theta

      eps1 = 0.0671 * eps0
      gamma1 = 20.20 + 146.4*theta + 316.0*theta*theta
      eps2 = 3.52 + 7.52*theta
      gamma2 = 39.8*gamma1

      eps_swd_l91dd = eps0 - f * ((eps0 - eps1)/CMPLX(f,gamma1)
     &              + (eps1 - eps2)/CMPLX(f,gamma2))

      RETURN
      END
