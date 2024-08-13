!----Add76_solar.f    v7.6 code for integrating over Solar fluxes, Rayleigh, YPar
!      covers Solar-J range 1:27 bins (i.e., RRTMG-SW)
!      uses both photons and Watts to do weighting
!      starts with pratmo full 77+SR bins and then generates the Solar-J bins (18+9)
! not used now     real*8, parameter :: H = 6.623e-34, C = 2.998e8

      implicit none
      integer, parameter :: NC_ = 199
      integer, parameter :: NB_ = 99
      integer, parameter :: NX_ = 6785
      integer, parameter :: NY_ = 40000

      real*8   SRB(15,8)
      real*8  WX(NX_),FX(NX_)
      real*8  WY(NY_),FY(NY_)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
      real*8, dimension(NC_) ::
     &        FFBIN,WWBIN,RRBIN,YYBIN,AABIN,BBBIN,CCBIN,DDBIN
      real*8, dimension(NB_) :: FBIN,WBIN,RBIN,YBIN,ABIN,BBIN,CBIN,DBIN
!       ABIN, BBIN, CBIN, DBIN are for added X-sections or Refractive Index averaging
      integer  I, J,J1,J2,  NC1,NC2,NC3,NC4, NB3,NB4
      real*8  WW, WNM, RAYLAY,YPAR, W11,W22
      character*80 TITLE
      character*6  TITLNEW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! Reads in the full set of wavelength bins needs to map the S-R bands
!     and not-adjacent bins into the Fast-JX 18 bins plus the Solar-J bins

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
        NB3 = IJX(NC3)
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


      open (3, file='SolarF_photons.dat', status='OLD')
        read(3,'(a)') TITLE
!            write(6,'(a)') TITLE
        read(3,*)
       do I = 1,NY_
        read(3,'(f10.4,e10.3)') WNM,FY(I)
        WY(I) = 1.d-3*WNM
       enddo
      close (3)
!!!!!!!!!!!!!!!!!!!!!!! finished setup !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!  NB, the  0.05 nm bins in v76 are NOT split but fall into ONE of the 76(+ for SJ) pratmo bins
!!!!!!!!!     Here the microbin 635.10 nm falls into the 635.10-778.00 nm pratmo bin (#76)
!!!!!!!!!     In some earlier versions it was triggered to fall into pratmo #75 = 560.10-635.10 nm
!!!!!!!!!     Differences are at most 5-7% and the allocation of solar flux is consistent
!!!!!!!!!  Thus some earlier versions of X-sections are slightly different.  Not all have been replaced
!!!! integration at high-res for photons and then later do Watts
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
        WW = WY(I)*1.d3

        FFBIN(J) = FFBIN(J) + FY(I)
        WWBIN(J) = WWBIN(J) + FY(I)*WY(I)
        RRBIN(J) = RRBIN(J) + FY(I)*RAYLAY(WY(I))
        YYBIN(J) = YYBIN(J) + FY(I)*YPAR(WY(I))
        I = I+1

       enddo
        WWBIN(J) = WWBIN(J)/FFBIN(J)
        RRBIN(J) = RRBIN(J)/FFBIN(J)
        YYBIN(J) = YYBIN(J)/FFBIN(J)
        AABIN(J) = AABIN(J)/FFBIN(J) !not used here
        BBBIN(J) = BBBIN(J)/FFBIN(J)
        CCBIN(J) = CCBIN(J)/FFBIN(J)
        DDBIN(J) = DDBIN(J)/FFBIN(J)
       if (I .ge. NY_) goto 2
      enddo
    2 continue

!!!! NC4 'pratmo' bins: Photon weighted values for w, Rayleigh, PAR
      write(6,'(a,3i5)') 'pratm bin#   w-eff  phot      X-Ray     Y-PAR'
      do J = 1,NC4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Phot ',J,
     &   WWBIN(J),FFBIN(J),RRBIN(J),YYBIN(J)
      enddo

!!!! Second integration from NC4 bins to NB4 bins for Fast-JX + Solar-J
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
          DBIN(J) = DBIN(J) + DDBIN(I)*DDBIN(J)
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
      write(6,'(2a)') ' S-J bin#    w-eff  phot/Watt X-Ray     Y-PAR'
!    &  ,  '     Liq:  nr + ni       Ice:  nr + ni'
      do J = 1,NB4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Phot ',J,
     &   WBIN(J),FBIN(J),RBIN(J),YBIN(J)
!     & ,ABIN(J),BBIN(J),CBIN(J),DBIN(J)
      enddo


        TITLNEW = 'w-eff '
        WBIN(1:27) = WBIN(1:27)*1000.
        write(6,'(a)') TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '01:06',WBIN(1:6),TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '07:12',WBIN(7:12),TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '13:18',WBIN(13:18),TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '19:24',WBIN(19:24),TITLNEW
        write(6,'(a5,3f10.0,31x,a6)')
     &  '25:27',WBIN(25:27),TITLNEW

        TITLNEW = 'solflx'
        FBIN(19:27) = 1.0
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',FBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',FBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',FBIN(13:18),TITLNEW
        write(6,'(a5,   6f10.1,1x,a6)')
     &  '    d',FBIN(19:24),TITLNEW
        write(6,'(a5,   3f10.1,31x,a6)')
     &  '    e',FBIN(25:27),TITLNEW

        TITLNEW = 'X-Rayl'
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    a',RBIN(1:6),TITLNEW,' phot'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    b',RBIN(7:12),TITLNEW,' phot'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    c',RBIN(13:18),TITLNEW,' phot'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    d',RBIN(19:24),TITLNEW,' phot'
        write(6,'(a5,1p,3e10.3,31x,a6,a5)')
     &  '    e',RBIN(25:27),TITLNEW,' phot'

        TITLNEW = 'Y-PAR '
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',YBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',YBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',YBIN(13:18),TITLNEW
        write(6,'(a5,   6F10.1,1x,a6)')
     &  '    d',YBIN(19:24)
        write(6,'(a5,   3F10.1,31x,a6)')
     &  '    e',YBIN(25:27)



!!!! WATTS: Repeat integration over wavelengths for Watts/m2:  FX() and WX()

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
        WW = WX(I)*1.d3
        FFBIN(J) = FFBIN(J) + FX(I)
        WWBIN(J) = WWBIN(J) + FX(I)*WX(I)
        RRBIN(J) = RRBIN(J) + FX(I)*RAYLAY(WX(I))
        YYBIN(J) = YYBIN(J) + FY(I)*YPAR(Wx(I))
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

      do J = 1,NC4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Watts',J,
     &   WWBIN(J),FFBIN(J),RRBIN(J),YYBIN(J)
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
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)         *SRB(I,J)
          WBIN(J) = WBIN(J) + FFBIN(I)*WWBIN(I)*SRB(I,J)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)*SRB(I,J)
          YBIN(J) = YBIN(J) + FFBIN(I)*YYBIN(I)*SRB(I,J)
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

!!!! NB4 bins: Watt-weighted values for w, Rayleigh, X-O3
      write(6,'(a,3i5)') 'pratm bin#   w-eff  Watts     X-Ray    Y-PAR'
      do J = 1,NB4
        write(6,'(a5,i4,f10.4,1p,7e10.3)') 'Watts',J,
     &   WBIN(J),FBIN(J),RBIN(J),YBIN(J)
!     &  ,ABIN(J),BBIN(J),CBIN(J),DBIN(J)
      enddo

        TITLNEW = 'w-eff '
        write(6,'(a)') TITLNEW
        WBIN(1:27) = WBIN(1:27)*1000.
        write(6,'(a5,6f10.0,1x,a6)')
     &  '    a',WBIN(1:6),TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '    b',WBIN(7:12),TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '    c',WBIN(13:18),TITLNEW
        write(6,'(a5,6f10.0,1x,a6)')
     &  '    d',WBIN(19:24)
        write(6,'(a5,6f10.0,1x,a6)')
     &  '    e',WBIN(25:27)

        TITLNEW = 'solarW'
        write(6,'(a)') TITLNEW
        write(6,'(a5,   6f10.4,1x,a6)')
     &  '    a',FBIN(1:6),TITLNEW
        write(6,'(a5,   6f10.4,1x,a6)')
     &  '    b',FBIN(7:12),TITLNEW
        write(6,'(a5,   6f10.3,1x,a6)')
     &  '    c',FBIN(13:18),TITLNEW
        write(6,'(a5,   6f10.3,1x,a6)')
     &  '    d',FBIN(19:24)
        write(6,'(a5,   6f10.3,1x,a6)')
     &  '    e',FBIN(25:27)

        TITLNEW = 'X-Rayl'
        write(6,'(2a)') TITLNEW,' watt'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    a',RBIN(1:6),TITLNEW,' watt'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    b',RBIN(7:12),TITLNEW,' watt'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    c',RBIN(13:18),TITLNEW,' watt'
        write(6,'(a5,1p,6e10.3,1x,a6,a5)')
     &  '    d',RBIN(19:24),TITLNEW,' watt'
        write(6,'(a5,1p,3e10.3,31x,a6,a5)')
     &  '    e',RBIN(25:27),TITLNEW,' watt'

        TITLNEW = 'Y-PAR '
!       skip, Y-PAR is only photon wtd

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



      function YPAR(WAVE)
c---Photosynthetically Active Radiation: action spectrum (quantum): Y-PAR
c---traced from:
c      McCree, Keith J. (1972a). "The action spectrum, absorptance and
c        quantum yield of photosynthesis in crop plants"
c        Agricultural and Forest Meteorology 9:191-216.
c      McCree, Keith J. (1972b). Agric. & Forest Meteorology 10:443-453.
c---PAR in PAR is normally quantified as µmol photons/m2/s =? µE/m2/s
c        photosynthetic photon flux (area) density, or PPFD.
      implicit none
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
        IWW = max(1, min(17, IWW))
        FWW = WWI - float(IWW)
        FWW = max(0.d0, min(1.d0, FWW))
      YPAR = Y(IWW) + (Y(IWW+1)-Y(IWW))*FWW
      return
      end
