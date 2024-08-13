c-------Add76_xO3.f    code for integrating O3 & q1D cross sections
!      covers Solar-J range 1:27 bins (i.e., RRTMG-SW)
!      starts with pratmo full 76+SR added bins and then generates the Fast-J bins(18)
!      also puts xO3 >778 nm (bin 19) into Fast-J bin 18

      implicit none
      integer, parameter :: NC_ = 199
      integer, parameter :: NB_ = 99
      integer, parameter :: NX_ = 6785
      integer, parameter :: NY_ = 40000

      real*8  SRB(15,8)
      real*8  WX(NX_),FX(NX_)
      real*8  WY(NY_),FY(NY_)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
! ABIN, BBIN, CBIN, DBIN are for wt'd sums of x-sections or Refractive Index or ...
      real*8, dimension(NC_) ::
     &        FFBIN,WWBIN,RRBIN,YYBIN,AABIN,BBBIN,CCBIN,DDBIN
      real*8, dimension(NB_) ::
     &        FBIN,WBIN,RBIN,YBIN,ABIN,BBIN,CBIN,DBIN
      integer  I, J,J1,J2,  NC1,NC2,NC3,NC4, NB3,NB4
      real*8  WW, WNM, RAYLAY,YPAR, W11,W22
      character*80 TITLE, TITLTBL
      character*6  TITLNEW,TSPEC(2)
      real*8   XNEW,QNEW,XT, RR18X,AA18X,CC18X
      integer  IT1,IT2,IT3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! Reads in the full set of wavelength bins needs to map the S-R bands
!     and not-adjacent bins into the Fast-JX 18 bins plus the Solar-J bins

      open (1, file='SolarJ_bins.dat', status='OLD')
        SRB(:,:) = 0.d0
        read(1,'(a)') TITLE
        read(1,'(a)') TITLE
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
        read(3,*)
       do I = 1,NY_
        read(3,'(f10.4,e10.3)') WNM,FY(I)
        WY(I) = 1.d-3*WNM
       enddo
      close (3)
!!!!!!!!!!!!!!!!!!!!!!! finished setup !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!! iniitalize the subroutine from external data tables if need be
       WW = 0.0
       XT = 200.
      call X_O3 (0, WW,XT,XNEW, TSPEC(1),TITLTBL)
      call Q_O3 (0, WW,XT,XNEW, TSPEC(2),TITLTBL)
       IT1 = 218
       IT2 = 258
       IT3 = 298

!!!! Begin integration of XO3 and qO1D over wavelengths
!!!!!!!!!  NB, the  0.05 nm bins in v76 are NOT split but fall into ONE of the 76 pratmo bins
!!!!!!!!!     Here the microbin 635.10 nm falls into the 635.10-778.00 nm pratmo bin (#76)
!!!!!!!!!     In some earlier versions it was triggered to fall into pratmo #75 = 560.10-635.10 nm
!!!!!!!!!     Differences are at most 5-7% and the allocation of solar flux is consistent
!!!!!!!!!  Thus some earlier versions of X-sections are slightly different.  Not all have been replaced
!!!!  FFBIN = solar flux sum,
!!!!  RRBIN,YYBIN,AABIN,BBBIN,CCBIN,DDBIN= sum for xO3, qO3 at 3 temperatures
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
!   allow wavelength integration out to IR bands to get last Chappuis bit inot bin 18
      do J = 1,NC4
       do while (WY(I) .lt. WCBIN(J+1) .and. I .lt. NY_)
        WW = WY(I)*1.d3
         XT = IT1
         call X_O3 (1, WW,XT,XNEW, TSPEC(1),TITLTBL)
         call Q_O3 (1, WW,XT,QNEW, TSPEC(2),TITLTBL)
          FFBIN(J) = FFBIN(J) + FY(I)
          RRBIN(J) = RRBIN(J) + FY(I)*XNEW
          YYBIN(J) = YYBIN(J) + FY(I)*XNEW*QNEW
         XT = IT2
         call X_O3 (1, WW,XT,XNEW, TSPEC(1),TITLTBL)
         call Q_O3 (1, WW,XT,QNEW, TSPEC(2),TITLTBL)
          AABIN(J) = AABIN(J) + FY(I)*XNEW
          BBBIN(J) = BBBIN(J) + FY(I)*XNEW*QNEW
         XT = IT3
         call X_O3 (1, WW,XT,XNEW, TSPEC(1),TITLTBL)
         call Q_O3 (1, WW,XT,QNEW, TSPEC(2),TITLTBL)
          CCBIN(J) = CCBIN(J) + FY(I)*XNEW
          DDBIN(J) = DDBIN(J) + FY(I)*XNEW*QNEW
        I = I+1
       enddo
       if (I .ge. NY_) goto 2
      enddo
    2 continue

      do J = 1,NC4
        YYBIN(J) = YYBIN(J)/max(RRBIN(J),1.d-40)
        BBBIN(J) = BBBIN(J)/max(AABIN(J),1.d-40)
        DDBIN(J) = DDBIN(J)/max(CCBIN(J),1.d-40)
      enddo
      do J = 1,NC4
        RRBIN(J) = RRBIN(J)/max(FFBIN(J),1.d-40)
        AABIN(J) = AABIN(J)/max(FFBIN(J),1.d-40)
        CCBIN(J) = CCBIN(J)/max(FFBIN(J),1.d-40)
      enddo


!!!! NC4 bins for pratmo, photon weighted X's
!      write(6,'(a)') TITLNEW
      write(6,'(2a)') 'pratmbin#    solflx XO3: 218K      258K      ',
     &    '298K q1D: 218K      258K      298K '
      do J = 1,NC4
        write(6,'(a5,i4,1p,7e10.3)') 'Phot ',J,FFBIN(J),
     &   RRBIN(J),AABIN(J),CCBIN(J),YYBIN(J),BBBIN(J),DDBIN(J)
      enddo

!!!! Second integration from NC4 bins to NB4 bins (Fast-JX + Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        RBIN(:) = 0.d0
        YBIN(:) = 0.d0
        ABIN(:) = 0.d0
        BBIN(:) = 0.d0
        CBIN(:) = 0.d0
        DBIN(:) = 0.d0
       do I=16,NC4
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)
          YBIN(J) = YBIN(J) + FFBIN(I)*RRBIN(I)*YYBIN(I)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)
          BBIN(J) = BBIN(J) + FFBIN(I)*AABIN(I)*BBBIN(I)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)
          DBIN(J) = DBIN(J) + FFBIN(I)*CCBIN(I)*DDBIN(I)
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)                  *SRB(I,J)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)         *SRB(I,J)
          YBIN(J) = YBIN(J) + FFBIN(I)*RRBIN(I)*YYBIN(I)*SRB(I,J)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)         *SRB(I,J)
          BBIN(J) = BBIN(J) + FFBIN(I)*AABIN(I)*BBBIN(I)*SRB(I,J)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)         *SRB(I,J)
          DBIN(J) = DBIN(J) + FFBIN(I)*CCBIN(I)*DDBIN(I)*SRB(I,J)
        enddo
       enddo

       do J=1,NB4
         YBIN(J) = YBIN(J)/max(RBIN(J),1.d-40)
         BBIN(J) = BBIN(J)/max(ABIN(J),1.d-40)
         DBIN(J) = DBIN(J)/max(CBIN(J),1.d-40)
       enddo
       do J=1,NB4
         RBIN(J) = RBIN(J)/max(FBIN(J),1.d-40)
         ABIN(J) = ABIN(J)/max(FBIN(J),1.d-40)
         CBIN(J) = CBIN(J)/max(FBIN(J),1.d-40)
       enddo

!!!! NB4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
!      write(6,'(a)') TITLNEW
      write(6,'(2a)') 'SJ   bin#    solflx XO3: 218K      258K      ',
     &    '298K q1D: 218K      258K      298K '
      do J = 1,NB4
        write(6,'(a5,i4,1p,7e10.3)') 'Phot ',J,
     &   FBIN(J),RBIN(J),ABIN(J),CBIN(J),YBIN(J),
     &   BBIN(J),DBIN(J)
      enddo

! fast-J v76 data for 'FJX_spec.dat'
! if XO3 slips over into SJ bin#19 consolidate BINS #18 & #19
!     and increase the effective XO3 in SJ bin #18 (applies to Chappuis)

      write(6,'(a)') 'NB: Xs in bin19 are consolidated in bin18 below'
       RR18X = (RBIN(18)*FBIN(18)+RBIN(19)*FBIN(19))/FBIN(18)
       AA18X = (aBIN(18)*FBIN(18)+ABIN(19)*FBIN(19))/FBIN(18)
       CC18X = (cBIN(18)*FBIN(18)+CBIN(19)*FBIN(19))/FBIN(18)
       RBIN(18) = RR18X
       ABIN(18) = AA18X
       CBIN(18) = CC18X


        write(6,'(a)') TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'a',RBIN(1:6),TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'b',RBIN(7:12),TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'c',RBIN(13:18),TSPEC(1)

        write(6,'(a)') TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'a',ABIN(1:6),TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'b',ABIN(7:12),TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'c',ABIN(13:18),TSPEC(1)

        write(6,'(a)') TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT3,'a',CBIN(1:6),TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT3,'b',CBIN(7:12),TSPEC(1)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT3,'c',CBIN(13:18),TSPEC(1)



        write(6,'(a)') TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'a',YBIN(1:6),TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'b',YBIN(7:12),TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'c',YBIN(13:18),TSPEC(2)

        write(6,'(a)') TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'a',BBIN(1:6),TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'b',BBIN(7:12),TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'c',BBIN(13:18),TSPEC(2)

        write(6,'(a)') TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT3,'a',DBIN(1:6),TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT3,'b',DBIN(7:12),TSPEC(2)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT3,'c',DBIN(13:18),TSPEC(2)

! Can redo this using Wats/m2 here, see Add76_solar.f, if better
!    heating rates in Chappuis band is wanted bins 17&18

      stop
      end




c-----------------------------------------------------------------------
      subroutine X_O3 (INIT, WW,XT,XNEW, TITLNEW,TITLTBL)
c-----------------------------------------------------------------------
c   WW = wavelength (nm) to calc Xsection for
c   XT = temerature (K) for interpolation
c   XNEW = cross section (cm2) as a function of WW and XT (and XP, XM)
c   INIT = initialization:
c     if INIT.eq.0 .then reads in any tables and sets Xsect name to TITLNEW
c-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: INIT
      real*8, intent(in) :: WW,XT
      real*8, intent(out) :: XNEW
      character*6, intent(out) :: TITLNEW
      character*80,intent(out) :: TITLTBL

      character*80 FTBL,FORMW
      real*8 W1(999),W2(999),XT1(999),XT2(999)
      real*8 XXT,T1,T2,TFACT,WW1
      integer NW,I,IW

      save NW,W1,W2,XT1,XT2,T1,T2

!---initialize O3 tables, be sure to "save" them for recall.
      if (INIT .eq. 0) then

        TITLNEW = 'xO3   '
        FTBL = 'XO3_JPL11Y.dat'
        open (3, file=FTBL, status='OLD')
          read(3,'(a80)') TITLTBL
            write(6,'(2a/a)') ' openfile=',FTBL, TITLTBL
          read(3,'(i4,1x,a)') NW,FORMW
         do I=1,NW
          read(3,FORMW) W1(I),W2(I),XT2(I),XT1(I)
!          write(6,'(i5,2f9.3,2f12.6)') I,W1(I),W2(I),XT2(I),XT1(I)
         enddo
         close(3)
          T2 = 295.d0
          T1 = 218.d0

      else

!---interpolate X-section vs. T, but use single mean value for entire wavelength bins
        XXT = min(T2, max(T1, XT))
        TFACT = (XXT - T1)/(T2 - T1)
!  note that W2(I) = W1(I+1) -- bins do not miss any wavelengths
!  but can cause problems when a point is on the boudary edge, so trick it:
       WW1 = WW - 0.01d0
        IW = 1
       do I=1,NW-1
        if (WW1 .gt. W2(I)) IW = I+1
       enddo
         XNEW = (XT1(IW) + TFACT*(XT2(IW)-XT1(IW)))
       WW1 = WW + 0.01d0
        IW = 1
       do I=1,NW-1
        if (WW1 .gt. W2(I)) IW = I+1
       enddo
        XNEW = 0.5d-20*(XT1(IW) + TFACT*(XT2(IW)-XT1(IW)) +XNEW)

      endif
      return
      end



c-----------------------------------------------------------------------
      subroutine Q_O3 (INIT, WW,XT,QNEW, TITLNEW,TITLTBL)
c-----------------------------------------------------------------------
      implicit none
      integer, intent(inout) :: INIT
      real*8, intent(in) :: WW,XT
      real*8, intent(out) :: QNEW
      character*6, intent(out) :: TITLNEW
      character*80,intent(out) :: TITLTBL

      real*8  T,QO1D,Q1,Q2,Q1Q2,EX1,EX2,EX3

!---quantum yield for O3 + hv => O(1D),
!---parametric fit for range w = 306-328 nm, T = 200-320K
!---JPL_2010 standard tables:
!     Table 4A-5. Parameters for the Calculation of O(1D) Quantum Yields.
      real*8, parameter::  X1 = 304.225d0
      real*8, parameter::  X2 = 314.957d0
      real*8, parameter::  X3 = 310.737d0
      real*8, parameter::  W1 = 5.576d0
      real*8, parameter::  W2 = 6.601d0
      real*8, parameter::  W3 = 2.187d0
      real*8, parameter::  A0 = 0.90d0
      real*8, parameter::  A1 = 0.8036d0
      real*8, parameter::  A2 = 8.9061d0
      real*8, parameter::  A3 = 0.1192d0
      real*8, parameter::  A4 = 0.0765d0
      real*8, parameter::  A4x = 0.08d0
      real*8, parameter::  V1 = 0.0d0
      real*8, parameter::  V2 = 825.518d0
      real*8, parameter::  RG = 0.695d0


      if (INIT .eq. 0) then
        TITLNEW = 'qO1D  '
        TITLTBL='Table 4A-5. Parameters for the Calculation of O(1D)'
        write(6,'(2a)') ' species:',TITLNEW
        write(6,'(a80)') TITLTBL

      else

c--set QO1D = 0.0 for W.gt.340.,  = 0.48 at 193 nm, = 0.90 at 225 nm
       if (WW .lt. 306.d0) then
          QO1D = A0
       elseif (wW .gt. 328.d0) then
          QO1D = A4x
       else
         T = min (320.d0, max(200.d0, XT))
         Q1 = exp(-V1/(RG*T))
         Q2 = exp(-V2/(RG*T))
         Q1Q2 = Q1/(Q1+Q2)
         EX1 = exp(-((X1-WW)/W1)**4)
         EX2 = exp(-((X2-WW)/W2)**2) * (T/300.d0)**2
         EX3 = exp(-((X3-WW)/W3)**2) * (T/300.d0)**1.5
         QO1D = Q1Q2*A1*EX1 + (1.d0-Q1Q2)*A2*EX2 + A3*EX3 + A4
       endif
       if (WW .lt. 220.d0) then
         QO1D = max(0.48d0, 0.48d0 + 0.42d0*(WW-190.d0)/30.d0)
       endif
       if (WW .gt. 340.d0) then
         QO1D = 0.d0
       endif
       QNEW = QO1D

      endif
      return
      end
