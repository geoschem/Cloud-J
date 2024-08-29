c-------Add76_xNO3.f   with user supplied subroutine that supplies X-section x q-yields
!   X-sections - generates only pratmo 76 bins abnd then Fast-J 18 bins

      implicit none
      integer, parameter :: NC_ = 199
      integer, parameter :: NB_ = 99
      integer, parameter :: NY_ = 40000

      real*8   SRB(15,8)
      real*8  WY(NY_),FY(NY_)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
      real*8, dimension(NC_) :: FFBIN,AABIN,BBBIN,CCBIN
      real*8, dimension(NB_) :: FBIN,ABIN,BBIN,CBIN
      integer  I, J,J1,J2,K,  NC1,NC2,NC3,NC4, NB3,NB4, NB76
      real*8  WW, WNM, RAYLAY,YPAR, W11,W22
      character*80 TITLE
      Character*80 TITLTBL, TITL76, TITL76a, TITL76b,TITL76c
      character*6  TITLNEW,TSPEC(2)
      real*8   XXWTa,XXWTb,XT,  TOTA, TOTB, QYLDA, JCLR
      integer :: TTT(2)
      character*1 :: ISX,ISP,ISXP


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
!>>>>  for all std x-sections, truncate integration at NC3, do not go into SJ 1
          if (NC4 .gt. NC_) stop
        read(1,'(5x,i5,f8.3)') (IJX(I),WCBIN(I), I=1,NC4+1)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,'(2x,15f5.1)') ((SRB(I,J),I=1,15),J=1,8)
!  I tracks the 76 or 85 or NC4 (w/Solar-J) bins from the high-res pratmo wavelengths
!  J or IJX(I) tracks the 1:18 (Cloud-J) ro 1:27 (Solar-J) bins
! convert all to microns
        do I = 1,NC4+1
          WCBIN(I) = WCBIN(I)*1.d-3
        enddo
        NB4 = IJX(NC4)
        NB3 = IJX(NC3)
          if (NB4 .gt. NB_) stop
      close (1)

! skip Watts for std X-section: open (2, file='SolarF_watts.dat', status='OLD')

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


!!!!!!!!!!!!! iniitalize the subroutine from external data tables if need be
       call X_NO3 (0, WW,XT, XXWTa,XXWTb,TITLNEW)
! flag for X sections, 'x' = strat only, 'p' = pressure interp (cannot be strat)
         ISXP = ' '
         TTT(1) = 190
         TTT(2) = 298

!!!!!!!!!!  major temperature-or-density loop for X-sections
      do K = 1,2
         XT = TTT(K)

!!!!!!!!!   flux-weighted means over the 76-pratmo bins
!!!!!!!!!  NB, the  0.05 nm bins in v76 are NOT split but fall into ONE of the 76 pratmo bins
!!!!!!!!!     Here the microbin 635.10 nm falls into the 635.10-778.00 nm pratmo bin (#76)
!!!!!!!!!     In some earlier versions it was triggered to fall into pratmo #75 = 560.10-635.10 nm
!!!!!!!!!     Differences are at most 5-7% and the allocation of solar flux is consistent
!!!!!!!!!  Thus some earlier versions of X-sections are slightly different.  Not all have been replaced
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
        BBBIN(:) = 0.d0
         J = 1
         I = 1
       do while (WY(I) .lt. WCBIN(J))
         I = I+1
       enddo
       do J = 1,NC3
        do while (WY(I) .lt. WCBIN(J+1) .and. I .lt. NY_)
           WW = WY(I)*1.d3

           call X_NO3 (1, WW,XT, XXWTa,XXWTb,TITLNEW)

           FFBIN(J) = FFBIN(J) + FY(I)
           AABIN(J) = AABIN(J) + FY(I)*XXWTa
           BBBIN(J) = BBBIN(J) + FY(I)*XXWTb
         I = I+1
        enddo
        if (I .ge. NY_) goto 2
       enddo
    2  continue
       do J=1,NC3
         AABIN(J) = AABIN(J)/FFBIN(J)
         BBBIN(J) = BBBIN(J)/FFBIN(J)
       enddo
       NB76 = NC3

!!!!!! Second integration from NC3 bins to NB3 bins (Fast-JX + not Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        ABIN(:) = 0.d0
        BBIN(:) = 0.d0
       do I=16,NC3
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)                  *SRB(I,J)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)         *SRB(I,J)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)         *SRB(I,J)
        enddo
       enddo
         TOTA = 0.
         TOTB = 0.
       do J=1,NB3
         TOTA = TOTA + ABIN(J)
         TOTB = TOTB + BBIN(J)
         ABIN(J) = ABIN(J)/FBIN(J)
         BBIN(J) = BBIN(J)/FBIN(J)
       enddo
         JCLR = TOTA+TOTB
         QYLDA = TOTA/(TOTA+TOTB)

! fast-J v76 data for 'FJX_spec.dat'

        TITL76a  = 'NO3a  = NO2+O   JPL2010'
        TITLNEW = 'NO3a'
        write(6,'(a)') TITL76a
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',ABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',ABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',ABIN(13:18),TITLNEW
        write(6,*)

        TITL76b  = 'NO3b  = NO+O2'
        TITLNEW = 'NO3b'
        write(6,'(a)') TITL76b
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',BBIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',BBIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',BBIN(13:18),TITLNEW
        write(6,*)

!  In case we just want the single X-section and fix the split. (likely)
        CBIN(1:18) = ABIN(1:18)+BBIN(1:18)
        TITL76c = 'NO3   q(NO2+O)=0.882, q(NO+O2)=0.112'
        TITLNEW = 'NO3'
        write(6,'(a)') TITL76c
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',CBIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',CBIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',CBIN(13:18),TITLNEW
        write(6,*)
        write(6,'(a30,f15.4)') 'clear sky total J (/sec): ',JCLR
        write(6,'(a30,f15.6)') 'q-yield path a: ',QYLDA
        write(6,*)
! 'pratmo' 76-bin tables
        write(6,'(a40,f5.0/(1p,8e10.3))') TITL76a,XT,(AABIN(I),I=1,NB76)
        write(6,*)
        write(6,'(a40,f5.0/(1p,8e10.3))') TITL76b,XT,(BBBIN(I),I=1,NB76)
        write(6,*)
        CCBIN(1:76) = AABIN(1:76)+BBBIN(1:76)
        write(6,'(a40,f5.0/(1p,8e10.3))') TITL76C,XT,(ccBIN(I),I=1,NB76)
        write(6,*)

!!!  end of loop over T
      enddo

      stop
      end




      subroutine X_NO3 (INIT, WW,TT, XXWTa,XXWTb,TITLNEW)
      implicit none
      real*8, intent(in) :: WW, TT
      integer, intent(in) :: INIT
      real*8, intent(out) :: XXWTa,XXWTb
      character*6, intent(out) :: TITLNEW

      real*8, parameter :: T1=-1096.4d0, T2=-529.5d0
     &                    ,T1_298=T1/298.d0,T2_298=T2/298.d0
      character*80 FTBL,TABLE,FORMW,FORMB
      real*8, dimension(999) :: W,XW,WQY,QYb298,QYb190  ! NO3b (NO+O2)
     &                         ,QYa298,QYa190           ! NO3a (NO2+O)
      real*8  TTL,XBFACT,XTFACT,XXW,FW,FWQ, XWI,XWIP1,TFACT
      real*8  QYa,QYb,QYa230,QYb230
      integer NW,NB,NQ, N, I,IW,IWQ

      save  W,XW,NW, WQY,QYb298,QYb190,QYa298,QYa190,NQ

      XTFACT(TTL) = (1.d0-dexp(T1/TTL)-2.d0*dexp(T2/TTL)) /
     &              (1.d0-dexp(T1_298)-2.d0*dexp(T2_298))

      TITLNEW = 'xNO3ab'
      FTBL = 'Add76_XNO3_JPL10.dat'

      if (INIT .eq. 0) then

        open (3, file=FTBL, status='OLD')
          read(3,'(a)') TABLE
             write(6,'(2a/a)') ' openfile=',FTBL, TABLE
          read(3,'(i4,1x,a)') NW
        do N = 1,NW
          read(3,*) W(N),XW(N)
        enddo
          read(3,*)
          read(3,*)
          read(3,*)
          read(3,'(a)') TABLE
             write(6,'(a)') TABLE
          read(3,*)
          read(3,*)
          read(3,'(i4,1x,a)') NQ,FORMB
        do N = 1,NQ
          read(3,FORMB)
     &  WQY(N),QYb298(N),QYb230,QYb190(N), QYa298(N),QYa230,QYa190(N)
!                    NO3b = (NO+O2)            NO3a = (NO2+O)
        enddo
        close(3)

      else

c---interpolate X-section vs. Wavelength, WW should be in nm!
        IW = 1
       do I = 2,NW-1
        if (WW .ge. W(I)) IW = I
       enddo

        FW = (WW - W(IW))/(W(IW+1) - W(IW))
        FW = min(1.d0, max(0.d0, FW))
        XXW = XW(IW) + FW*(XW(IW+1)-XW(IW))

c---interpolate QY. vs. Wavelength
        IWQ = 1
       do I = 2,NQ-1
        if (WW .ge. WQY(I)) IWQ = I
       enddo
        FWQ = (WW - WQY(IWQ))/(WQY(IWQ+1) - WQY(IWQ))
        FWQ = min(1.d0, max(0.d0, FWQ))
       if (TT .gt. 244.d0) then
        QYa = 1.d-3 * (QYa298(IWQ) + FWQ*(QYa298(IWQ+1)-QYa298(IWQ)) )
        QYb = 1.d-3 * (QYb298(IWQ) + FWQ*(QYb298(IWQ+1)-QYb298(IWQ)) )
       else
        QYa = 1.d-3 * (QYa190(IWQ) + FWQ*(QYa190(IWQ+1)-QYa190(IWQ)) )
        QYb = 1.d-3 * (QYb190(IWQ) + FWQ*(QYb190(IWQ+1)-QYb190(IWQ)) )
       endif
c---T dependence: set for NO3
       XXWTa = XTFACT(TT) * XXW * 1.d-20 * QYa
       XXWTb = XTFACT(TT) * XXW * 1.d-20 * QYb

      endif
      return
      end
