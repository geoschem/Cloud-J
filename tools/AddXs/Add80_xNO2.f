!-------Add80_XNO2.f  -- calculates both X and q*X for NO2 to implement NO2 as absorber
! NB. v80 results differ slightly from previous/recent values because of .le. vs .lt.
!     J=67 (345.00-355.05) now includes I=4900:5100 (345.00:355.00) vs. I=4901:5101
!     Cross sections differ differ by <0.5% but are consistent with Solar fluxes.

      implicit none
      integer, parameter :: NC_ = 199
      integer, parameter :: NB_ = 99
      integer, parameter :: NY_ = 40000
      real*8   SRB(15,8)
      real*8  WY(NY_),FY(NY_)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
      real*8, dimension(NC_) :: FFBIN,AABIN,BBBIN,CCBIN,DDBIN
      real*8, dimension(NB_) :: FBIN,ABIN,BBIN,CBIN,DBIN
      integer  I, J,J1,J2,K,  NC1,NC2,NC3,NC4, NB3,NB4, NB76
      real*8  WW, XT, WNM
      character*80 TITLE, TITL76
      character*6  TITLNEW
      integer :: TTT(4)
      character*1 :: ISX,ISP,ISXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!   setup NO2 data  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character*6 TITLTBL(4)
      character*80 TABLE,FORMW
      real*8 WNMX(99),XNO2A(98),XNO2B(98), TTX(2)
      real*8 WNMQ(30),QNO2A(30),QNO2B(30), TTQ(2)
      integer N, NX, NQ
      real*8 XNO2,XNO2e, QNO2,QNO2e, W2, XSCALE

C!!!!!!!!!!!!!!!!!! read in NO2 data tables !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       TITLNEW = 'NO2   '
       ISXP = ' '
          write(6,'(a9,a6)')  'species: ',TITLNEW
       open (3, file='XNO2_JPL15.dat', status='OLD')
          read(3,'(a80)') TABLE
            write(6,'(a80)') TABLE
          read(3,'(2f5.0,e10.3)') TTX, XSCALE
            write(6,'(2f10.1,1p,e10.2)') TTX, XSCALE
          read(3,'(i5,1x,a)') NX,FORMW
            write(6,'(i5,a80)') NX,FORMW
        do N=1,NX
         read(3,FORMW) WNMX(N),W2,XNO2A(N),XNO2B(N)
           XNO2A(N) = XNO2A(N)*XSCALE
           XNO2B(N) = XNO2B(N)*XSCALE
        enddo
          read(3,'(a80)') TABLE
            write(6,'(a80)') TABLE
          read(3,'(2f5.0,e10.3)') TTQ
            write(6,'(2f10.1,1p,e10.2)') TTQ
          read(3,'(i5,1x,a)') NQ,FORMW
            write(6,'(i5,a70)') NQ,FORMW
        do N=1,NQ
         read(3,FORMW) WNMQ(N),QNO2A(N),QNO2B(N)
        enddo
       close(3)
       TITLTBL(1) = 'NO2   '
       TITLTBL(2) = 'NO2e  '
       TITLTBL(3) = 'xNO2  '
       TITLTBL(4) = 'xNO2e '
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
!!!!!!!!!!!!!!!!!!!!!!! finished setup !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!! iniitalize cross section & quantum yield data !!!!!!!!!!!!
!  calcualte 2 tables, at 200K & 300K
        TTT(1) = 200
        TTT(2) = 300
      do K = 1,2
        XT = TTT(K)

!---now ready to do any flux-weighted means over the 76-pratmo bins
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
        BBBIN(:) = 0.d0
        CCBIN(:) = 0.d0
        DDBIN(:) = 0.d0
         J = 1
         I = 1
! primary high-resolution wavelength loop - generate input for pratmo reference J's
       do while (WY(I) .lt. WCBIN(J))
         I = I+1
       enddo
       do J = 1,NC3
        do while (WY(I) .lt. WCBIN(J+1) .and. I .lt. NY_)
           WW = WY(I)*1.d3

c---      XNO2  is linear, no extrapol beyond range 220 to 294K
c---      XNO2e does log interpolation and extrapolates beyond 220K & 294K
          call X_NO2 (WNMX,XNO2A,XNO2B,TTX,NX, WW,XT,XNO2,XNO2e)

c---      QNO2 is linear, limited to range Q248 to Q298
c---      QNO2e does log intepolation and extrapolates beyond 248K & 298K
          call Q_NO2 (WNMQ,QNO2A,QNO2B,TTQ,NQ, WW,XT,QNO2,QNO2e)

           FFBIN(J) = FFBIN(J) + FY(I)
           AABIN(J) = AABIN(J) + FY(I)*XNO2 *QNO2
           BBBIN(J) = BBBIN(J) + FY(I)*XNO2e*QNO2e
           CCBIN(J) = CCBIN(J) + FY(I)*XNO2
           DDBIN(J) = DDBIN(J) + FY(I)*XNO2e
         I = I+1
        enddo
        if (I .ge. NY_) goto 2
       enddo
    2   continue
       do J=1,NC3
         AABIN(J) = AABIN(J)/FFBIN(J)
         BBBIN(J) = BBBIN(J)/FFBIN(J)
         CCBIN(J) = CCBIN(J)/FFBIN(J)
         DDBIN(J) = DDBIN(J)/FFBIN(J)
       enddo
       NB76 = NC3

!!!!!! Second integration from NC3 bins to NB3 bins (Fast-JX + not Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        ABIN(:) = 0.d0
        BBIN(:) = 0.d0
        CBIN(:) = 0.d0
        DBIN(:) = 0.d0
       do I=16,NC3
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)
          DBIN(J) = DBIN(J) + FFBIN(I)*DDBIN(I)
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)                  *SRB(I,J)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)         *SRB(I,J)
          BBIN(J) = BBIN(J) + FFBIN(I)*BBBIN(I)         *SRB(I,J)
          CBIN(J) = CBIN(J) + FFBIN(I)*CCBIN(I)         *SRB(I,J)
          DBIN(J) = DBIN(J) + FFBIN(I)*DDBIN(I)         *SRB(I,J)
        enddo
       enddo
       do J=1,NB3
         ABIN(J) = ABIN(J)/FBIN(J)
         BBIN(J) = BBIN(J)/FBIN(J)
         CBIN(J) = CBIN(J)/FBIN(J)
         DBIN(J) = DBIN(J)/FBIN(J)
       enddo

! fast-J v80 data for 'FJX_spec.dat'

        TITLNEW = TITLTBL(1)
        write(6,'(a6,a44)') TITLNEW,'  NO2+hv=NO+O   JPL15'
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',ABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'b',ABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'c',ABIN(13:18),TITLNEW

        TITLNEW = TITLTBL(2)
        write(6,'(a6,a44)') TITLNEW,'  NO2+hv=NO+O   JPL15'
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',BBIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'b',BBIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'c',BBIN(13:18),TITLNEW

        TITLNEW = TITLTBL(3)
        write(6,'(a6,a44)') TITLNEW,'  NO2+hv=absorp JPL15'
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',CBIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'b',CBIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'c',CBIN(13:18),TITLNEW

        TITLNEW = TITLTBL(4)
        write(6,'(a6,a44)') TITLNEW,'  NO2+hv=absorp JPL15'
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',DBIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'b',DBIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &    ' ',TTT(K),'c',DBIN(13:18),TITLNEW

        TITLNEW = 'Sphot'
        write(6,'(a15)') 'SolarFlux new80'
        write(6,'(a1,a4,1p,6e10.3,1x,a6)')
     &   ISXP,'a',FBIN(1:6),TITLNEW
        write(6,'(a1,a4,1p,6e10.3,1x,a6)')
     &    ' ','b',FBIN(7:12),TITLNEW
        write(6,'(a1,a4,1p,6e10.3,1x,a6)')
     &    ' ','c',FBIN(13:18),TITLNEW

! 'pratmo' 76-bin tables
      write(6,'(a6,f5.0/(1p,8e10.3))')TITLTBL(1),XT,(AABIN(I),I=1,NB76)
      write(6,'(a6,f5.0/(1p,8e10.3))')TITLTBL(2),XT,(BBBIN(I),I=1,NB76)
      write(6,'(a6,f5.0/(1p,8e10.3))')TITLTBL(3),XT,(CCBIN(I),I=1,NB76)
      write(6,'(a6,f5.0/(1p,8e10.3))')TITLTBL(4),XT,(DDBIN(I),I=1,NB76)
        write(6,*)

!!!  end of K-loop over T
      enddo

      stop
      end


!!!!!!!!!!!!!!!!!  NO2 subroutines  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine X_NO2(WX,XNO2A,XNO2B,TTX,NX, W,XT,XNO2,XNO2X)
      implicit none
      real*8,  intent(in) :: WX(NX),XNO2A(NX),XNO2B(NX),TTX(2)
      integer, intent(in) :: NX
      real*8,  intent(in) :: W,XT
      real*8, intent(out) :: XNO2,XNO2X
      integer IW,I
      real*8 XA,XB,FT,FTX
c---for NO2 these 2 TTX are 220 & 294 K
        FTX = (XT - TTX(1)) / (TTX(2) - TTX(1))
        FT = max (0.d0, min (1.d0, FTX))
c---assume WX and WW are both in nm (or microns)
c---for X's assume zero value outside range of tables
      if (W .lt. WX(1) .or. W .gt. WX(Nx)) then
        XA = 0.d0
        XB = 0.d0
        goto 1
      endif
      do I=1,NX-1
       if (W .gt. WX(I)) IW = I
      enddo
c---for X's the table value is mean of the inteval IW:IW+1
        XA = XNO2A(IW)
        XB = XNO2B(IW)
    1 continue
      XNO2 = XA + FT*(XB-XA)
        XNO2X = XNO2
        if (XA .gt. 0.d0) then
      XNO2X = XA * (XB/XA)**FTX
        endif
    2 continue
      return
      end

      subroutine Q_NO2(WQ,QNO2A,QNO2B,TTQ,NQ, W,XT,QNO2,QNO2X)
      implicit none
      real*8,  intent(in) :: WQ(NQ),QNO2A(NQ),QNO2B(NQ),TTQ(2)
      integer, intent(in) :: NQ
      real*8,  intent(in) :: W,XT
      real*8, intent(out) :: QNO2,QNO2X
      integer IW,I
      real*8 QA,QB,FT,FTX,FW
c---for q-NO2 these 2 TTQ are 248 & 298 K
        FTX = (XT - TTQ(1)) / (TTQ(2) - TTQ(1))
        FT = max (0.d0, min (1.d0, FTX))
c---assume WX and WW are both in nm (or microns)
c---for Q's extrapolate first and last values
      if (W .le. WQ(1)) then
        QA = QNO2A(1)
        QB = QNO2B(1)
        goto 1
      endif
      if (W .gt. WQ(NQ)) then
        QA = QNO2A(NQ)
        QB = QNO2B(NQ)
        goto 1
      endif
      do I=1,NQ-1
       if (W .gt. WQ(I)) IW = I
      enddo
c---interpolate in wavelength between tabulated values
        FW = (W - WQ(IW))/(WQ(IW+1) - WQ(IW))
        QA = QNO2A(IW) + FW*(QNO2A(IW+1)-QNO2A(IW))
        QB = QNO2B(IW) + FW*(QNO2B(IW+1)-QNO2B(IW))
    1 continue
      QNO2 = QA + FT*(QB-QA)
        QNO2X = QNO2
        if (QA .gt. 0.d0) then
      QNO2X = QA * (QB/QA)**FTX
        endif
    2 continue
      return
      end
