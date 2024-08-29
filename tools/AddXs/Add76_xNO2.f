c-------Add76_xNO2.f   with user supplied subroutine that supplies X-section x q-yields
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
      real*8  WW, XT, WNM
      character*80 TITLE, TITLTBL,TITL76
      character*6  TITLNEW
      real*8 XNO2, XNO2X, QNO2, QNO2X, TOTA, TOTB
      integer :: TTT(4)
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


!!!!!!!!!!!!! iniitalize
!not extrapolated, uses Q-yld at 248K min
        TTT(1) = 220
        TTT(2) = 294
!        TITLTBL(2) = 'NO2   '
!std log extrapolation
        TTT(3) = 200
        TTT(4) = 300
!        TITLTBL(4) = 'NO2ex '
        ISXP = ' '

      do K = 1,4
        XT = TTT(K)

!!!!!!!!!  flux-weighted means over the 76-pratmo bins
!!!!!!!!!  some earlier versions of X-sections are slightly (<1%) different.  Not all have been replaced
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
        BBBIN(:) = 0.d0
         J = 1
         I = 1
! primary high-resolution wavelength loop - generate input for pratmo reference J's
       do while (WY(I) .lt. WCBIN(J))
         I = I+1
       enddo
       do J = 1,NC3
        do while (WY(I) .lt. WCBIN(J+1) .and. I .lt. NY_)
           WW = WY(I)*1.d3

           call X_NO2 (WW, XT, XNO2,XNO2X)
           call Q_NO2 (WW, XT, QNO2,QNO2X)

           FFBIN(J) = FFBIN(J) + FY(I)
           AABIN(J) = AABIN(J) + FY(I)*XNO2 *QNO2
           BBBIN(J) = BBBIN(J) + FY(I)*XNO2X*QNO2X
         I = I+1
        enddo
        if (I .ge. NY_) goto 2
       enddo
    2   continue
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

       write(6,'(a20,i5,2f10.6)') 'TTT, J-NO2, J-NO2ex',TTT(K),TOTA,TOTB

! fast-J v76 data for 'FJX_spec.dat'
        TITLNEW = 'NO2   '
        TITL76  = 'NO2 bounded interp   JPL2010'
        write(6,'(a)') TITL76
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',ABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',ABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',ABIN(13:18),TITLNEW
        write(6,*)

! 'pratmo' 76-bin tables
        write(6,'(a40,f5.0/(1p,8e10.3))') TITL76,XT,(AABIN(I),I=1,NB76)
        write(6,*)
        TITL76  = 'NO2 extrap    JPL2010'

        TITLNEW = 'NO2ex '
        write(6,'(a)') TITL76
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',BBIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',BBIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',BBIN(13:18),TITLNEW
        write(6,*)

! 'pratmo' 76-bin tables
        write(6,'(a40,f5.0/(1p,8e10.3))') TITL76,XT,(BBBIN(I),I=1,NB76)
        write(6,*)
        write(6,*)

!!!  end of loop over T
      enddo

      stop
      end



      subroutine X_NO2 (WW, TT, XNO2,XNO2X)
c---JPL_2010 standard tables.
c---interpolates NO2 cross section vs. wavelength and temperature.
c---    XNO2 is linear, limited to range 220 to 294K
c---    XNO2X does log interpolation and extrapolates beyond 220K & 294K
c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: XNO2,XNO2X
      integer IW,I
      real*8 Q1,Q2,FT,FTX

C----JPL 2010 Table 4C-2.
      real*4, parameter, dimension(43)  :: WNM =
     &[240.964, 243.902, 246.914, 250.000, 253.165, 256.410, 259.740,
     & 263.158, 266.667, 270.270, 273.973, 277.778, 281.690, 285.714,
     & 289.855, 294.118, 298.507, 303.030, 307.692, 312.5, 317.5, 322.5,
     & 327.5, 332.5, 337.5, 342.5, 347.5, 352.5, 357.5, 362.5, 367.5,
     & 372.5, 377.5, 382.5, 387.5, 392.5, 397.5, 402.5, 407.5, 412.5,
     & 417.5, 422.5, 427.5]
      real*4, parameter, dimension(42)  :: X220K =
     &[4.14, 0.961, 0.859, 0.191, 0.496, 0.872, 1.26, 1.77, 2.36, 3.03,
     & 3.94, 5.16, 6.29, 7.72, 9.64, 11.6, 13.2, 16.0, 18.5, 20.8, 24.2,
     & 27.2, 29.4, 33.0, 37.0, 38.6, 43.5, 47.7, 49.2, 53.7, 55.2, 58.4,
     & 58.5, 59.2, 62.4, 58.5, 64.0, 57.0, 61.8, 58.3, 59.3, 56.0]
      real*4, parameter, dimension(42)  :: X294K =
     &[5.77,  2.79,  1.62,  0.998, 1.05,  1.28, 1.58, 2.05, 2.64, 3.24,
     & 4.07, 5.21, 6.23, 7.59, 9.51, 11.5, 13.2, 16.1, 18.8, 21.6, 25.3,
     & 28.7, 31.7, 35.8, 40.2, 41.8, 46.2, 49.7, 50.9, 54.9, 56.1, 59.0,
     & 59.3, 60.1, 63.0, 59.7, 64.4, 58.2, 62.4, 59.1, 59.9, 57.0]

        XNO2 = 0.d0
        XNO2X = 0.d0
      if (WW .lt. WNM(1) .or. WW .gt. WNM(43)) goto 2
      do I=1,42
       if (WW .gt. WNM(I)) IW = I
      enddo
        Q1 = X220K(IW)*1.d-20
        Q2 = X294K(IW)*1.d-20
        FTX = (TT - 220.d0) / (294.d0 - 220.d0)
        FT = max (0.d0, min (1.d0, FTX))
        XNO2 = Q1 + FT*(Q2-Q1)
        XNO2X = XNO2
      if (Q1 .gt. 0.d0) then
        XNO2X = Q1 * (Q2/Q1)**FTX
      endif
    2 continue
      return
      end


      subroutine Q_NO2 (WW, TT, QNO2,QNO2X)
c---JPL_2010 standard tables.
c---interpolates NO2 quantum yield vs. wavelength and temperature.
c---    QNO2 is linear, limited to range Q248 to Q298
c---    QNO2X does log interpolation and extrapolates beyond 248K & 298K
c---JPL 2006 recommendations
      implicit none
      real*8, intent(in) :: WW, TT
      real*8, intent(out):: QNO2,QNO2X
      integer IW
      real*8 FW,Q1,Q2,FT,FTX

C----------JPL 2010 Table 4C-3.
      real*4, parameter, dimension(26)  :: WNM =
     & [398., 399., 400., 401., 402., 403., 404., 405., 406., 407.,
     &  408., 409., 410., 411., 412., 413., 414., 415., 416., 417.,
     &  418., 419., 420., 421., 422., 423.]
      real*4, parameter, dimension(26)  :: Q298 =
     & [1.00, 0.95, 0.88, 0.75, 0.62, 0.53, 0.44, 0.37, 0.30, 0.26,
     &  0.22, 0.18, 0.15, 0.13, 0.11, 0.09, 0.08, 0.06, 0.05, 0.04,
     &  0.03, 0.02, 0.02, .015, 0.01, 0.00]
      real*4, parameter, dimension(26)  :: Q248 =
     & [1.00, 0.94, 0.86, 0.69, 0.56, 0.44, 0.34, 0.28, 0.22, 0.18,
     &  0.14, 0.12, 0.10, 0.08, 0.07, 0.06, 0.04, 0.03, 0.02, 0.02,
     &  0.02, 0.01, 0.01, 0.01, 0.01, 0.00]

      IW = (WW - 397.d0)
      IW = max(1, min(25, IW))
      FW = WW - WNM(IW)
      FW = max(0.d0, min(1.d0, FW))
      Q2 = Q298(IW) + FW*(Q298(IW+1)-Q298(IW))
      Q1 = Q248(IW) + FW*(Q248(IW+1)-Q248(IW))
      FTX = (TT - 248.d0) / (298.d0 - 248.d0)
      FT = max (0.d0, min (1.d0, FTX))
      QNO2 = Q1 + FT*(Q2-Q1)
      QNO2X = QNO2
      if (Q1 .gt. 0.d0) then
        QNO2X = Q1 * (Q2/Q1)**FTX
        QNO2X = min (1.d0, QNO2X)
      endif
      return
      end
