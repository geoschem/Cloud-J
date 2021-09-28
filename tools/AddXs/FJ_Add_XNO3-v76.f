c-------FJ_Add_XNO3-v76.f    modern code for integrating X-NO3 cross sections
! uses both photons and Watts to do weighting  Example = water liqu & ice
! starts with pratmo full 77+SR bins and then generates the Solar-J bins (18+9)

      implicit none
      integer, parameter :: NC_ = 99
      integer, parameter :: NB_ = 27
      integer, parameter :: NX_ = 6785
      integer, parameter :: NY_ = 40000


      real*8   SRB(15,8)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
      integer  IBINX(NX_),IBINY(NY_)

      real*8, dimension(NC_) ::
     &          FFBIN,WWBIN,RRBIN,AABIN,BBBIN,CCBIN,DDBIN

      real*8  W1(NB_),W2(NB_),WNM, RAYLAY,YPAR
      real    WL,RNW,CNW,RNI,CNI,A1,A2
      real*8  WWL

      real*8  WX(NX_),FX(NX_)
      real*8  WY(NY_),FY(NY_)

      real*8, dimension(NB_) :: FBIN,WBIN,RBIN,ABIN,BBIN,CBIN,DBIN
      real*8, dimension(NB_) ::
     &           FBINw,WBINw,RBINw,ABINw,BBINw,CBINw,DBINw
      character*78 TITLE,TITLNEW
      real*8       TT, XXWTa, XXWTb, QYa, QYb
      integer  I, J,J1,J2, INIT, IT1,IT2,It3,IT4,IT5,IT6

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
!---now assign bin #(I=1:NC4) to each ASTM2 microbin J (1:1697)
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
!!!!!!!!!!!!!!!!!!!!!!! finished setup !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       WWL = 0.0
       INIT = 0
       call X_NO3 (WWL,TT, XXWTa,XXWTb,INIT,TITLNEW)
       INIT = 1
       IT1 = 190
       IT2 = 298

!!!! integration of refractive index at high-res for photons NY_, WY(),FY()
        FFBIN(:) = 0.d0
        WWBIN(:) = 0.d0
        RRBIN(:) = 0.d0
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
        WWL = WY(I)
          TT = IT1
          call X_NO3 (WWL,TT, XXWTa,XXWTb,INIT,TITLNEW)
        FFBIN(J) = FFBIN(J) + FY(I)
        WWBIN(J) = WWBIN(J) + FY(I)*WY(I)
        RRBIN(J) = RRBIN(J) + FY(I)*RAYLAY(WY(I))
        AABIN(J) = AABIN(J) + FY(I)*XXWTa
        BBBIN(J) = BBBIN(J) + FY(I)*XXWTb
          TT = IT2
          call X_NO3 (WWL,TT, XXWTa,XXWTb,INIT,TITLNEW)
        CCBIN(J) = CCBIN(J) + FY(I)*XXWTa
        DDBIN(J) = DDBIN(J) + FY(I)*XXWTb
        I = I+1
       enddo
        WWBIN(J) = WWBIN(J)/FFBIN(J)
        RRBIN(J) = RRBIN(J)/FFBIN(J)
        AABIN(J) = AABIN(J)/FFBIN(J)
        BBBIN(J) = BBBIN(J)/FFBIN(J)
        CCBIN(J) = CCBIN(J)/FFBIN(J)
        DDBIN(J) = DDBIN(J)/FFBIN(J)
       if (I .ge. NY_) goto 2
      enddo
    2 continue

!!!! NC4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
!      write(6,'(a)') TITLNEW
      write(6,'(2a)') 'pratmbin#    solflx    w-eff  X-Raylay    ',
     &    'XNO3a & XNO3b @ 190K & 298K'
      do J = 1,NC4
        write(6,'(a5,i4,1p,e10.3,0p,f10.4,1p,5e10.3)') 'Phot ',J,
     &   FFBIN(J),WWBIN(J),RRBIN(J),AABIN(J),BBBIN(J),CCBIN(J),DDBIN(J)
      enddo

!!!! Second integration from NC4 bins to NB4 bins (Fast-JX + Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        WBIN(:) = 0.d0
        RBIN(:) = 0.d0
        ABIN(:) = 0.d0
        BBIN(:) = 0.d0
        CBIN(:) = 0.d0
        DBIN(:) = 0.d0
       do I=16,NC4
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          WBIN(J) = WBIN(J) + FFBIN(I)*WWBIN(I)
          RBIN(J) = RBIN(J) + FFBIN(I)*RRBIN(I)
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
          ABIN(J) = ABIN(J)/FBIN(J)
          BBIN(J) = BBIN(J)/FBIN(J)
          CBIN(J) = CBIN(J)/FBIN(J)
          DBIN(J) = DBIN(J)/FBIN(J)
        endif
       enddo

!!!! NB4 bins: Photon weighted values for w, Rayleigh, liq-water, ice-water
!      write(6,'(a)') TITLNEW
      write(6,'(2a)') 'pratmbin#    solflx    w-eff  X-Raylay    ',
     &    'XNO3a & XNO3b @ 190K & 298K'
      do J = 1,NB4
        write(6,'(a5,i4,1p,e10.3,0p,f10.4,1p,5e10.3)') 'Phot ',J,
     &   FBIN(J),WBIN(J),RBIN(J),ABIN(J),BBIN(J),CBIN(J),DBIN(J)
      enddo

! fast-J v73 data for 'FJX_spec.dat'

        TITLNEW = 'solflx'
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',FBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',FBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',FBIN(13:18),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    d',FBIN(19:21)

        TITLNEW = 'X-rayl'
        write(6,'(a)') TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    a',RBIN(1:6),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    b',RBIN(7:12),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    c',RBIN(13:18),TITLNEW
        write(6,'(a5,1p,6e10.3,1x,a6)')
     &  '    d',RBIN(19:21)

        TITLNEW = 'w-eff '
        write(6,'(a)') TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    a',WBIN(1:6),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    b',WBIN(7:12),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    c',WBIN(13:18),TITLNEW
        write(6,'(a5,6f10.3,1x,a6)')
     &  '    d',WBIN(19:21)

        TITLNEW = 'NO3a  = NO2+O'
        write(6,'(a)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'a',ABIN(1:6),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'b',ABIN(7:12),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'c',ABIN(13:18),TITLNEW(1:6)
        write(6,'(a)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'a',CBIN(1:6),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'b',CBIN(7:12),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'c',CBIN(13:18),TITLNEW(1:6)

        TITLNEW = 'NO3b  = NO+O2'
        write(6,'(a)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'a',BBIN(1:6),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'b',BBIN(7:12),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'c',BBIN(13:18),TITLNEW(1:6)
        write(6,'(a)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'a',DBIN(1:6),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'b',DBIN(7:12),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'c',DBIN(13:18),TITLNEW(1:6)


        do I = 1,18
        ABIN(I) = ABIN(I)+BBIN(I)
        CBIN(I) = CBIN(I)+DBIN(I)
        enddo
        TITLNEW = 'NO3   q(NO2+O)=0.882, q(NO+O2)=0.112'
        write(6,'(a)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'a',ABIN(1:6),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'b',ABIN(7:12),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT1,'c',ABIN(13:18),TITLNEW(1:6)
        write(6,'(a)') TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'a',CBIN(1:6),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'b',CBIN(7:12),TITLNEW(1:6)
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',IT2,'c',CBIN(13:18),TITLNEW(1:6)

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


      subroutine X_NO3 (WWX,TT, XXWTa,XXWTb,INIT,TITLNEW)
      implicit none
      real*8, intent(in) :: WWX, TT
      integer, intent(inout) :: INIT
      real*8, intent(out) :: XXWTa,XXWTb
      character*6, intent(out) :: TITLNEW

      real*8, parameter :: T1=-1096.4d0, T2=-529.5d0
     &                    ,T1_298=T1/298.d0,T2_298=T2/298.d0
      character*80 FTBL,TABLE,FORMW,FORMB
      real*8, dimension(999) :: W,XW,WQY,QYb298,QYb190  ! NO3b (NO+O2)
     &                         ,QYa298,QYa190           ! NO3a (NO2+O)
      real*8  WW,TTL,XBFACT,XTFACT,XXW,FW,FWQ, XWI,XWIP1,TFACT
      real*8  QYa,QYb,QYa230,QYb230
      integer NW,NB,NQ, N, I,IW,IWQ

      save  W,XW,NW, WQY,QYb298,QYb190,QYa298,QYa190,NQ

      XTFACT(TTL) = (1.d0-dexp(T1/TTL)-2.d0*dexp(T2/TTL)) /
     &              (1.d0-dexp(T1_298)-2.d0*dexp(T2_298))

      TITLNEW = 'xNO3ab'
      FTBL = 'XNO3_JPL10tbl.dat'
      if (INIT .eq. 0) then
        open (3, file=FTBL, status='OLD')
          read(3,'(a)') TABLE
           write(6,'(2a/a)') ' openfile=',FTBL, TABLE
          read(3,'(i4,1x,a)') NW
        do N = 1,NW
          read(3,*) W(N),XW(N)
!           write(6,*) N,W(N),XW(N)
         enddo

!        write(6,*) '  finish XNO3'
        do N = 1,9
          read(3,*)
        enddo
          read(3,'(i4,1x,a)') NQ
        do N = 1,NQ
          read(3,*) WQY(N),QYb298(N),QYb230,QYb190(N),
     &           QYa298(N),QYa230,QYa190(N)
!            write(6,*) N,WQY(N), QYb298(N),QYa190(N)
        enddo
        close(3)
        INIT = 1
      else

c---interpolate X-section vs. Wavelength, convert to nm!
        WW = 1.d3*WWX

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
