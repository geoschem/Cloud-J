!-------Add80_xH2O.f    calculate H2O X-section for v8.0 code
!  same format as modern v7.6

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
      real*8  WW, WNM, W11,W22
      character*80 TITLE,TITLTBL, TITL76
      character*6  TITLNEW,TSPEC(2)
      real*8   XNEW, XT
      integer :: TTT(2)
      character*1 :: ISXP


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
       call X_H2O (0, WW,XT,XNEW, TITLNEW,TITLTBL)

! flag for X sections, 'x' = strat only, 'p' = pressure interp (cannot be strat)
         ISXP = ' '
         TTT(1) = 300

!!!!!!!!!!  major temperature-or-density loop for X-sections
      do K = 1,1
         XT = TTT(K)

!!!!!!!!!   flux-weighted means over the 76-pratmo bins
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

           call X_H2O (1, WW,XT,XNEW, TITLNEW,TITLTBL)

           FFBIN(J) = FFBIN(J) + FY(I)
           AABIN(J) = AABIN(J) + FY(I)*XNEW
         I = I+1
        enddo
        if (I .ge. NY_) goto 2
       enddo
    2  continue
       do J=1,NC3
         AABIN(J) = AABIN(J)/FFBIN(J)
       enddo
       NB76 = NC3

!!!!!! Second integration from NC3 bins to NB3 bins (Fast-JX + not Solar-J)
!---combine fast-JX bins: non-SR bands (16:NB) are assigned a single JX bin
!---    SR bands are split (by Opacity Distrib Fn) into a range of JX bins
        FBIN(:) = 0.d0
        ABIN(:) = 0.d0
       do I=16,NC3
        J = IJX(I)
          FBIN(J) = FBIN(J) + FFBIN(I)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)
       enddo
       do I=1,15
        do J=1,8
          FBIN(J) = FBIN(J) + FFBIN(I)                  *SRB(I,J)
          ABIN(J) = ABIN(J) + FFBIN(I)*AABIN(I)         *SRB(I,J)
        enddo
       enddo
       do J=1,NB3
         ABIN(J) = ABIN(J)/FBIN(J)
       enddo

! fast-J v76-80 data for 'FJX_spec.dat'

        TITL76  = 'H2O     Pei2019'
        TITLNEW = 'H2O   '
        write(6,'(a)') TITL76
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'a',ABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',ABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',ABIN(13:18),TITLNEW
        write(6,*)

!!!  end of loop over T
      enddo

      stop
      end



c-----------------------------------------------------------------------
      subroutine X_H2O (INIT, WW,XT,XNEW, TITLNEW,TITLTBL)
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
      real*8 W1(999),XT1(999)
      real*8 T1,WW1
      integer NW,I,IW

      save NW,W1,XT1,T1

!---initialize Xsect tables, be sure to "save" them for recall.
      if (INIT .eq. 0) then

          TITLNEW = 'H2O   '
          FTBL = 'Add80_xH2O_Pei.dat'
        open (3, file=FTBL, status='OLD')
          read(3,'(a80)') TITLTBL
            write(6,'(2a/a)') ' openfile=',FTBL, TITLTBL
          read(3,'(i4,1x,a)') NW,FORMW
         do I=1,NW
          read(3,FORMW) W1(I),XT1(I)
         enddo
         close(3)
          T1 = 298.d0

      else

!---NOTE T-interpolation removed
!   BUT  this causes problems when a point is on the boudary edge, so trick it:
       WW1 = WW - 0.01d0
        IW = 1
       do I=1,NW-1
        if (WW1 .gt. W1(I)) IW = I+1
       enddo
         XNEW = XT1(IW)
       WW1 = WW + 0.01d0
        IW = 1
       do I=1,NW-1
        if (WW1 .gt. W1(I)) IW = I+1
       enddo

!       write(6,'(a,i5,f8.3,1p,3e10.3)')'XH2O:',IW,WW1,XNEW,XT1(IW)
       XNEW = 0.5d-25*(XT1(IW) + XNEW)

      endif
      return
      end
