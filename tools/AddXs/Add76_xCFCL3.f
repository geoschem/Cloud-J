c-------Add76_xCFCl3.f   with user supplied subroutine that supplies X-section x Q-yield
!   X-sections - generates only pratmo 76 bins abnd then Fast-J 18 bins

      implicit none
      integer, parameter :: NC_ = 199
      integer, parameter :: NB_ = 99
      integer, parameter :: NY_ = 40000

      real*8   SRB(15,8)
      real*8  WY(NY_),FY(NY_)
      real*8, dimension(NC_) :: WCBIN
      integer,dimension(NC_) :: IJX
      real*8, dimension(NC_) :: FFBIN,AABIN
      real*8, dimension(NB_) :: FBIN,ABIN
      integer  I, J,J1,J2,K,  NC1,NC2,NC3,NC4, NB3,NB4, NB76
      real*8  WW, WNM, RAYLAY,YPAR, W11,W22
      character*80 TITLE
      Character*90 TITLTBL, TITL76
      character*6  TITLNEW,TSPEC(2)
      real*8   XNEW,QNEW,XT,XP,XM, MM(3)
      integer :: TTT(3),PPP(3)
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

      call XCFCL3 (0, WW,XT,XP,XM, XNEW,
     &     MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)

! flag for X sections, 'x' = strat only, 'p' = pressure interp (cannot be strat)
            ISXP = ' '
         if (ISX .eq. 'x') then
            ISXP = ISX
         else if (ISP .eq. 'p') then
            ISXP = ISP
         endif


!!!!!!!!!!  major temperature-or-density loop for X-sections
      do K = 1,3
      if (TTT(K) .gt. 0) then
         XT = TTT(K)
         XP = PPP(K)
         XM = MM(K)

!!!!!!!!!   flux-weighted means over the 76-pratmo bins
!!!!!!!!!  NB, the  0.05 nm bins in v76 are NOT split but fall into ONE of the 76 pratmo bins
!!!!!!!!!     Here the microbin 635.10 nm falls into the 635.10-778.00 nm pratmo bin (#76)
!!!!!!!!!     In some earlier versions it was triggered to fall into pratmo #75 = 560.10-635.10 nm
!!!!!!!!!     Differences are at most 5-7% and the allocation of solar flux is consistent
!!!!!!!!!  Thus some earlier versions of X-sections are slightly different.  Not all have been replaced
        FFBIN(:) = 0.d0
        AABIN(:) = 0.d0
         J = 1
         I = 1
       do while (WY(I) .lt. WCBIN(J))
         I = I+1
       enddo
       do J = 1,NC3
        do while (WY(I) .lt. WCBIN(J+1) .and. I .lt. NY_)
           WW = WY(I)*1.d3

           call XCFCL3 (1, WW,XT,XP,XM, XNEW,
     &        MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)

           FFBIN(J) = FFBIN(J) + FY(I)
           AABIN(J) = AABIN(J) + FY(I)*XNEW
         I = I+1
        enddo
        if (I .ge. NY_) goto 2
       enddo
    2  continue
       do J=1,NC3
        if (FFBIN(J) .gt. 0.d0) then
          NB76 = J
          AABIN(J) = AABIN(J)/FFBIN(J)
        endif
       enddo

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
         ABIN(J) = ABIN(J)/max(FBIN(J),1.d-40)
       enddo

! fast-J v76 data for 'FJX_spec.dat'
        write(6,'(a90)') TITLTBL
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ISXP,TTT(K),'a',ABIN(1:6),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &  ' ',TTT(K),'b',ABIN(7:12),TITLNEW
        write(6,'(a1,i3,a1,1p,6e10.3,1x,a6)')
     &   ' ',TTT(K),'c',ABIN(13:18),TITLNEW
        write(6,*)

! 'pratmo' 76-bin tables
        TITL76 = TITLTBL
        write(6,'(a20,f5.0/(1p,8e10.3))') TITL76,XT,(AABIN(I),I=1,NB76)
        write(6,*)

!!!  end of loop over T & p
      endif
      enddo

      stop
      end



c-------------other species can be modeled after sub XCFCL3
c-----------------------------------------------------------------------
      subroutine XCFCL3 (INIT, WW,XT,XP,XM,XNEW,
     &     MM,TTT,PPP,ISX,ISP,TITLNEW,TITLTBL)
c-----------------------------------------------------------------------
c   WW = wavelength (nm) to calc Xsection for
c   XT = temerature (K) for interpolation
c   XP = pressure (hPa) can be used in Stern-Volmer formula if need be
c   XM = air density (#/cm3), ditto
c   XNEW = cross section (cm2) as a function of WW and XT (and XP, XM)
c   INIT = initialization:
c     if INIT.eq.0 .then reads in any tables and sets Xsect name to TITLNEW
c-----------------------------------------------------------------------
      implicit none
      save W,XW,NW,XB,NB

      integer, intent(in) :: INIT
      real*8, intent(in) :: WW,XT,XP,XM
      real*8, intent(out) :: XNEW
      real*8, intent(out) :: MM(3)
      integer, intent(out) :: TTT(3),PPP(3)
      character*1, intent(out) :: ISX,ISP
      character*6, intent(out) :: TITLNEW
      character*90, intent(out) :: TITLTBL

      character*80 FTBL,TABLE,FORMW,FORMB
      real*8 W(999), XW(999),XB(99), WWL,TTL,XBFACT,XTFACT,XXW,FW
      integer NW,NB,  N, I,IW

      FTBL = 'Add76_xCFCl3_JPL11.dat'
      TITLNEW = 'CFCL3 '

!---include the Temperatures (p, M) that you want to interpolate to and give tables.
!---general format for generating all Xsections
      if (INIT .eq. 0) then

          write(6,'(2a)') ' species:',TITLNEW
        open (3, file=FTBL, status='OLD')
            write(6,'(2a/a)') ' openfile=',FTBL
          read(3,'(a80)') TABLE
            write(6,'(a80)') TABLE
          read(3,'(4x,a1)') ISX
          read(3,'(4x,a1)') ISP
            write(6,'(a4,a2,4x,a4,a2)') 'ISX=',ISX,'ISP=',ISP
         if (ISP .eq. 'p') then
! if doing pressure interpolation for q-yield, specify p (hPa, int) and density (e19)
          read(3,'(5x,3i5)') PPP
               write(6,'(a4,4i4)') 'ppp=',PPP
          read(3,'(5x,3f5.0)') MM
             MM(:) = MM(:)*1.e19
         endif
! always read Temperature, can be up to 3, in increasing order if data tables exist
         read(3,'(5x,3i5)') TTT
             write(6,'(a4,4i4)') 'TTT=',TTT
         read(3,'(a90)') TITLTBL
!---this is all table specific formats, may need to add columns for T's
         read(3,'(i4,1x,a)') NW,FORMW
        do N=1,NW
         read(3,FORMW) W(N),XW(N)
        enddo
         read(3,'(i4,1x,a)') NB,FORMB
        do N=1,NB
         read(3,FORMB) XB(N)
        enddo
        close(3)

      else

!-- this section may need to be cutomizzed for the type of tabulated x-sections
!-- CFCl3 has some odd polynomial functions for T & wavelength

!---interpolate X-section vs. Wavelength
        IW = 1
      do I=2,NW-1
        if (WW .gt. W(I)) IW = I
      enddo
        FW = (WW - W(Iw))/(W(IW+1) - W(IW))
        FW = min(1.d0, max(0.d0, FW))
      XXW = XW(IW) + FW*(XW(IW+1)-XW(IW))
c---NB CFCl3 Xsections scaling valid for 220-296 K and 200-238 nm only
c---apply min-max range for T dependence: set for CFCl3
        WWL = min(238.d0, max(200.d0, WW))
        XBFACT = XB(1) + XB(2)*(WWL-200.d0) + XB(3)*(WWL-200.d0)**2
        TTL = min(296.d0, max(220.d0, XT))
        XTFACT = exp((TTL-296.d0)*XBFACT)
      XNEW = XXW  * XTFACT * 1.e-20

      endif

      return
      end
