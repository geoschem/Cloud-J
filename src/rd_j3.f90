      implicit none

      character*40  ICTfile
      character*72  TITLE1
      character*9   TJ(10)

      integer       I,J,L,N

      real*4        SZA(20), XJ(57,10,20), Z(58)

      real*4     H(57),X1(57),X2(57),X3(57)

      SZA(:) = 0.0

      ICTfile = 'J3.out'
      open (11,file=ICTFILE,status='old',err=9)

        read (11,*,err=9)  TITLE1
        read (11,*,err=9)  TITLE1
        do I=1,314
        read(11,*)
        enddo

      do J = 1,18
        read(11,'(28x,f6.2)') SZA(J)
        do I=1,306
        read(11,*)
        enddo
        read(11,'(67x,f7.2)') (H(L), L=1,57)
        read(11,*)
        read(11,*)
        read(11,'(13x,e9.2,54x,2e9.2)') (X1(L),X2(L),X3(L), L=1,57)

        do L=1,57
         XJ(L,1,J) = H(L)
         XJ(L,2,J) = X1(L)
         XJ(L,3,J) = X2(L)
         XJ(L,4,J) = X3(L)
        enddo

      enddo

      write(6,'(a)') ' 1:18 heating rates'
      write(6,'(5x,18f9.1)') (SZA(J), J=1,18)
      do L=1,57
      write(6,'(i5,18f9.2)') L,(XJ(L,1,J), J=1,18)
      enddo

      write(6,'(a)') ' J-O3'
      write(6,'(5x,18f9.1)') (SZA(J), J=1,18)
      do L=1,57
      write(6,'(i5,1p,18e9.2)') L,(XJ(L,2,J), J=1,18)
      enddo

      write(6,'(a)') ' J-NO2'
      write(6,'(5x,18f9.1)') (SZA(J), J=1,18)
      do L=1,57
      write(6,'(i5,1p,18e9.2)') L,(XJ(L,3,J), J=1,18)
      enddo

    9 stop
      end
