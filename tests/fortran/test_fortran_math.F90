      program test_fortran_math
      USE CLDJ_CMN_MOD
      USE CLDJ_FJX_SUB_MOD, ONLY : X_INTERP
      implicit none
      
      real*8 :: temp, t1, x1, t2, x2, t3, x3, res
      integer :: lqq_in

      ! Read inputs from stdin line-by-line
      ! Format: temp, t1, x1, t2, x2, t3, x3, lqq_in
      do while (.true.)
         read(*,*,end=99) temp, t1, x1, t2, x2, t3, x3, lqq_in
         call X_INTERP(temp, res, t1, x1, t2, x2, t3, x3, lqq_in)
         write(*,'(1p,e24.15)') res
      enddo

   99 continue
      end program test_fortran_math
