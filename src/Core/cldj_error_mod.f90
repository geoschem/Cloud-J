! INTERFACE:

      MODULE CLDJ_ERROR_MOD

      IMPLICIT NONE

      PUBLIC  :: SAFE_DIV

      CONTAINS

!-----------------------------------------------------------------------
      function SAFE_DIV( numer, denom, alt_nan, alt_overflow, alt_underflow ) &
       result( quot )
      ! This funtion performs "safe division", that is to prevent overflow,
      ! underlow, NaN, or infinity errors.  An alternate value is returned
      ! if the division cannot be performed.

      ! This function is based on the equivalent function in GEOS-Chem
      ! (https://github.com/geoschem/geos-chem). 
!-----------------------------------------------------------------------

      ! Numerator and denominator
      real*8, intent(in) :: numer
      real*8, intent(in) :: denom

      ! Alternate value to be returned if the division is either NAN (0/0) or
      ! leads to overflow (i.e., a too large number)
      real*8, intent(in) :: alt_nan

      ! Alternate value to be returned if the division leads to overflow
      ! (default is ALT_NAN)
      real*8, optional, intent(in) :: alt_overflow
 
      ! Alternate value to be returned if the division leads to underflow
      ! (default is 0, but you could use TINY() if you want a non-zero result).
      real*8, optional, intent(in) :: alt_underflow

      ! Return value is output from division
      real*8 :: quot

      ! NaN
      if ( numer == 0 .and. denom == 0 ) THEN
         quot = alt_nan

      ! Overflow
      else if ( EXPONENT(numer) - EXPONENT(denom) >= MAXEXPONENT(numer) &
                .OR. Denom == 0 ) then
         quot = alt_nan
         if ( PRESENT(alt_overflow) ) quot = alt_overflow

      ! Underflow
      else if ( EXPONENT(numer) - EXPONENT(denom) <= MINEXPONENT(numer) ) then

         quot = 0D0
         if ( PRESENT(alt_underflow) ) quot = alt_underflow
  
      else
  
         ! No problem
         quot = numer / denom
  
      endif

      end function SAFE_DIV

      END MODULE CLDJ_ERROR_MOD
