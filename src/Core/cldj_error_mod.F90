MODULE CLDJ_ERROR_MOD

  IMPLICIT NONE

  PUBLIC  :: CLOUDJ_ERROR
  PUBLIC  :: CLOUDJ_ERROR_STOP
  PUBLIC  :: SAFE_DIV

  INTEGER, PUBLIC, PARAMETER :: CLDJ_SUCCESS =  0   ! Routine returns success
  INTEGER, PUBLIC, PARAMETER :: CLDJ_FAILURE = -1   ! Routine returns failure

CONTAINS

  !-----------------------------------------------------------------------
  subroutine CLOUDJ_ERROR_STOP( errmsg, loc )

#if defined( MAPL_ESMF )
#include "MAPL_Generic.h"
    USE MAPLBase_Mod
#elif defined( MODEL_CESM )
    USE CAM_ABORTUTILS,     ONLY : ENDRUN
#endif

    character(len=*), intent(in) :: errmsg
    character(len=*), intent(in) :: loc

    character(len=1023) :: msg
#if defined( MAPL_ESMF )
    integer :: rc
#endif

    ! Define message
    msg = 'CLOUDJ_ERROR_STOP: '//TRIM(errmsg)//' --> LOCATION: '//TRIM(loc)

    ! Stop run
#if defined( MAPL_ESMF )
    RC = CLDJ_FAILURE
    _ASSERT(RC==CLDJ_SUCCESS,TRIM(msg))
#elif defined( MODEL_CESM )
    call endrun(msg)
#elif defined( MODEL_GCCLASSIC )
    WRITE(6,*) TRIM(msg)
    call exit( 99999 )
#else
    WRITE(*,*) TRIM(msg)
    stop
#endif

  end subroutine CLOUDJ_ERROR_STOP

  !-----------------------------------------------------------------------
  subroutine CLOUDJ_ERROR( errmsg, loc, rc )

#if defined( MODEL_CESM )
    USE CAM_ABORTUTILS, ONLY : ENDRUN
#elif MAPL_ESMF
#include "MAPL_Generic.h"
    USE ESMF
    USE MAPLBase_Mod
#endif

    CHARACTER(LEN=*), INTENT(IN   )  :: errmsg  ! Message to display
    CHARACTER(LEN=*), INTENT(IN   )  :: loc     ! Location of error
    INTEGER,          INTENT(INOUT)  :: rc      ! Error code

    CHARACTER(LEN=1023) :: msg
#if defined( MAPL_ESMF )
    INTEGER             :: localPET, status
    CHARACTER(4)        :: localPETchar
    TYPE(ESMF_VM)       :: vm
#endif

#ifdef MAPL_ESMF
    ! Get current thread number
    CALL ESMF_VMGetCurrent(vm, RC=STATUS)
    CALL ESMF_VmGet( vm, localPET=localPET, __RC__ )
    WRITE(localPETchar,'(I4.4)') localPET
    msg = 'CLOUDJ_ERROR ['//TRIM(localPETchar)//']: '//TRIM(errmsg) &
         //' --> LOCATION: ' // TRIM(loc)
    WRITE(*,*) TRIM(msg)
#elif MODEL_CESM
    CALL ENDRUN('CLOUDJ_ERROR: '//TRIM(errmsg)//' --> LOCATION: '//TRIM(loc))
#else
    msg = 'CLOUDJ_ERROR: '//TRIM(errmsg)// ' --> LOCATION: '//TRIM(loc)
    WRITE(6,*) TRIM(msg)
    CALL Flush( 6 )
#endif

    ! Return with failure, but preserve existing error code
    IF ( rc == CLDJ_SUCCESS ) THEN
       rc = CLDJ_FAILURE
    ENDIF

  end subroutine CLOUDJ_ERROR

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
