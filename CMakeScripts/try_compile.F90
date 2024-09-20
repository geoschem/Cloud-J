program try_compile
    implicit none

    character (len = *), parameter :: FILE_NAME = "simple_xy.nc"
    integer, parameter :: NDIMS = 2
    integer, parameter :: NX = 6, NY = 12
    integer :: ncid, varid, dimids(NDIMS)
    integer :: x_dimid, y_dimid
    integer :: data_out(NY, NX)
    integer :: x, y
#ifndef NO_OMP
    integer nthreads, tid, OMP_GET_THREAD_NUM

    write(*,*) 'About to start multiple threads'
    ! Try OpenMP
    !$OMP PARALLEL PRIVATE(nthreads, tid)
    TID = OMP_GET_THREAD_NUM()
    write(*,*) 'Hello from thread ', tid
    !$OMP END PARALLEL
#endif

  contains

end program
