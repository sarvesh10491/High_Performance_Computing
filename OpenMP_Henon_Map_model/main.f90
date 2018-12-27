!   Script to find Fractal dimensions and C parameter of the uncertainty
!   calculation  
!   Main program for the Assignent 5
    program test_naive
    use henondim 
    use precision
    use lsq
    use, intrinsic:: iso_fortran_env
    implicit none
    
    integer:: n  ! number of points on a side
    integer:: ierr,k
    real(DP)::start,finish
    real(DP)::noise(10)
    real(DP)::Neps(10)
    real(DP)::param(2)
    

    print *,"Program execution started\n"
    n=4096
    
    noise(:)=(/((2**k),k=12,21)/)
    noise=1/noise
   
    call basin_alg(n,noise,Neps,ierr)
    
    write(*,*) "Count vector",Neps

    call linfit("p",size(Neps),noise,Neps,param,ierr)

    write(*,*) "Parameter C and P are:",param
    Write(*,*) "Thus the fract dimensions:",param(2)+1
end program test_naive
