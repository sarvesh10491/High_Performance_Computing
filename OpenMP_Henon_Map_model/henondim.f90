module henondim
!  HENONDIM - Implement Algorithm D to determine the dimension of the
!  basin of
!  infinity of the Henon map.
!
!  SYNOPSIS
!    use henondim
!
!  DESCRIPTION
!  WE first find the basin vector on original points by interating the grid
!  formed by equal distribution from -3 to 3 for n=4096 for maxiter which in
!  our case is 100, for doing this operation we have defined a subroutine
!  "hennon_map".Then we are checking that adding the uncertainty into the
!  original grid's  horizontal coordinates will result into change in the basin
!  vector or not.Since each computation is just the change in noise value we are
!  multithreading based on the noise values thus for checking 20 noise values
!  from +/-2^-12 to +/- 2^-21 we are runnig a do loop and parallelizing using
!  openMP.For grid size n=4096 we are getting 3 times speed up against the
!  serial program.    
!    
!
!  PUBLIC ROUTINES DEFINED
!    Brief list of procedures that can be called from outside this
!    module.
!    1] henon_map()
!    2] basin_algo()
!
!  REQUIRED DEPENDENCIES
!    'precision' - to define KINDs for single and double floating-point.
!    And any other modules that you define.
!   
!   'lsq':Interface module for LAPACK QR routines, and a special-purpose
!    subroutine for simple Y versus X linear regressions.
!
!  REVISION HISTORY
!   10/28/18 - Final draft.
!
!  PROGRAMMER
!   Sarvesh Patil, sspati13@asu.edu
!
!  Collaborator
!   Shubham Nandanwankar, snandanw@asu.edu
!------------------------------------------------------------------------------
!  External Modules and Global variable definations
    use precision
    use lsq
    implicit none
    real(DP), parameter:: LOCKOUT=100.0
    real(DP), parameter:: BOXMIN=-3.0, BOXMAX=3.0 
    contains
!-------------------------------------------------------------------------------
!  Iterate the Henon map for one or more initial conditions up to
!  MAXITER times or until the orbit of one of the points is further than
!  LOCKOUT units away from the origin.  (It suffices to check whether
!  the
!  |X_n| > LOCKOUT on the nth iteration.)
!  
!  To limit the computation WHERE block is used, this block will limit the
!  computations to those grid vector points whose distance from the original
!  state is less than lockout.

!  Argument declarations go here.  Explicitly declare all arguments and
!  give
!  each of them an INTENT.
!  prev:   Read only grid value one copy of the original
!  n:      n**2 is the size of dimension in either X or Y direction
!  maxiter:Maximum no. of iterations to find the basin vector
!  param_a:a,b values in henon map calculation at 1,2 index resp
!  lockout:Boundeness check value
!  basin:  Writable vector which will store the calculated basin vector
!  noise:  A integer containing epsilon value both +/-  

subroutine henon_map(prev,n,maxiter,param_a,lockout,basin,noise,ierr_h)
   
    integer, intent(in)::n
    real(DP), intent(in):: prev(n*n,2)
    real(DP), intent(in):: maxiter
    real(DP), intent(in):: param_a(2)
    real(DP), intent(in):: lockout
    real(DP), intent(inout)::basin(n*n,1)
    integer, intent(out):: ierr_h
    real(DP), intent(in)::noise    

    integer::myid,count_thread
    integer::thread_block,num_threads
    integer, external::OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
    
    real(DP), dimension(:,:), allocatable::X_MESH,prev_X,v,d    
    integer::j
    
    allocate(X_MESH(n*n,2),v(n*n,2),prev_X(n*n,1))
    X_MESH=prev
!    d=0.0
    v=0.0
    X_MESH(:,1)=X_MESH(:,1)+noise
!    num_threads=OMP_GET_NUM_THREADS()
!    thread_block=(n*n)/num_threads
!    count_thread=0
   do j=1,(maxiter)
 !       myid=OMP_GET_THREAD_NUM() 
       WHERE (basin(:,1)<lockout)
            prev_X(:,1)=X_MESH(:,1)
            X_MESH(:,1)=param_a(1)-X_MESH(:,1)**2+param_a(2)*X_MESH(:,2)
            X_MESH(:,2)=prev_X(:,1)
            v(:,1)=(X_MESH(:,1)-prev(:,1))**2!,[thread_block])
            v(:,2)=(X_MESH(:,2)-prev(:,2))**2!,[thread_block])
            basin(:,1)=SQRT(v(:,1)+v(:,2))!,[thread_block])
        ENDWHERE
    enddo
!    write(*,*) "Dist vect is",d
    WHERE(basin(:,1)>lockout)
        basin(:,1)=1.0
    ELSEWHERE
        basin(:,1)=0.0
    ENDWHERE
!    basin=d
    ierr_h=0
    deallocate(X_MESH,v,prev_X)
    return
    end subroutine henon_map
!--------------------------------------------------------------------------------------------
!Description:
!   This is just a wrapper around henon_map to ease the multithreading memory
!   operations.Doing so will make this subroutine similar to that of a thread
!   function while using posix thread.This will also help in reduing the stack
!   requirement of each thread.
!Argument:
!   X_MESH:      Read only grid value one copy of the original
!   n:           n**2 is the size of dimension in either X or Y direction
!   maxiter:     Maximum no. of iterations to find the basin vector
!   param_a:     a,b values in henon map calculation at 1,2 index resp
!   lockout:     Boundeness check value
!   basin_arg:   Basin vector calculated when no noise is added
!   noise:       A integer containing epsilon value both +/-                      
!   Neps:        Matrix which stores no. of points differering from orginal
!                basin vector for both +/- noise
!   m:           Index to write values properly in Neps
!   ierr_w:      error variable
subroutine wrap_henon_map(X_MESH,n,maxiter,param_a,lockout,basin_org,noise,Neps,m,ierr_w)
    
    integer,intent(in)::n
    real(DP),intent(in)::X_MESH(n*n,1)
    real(DP),intent(in)::maxiter
    real(DP),intent(in)::param_a(2)
    real(DP),intent(in)::lockout
    real(DP),intent(inout)::Neps(n)
    real(DP),intent(in)::noise
    real(DP),intent(in)::basin_org(n*n,1)
    integer, intent(out):: ierr_w
    integer,intent(in)::m

!    real(DP), dimension(:,:), allocatable::basin_eps
    real(DP)::basin_eps(n*n,2)
    integer::ierr_h,z,pos,neg

    write(*,*) "Inside wrapper with noise:",noise
    pos=0
    neg=0
!    if(noise.eq.0.d0) then
!        call henon_map(X_MESH,n,maxiter,param_a,lockout,basin_vect(:,1),0.d0,ierr_h)
!        return
!    else 
!    allocate(basin_eps(n*n,2))
    call henon_map(X_MESH,n,maxiter,param_a,lockout,basin_eps(:,1),noise,ierr_h)            
    call henon_map(X_MESH,n,maxiter,param_a,lockout,basin_eps(:,2),-noise,ierr_h)
    Neps(m)=count((basin_org(:,1).ne.basin_eps(:,1)).OR.(basin_org(:,1).ne.basin_eps(:,2)))
!        do z=1,n*n
!            pos=(basin_org(z,1)/=basin_eps(z,1))
!            neg=(basin_org(z,1)/=basin_eps(z,2))
!            if(pos.OR.neg) then
!                Neps(m)=Neps(m)+1.0
!        endif
 
   write(*,*)"count for noise:",noise," is:",Neps(m)
!    end if
    ierr_w=0
!    deallocate(basin_eps)
    return    

end subroutine wrap_henon_map
!-------------------------------------------------------------------------------
subroutine basin_alg(n,noise,Neps,ierr)
! BASIN_ALG - implements Algorithm D.
!   n:    grid size received from the main program
!   ierr: error indicator
    integer, intent(in)::n
    integer , intent(out)::ierr  
    real(DP), intent(in)::noise(10)
    real(DP), intent(out)::Neps(10)

    real(DP):: h
    integer :: i=0,m,z
    integer :: size_x_loc,x,y
    integer :: ierr_h
    real(DP), dimension(:,:), allocatable:: X_MESH,basin_eps,basin_org
    real(DP), dimension(:),allocatable:: x_loc
    real(DP):: count_b
    real(DP):: param_a(2)
!    real(DP):: basin_org(n*n,1)
!    real(DP)::x_loc(n)
    real(DP):: maxiter
    real(DP)::lockout
    logical :: pos,neg

    allocate(X_MESH(n*n,2),basin_org(n*n,1),x_loc(n))
    print *,"Inside subroutine call \n"
    param_a(1)=2.12
    param_a(2)=-0.3
    maxiter=100
    lockout=100.0
    
    write(*,*) "noise vector",noise
    
    h=(BOXMAX-BOXMIN)/(n-1)
    x_loc=BOXMIN+h*(/(i,i=0,n-1)/)
    size_x_loc=size(x_loc)
   

    do x=1,n
        X_MESH((x-1)*n+1:x*n,1)=x_loc(x)
        do y=1,n
            X_MESH(n*(x-1)+y,2)=x_loc(y)
        enddo
    enddo
   
    call henon_map(X_MESH,n,maxiter,param_a,lockout,basin_org,0.d0,ierr_h)
    
!$OMP PARALLEL DO       
    do m=1,10
!        count_b=0
!        allocate(basin_pos(n*n,2))
        write(*,*) "here here :"
        call wrap_henon_map(X_MESH,n,maxiter,param_a,lockout,basin_org(:,1),noise(m),Neps,m,ierr_h)            
!        call henon_map(X_MESH,n,maxiter,param_a,lockout,basin_eps(:,2),-noise(m),ierr_h)
!        Neps(m)=count((basin_org(:,1).ne.basin_eps(:,1)).OR.(basin_org(:,1).ne.basin_eps(:,2)))
!        count_b=0
!        do z=1,n*n
!            pos=(basin_org(z,1)/=basin_eps(z,1))
!            neg=(basin_org(z,1)/=basin_eps(z,2))
!            if(pos.OR.neg) then
!                count_b=count_b+1.0
!            endif
!        end do
!    Neps(m)=count_b
    end do
!$OMP END PARALLEL DO

    ierr=0
    deallocate(X_MESH,basin_org,x_loc)!basin_eps_pos,basin_eps_neg)
    return
    end subroutine basin_alg
!-------------------------------------------------------------------------------
    end module henondim

