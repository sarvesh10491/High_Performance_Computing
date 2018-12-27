!  JACOBIMPI - Example program for the MPI version of the Jacobi Iteration calculations.  
!
!  DESCRIPTION
!   MPI based program to calculate Jacobi iteration for Poisson’s equation on the unit square. 
!   No. of processors are set as 16 in 4x4 grid formation and based on that
!   partitions are created for each processor on the grid. Each process
!   calculates the partitioned data based on its rank thus calulates the
!   true solution for equation. MPI_REDUCE call then adds all the
!   delta-c value array local to each processor on the root
!   procesor with rank(procnum)=0   
!   A row/column padding is used on each border of grid which serve as ghost points in calculations to ease the computation complexity
!   Processor ranking is done in 4x4 grid in row-wise manner starting from bottom left to end at top right fashion
!   For each block addition row/column buffer is added to each blocks compute matrix to use it when doin message passing across boundary
!
!  REVISION HISTORY
!   12/6/2018
!
!  PROGRAMMER
!  Sarvesh Patil, sspati13@asu.edu
!  ASU ID : 1213353386

!############################################################################

use, intrinsic::iso_fortran_env  ! Fortran 2003 and later
use mpi

implicit none
integer, parameter:: ROOT=0
integer, parameter:: M=255  ! Upper limit for datapoints
integer, parameter:: INTERNAL_ERROR=1  ! shell error code; must be nonzero
integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision

!  Local variables
integer:: procnum,init_x,init_y,bound,ierr,npes,Maxiter
real(DOUBLE)::fmat(256,256) ! f 
Maxiter=5000    ! iteration limitvalue matrix
bound=0

call mpi_init(ierr)  ! should be first executable statement
call mpi_comm_rank(MPI_COMM_WORLD, procnum, ierr)  ! processor rank
call mpi_comm_size(MPI_COMM_WORLD, npes, ierr) ! #processes


! Evaluate function values for root
if(procnum.eq.0) then
    call Eval_f(M, fmat)
end if

call mpi_bcast(fmat, size(fmat), MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)

call create_grid(procnum, init_x, init_y)   ! Function to allocate index for each processor block in grid

call Eval_jacobi(procnum,init_x,init_y,bound,M,Maxiter,fmat)    ! Function to evaluate Jacobi iterations

call mpi_finalize(ierr)

stop

911 continue
    call mpi_abort(MPI_COMM_WORLD, INTERNAL_ERROR, ierr)
end



!################################################################

subroutine create_grid(procnum, init_x, init_y)
    use, intrinsic::iso_fortran_env
    implicit none

    integer, intent(in):: procnum
    integer, intent(inout):: init_x
    integer, intent(inout):: init_y

    ! 1st row
    if(procnum.le.3) then
        init_y=1;

        if(procnum.eq.0) then
            init_x=1;
        endif
        if(procnum.eq.1) then
            init_x=65
        endif
        if(procnum.eq.2) then
            init_x=129
        endif
        if(procnum.eq.3) then
            init_x=193
        endif

    ! 2nd row
    else if(procnum.le.7) then
        init_y=65;

        if(procnum.eq.4) then
            init_x=1;
        endif
        if(procnum.eq.5) then
            init_x=65
        endif
        if(procnum.eq.6) then
            init_x=129
        endif
        if(procnum.eq.7) then
            init_x=193
        endif

    ! 3rd row
    else if(procnum.le.11) then
        init_y=129;

        if(procnum.eq.8) then
            init_x=1;
        endif
        if(procnum.eq.9) then
            init_x=65
        endif
        if(procnum.eq.10) then
            init_x=129
        endif
        if(procnum.eq.11) then
            init_x=193
        endif

    ! 4th row
    else
        init_y=193;

        if(procnum.eq.12) then
        init_x=1;
        endif
        if(procnum.eq.13) then
        init_x=65
        endif
        if(procnum.eq.14) then
        init_x=129
        endif
        if(procnum.eq.15) then
        init_x=193
        endif

    endif    

    return
end subroutine create_grid

!----------------------------------------------------------------------------

subroutine Eval_f(m, cur_f)
    use, intrinsic::iso_fortran_env
    implicit none

    integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision
    integer, intent(in):: m
    real(DOUBLE), intent(out)::cur_f(m+1,m+1)   ! matrix to store current f function values
    real(DOUBLE)::x
    real(DOUBLE)::y
    integer::j
    integer::k

    cur_f = 0
    do j=2,m
        x=(j-1)
        x=x/m
        do k=2,m
            y=(k-1)
            y=y/m
            ! f(x,y) = (e^(x+y))[(x^2 + 3x)(y^2 - y) + (y^2 + 3y)(x^2 - x)]
            cur_f(j,k)=((exp(x+y)) * ((((x**2)+3*x) * ((y**2)-y)) + (((y**2)+3*y) * ((x**2)-x))))
        enddo
    enddo

    return
end subroutine Eval_f

!----------------------------------------------------------------------------

subroutine update_gridborder(procnum, old_x, bound)
    use, intrinsic::iso_fortran_env
    implicit none

    integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision
    integer, intent(in):: procnum, bound
    real(DOUBLE), intent(inout)::old_x(66,66)

    ! 1st row
    if(procnum.lt.4) then
        old_x(:,2)=bound
    endif

    ! last row
    if(procnum.gt.11) then
        old_x(:,65)=bound
    endif

    ! 1st column
    if((procnum.eq.0).OR.(procnum.eq.4).OR.(procnum.eq.8).OR.(procnum.eq.12)) then
        old_x(2,:)=bound
    endif

    ! last column
    if((procnum.eq.3).OR.(procnum.eq.7).OR.(procnum.eq.11).OR.(procnum.eq.15)) then
    old_x(65,:)=bound
    endif

    return
end subroutine update_gridborder

!----------------------------------------------------------------------------

subroutine set_looper(procnum, init_x, init_y, final_x, final_y)
    use, intrinsic::iso_fortran_env
    implicit none

    integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision
    integer, intent(in):: procnum
    integer, intent(inout)::init_x, init_y, final_x, final_y

    ! 1st row i.e. bottom border of grid
    if((procnum.eq.0).OR.(procnum.eq.1).OR.(procnum.eq.2).OR.(procnum.eq.3)) then
        init_y=3

        if(procnum.eq.0) then
            init_x=3
            final_x=65
            final_y=65
        else
            init_x=2
            if(procnum.eq.3) then
                final_x=64
                final_y=65
            else
                final_y=65
                final_x=65
            endif
        endif

    ! Left border of grid excluding 0th processor
    else if((procnum.eq.4).OR.(procnum.eq.8).OR.(procnum.eq.12)) then
        init_x=3
        init_y=2
        final_x=65

        if(procnum.eq.12) then
            final_y=64
        else
            final_y=65
        endif

    ! Top border of grid excluding 12th processor
    else if((procnum.eq.13).OR.(procnum.eq.14).OR.(procnum.eq.15)) then
        init_y=2
        init_x=2
        final_y=64

        if(procnum.eq.15) then
            final_x=64
        else
            final_x=65
        endif

    ! Right border of grid excluding 3rd & 15th processors
    else if((procnum.eq.11).OR.(procnum.eq.7)) then
        init_y=2
        init_x=2
        final_x=64
        final_y=65

    else
        init_y=2
        init_x=2
        final_x=65
        final_y=65
    endif

    return
end subroutine set_looper

!----------------------------------------------------------------------------


subroutine Eval_trueU(j_init, k_init, cur_u, bound)
    use, intrinsic::iso_fortran_env
    implicit none

    integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision
    integer, intent(in)::bound
    integer, intent(in):: j_init
    integer, intent(in):: k_init
    real(DOUBLE), intent(out)::cur_u(64,64) ! matrix for storing currently computed U values
    real(DOUBLE)::x
    real(DOUBLE)::y
    integer::j
    integer::k

    cur_u = bound

    do j=1,64
        x=(j_init+j-2)
        x=x/255
        if((x.gt.0).AND.(x.lt.1)) then
            do k=1,64
            y=(k_init+k-2)
            y=y/255
                if((y.gt.0).AND.(y.lt.1)) then
                    ! u(x,y) = (e^(x+y))*(x^2 - x)*(y^2 - y)
                    cur_u(j,k) = ((exp(x+y)) * ((x**2)-x) * ((y**2)-y))
                endif
            enddo
        endif
    enddo

    return
end subroutine Eval_trueU

!----------------------------------------------------------------------------

subroutine Eval_deltaC(temp_u, new_u, true_u, del_c)
    use, intrinsic::iso_fortran_env
    implicit none

    integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision
    real(DOUBLE),intent(in):: new_u(66,66)
    real(DOUBLE),intent(in):: temp_u(66,66)
    real(DOUBLE),intent(in):: true_u(64,64)
    real(DOUBLE),intent(out)::del_c(2)
    real(DOUBLE)::absval
    real(DOUBLE)::c
    real(DOUBLE)::delta

    integer::j
    integer::k
    absval = 0
    c = 0
    delta = 0

    ! Eval delta value
    do j=2,65
        do k=2,65
            absval = sqrt((new_u(j,k) - temp_u(j,k))**2)
            if(absval.gt.delta) then
                delta = absval
            endif
        enddo
    enddo

    ! Eval C value
    do j=2,65
        do k=2,65
            absval = sqrt((new_u(j,k) - true_u((j-1),(k-1)))**2)
            if(absval.gt.c) then
                c = absval
            endif
        enddo
    enddo

    del_c(1) = delta
    del_c(2) = c

    return
end subroutine Eval_deltaC

!----------------------------------------------------------------------------

subroutine Eval_jacobi(procnum, begX, begY, bound, M, Maxiter, fmat)
    use mpi
    integer, parameter:: DOUBLE=kind(1.0d0)

    integer, intent(in)::procnum, begX, begY, bound, M, Maxiter
    real(DOUBLE), intent(in)::fmat(256,256)

    integer:: status(MPI_STATUS_SIZE)
    integer::init_x,init_y,final_x,final_y,x,y,n, ierr
    integer::send_to_north,rec_from_south,send_to_south,rec_from_north,send_to_east,rec_from_west,send_to_west,rec_from_east ! Tags for direction

    ! biffers to send over grid block boundaries
    real(DOUBLE):: send_NS_buff(66)
    real(DOUBLE):: rec_NS_buff(66)
    real(DOUBLE):: send_NS_buff1(66)
    real(DOUBLE):: rec_NS_buff1(66)
    
    real(DOUBLE):: del_c(2),del_c_red(2)    ! vector to hold delta & C values

    real(DOUBLE), dimension(:,:), allocatable:: new_x, old_x, temp_x, true_u
    allocate(new_x(66,66),old_x(66,66),temp_x(66,66),true_u(64,64))


    ! Assigning tags for msg passing
    send_to_north = 0
    send_to_south = 1
    send_to_east = 2
    send_to_west = 3

    rec_from_south = 0
    rec_from_north = 1
    rec_from_west = 2
    rec_from_east = 3

    old_x = 0
    call update_gridborder(procnum, old_x, bound)   ! Update the grid boundary processors
    new_x = old_x
    temp_x = old_x

    call set_looper(procnum, init_x, init_y, final_x, final_y)  ! Initialize loop indices

    call Eval_trueU(begX,begY,true_u,bound) ! Evaluate true U values for the iteration

    ! Do the computations in each block of grid & pass results
    do n=1,Maxiter
        do x=init_x,final_y
            do y=init_y,final_y
                new_x(x,y) = (0.25 * ( old_x(x+1,y) + old_x(x-1,y) + old_x(x,y+1) + old_x(x,y-1) - ( (fmat(begX+x-init_x, begY+y-init_y) )/( M**2 ) ) ))
            enddo
        enddo

        if(mod(n,10).eq.0) then
            call Eval_deltaC(temp_x,new_x,true_u,del_c)
            call mpi_barrier(MPI_COMM_WORLD, ierr)  ! wait for everyone to finish
            call mpi_allreduce(del_c, del_c_red, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)
            
            if(del_c_red(1).lt.0.0001) then
                goto 108
            else
                temp_x=new_x
            endif
        endif


        if(procnum.lt.12) then
            call MPI_Sendrecv( send_NS_buff, 66, MPI_DOUBLE,(procnum+4),send_to_north, rec_NS_buff,66, MPI_DOUBLE,(procnum+4),rec_from_north, MPI_COMM_WORLD, status, ierr)
            new_x(:,1)=rec_NS_buff
        endif

        if(procnum.gt.3) then
            send_NS_buff1=new_x(:,2)
            call MPI_Sendrecv( send_NS_buff1, 66,MPI_DOUBLE,(procnum-4),send_to_south, rec_NS_buff1,66,MPI_DOUBLE,(procnum-4), rec_from_south, MPI_COMM_WORLD, status, ierr)
            new_x(:,66)=rec_NS_buff1
        endif

        old_x=new_x

    enddo

    108 continue
    if(procnum.eq.ROOT) then
        write(*,*) "********************************************"
        write(*,*) "Breaking out after ", n," iterations."
        write(*,*) "Delta value : ", del_c_red(1)
        write(*,*) "C value : ", del_c_red(2)
        write(*,*) "********************************************"
    endif
    
    return

end subroutine Eval_jacobi


