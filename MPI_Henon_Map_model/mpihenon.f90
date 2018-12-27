!  MPIHENON - Example main program for the MPI version of the Henon basin
!  uncertainty calculation.  All values are taken from a file named
!  on the command line (or 'henon.txt' otherwise).
!
!  DESCRIPTION
!   MPI based program to calculate no. of uncertain points when noise is induced
!   into the data. No. of processors are read dynamocally and based on that
!   partitions are created for each processor on the grid data. Each process
!   calculates the partitioned data based on its rank thus calulates the total
!   no. of uncertain points into uncert array. MPI_REDUCE call then adds all the
!   uncert array local to each processor into uncertsum array on the root
!   procesor with rank(me)=0
!
!  PUBLIC ROUTINES DEFINED
!    Brief list of procedures that can be called from outside this
!    module.
!    1] basin_dim()
!
!  REQUIRED DEPENDENCIES
!   'precision' - to define KINDs for single and double floating-point.
!   And any other modules that you define.
!   'henondim' - to provide basin vector calculation mechanism and to provide
!   function call to calculate no. of uncertain points    
!
!  REVISION HISTORY
!   11/14/2018
!
!  PROGRAMMER
!  Sarvesh Patil, sspati13@asu.edu
!
!  Collaborator
!  Shubham Nandanwankar, snandanw@asu.edu

    use, intrinsic::iso_fortran_env  ! Fortran 2003 and later
    use mpi
    use henondim
    use precision
    ! use henondim or whatever module contains your actual Henon basin code
    implicit none
    integer, parameter:: ROOT=0
    integer, parameter:: NEPS=10
    integer, parameter:: INTERNAL_ERROR=1  ! shell error code; must be nonzero
    integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision

!  Local variables

    integer:: j, ierr, nargs
    real(DP),parameter:: HALF=0.5

!  These are the epsilon values requested in previous assignments, and you
!  will use the same values for this program.

    real(DP)::epsilon(NEPS)=[(HALF**j, j=12,21)]

!  Map parameters, as a 2-vector

    real(DP)::henonparams(2)  !  parameters a and b

!  Other variables

!  The first four values determine the X and Y limits of the grid,
!  and the second two values are the number of grid points.  XNUM and YNUM
!  are intended to be integer valued, but they are declared as double
!  precision so that they can be sent as a single message.  (As integer
!  values, XNUM and YNUM would have to be sent either as a separate
!  message, or all these values would have to be wrapped as a user-defined
!  type and registered with the MPI library.)  

    real(DP):: xextent(6)  ! xmin, xmax, ymin, ymax, xnum, ynum

!  Other variables

    character(200):: errmsg='none'  ! sufficiently long for typical messages
    character(80):: filename='henon.txt'  ! default input file name
    integer:: uncert(NEPS)  ! number of uncertain points for each epsilon
    integer:: uncertsum(NEPS)  ! total number of uncertain points
    integer:: me  ! my MPI rank (0 = root process)
    integer:: npes  ! Number of Processing ElementS (> 0), i.e., the
    namelist/basinparams/ henonparams, xextent

!  Initialize.

    call mpi_init(ierr)  ! should be first executable statement
    call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)  ! my rank
    call mpi_comm_size(MPI_COMM_WORLD, npes, ierr) ! number of processes

!  Start the Henon basin computation.

    uncert=0

!  Read data on the root process and broadcast.  We MUST NOT allow the Fortran
!  runtime to kill the root process if the input data file cannot be opened or
!  if there is an i/o error.  Otherwise, the root process exits without
!  calling MPI_FINALIZE, the other processes hang while waiting from data 
!  from the now-dead root process.
!  General rule: !  ALL i/o statements in an MPI Fortran program
!  should include the IOSTAT= and IOMSG= specifiers, and an appropriate
!  branch taken to display the error message and kill the entire MPI program.
!  Similarly, ALL exceptions in a C++ program must be caught in a global
!  catch{} block that takes analogous action.

    if(me.eq.ROOT) then
       nargs=command_argument_count()
       if(nargs.gt.0) call get_command_argument(1,filename)
       open(unit=4, file=filename, status='old', iostat=ierr, iomsg=errmsg)
       if(ierr.ne.0) goto 911
       read(4,basinparams,iostat=ierr,iomsg=errmsg)
       if(ierr.ne.0) goto 911
       close(4,iostat=ierr,iomsg=errmsg)  ! close the data file
       if(ierr.ne.0) goto 911
    endif

!  Broadcast the map parameters to everyone.
!  All processes must call MPI_BCAST to receive the broadcast messages.

    call mpi_bcast(henonparams, 2, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)
!    call mpi_bcast(n,1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)
!  Likewise, broadcast the grid extent and number of grid points in 
!  each direction (as floating-point numbers).

    call mpi_bcast(xextent, size(xextent), MPI_DOUBLE, ROOT, MPI_COMM_WORLD,ierr)

!  Divide up the grid according to the number of processes and calculate
!  the number of uncertain points for each epsilon on our portion of the grid.
!  Include whatever other parameters that you need in this call.
!  You should include appropriate code to handle unexpected errors by
!  writing a message to ERROR_UNIT and calling MPI_ABORT.

    call basin_dim(xextent,me,npes,henonparams,epsilon,uncert)

!  Sum up the total number of uncertain points from each process for each
!  epsilon and place the results on the root.  MPI computes a global sum
!  by adding corresponding elements of UNCERT from each process.
!  Other global operations include MPI_PROD (product), MPI_MAX (maximum of each
!  element), and MPI_MIN (minimum of each element).  It's also possible to
!  define your own, but we won't cover that topic here.

    call mpi_reduce(uncert, uncertsum, NEPS, MPI_INTEGER, MPI_SUM, ROOT,MPI_COMM_WORLD, ierr)

!  Let the root process report the results on the standard output.
!  UNCERTSUM is undefined on the other processes.

    if(me.eq.ROOT) then
       do j=1, NEPS
         write(OUTPUT_UNIT,110,iostat=ierr,iomsg=errmsg) epsilon(j),uncertsum(j)
       enddo
    endif
110 format(es13.6, 2x, i7)

!  Clean up and exit.

    call mpi_barrier(MPI_COMM_WORLD, ierr)  ! wait for everyone to finish
    call mpi_finalize(ierr)  ! should be last executable statement before STOP
    stop
!---------   End of normal computation

!  Exception handling:
!  Abort on i/o errors, because if we can't read the data on the root process, 
!  then we can't compute a useful result.
!  Advice: when composing error messages, it is helpful to include the rank
!  of the MPI process, so you will know who is complaining.
!  Always write error messages to the standard error unit, not to the standard
!  output, as MPI launchers usually create separate error log files.
!  In such cases, your error messages won't get mixed in with other
!  computational output.  Error log files should be empty upon successful
!  completion or (if you desire) contain only debugging messages.

911 continue
    if(me.eq.ROOT) write(ERROR_UNIT,999) me, errmsg
999 format('MPI process ', i0, ': ', a)
    call mpi_abort(MPI_COMM_WORLD, INTERNAL_ERROR, ierr)
    end

!---------CUT HERE---------------------------------------------------------
   subroutine basin_dim(xextent,me,npes,params,epsilon_b,uncert)
!  BASIN_DIM - dummy subroutine that should be deleted here and replaced with
!  a suitably modified routine in module henondim.
!  This version simply prints each process's copy of the inputs to
!  the standard output.  The outputs are likely to be jumbled when
!  run on multiple processors. This subroutine uses henondim module and
!  precision module from assignment 5(OpenMP), imported modules provide
!  functions to calculate function call to calculate uncertain no. of points and custom precision
!  floating point respectively. Each process then calculates stores results i.e.
!  no. of uncertain points in uncert array  and returns to the main script where
!  these values are summed up based on the epsilon values   

    use, intrinsic::iso_fortran_env
    use henondim
    use precision
    implicit none
    real(DP), intent(in)::xextent(6)
    real(DP),intent(in):: params(2)
    real(DP),intent(in)::epsilon_b(10)
    integer, intent(in)::me, npes
    integer, intent(inout)::uncert(10)
    
    integer::x,y,index_node,start_index,end_index,i,ierr_mpi
    real(DP)::h,maxiter_ext,lockout_ext
    real(DP),dimension(:,:),allocatable::X_MESH,basin_org
    real(DP),dimension(:),allocatable::x_loc
    integer::n,m
    
!   Maxiter, lockout and grid length assignment
    maxiter_ext=100.d0
    lockout_ext=100.d0
    n=int(xextent(5))

!   index_node is length of the vector each processor does operations on    
    index_node=(n*n)/npes

    allocate(X_MESH(n*n,2),x_loc(n),basin_org(index_node,1))

!   calculating n equally distributed points based on the values given in the
!   xextent array    
    h=(xextent(2)-xextent(1))/(n-1)
    x_loc=xextent(1)+h*(/(i,i=0,n-1)/)    

!   start and end index for each processor based on the rank value 
    start_index=(me*index_node)+1
    end_index=(me*index_node)+index_node
!    write(*,*) "xentent value is:",xextent,"n is:",n

!   grid calculation i.e. (n*n,2) array with each x,y coordinate pair    
    do x=1,n
        X_MESH((x-1)*n+1:x*n,1)=x_loc(x)
        do y=1,n
             X_MESH(n*(x-1)+y,2)=x_loc(y)
        enddo 
    enddo

!   function call to calculate original(without any noise) basin vector
    call henon_map(X_MESH(start_index:end_index,1:2),index_node,maxiter_ext,params,lockout_ext,basin_org,0.d0,ierr_mpi)
!    write(*,*) "size basin_org:",size(basin_org)

!   function call to calulate no. of uncertain points per procesor
    do m=1,10
        call wrap_henon_map(X_MESH(start_index:end_index,:),index_node,maxiter_ext,params,lockout_ext,basin_org(:,1),epsilon_b(m),uncert,m,ierr_mpi)
    enddo

    write(*,*) "Uncertain points for Processor:",me,", are:",uncert
    write(OUTPUT_UNIT, 100) me, npes
    write(OUTPUT_UNIT, 101) me, params
    write(OUTPUT_UNIT, 102) me, xextent
100 format('process ', i0, ': process count: ', i0)
101 format('process ', i0, ': parameters: ',2f8.2)
102 format('process ', i0, ': xextent: ', 6f8.2)

!  Place dummy values in the first 3 elements of UNCERT

!   uncert=[me, 2*me, 3*me]  ! depends on rank and number of processors
   deallocate(X_MESH,x_loc,basin_org) 
   return
   end subroutine basin_dim

