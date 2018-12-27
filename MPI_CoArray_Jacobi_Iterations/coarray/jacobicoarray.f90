!  JACOBICOARRAY - Example program for the COARRAY version of the Jacobi Iteration calculations.  
!
!  DESCRIPTION
!   COARRAY based program to calculate Jacobi iteration for Poisson’s equation on the unit square. 
!   No. of processors are set as 16 in 4x4 grid formation in ribbon sets formation and based on that
!   partitions are created for each processor on the grid. 
!   Each image calculates the partitioned data based on its image number thus calulates the
!   true solution for equation.
!
!  REVISION HISTORY
!   12/6/2018
!
!  PROGRAMMER
!  Sarvesh Patil, sspati13@asu.edu
!  ASU ID : 1213353386

!############################################################################

program jacobicoarray
    use, intrinsic::iso_fortran_env  ! Fortran 2003 and later
    implicit none
    
    real u(0:255,0:17,0:1)[*]   ! Creating a coarray for 16 processors with 3rd dimesion to be used for storing two compute values
    real fmat(0:255,0:255),trueU(0:255,0:255) ! Grid for computing f & true U values
    real diff(16)   ! vector to store differences for delta & c calculations
    integer num_imgs, cur_img, old, new, Maxiter, bound, m, k, p
    real x,y,i,j,c,delta

    m=255   ! Upper limit of datapoints
    Maxiter=5000    ! Iterations limit
    bound=0 ! boundary
    u=0
    num_imgs = NUM_IMAGES() ! To get total no. of processor images
    cur_img = THIS_IMAGE()  ! to get current image number


    !write(*,*) "Total images",num_imgs
    !write(*,*) "Current image",cur_img

    !=======================================================
    ! Eval True solution of equation
    ! u(x,y) = (e^(x+y))*(x^2 - x)*(y^2 - y)
    trueU = bound
    call Eval_trueU(trueU, m)
    sync all
    !write(*,*),trueU


    ! Eval f matrix
    ! f(x,y) = (e^(x+y))[(x^2 + 3x)(y^2 - y) + (y^2 + 3y)(x^2 - x)]
    fmat = bound
    call Eval_f(fmat, m)
    

    !========================================================
    ! Update 1st and last column in grid on end processors
    if(cur_img==1) then
        u(:,1,:)=bound
    end if
    if(cur_img==16) then
        u(:,16,:)=bound
    end if

    ! Update 1st and last rows in grid for all processors
    u(0,:,:)=bound
    u(255,:,:)=bound

    !=========================================================
    new = 1
    old = 1-new
    ! Loop to evaluate u values & update coarray
    do k=1, Maxiter
        ! For all images except 1st image update copy
        if (cur_img.gt.1) then
            u(:,17,old)[cur_img-1] = u(:,1,old)
        endif
        ! For all images except last image update copy
        if (cur_img.lt.16) then
             u(:,0,old)[cur_img+1] = u(:,16,old)
        endif
        sync all

        ! Update u value on all images
        do j=1,16
            if(j.ne.cur_img) then
                do i=1, 254
                    u(i,j,new) = 0.25 * ( u(i+1,j,old) + u(i-1,j,old) + u(i,j+1,old) + u(i,j-1,old) - ( (1/(m**2)) * fmat(i,((cur_img-1)*16)+j-1) ) )
                enddo
            endif
        enddo


        ! Eval delta on each of 10th iteration
        if(mod(k,10).eq.0) then
            ! 1st image
            if(cur_img.eq.1) then
                diff(1)=maxval(SQRT(((u(1:254,2:16,1)-u(1:254,2:16,0)) ** 2)))
            end if

            ! All interior images
            do p=2,16
                if(cur_img.eq.p) then
                    diff(p)=maxval(SQRT((u(1:254,1:16,1)-u(1:254,1:16,0)) ** 2))
                end if
            end do

            ! Last image
            if(cur_img.eq.16) then
                diff(16)=maxval(SQRT((u(1:254,1:15,1)-u(1:254,1:15,0)) ** 2))
            end if

            if(delta.lt.maxval(diff)) then
                delta=maxval(diff)
            endif
        endif


        ! Copy old value to new value for the next iteration
        new = old
        old = 1-new
     enddo


    !=========================================================================

    ! Eval c values
    ! 1st image
    if(cur_img.eq.1) then
        diff(1)=maxval(SQRT(((u(1:254,1:16,1)-trueU(1:254,0:15)) ** 2)))
    end if

    ! All interior images
    do k=2,16
        if(cur_img.eq.k) then
            diff(k)=maxval(SQRT(((u(1:254,1:16,1)-trueU(1:254,(k-1)*16:((k-1)*16)+15)) **2)))
        end if
    end do

    ! Last image
    if(cur_img.eq.16) then
        diff(16)=maxval(SQRT(((u(1:254,1:16,1)-trueU(1:254,240:255)) ** 2)))
    end if

    c = maxval(diff)



    if(cur_img.eq.1) then
        write(*,*) "********************************************"
        write(*,*) "Delta value : ", delta
        write(*,*) "C value : ", c
        write(*,*) "********************************************"
    end if

end program jacobicoarray

!##########################################################################3

subroutine Eval_trueU(trueU, m)
    use, intrinsic::iso_fortran_env
    implicit none

    real, intent(inout)::trueU(0:255,0:255)
    integer, intent(in):: m
    real x,y,i,j
    
    do i=1,m-1
        x=i/m
        do j=1,m-1
            y=j/m
            trueU(i,j)=((exp(x+y))*((x**2)-x)*((y**2)-y))
        enddo
    enddo

    return
end subroutine Eval_trueU

!----------------------------------------------------------------------------

subroutine Eval_f(fmat, m)
    use, intrinsic::iso_fortran_env
    implicit none

    real, intent(inout)::fmat(0:255,0:255)
    integer, intent(in):: m
    real x,y,i,j
    
    do i=1,m-1
        x=i/m
        do j=1,m-1
            y=j/m
            fmat(i,j)=((exp(x+y))*((((x**2)+3*x)*((y**2)-y))+(((y**2)+3*y)*((x**2)-x))))
        enddo
    enddo

    return
end subroutine Eval_f