! Name : Sarvesh Patil
! ASU ID: 1213353386
! In collaboration with: Shubham Nandanwankar (1213350370)
! File : precision.f90
! DESCRIPTION : 
!  precision setting for flating point numbers
! REVISION DRAFT:
!  11/14/2018
! ====================================================================================================

   module precision
   integer, parameter:: SP=kind(1.0)  ! single precision
   integer, parameter:: DP=kind(1.0d0)  ! double precision
   end module precision

! ====================================================================================================
! End of precision.f90
! ====================================================================================================
