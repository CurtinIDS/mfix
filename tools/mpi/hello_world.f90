!======================================================================
!  This is a simple "Hello world" Fortran example  
!  to verify the MPI installation
!  Compile with: 
!
!  mpif90 -o hello_world.exe hello_world.f90
!
!  Run on 4 processors with: 
!
!  mpirun -np 4 hello_world.exe
!
!  It should display something like this (the order may be different)
!  and will most likely vary everytime the program is executed :
!
! Hello world from rank            1  of            4
! Hello world from rank            2  of            4
! Hello world from rank            0  of            4
! Hello world from rank            3  of            4
!
! To change the number of processors, modify the number after -np
! i.e., to run with 8 processors:
!
!  mpirun -np 8 hello_world.exe
!
!======================================================================

   program hello
   include 'mpif.h'
   integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)
  
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
   print*, 'Hello world from rank ',rank,' of ',size
   call MPI_FINALIZE(ierror)
   end


