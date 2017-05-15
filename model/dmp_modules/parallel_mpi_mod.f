        module parallel_mpi

!       A module to carry out init, finalize and check for any parallel errors

        use geometry
        use compar
        implicit none

        contains

        subroutine parallel_init()
          implicit none
          integer :: ierr

          numPEs = 1
          myPE = 0

#ifdef MPI
          call MPI_Init(ierr)
          call MPI_Check( 'parallel_init:MPI_Init ', ierr)

          call MPI_COMM_SIZE( MPI_COMM_WORLD, numPEs, ierr )
          call MPI_Check( 'parallel_init:MPI_Comm_size ', ierr )

          call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
          call MPI_Check( 'parallel_init:MPI_Comm_size ', ierr )
#endif
          return
        end subroutine parallel_init

        subroutine parallel_fin()
          implicit none

#ifdef MPI
          integer :: ierr

          call MPI_Finalize(ierr)
          call MPI_Check( 'parallel_init:MPI_Finalize ', ierr)
#endif

          return
        end subroutine parallel_fin

        subroutine MPI_Check( msg, ierr )
#ifdef MPI
          use mpi, only: MPI_SUCCESS ! ignore-depcomp
#endif
          implicit none
          character(len=*),intent(in) :: msg
          integer, intent(in) :: ierr

#ifdef MPI
          character(len=512) :: errmsg
          integer :: resultlen, ierror

          if (ierr .ne. MPI_SUCCESS ) then
             call MPI_Error_string( ierr, errmsg, resultlen, ierror )
             print*, 'Error: ', msg
             print*, errmsg(1:resultlen)
             stop '** ERROR ** '
          endif
#endif

          return
        end subroutine MPI_Check


        end module parallel_mpi

