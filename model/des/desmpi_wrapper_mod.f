!------------------------------------------------------------------------
! Module           : desmpi_wrapper
! Purpose          : Contains wrapper class for mpi communications- send,recv
!                    and scatter, gather
! Author           : Pradeep.G
!------------------------------------------------------------------------
#include "version.inc"

module desmpi_wrapper

      use parallel_mpi
      use mpi_utility
      use compar

      interface des_mpi_irecv
         module procedure des_mpi_irecv_db
      end interface

      interface des_mpi_isend
         module procedure des_mpi_isend_db
      end interface

      interface des_mpi_scatterv
         module procedure des_mpi_scatterv_i
         module procedure des_mpi_scatterv_db
      end interface

      interface des_mpi_gatherv
         module procedure des_mpi_gatherv_i
         module procedure des_mpi_gatherv_db
      end interface

      contains

!------------------------------------------------------------------------
! Subroutine       : des_mpi_barrier
! Purpose          : Wrapper class for barrier
!
!------------------------------------------------------------------------
      subroutine des_mpi_barrier
      implicit none
! local variables
      character(len=80), parameter :: name = 'des_mpi_barrier'
#ifdef MPI
      integer lerr
      call mpi_barrier(mpi_comm_world, lerr)
      call mpi_check( name //':mpi_barrier ', lerr )
#endif
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_irecv_db
! Purpose          : Wrapper class for mpi_irecv
!------------------------------------------------------------------------
      subroutine des_mpi_irecv_db(precvbuf,precvcnt,ptoproc,ptag,preq,perr)
      implicit none
! dummy variables
      double precision, dimension(:) :: precvbuf
      integer :: precvcnt,ptoproc,ptag,preq,perr
#ifdef MPI
      call mpi_irecv(precvbuf,precvcnt,mpi_double_precision,ptoproc,ptag,mpi_comm_world,preq,perr)
#endif
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_isend_db
! Purpose          : Wrapper class for mpi_isend
!------------------------------------------------------------------------
      subroutine des_mpi_isend_db(psendbuf,psendcnt,ptoproc,ptag,preq,perr)
      implicit none
! dummy variables
      double precision, dimension(:) :: psendbuf
      integer :: psendcnt,ptoproc,ptag,preq,perr
#ifdef MPI
      call mpi_isend(psendbuf,psendcnt,mpi_double_precision,ptoproc,ptag,mpi_comm_world,preq,perr)
#endif
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_wait
! Purpose          : Wrapper class for mpi_wait
!------------------------------------------------------------------------
      subroutine des_mpi_wait(preq,perr)
      implicit none
! dummy variables
      integer :: preq,perr
#ifdef MPI
! local variables
      integer :: lmpi_status(mpi_status_size)
      call mpi_wait(preq,lmpi_status,perr)
#endif
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_scatterv_db
! Purpose          : Wrapper class for mpi_scatterv
!------------------------------------------------------------------------
      subroutine des_mpi_scatterv_db(prootbuf,pscattercnts,pdispls, &
                           pprocbuf,precvcnt,proot,perr )
      implicit none
! dummy variables
      double precision, dimension(:):: prootbuf,pprocbuf
      integer, dimension (:) :: pdispls,pscattercnts
      integer :: precvcnt,proot,perr

#ifdef MPI
      call mpi_scatterv(prootbuf,pscattercnts,pdispls,mpi_double_precision, &
                        pprocbuf,precvcnt,mpi_double_precision,proot,mpi_comm_world,perr )
#else
      pprocbuf = prootbuf
#endif
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_scatterv_i
! Purpose          : Wrapper class for mpi_scatterv
!------------------------------------------------------------------------
      subroutine des_mpi_scatterv_i(prootbuf,pscattercnts,pdispls, &
                           pprocbuf, precvcnt,proot,perr )
      implicit none
! dummy variables
      integer, dimension(:) :: prootbuf,pprocbuf
      integer, dimension (:) :: pdispls,pscattercnts
      integer :: precvcnt,proot,perr

#ifdef MPI
      call mpi_scatterv(prootbuf,pscattercnts,pdispls,mpi_integer, &
                        pprocbuf,precvcnt,mpi_integer,proot,mpi_comm_world,perr )
#else
      pprocbuf = prootbuf
#endif
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_gatherv_db
! Purpose          : Wrapper class for mpi_gatherv
!------------------------------------------------------------------------
      subroutine des_mpi_gatherv_db(psendbuf,psendcnt,precvbuf, &
                                    precvcnts, pdispls,proot,perr )
      implicit none
! dummy variables
      double precision, dimension(:) :: psendbuf,precvbuf
      integer, dimension (:) :: pdispls,precvcnts
      integer :: psendcnt,proot,perr

#ifdef MPI
      call mpi_gatherv(psendbuf,psendcnt,mpi_double_precision,precvbuf,precvcnts, &
                       pdispls,mpi_double_precision,proot,mpi_comm_world,perr )
#else
      precvbuf = psendbuf
#endif
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_gatherv_i
! Purpose          : Wrapper class for mpi_gatherv
!------------------------------------------------------------------------
      subroutine des_mpi_gatherv_i(psendbuf,psendcnt,precvbuf, &
                                    precvcnts, pdispls,proot,perr )
      implicit none
! dummy variables
      integer, dimension(:) :: psendbuf,precvbuf
      integer, dimension (:) :: pdispls,precvcnts
      integer :: psendcnt,proot,perr

#ifdef MPI
      call mpi_gatherv(psendbuf,psendcnt,mpi_integer,precvbuf,precvcnts, &
                       pdispls,mpi_integer,proot,mpi_comm_world,perr )
#else
      precvbuf = psendbuf
#endif
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_stop
! Purpose          : Wrapper class for mpi_abort
!------------------------------------------------------------------------
      subroutine des_mpi_stop(myid)

      USE funits
      implicit none

      integer, optional, intent(in) :: myid

#ifdef MPI

      INTEGER :: mylid, ERRORCODE

      if (.not. present(myid)) then
         mylid = myPE
      else
         mylid = myid
      endif

      write(*,100) mylid
      write(UNIT_LOG,100) mylid

 100  format(/,'*****************',&
      '********************************************',/, &
      '(PE ',I2,') : A fatal error occurred in des routines',/,9X, &
      '*.LOG file may contain other error messages ',/,'*****************', &
      '********************************************',/)

!     call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_ABORT(MPI_COMM_WORLD, ERRORCODE, mpierr)
      write(*,"('(PE ',I2,') : MPI_ABORT return = ',I2)") &
      mylid,mpierr

      call MPI_Finalize(mpierr)

      ERROR_STOP 'MPI terminated from des_mpi_stop'
#else
      ERROR_STOP 'terminated from des_mpi_stop'
#endif
      end subroutine

      end module
