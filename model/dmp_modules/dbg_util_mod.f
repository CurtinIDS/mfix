!       Debug printout utilities written by Aytekin Gel from Aeolus Research Inc.
!       following the convention in ORNL's mpi_utility template

    module dbg_util

    use compar
    use geometry
    use parallel_mpi
    use indices

    implicit none

!       Object-oriented approach to direct to the correct procedure
!       depending on the argument type. i stands for integer, r for real
!       and d for double precision. 0 for scalar, 1 for vector, 2 for
!       2-D array and similarly 3.

      interface dbgprn
        module  procedure  dbgprn_0i,dbgprn_1i, &
                           dbgprn_0r,dbgprn_1r
!               , , dbgprn_2i, dbgprn_3i, &
!                          dbgprn_0r, dbgprn_1r, dbgprn_2r, dbgprn_3r, &
!                          bcast_0d, bcast_1d, bcast_2d, bcast_3d
      end interface

      interface prnfield
         module  procedure prnfield_1d,prnfield_1r,prnfield_2d
      end interface

      contains

      subroutine dbgprn_0i( buffer, ncount, myid )
      integer, intent(inout) :: buffer
      integer, intent(in) :: ncount
      integer, optional, intent(in) :: myid
      integer :: dbgout = 75

      open(convert='big_endian',unit=dbgout, FILE='dbg'//fbname//'.dat', STATUS='UNKNOWN')
      write(dbgout,"('(PE ',I3,') :')") myPE
      write(dbgout,"(10X,'buffer = ',I6)") buffer
      close(dbgout)
      return
      end subroutine dbgprn_0i

      subroutine dbgprn_1i( buffer, ncount, myid )
      integer, intent(inout), dimension(:) :: buffer
      integer, intent(in) :: ncount
      integer, optional, intent(in) :: myid
      integer :: dbgout = 75
      integer :: i

      open(convert='big_endian',unit=dbgout, FILE='dbg'//fbname//'.dat', STATUS='UNKNOWN')
      write(dbgout,"('(PE ',I3,') :')") myPE
      do i=1,ncount
        write(dbgout,"(10X,'buf(',I3,')= ',I6)") i,buffer(i)
      end do
      close(dbgout)
      return
      end subroutine dbgprn_1i

      subroutine dbgprn_0r( buffer, ncount, myid )
      real, intent(inout) :: buffer
      integer, intent(in) :: ncount
      integer, optional, intent(in) :: myid
      integer :: dbgout = 75

      open(convert='big_endian',unit=dbgout, FILE='dbg'//fbname//'.dat', STATUS='UNKNOWN')
      write(dbgout,"('(PE ',I3,') :')") myPE
      write(dbgout,"(10X,'buffer = ',E14.6)") buffer
      close(dbgout)
      return
      end subroutine dbgprn_0r

      subroutine dbgprn_1r( buffer, ncount, myid )
      real, intent(inout), dimension(:) :: buffer
      integer, intent(in) :: ncount
      integer, optional, intent(in) :: myid
      integer :: dbgout = 75
      integer :: i

      open(convert='big_endian',unit=dbgout, FILE='dbg'//fbname//'.dat', STATUS='UNKNOWN')
      write(dbgout,"('(PE ',I3,') :')") myPE
      do i=1,ncount
        write(dbgout,"(10X,'buf(',I3,')= ',E14.6)") i,buffer(i)
      end do
      close(dbgout)
      return
      end subroutine dbgprn_1r

      subroutine prnfield_1d (gbuf,varname,flagl)

      use functions
      implicit none

      double precision, intent(in), dimension(:) :: gbuf
      character(len=3), intent(in)   :: flagl
      character(len=*), intent(in)   :: varname
      integer :: ldbg = 35
      integer :: i,j,k
!      integer, optional, intent(in) :: mroot, idebug

       OPEN(CONVERT='BIG_ENDIAN',unit=ldbg,file=flagl//fbname//'.LOG',status='UNKNOWN')
       write(ldbg,"('Dumping variable : ',A10)") varname
       DO K = kstart3, kend3                               !//AIKEPARDBG
         write(ldbg,"('K = ',I5)") K                !//AIKEPARDBG
         write(ldbg,"(12X,14(I3,11X))") (I,i=Istart3,Iend3)  !//AIKEPARDBG
          DO J = jstart3, Jend3                            !//AIKEPARDBG
            write(ldbg,"(I3,')')",ADVANCE="NO") J               !//AIKEPARDBG
            DO I = istart3, Iend3                          !//AIKEPARDBG
              write(ldbg,"(2X,E12.4)",ADVANCE="NO") gbuf(FUNIJK(I,J,K)) !//AIKEPARDBG
            END DO                                       !//AIKEPARDBG
            write(ldbg,"(/)")                        !//AIKEPARDBG
          END DO                                         !//AIKEPARDBG
       END DO                                            !//AIKEPARDBG
       close(35)
      end subroutine prnfield_1d


      subroutine prnfield_1r (gbuf,varname,flagl)

      use functions
      implicit none

      real, intent(in), dimension(:) :: gbuf
      character(len=3), intent(in)   :: flagl
      character(len=*), intent(in)   :: varname
      integer :: ldbg = 35
      integer :: i,j,k
!      integer, optional, intent(in) :: mroot, idebug

       OPEN(CONVERT='BIG_ENDIAN',unit=ldbg,file=flagl//fbname//'.LOG',status='UNKNOWN')
       write(ldbg,"('Dumping variable : ',A10)") varname
       DO K = kstart3, kend3                               !//AIKEPARDBG
         write(ldbg,"('K = ',I5)") K                !//AIKEPARDBG
         write(ldbg,"(12X,14(I3,11X))") (I,i=Istart3,Iend3)  !//AIKEPARDBG
          DO J = jstart3, Jend3                            !//AIKEPARDBG
            write(ldbg,"(I3,')')",ADVANCE="NO") J               !//AIKEPARDBG
            DO I = istart3, Iend3                          !//AIKEPARDBG
              write(ldbg,"(2X,E12.4)",ADVANCE="NO") gbuf(FUNIJK(I,J,K)) !//AIKEPARDBG
            END DO                                       !//AIKEPARDBG
            write(ldbg,"(/)")                        !//AIKEPARDBG
          END DO                                         !//AIKEPARDBG
       END DO                                            !//AIKEPARDBG
       close(35)
      end subroutine prnfield_1r

      subroutine prnfield_2d (gbuf,varname,flagl)

      use functions
      implicit none

      double precision, intent(in), dimension(:,:) :: gbuf
      character(len=3), intent(in)   :: flagl
      character(len=*), intent(in)   :: varname
      integer :: ldbg = 35
      integer :: i,j,k
!      integer, optional, intent(in) :: mroot, idebug

       OPEN(CONVERT='BIG_ENDIAN',unit=ldbg,file=flagl//fbname//'.LOG',status='UNKNOWN')
       write(ldbg,"('Dumping variable : ',A10)") varname
       DO K = kstart3, kend3                               !//AIKEPARDBG
         write(ldbg,"('K = ',I5)") K                !//AIKEPARDBG
         write(ldbg,"(12X,14(I3,11X))") (I,i=Istart3,Iend3)  !//AIKEPARDBG
          DO J = jstart3, Jend3                            !//AIKEPARDBG
            write(ldbg,"(I3,')')",ADVANCE="NO") J               !//AIKEPARDBG
            DO I = istart3, Iend3                          !//AIKEPARDBG
              write(ldbg,"(2X,E12.4)",ADVANCE="NO") gbuf(FUNIJK(I,J,K),1) !//AIKEPARDBG
            END DO                                       !//AIKEPARDBG
            write(ldbg,"(/)")                        !//AIKEPARDBG
          END DO                                         !//AIKEPARDBG
       END DO                                            !//AIKEPARDBG
       close(35)
      end subroutine prnfield_2d

    end module dbg_util
