        module debug
!//BUGFIX 0904 added funits module here for declaration of UNIT_LOG
        USE funits
!       USE dbg_util
        implicit none

        integer :: idebug = 0

        interface assert
        module procedure assert_i, assert_d, assert_i2, assert_d2
        end interface
!//BUG 0904 unit_log declaration shouldn't be here, it causes conflict when
!//BUG both MPI_UTILITY (via DEBUG.mod) and FUNITS are USEd in same routine
!//BUG  integer :: unit_log = 13
        interface write_debug
        module procedure write_debug_0i, write_debug_0d, write_debug_0, &
        write_debug_1i, write_debug_1d, write_debug_0l
        end interface


        contains

!       Some debugging means

        subroutine debug_init(myPE)
        integer, intent(in) :: myPE

        character(LEN=255) :: filename

        write(filename,'("debug",f6.3)') dble(myPE)/dble(1000)
        print*,'filename ', filename

        open(unit_log, file=filename, access='sequential',form='formatted',convert='big_endian')
        rewind(unit_log)

        return
        end subroutine debug_init


        subroutine assert_i( lcond, msg, value )
        logical, intent(in) :: lcond
        character(len=*), intent(in) :: msg
        integer, intent(in) :: value

        if (.not. lcond) then
           print*,'Assertion error: ', msg, value
           stop '** ERROR ** '
        endif

        return
        end subroutine assert_i


        subroutine assert_i2( lcond, msg, value, value2 )
        logical, intent(in) :: lcond
        character(len=*), intent(in) :: msg
        integer, intent(in) :: value, value2

        if (.not. lcond) then
           print*,'Assertion error: ', msg, value, value2
           stop '** ERROR ** '
        endif

        return
        end subroutine assert_i2




        subroutine assert_d( lcond, msg, value )
        logical, intent(in) :: lcond
        character(len=*), intent(in) :: msg
        double precision, intent(in) :: value

        if (.not. lcond) then
           print*,'Assertion error: ', msg, value
           stop '** ERROR ** '
        endif

        return
        end subroutine assert_d

        subroutine assert_d2( lcond, msg, value, value2 )
        logical, intent(in) :: lcond
        character(len=*), intent(in) :: msg
        double precision, intent(in) :: value, value2

        if (.not. lcond) then
           print*,'Assertion error: ', msg, value, value2
           stop '** ERROR ** '
        endif

        return
        end subroutine assert_d2


        subroutine write_debug_0( name, msg )
        character(len=*), intent(in) :: name, msg

        character(len=80) :: line(1)

        line(1) = msg
        call write_error( name, line, 1 )

        return
        end subroutine write_debug_0

        subroutine write_debug_1i( name, msg, x )
        character(len=*), intent(in) :: name, msg
        integer, intent(in), dimension(:) :: x

!       ---------------
!       local variables
!       ---------------
        character(len=80) :: line(1+size(x))
        integer :: i, ip

        line(1) = " "
        line(1) = msg

        ip = 2
        do i=lbound(x,1),ubound(x,1)
          line(ip) = " "
          write(line(ip), 9001) i, x(i)
 9001     format('i = ', i7,' value = ', i9 )

          ip = ip + 1
        enddo

        call  write_error( name, line, 1+size(x) )
        return
        end subroutine write_debug_1i


        subroutine write_debug_1d( name, msg, x )
        character(len=*), intent(in) :: name, msg
        double precision, intent(in), dimension(:) :: x

!       ---------------
!       local variables
!       ---------------
        character(len=80) :: line(1+size(x))
        integer :: i, ip

        line(1) = " "
        line(1) = msg

        ip = 2
        do i=lbound(x,1),ubound(x,1)
          line(ip) = " "
          write(line(ip), 9001) i, x(i)
 9001     format('i = ', i7,' value = ', 1pd30.10 )

          ip = ip + 1
        enddo

        call  write_error( name, line, 1+size(x) )
        return
        end subroutine write_debug_1d



        subroutine write_debug_0i(name, msg, x1, x2, x3, x4 )
        character(len=*), intent(in) :: name, msg
        integer, intent(in) :: x1
        integer, intent(in), optional :: x2,x3,x4

        character(len=80) :: line(1)
        integer :: narg

        narg = 1
        if (present(x2)) then
           narg = narg + 1
        endif
        if (present(x3)) then
           narg = narg + 1
        endif
        if (present(x4)) then
           narg = narg + 1
        endif

        select case ( narg )
        case (1)
           write(line(1),*) msg, x1
        case (2)
           write(line(1),*) msg, x1, x2
        case (3)
           write(line(1),*) msg, x1, x2, x3
        case (4)
           write(line(1),*) msg, x1, x2, x3, x4
        case default
           write(line(1),*) msg
        end select

        call write_error( name, line, 1 )

        return
        end subroutine write_debug_0i



        subroutine write_debug_0d(name, msg, x1, x2, x3, x4 )
        character(len=*), intent(in) :: name, msg
        double precision, intent(in) :: x1
        double precision, intent(in), optional :: x2,x3,x4

        character(len=80) :: line(1)
        integer :: narg

        narg = 1
        if (present(x2)) then
           narg = narg + 1
        endif
        if (present(x3)) then
           narg = narg + 1
        endif
        if (present(x4)) then
           narg = narg + 1
        endif

        select case ( narg )
        case (1)
           write(line(1),*) msg, x1
        case (2)
           write(line(1),*) msg, x1, x2
        case (3)
           write(line(1),*) msg, x1, x2, x3
        case (4)
           write(line(1),*) msg, x1, x2, x3, x4
        case default
           write(line(1),*) msg
        end select

        call write_error( name, line, 1 )

        return
        end subroutine write_debug_0d

        subroutine write_debug_0l(name, msg, x1, x2, x3, x4 )
        character(len=*), intent(in) :: name, msg
         logical, intent(in) :: x1
         logical, intent(in), optional :: x2,x3,x4

        character(len=80) :: line(1)
        integer :: narg

        narg = 1
        if (present(x2)) then
           narg = narg + 1
        endif
        if (present(x3)) then
           narg = narg + 1
        endif
        if (present(x4)) then
           narg = narg + 1
        endif

        select case ( narg )
        case (1)
           write(line(1),*) msg, x1
        case (2)
           write(line(1),*) msg, x1, x2
        case (3)
           write(line(1),*) msg, x1, x2, x3
        case (4)
           write(line(1),*) msg, x1, x2, x3, x4
        case default
           write(line(1),*) msg
        end select

        call write_error( name, line, 1 )

        return
        end subroutine write_debug_0l

!//     --------------------------------------
!//S    should be linked with mfix write_error
!//     --------------------------------------

        subroutine write_error( name, line, lmax )
        integer, intent(in) :: lmax
        character(len=*) name,line(*)

        integer :: L


      WRITE (UNIT_LOG, 1000) NAME
      DO L = 1, LMAX
         WRITE (UNIT_LOG, 1010) LINE(L)
      END DO
      WRITE (UNIT_LOG, 1020)
      flush(UNIT_LOG)

      RETURN
 1000 FORMAT(1X,70('*'),/,/,1X,'From : ',A)
 1010 FORMAT(1X,A)
 1020 FORMAT(/,/,1X,70('*'))
      END SUBROUTINE WRITE_ERROR



        end module debug

