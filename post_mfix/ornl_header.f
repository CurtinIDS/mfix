      subroutine ornl_routines

      use geometry
      use run
      use usr_input

      implicit none

      integer :: selection , L, nt, nt2
      real    :: time_in_res

      real, allocatable :: arr(:)
      real(kind=8), allocatable :: time_series(:), time_series2(:)

      allocate (arr(ijkmax2))
      allocate (time_series(100000))
      allocate (time_series2(100000))
!
      call init_usr_input
!
 10   continue
      call header_zone(selection,0,2)
      if (selection .eq. 0) then
         deallocate(arr)
         deallocate(time_series)
         deallocate(time_series2)
         return
      end if
!
      call read_res1
      time_in_res = time

      DO L = 1, N_SPX
         REC_POINTER(L) = 4
         AT_EOF(L)      = .false.
      END DO
!
      if (selection .eq. 1) call do_zonal_avg(time_in_res,arr,time_series,nt)
      if (selection .eq. 2) then
         call do_zonal_avg(time_in_res,arr,time_series,nt)
         call do_zonal_avg(time_in_res,arr,time_series2,nt2)
      end if
!
      call header_low_pass(selection,0,2)
      if (selection .eq. 1) call do_low_pass(time_series,nt)
      if (selection .eq. 2) then
         call do_low_pass(time_series,nt)
         call do_low_pass(time_series2,nt2)
         call do_cross_corr(time_series, time_series2, nt)
      end if
!
      call header_daw(selection,0,5)
!
      if (selection .eq. 0) then
         deallocate(arr)
         deallocate(time_series)
         deallocate(time_series2)
         return
      end if
!
      if (selection .eq. 1) call do_simple_stats(time_series,nt)
      if (selection .eq. 2) call do_auto_corr(time_series,nt)
      if (selection .eq. 3) call do_power_spectrum(time_series,nt)
      if (selection .eq. 4) call do_t3_symm(time_series,nt)
      if (selection .eq. 5) call do_trevsgn(time_series,nt)


      goto 10
      deallocate(arr)
      deallocate(time_series)
      deallocate(time_series2)
      end
!
! *********************************************************************
! *************************** header_zone *****************************
! *********************************************************************
!
      SUBROUTINE header_zone(selection,min_sel,max_sel)
!
      IMPLICIT NONE
!
      integer :: selection , min_sel , max_sel

!             NUMBER OF INCORRECT INPUT SELECTIONS MADE
      INTEGER :: NERROR
!
      NERROR = -1
!
10    NERROR = NERROR + 1
      IF (NERROR.GT.10) THEN
         WRITE (*,*) ' HEADER_ORNL : TOO MANY INCORRECT INPUTS'
         selection = -1
         return
      END IF
      write (*,*) ' '
      WRITE (*,*)&
        ' *************************************************'
      WRITE (*,*)&
        '  0   - Exit ORNL routines'
      WRITE (*,*)&
        '  1   - do zonal average, for no averaging give the same value for minimum & maximum'
      WRITE (*,*)&
        '  2   - do zonal average for cross-correlation'
      WRITE (*,*)&
     &  ' *************************************************'
!
      CALL GET_SELECTION (SELECTION)

      if (selection .lt. min_sel) goto 10
      if (selection .gt. max_sel) goto 10
!
      RETURN
      END

! *********************************************************************
! *************************** header_low_pass *************************
! *********************************************************************
!
      SUBROUTINE header_low_pass(selection,min_sel,max_sel)
!
      IMPLICIT NONE
!
      integer :: selection , min_sel , max_sel

!             NUMBER OF INCORRECT INPUT SELECTIONS MADE
      INTEGER :: NERROR
!
      NERROR = -1
!
10    NERROR = NERROR + 1
      IF (NERROR.GT.10) THEN
         WRITE (*,*) ' HEADER_ORNL : TOO MANY INCORRECT INPUTS'
         selection = -1
         return
      END IF
      write (*,*) ' '
      WRITE (*,*)&
        ' *************************************************'
      WRITE (*,*)&
        '  0   - Skip Low-pass Filter'
      WRITE (*,*)&
        '  1   - Apply Low-pass Filter'
      WRITE (*,*)&
        '  2   - Apply Low-pass Filter for cross-correlation'
      WRITE (*,*)&
        ' *************************************************'
!
      CALL GET_SELECTION (SELECTION)

      if (selection .lt. min_sel) goto 10
      if (selection .gt. max_sel) goto 10
!
      RETURN
      END

! *********************************************************************
! *************************** header_daw ******************************
! *********************************************************************
!
      SUBROUTINE header_daw(selection,min_sel,max_sel)
!
      IMPLICIT NONE
!
      integer :: selection , min_sel , max_sel

!             NUMBER OF INCORRECT INPUT SELECTIONS MADE
      INTEGER :: NERROR
!
      NERROR = -1
!
10    NERROR = NERROR + 1
      IF (NERROR.GT.10) THEN
         WRITE (*,*) ' HEADER_ORNL : TOO MANY INCORRECT INPUTS'
         selection = -1
         return
      END IF
      write (*,*) ' '
      WRITE (*,*)&
        ' *************************************************'
      WRITE (*,*)&
        '  0   - Exit ORNL routines'
      WRITE (*,*)&
        '  1   - simple statistics on a 1D time series'
      WRITE (*,*)&
        '  2   - auto  correlation on a 1D time series'
      WRITE (*,*)&
        '  3   - power spectral density of a 1D time series'
       WRITE (*,*)&
        '  4   - cubed temporal asymmetry routine for 1D time series'
       WRITE (*,*)&
        '  5   - signed temporal asymmetry routine for 1D time series'
      WRITE (*,*)&
        ' *************************************************'
!
      CALL GET_SELECTION (SELECTION)

      if (selection .lt. min_sel) goto 10
      if (selection .gt. max_sel) goto 10
!
      RETURN
      END
