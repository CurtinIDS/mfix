!
! *********************************************************************
! ******************************* do_zonal_avg ************************
! *********************************************************************
!
      subroutine do_zonal_avg(time_in_res,arr,time_series,nt)

      use geometry
      use fldvar
      use run
      use indices
      use compar
      use usr_input
      use functions

      implicit none


      integer   :: i , j , k , ijk , L , spx_num , nt , nstep_1

      real      :: time_in_res , time_now , value, time_prev
      real      :: arr(*)
      real(kind=8)    :: time_series(*)

      logical :: ask_for_times , init_read_res

      ask_for_times = .true.
      init_read_res = .true.

      usr_allow_tavg = .false.        ! can not time average
      usr_must_i_avg = .true.
      usr_must_j_avg = .true.
      usr_must_k_avg = .true.
!
      call get_usr_input(time_in_res,ask_for_times,init_read_res)
      if (usr_done) return
!
      call usr_start_time(time_in_res)    ! goto time = usr_t1
      if (usr_status .ne. 0) goto 20
!
!
! determine a SPX file to use to check for EOF
!
      do L = 1,N_SPX
         if (read_SPX(L)) spx_num = L
      end do
!
      nt = 0     ! number of times in the series
      time_prev = 0.0d0

 15   continue   ! start of time loop
!
!
      call get_same_time (READ_SPX, REC_POINTER,&
                           AT_EOF, TIME_NOW, TIME_REAL, NSTEP_1)
      if(time_now.le.time_prev) goto 16
      time_prev = time_now
      if (at_eof(spx_num)) goto 20
      if (time_now .gt. usr_t2) goto 20
!
      call usr_set_array(ijkmax2,arr,usr_var_num,usr_m,usr_n)
      call spatial_averaging(arr,usr_i1,usr_i2,usr_j1,usr_j2, &
                usr_k1,usr_k2,usr_i_avg,usr_j_avg,usr_k_avg)

      IJK = FUNIJK(usr_i1,usr_j1,usr_k1)
      time_series(nt+1) = arr(ijk)
      nt = nt + 1
16    continue
!
      goto 15


 20   continue

      return
      end
