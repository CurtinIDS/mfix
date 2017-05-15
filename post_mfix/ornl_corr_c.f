!
! *********************************************************************
! ******************************* do_auto_corr ************************
! *********************************************************************
!
      subroutine do_auto_corr(time_series, nt)

      use geometry
      use fldvar
      use run
      use indices
      use compar
      use usr_input
      use functions

      implicit none


      integer :: minlag , maxlag
      integer :: i , j , k , ijk , L , spx_num , nt , nstep_1
      integer :: lagstep

      real(kind=8)    :: time_series(*)

      real(kind=8), allocatable    :: acf(:,:)
      integer, allocatable   :: lags(:)

      write (*,*) ' enter minlag , maxlag'
      read  (*,*) minlag , maxlag

      i = abs(maxlag-minlag) + 1
      allocate ( acf(2,i) )
      allocate ( lags(i) )

      call auto_correlation(time_series,1,nt,minlag,maxlag,lagstep,lags,acf)


      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing auto correlation for :'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing auto correlation for :'
         write (40,*)  ' '
      end if

      call usr_write_input(1)

      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' minlag :' , minlag
         write (*,*) ' maxlag :' , maxlag
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' minlag :' , minlag
         write (40,*) ' maxlag :' , maxlag
         write (40,*)  ' '
      end if

      if (minlag .le. maxlag) then
         lagstep = 1
      else
         lagstep = -1
      end if
      L = 0 ! initial write index
      do k = minlag,maxlag,lagstep
         L = L + 1
         if (usr_fname(1:1) .eq. '*') then
            write (*,*) acf(1,L) , acf(2,L)
         else
            write (40,*) acf(1,L) , acf(2,L)
         end if
      end do ! k

      deallocate (acf)
      deallocate (lags)

 !
      return
      end


!
! *********************************************************************
! ******************************* do_cross_corr ***********************
! *********************************************************************
!
      subroutine do_cross_corr(time_series, time_series2, nt)

      use geometry
      use fldvar
      use run
      use indices
      use compar
      use usr_input
      use functions

      implicit none


      integer :: minlag , maxlag , lagstep
      integer :: i , j , k , ijk , L , spx_num , nt , nstep_1

      real    :: time_series(*) , time_series2(*)

      real(kind=8), allocatable    :: ccf(:,:)
      integer, allocatable   :: lags(:)

      write (*,*) ' enter minlag , maxlag, lagstep'
      read  (*,*) minlag , maxlag, lagstep
!
      i = abs(maxlag-minlag) + 1
      allocate ( ccf(2,i) )
      allocate ( lags(i) )

      call cross_correlation(time_series,time_series2,1,nt, &
                                        minlag,maxlag,lagstep,lags,ccf)

!

      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing cross correlation for :'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing cross correlation for :'
         write (40,*)  ' '
      end if

      call usr_write_input(2)
      call usr_write_input(1)

      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' minlag :' , minlag
         write (*,*) ' maxlag :' , maxlag
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' minlag :' , minlag
         write (40,*) ' maxlag :' , maxlag
         write (40,*)  ' '
      end if

      if (minlag .le. maxlag) then
         lagstep = 1
      else
         lagstep = -1
      end if
      L = 0 ! initial write index
      do k = minlag,maxlag,lagstep
         L = L + 1
         if (usr_fname(1:1) .eq. '*') then
            write (*,*) ccf(1,L) , ccf(2,L)
         else
            write (40,*) ccf(1,L) , ccf(2,L)
         end if
      end do ! k

      deallocate (ccf)
      deallocate (lags)

!
      return
      end

!
! *********************************************************************
! ******************************** do_t3_symm *************************
! *********************************************************************
!
      subroutine do_t3_symm(time_series, nt)

      use geometry
      use fldvar
      use run
      use indices
      use compar
      use usr_input
      use functions

      implicit none


      integer :: minlag , maxlag
      integer :: i , j , k , ijk , L , spx_num , nt , nstep_1
      integer :: lagstep

      real(kind=8)    :: time_series(100000)

      real(kind=8), allocatable    :: tlag(:)
      real(kind=8), allocatable    :: tsym(:)

      write (*,*) ' enter minlag , maxlag'
      read  (*,*) minlag , maxlag

      if (minlag .le. maxlag) then
         lagstep =  1
      else
         lagstep = -1
      end if
!
      i = abs(maxlag-minlag) + 1
      allocate ( tlag(i) )
      allocate ( tsym(i) )

      call t3sym(time_series,1,nt,minlag,maxlag,lagstep,tlag,tsym)


      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing temporal asymmetry (t3sym) for :'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing temporal asymmetry (t3sym) for :'
         write (40,*)  ' '
      end if

      call usr_write_input(1)

      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' minlag :' , minlag
         write (*,*) ' minlag :' , maxlag
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' minlag :' , minlag
         write (40,*) ' minlag :' , maxlag
         write (40,*)  ' '
      end if

      L = 0 ! initial write index
      do k = minlag,maxlag,lagstep
         L = L + 1
         if (usr_fname(1:1) .eq. '*') then
            write (*,*) tlag(L) , tsym(L)
         else
            write (40,*) tlag(L) , tsym(L)
         end if
      end do ! k

      deallocate (tlag)
      deallocate (tsym)

      usr_allow_tavg = .true.
      usr_must_i_avg = .false.
      usr_must_j_avg = .false.
      usr_must_k_avg = .false.
 !
      return
      end
!
! *********************************************************************
! ****************************** do_trevsgn *************************
! *********************************************************************
!
      subroutine do_trevsgn(time_series, nt)

      use geometry
      use fldvar
      use run
      use indices
      use compar
      use usr_input
      use functions

      implicit none


      integer :: minlag , maxlag
      integer :: i , j , k , ijk , L , spx_num , nt , nstep_1
      integer :: lagstep

      real(kind=8)    :: time_series(*)

      real(kind=8), allocatable    :: tlag(:)
      real(kind=8), allocatable    :: tsym(:)

      write (*,*) ' enter minlag , maxlag'
      read  (*,*) minlag , maxlag

      if (minlag .le. maxlag) then
         lagstep =  1
      else
         lagstep = -1
      end if
!
      i = abs(maxlag-minlag) + 1
      allocate ( tlag(i) )
      allocate ( tsym(i) )

      call trevsgn(time_series,1,nt,minlag,maxlag,lagstep,tlag,tsym)


      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing temporal asymmetry (t3sym) for :'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing temporal asymmetry (t3sym) for :'
         write (40,*)  ' '
      end if

      call usr_write_input(1)

      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' minlag :' , minlag
         write (*,*) ' minlag :' , maxlag
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' minlag :' , minlag
         write (40,*) ' minlag :' , maxlag
         write (40,*)  ' '
      end if

      L = 0 ! initial write index
      do k = minlag,maxlag,lagstep
         L = L + 1
         if (usr_fname(1:1) .eq. '*') then
            write (*,*) tlag(L) , tsym(L)
         else
            write (40,*) tlag(L) , tsym(L)
         end if
      end do ! k

      deallocate (tlag)
      deallocate (tsym)

      usr_allow_tavg = .true.
      usr_must_i_avg = .false.
      usr_must_j_avg = .false.
      usr_must_k_avg = .false.
 !
      return
      end
