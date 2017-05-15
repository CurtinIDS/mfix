!
! *********************************************************************
! **************************** do_simple_stats ************************
! *********************************************************************
!
      subroutine do_simple_stats(time_series,nt)

      use geometry
      use fldvar
      use run
      use indices
      use compar
      use usr_input
      use functions

      implicit none

      real(kind=8)  :: vmin , vmax , vavg , variance , skw , krt
      real(kind=8)  :: dev  , AAD , tOrb  , tDev

      integer :: nOrb , L , spx_num , nt , nstep_1
      integer :: i , j , k , ijk

      real(kind=8)    :: time_series(*)

      call simple_statistics(time_series,1,nt,vmin,vmax,vavg, &
                      variance,skw,krt,dev,AAD,tOrb,nOrb,tDev)

      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing simple statistics for :'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing simple statistics for :'
         write (40,*)  ' '
      end if

      call usr_write_input(1)

      if (usr_fname(1:1) .eq. '*') then
         write (*,*)  ' min  = ' , vmin
         write (*,*)  ' max  = ' , vmax
         write (*,*)  ' avg  = ' , vavg
         write (*,*)  ' var  = ' , variance
         write (*,*)  ' skw  = ' , skw
         write (*,*)  ' krt  = ' , krt
         write (*,*)  ' dev  = ' , dev
         write (*,*)  ' AAD  = ' , AAD
         write (*,*)  ' tOrb = ' , tOrb
         write (*,*)  ' nOrb = ' , nOrb
      else
         write (40,*) ' min  = ' , vmin
         write (40,*) ' max  = ' , vmax
         write (40,*) ' avg  = ' , vavg
         write (40,*) ' var  = ' , variance
         write (40,*) ' skw  = ' , skw
         write (40,*) ' krt  = ' , krt
         write (40,*) ' dev  = ' , dev
         write (40,*) ' AAD  = ' , AAD
         write (40,*) ' tOrb = ' , tOrb
         write (40,*) ' nOrb = ' , nOrb
      end if

      usr_allow_tavg = .true.
      usr_must_i_avg = .false.
      usr_must_j_avg = .false.
      usr_must_k_avg = .false.
 !
      return
      end
