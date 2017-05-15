!
! *********************************************************************
! ******************************* do_low_pass ************************
! *********************************************************************
!
      subroutine do_low_pass(time_series,nt)

      use usr_input

      implicit none


      integer(kind=4) :: ibeg, iend, fo
      integer :: nt

      real(kind=8)  :: cf, sr
      real(kind=8)    :: time_series(*)

      real(kind=8), allocatable    :: scratch(:), fts(:)

      allocate ( scratch(nt))
      allocate ( fts(nt))

      ibeg = 1
      iend = nt

      write(*,*) 'Data sampling rate [samples/time unit]-suggested=100'
      read(*,*) sr
      write(*,*) 'filter cutoff frequency [cycles/time unit]-sugg.=30'
      read(*,*) cf
      write(*,*) 'filter order (number of passes)-sugg.=4'
      read(*,*) fo

!     write(*,*) time_series(1:nt)

      call lpf(time_series,scratch,ibeg,iend,sr,cf,fo,fts)

      time_series(1:nt) = fts(1:nt)
!     write(*,*) time_series(1:nt)

      deallocate (scratch)
      deallocate (fts)
 !
      return
      end
