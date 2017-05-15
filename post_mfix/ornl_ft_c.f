!
! *********************************************************************
! ******************************* do_power_spectrum ************************
! *********************************************************************
!

! Updated 2005-01-03 ceaf: changed suggestion comments for reflect default
!   choices of ol (suggest 50, rather than 0), nb (no upper bound, given
!   sufficient length of time_series, and hng (typically use 3-4).

      subroutine do_power_spectrum(time_series,nt)
      use usr_input

      implicit none


      integer(kind=4) :: nu, nb, nrm, ibeg, iend
      integer :: nt,k

      real(kind=8)  :: sr, ol, hng
      real(kind=8)  :: time_series(*)

      real(kind=8), allocatable    :: f(:), p(:)

      allocate (f(nt))
      allocate (p(nt))

      ibeg = 1
      iend = nt

      write(*,*) 'Data sampling rate [samples/time unit]'
      read(*,*) sr
      write(*,*) 'binary power of FFT block size (block=2**nu) '
      write(*,*) '(typically, this depends on the desired frequency '
      write(*,*) 'resolution, but 10-13 are useful values; frequency '
      write(*,*) 'resolution (PSD bin width) is  df = sr / 2**(nu-1))'
      read(*,*) nu
      write(*,*) 'number of averaging blocks (1 is minimum; nb <= fix(length(TS)/(2**nu)).)'
      read(*,*) nb
      write(*,*) 'percentage overlap in blocks [0,99] (typically, 50)'
      read(*,*) ol
      write(*,*) 'Hanning-window power (0=no windowing; use 2-5, typically 3-4)'
      read(*,*) hng
      write(*,*) 'normalization of power (0=none, 1=sum of power, 2=data variance)'
      read(*,*) nrm

      if(nb*(1-ol/100.)*(2**nu).gt.nt) then
         write(*,*) 'TS is not long enough for nb*(1-ol/100)*2**nu'
         write(*,*) '  nb =', nb
         write(*,*) '  2**nu =', 2**nu
         write(*,*) '  ol =', ol
         write(*,*) '  length(TS) =', nt
         return
      endif


      call psd(time_series, sr, nu, nb, ol, hng, nrm, f, p)

 !
      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing power spectral density:'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing power spectral density:'
         write (40,*)  ' '
      end if

      call usr_write_input(1)

      do k = 1, 2** (nu-1)
         if (usr_fname(1:1) .eq. '*') then
            write (*,*) f(k), p(k)
         else
            write (40,*) f(k), p(k)
         end if
      end do ! k
!
      deallocate (f)
      deallocate (p)
 !
      return
      end


