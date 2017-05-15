!***********************************************************************
!* ID => ornl_ft                                                       *
!* Function => Fourier-transform routines                              *
!* Release date => 2005-01-04                                          *
!***********************************************************************
!* Packager => C.E.A. Finney (Oak Ridge National Laboratory)           *
!* Contact information (as of 2005-01-04) :                            *
!*   Email => <finneyc@ornl.gov>                                       *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Routine                      Last modified                          *
!* psd                          2005-01-03 22:00Z                      *
!* fft                          2002-02-11 23:53Z                      *
!* bitrev                       2001-12-20 22:00Z                      *
!***********************************************************************
!* Changes:                                                            *
!*   2001-11-02 ceaf => updated psd with float hng and new nrm         *
!*              parameter                                              *
!*   2001-12-20 ceaf => changed order of arguments in bitrev           *
!*   2002-02-11 ceaf => fixed declarations in fft which caused compi-  *
!*              lation errors on some systems                          *
!*   2005-01-03 ceaf => changed default parameters in psd              *
!***********************************************************************

!***********************************************************************
!* ID => psd                                                           *
!* Function => estimates power spectral density of a time series       *
!* Date => 2001-03-22                                                  *
!* Last modified => 2005-01-03 22:00Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2005-01-03) :                            *
!*   Email => <finneyc@ornl.gov>                                       *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-04-01 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   sr => data sampling rate [samples / time unit] (this is a         *
!*         characteristic of the data, but if none specified, then use *
!*         2, so that reported frequencies will range from about 0 to  *
!*         one Nyquist frequency)                                      *
!*   nu => binary power of FFT block size (block=2**nu) (typically,    *
!*         this depends on the desired frequency resolution, but 10-13 *
!*         are useful values; frequency resolution (PSD bin width) is  *
!*         df = sr / 2**(nu-1))                                        *
!*   nb => number of averaging blocks (1 is minimum, 15 is a useful    *
!*         maximum; the constraint is nb <= fix(length(TS)/(2**nu)).)  *
!*   ol => percentage overlap in blocks [0,10]; typically, use 0.      *
!*   hng => Hanning-window power (0 -> no windowing, 2-5 -> useful)    *
!*   nrm => normalization of power (1=sum of power, 2=data variance)   *
!* Output parameters:                                                  *
!*   f => frequency bin centers                                        *
!*   p => power at each corresponding bin in f                         *
!***********************************************************************
!* Changes:                                                            *
!*   2001-11-02 ceaf => changed hng from integer(kind=4) to real*4 to allow  *
!*              exotic windowing powers; added nrm input parameter.    *
!*   2001-12-20 ceaf => changed all floats from real*4 to real(kind=8)       *
!*   2005-01-03 ceaf => changed default values for nb,ol and ranges in *
!*              variance calculation/normalization to improve accuracy *
!***********************************************************************

      subroutine psd(TS,sr,nu,nb,ol,hng,nrm,f,p)

      external fft

      integer(kind=4) NUMAX !................. maximal value of nu (parameter)
      parameter (NUMAX=16)
      integer(kind=4) NMAX !................... maximal value of n (parameter)
      parameter (NMAX=2**NUMAX)

      real(kind=8) TS(1:*) !.............................. time series (input)
      real(kind=8) sr !............................ data sampling rate (input)
      integer(kind=4) nu !............. binary power of FFT block size (input)
      integer(kind=4) nb !................. number of averaging blocks (input)
      real(kind=8) ol !...................... block overlap percentage (input)
      real(kind=8) hng !......................... Hanning-window power (input)
      integer(kind=4) nrm !... normalization (1=sum power, 2=variance) (input)
      real(kind=8) f(1:*) !.................... frequency bin centers (output)
      real(kind=8) p(1:*) !..... power at each corresponding bin in f (output)
      integer(kind=4) n !...................................... FFT block size
      integer(kind=4) n2 !................................ FFT half block size
      integer(kind=4) nol !......................... number of overlap records
      real(kind=8) yr(0:NMAX-1) !.............................. real component
      real(kind=8) yi(0:NMAX-1) !......................... imaginary component
      integer(kind=4) k,l !........................................... indices
      integer(kind=4) bi,ei !.......... pointers in TS for block begin and end
      real(kind=8) var !........................................ data variance
      real(kind=8) sum !........................... sum of block data or power
      real(kind=8) avg !...................................... block data mean

! --- Sanity and range checks ---
! !!! Calling program should verify that TS is long enough for nb*n. !!!
      if (nu.gt.NUMAX) then
       write(*,*) 'ornl.psd : nu > NUMAX.  Reduce nu or recompile psd.'
       return
      endif
      if (nu.lt.1) then
       write(*,*) 'ornl.psd: nu < 1.  Respecify valid nu.'
       return
      endif
      n = 2 ** nu
      n2 = 2 ** (nu-1)
      if (nb.lt.0) then
       nb = 1
       write(*,*) 'ornl.psd: nb < 0; set to 1.'
      endif
! 2005-01-03 ceaf: Removed - given short data sets in MFIX simulations,
!   all data may be be needed!
!     if (nb.gt.15) then
!      write(*,*) 'ornl.psd: nb > 15.  WARNING: abnormally high value.'
!     endif
      if (ol.lt.0.) then
       ol = 50.
       write(*,*) 'ornl.psd: ol < 0; set to 50.'
      endif
      if (ol.gt.50.) then
       write(*,*) 'ornl.psd: ol > 50.  WARNING: biased results ?'
      endif
      if (ol.gt.99.) then
       ol = 50.
       write(*,*) 'ornl.psd: ol > 99; set to 50.'
      endif
      nol = int(ol / 100. * n)
      if (hng.lt.0) then
       hng = 4
       write(*,*) 'ornl.psd: hng < 0; set to 4.'
      endif
      if (hng.eq.0) then
       write(*,*) 'ornl.psd: hng = 0.  WARNING: no windowing function !'
       write(*,*) 'ornl.psd: Lack of windowing violates DFT assumption.'
      endif
      if (hng.gt.10) then
       write(*,*) 'ornl.psd: hng > 10.  WARNING: biased results ?'
       write(*,*) 'ornl.psd: Too much frequency info may be damped out.'
      endif

! --- Initialize ---
      do k=1,n2
       f(k) = (sr / n / 2) + (k - 1) * sr / n
       p(k) = 0.
      enddo ! k

! --- Calculate PSD over nb blocks ---
      bi = 1
      ei = bi + n - 1
      do l=1,nb
! ... calculate block data mean ...
       sum = 0.
       do k=0,(n-1)
        sum = sum + TS(bi + k)
       enddo ! k
       avg = sum / n
! ... subtract block data mean ...
       do k=0,(n-1)
        yr(k) = TS(bi + k) - avg
        yi(k) = 0.
       enddo ! k
       call fft(yr, yi, nu, 1, hng) ! 1 is forward FFT
       do k=1,n2
        p(k) = p(k) + yr(k-1)**2 + yi(k-1)**2
       enddo ! k
       bi = bi + n - nol ! advance pointer for next block
       ei = bi + n - 1 ! ditto
      enddo ! l

! --- Average power over number of blocks ---
      do k=1,n2
       p(k) = p(k) / (nb * n2)
      enddo ! k

! --- Normalize if requested ---
! ... Normalize by sum of power ...
      if (nrm.eq.1) then
       ! find sum
       sum = 0.
       do k=1,n2
        sum = sum + p(k)
       enddo ! k
       ! normalize
       do k=1,n2
        p(k) = p(k) / sum
       enddo ! k
      endif
! ... Standardize by data variance ...
      if (nrm.eq.2) then
       ! find data mean
       sum = 0.
! 2005-01-03 ceaf: changed loop to 1,ei from nb*n to calculate
!   variance based only on the data used.
       do l=1,ei
        sum = sum + TS(l)
       enddo ! l
       avg = sum / n
       ! find data variance
       var = 0.
! 2005-01-03 ceaf: changed loop to 1,ei from nb*n to calculate
!   variance based only on the data used.
       do l=1,ei ! 2005-01-03: was nb*n
        var = var + (TS(l) - avg)**2
       enddo ! l
       var = var / (ei - 1)
       ! standardize
       if (var.gt.0) then
        do k=1,n2
         p(k) = p(k) / var
        enddo ! k
       else
        write(*,*) 'ornl.psd: WARNING: data have zero variance !'
       endif
      endif

      return
      end

!***********************************************************************
!* ID => fft                                                           *
!* Function => estimates Fast Fourier Transform of a time-series block *
!* Date => 2001-01-04                                                  *
!* Last modified => 2002-02-11 23:53Z                                  *
!***********************************************************************
!* Programmer => C.S. Daw (Oak Ridge National Laboratory)              *
!* Contact information (as of 2005-01-03) :                            *
!*   Email => <dawcs@ornl.gov>                                         *
!*   Telephone => 865-946-1341                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.S. Daw.                                     *
!* Placed in the PUBLIC DOMAIN 2001-01-07 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   yr => array with data real component, dimensioned 0:* (modified)  *
!*   yi => array with data imaginary component, of same dimension as yr*
!*         (modified)                                                  *
!*   nu => binary power of FFT block size (block=2**nu)                *
!*   dir => transform direction (>0 -> forward, <0 -> reverse)         *
!*   hng => Hanning-window power (0 -> no windowing, 2-5 -> useful)    *
!* Output parameters:                                                  *
!*   yr => real component of Fourier transform                         *
!*   yi => imaginary component of Fourier transform                    *
!***********************************************************************
!* Changes:                                                            *
!*   2001-11-02 ceaf => changed hng from integer(kind=4) to real*4 to allow  *
!*              exotic windowing powers                                *
!*   2001-12-20 ceaf => changed all floats to real(kind=8); added declaration*
!*              for bitrev; changed calls to bitrev for new order of   *
!*              input arguments                                        *
!*   2002-02-11 ceaf => removed external and declaration statement for *
!*              bitrev - caused compilation errors on some systems     *
!*   2006-03-01 sp => Changed real*4 to real(kind=8) to match the calling    *
!*              routine                                                *
!***********************************************************************

      subroutine fft(yr,yi,nu,dir,hng)

      real(kind=8) yr(0:*) !.................... real component (input/output)
      real(kind=8) yi(0:*) !............... imaginary component (input/output)
      integer(kind=4) nu !............. binary power of FFT block size (input)
      integer(kind=4) dir !....................... transform direction (input)
      real(kind=8) hng !......................... Hanning-window power (input)
      real(kind=8) pi
      integer(kind=4) n !...................................... FFT block size
      integer(kind=4) k,n2,nu1,l,i,a,b,kn2
      real(kind=8) hn,arg,s,c,tr,ti,yt

      pi = 3.141592654
      n = 2 ** nu

! --- Apply windowing function ---
      if (hng.lt.0) then
       hng = 4.
       write(*,*) 'ornl.fft: NOTICE: Hanning-window power set to 4.'
      endif
      if (hng.gt.0.) then
       do k=0,(n-1)
        hn = ((1 - cos(2 * pi * k / n)) / 2) ** hng
        yr(k) = yr(k) * hn
        yi(k) = yi(k) * hn
       enddo ! k
      endif

! --- Change imaginary-component signs for reverse FFT ---
      if (dir.lt.0) then
       do k=0,(n-1)
        yi(k) = -yi(k)
       enddo ! k
      endif

! --- FFT ---
      n2 = n / 2
      nu1 = nu - 1
      k = 0

      do l = 1, nu
 1     do i = 1,n2
        a = int(k / (2 ** nu1))
        call bitrev(a, nu, b)
        arg = 2 * pi * b / dble(n)
        c = cos(arg)
        s = sin(arg)
        kn2 = k + n2
        tr = yr(kn2) * c + yi(kn2) * s
        ti = yi(kn2) * c - yr(kn2) * s
        yr(kn2) = yr(k) - tr
        yi(kn2) = yi(k) - ti
        yr(k) = yr(k) + tr
        yi(k) = yi(k) + ti
        k = k + 1
       enddo ! i
       k = k + n2
       if (k.lt.(n-1)) goto 1
       k = 0
       nu1 = nu1 - 1
       n2 = n2 / 2
      enddo ! l

      do k=0,(n-1)
       a = k
       call bitrev(a, nu, b)
       i = b
       if (i.le.k) goto 2
       yt = yr(i)
       yr(i) = yr(k)
       yr(k) = yt
       yt = yi(i)
       yi(i) = yi(k)
       yi(k) = yt
 2     continue
      enddo ! k

! --- Scale components and reinvert yi for reverse FFT ---
      if (dir.lt.0) then
       do k=0,(n-1)
        yi(k) = -yi(k) / n
        yr(k) = yr(k) / n
       enddo ! k
      endif

      return
      end

!***********************************************************************
!* ID => bitrev                                                        *
!* Function => reverse bits in an integer (used by FFT)                *
!* Date => 2001-01-07                                                  *
!* Last modified => 2001-12-20 22:00Z                                  *
!***********************************************************************
!* Programmer => C.S. Daw (Oak Ridge National Laboratory)              *
!* Contact information (as of 2005-01-03) :                            *
!*   Email => <dawcs@ornl.gov>                                         *
!*   Telephone => 865-946-1341                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.S. Daw.                                     *
!* Placed in the PUBLIC DOMAIN 2001-01-07 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   a => integer                                                      *
!*   nu => number of bits                                              *
!* Output parameters:                                                  *
!*   b => integer                                                      *
!***********************************************************************
!* Changes:                                                            *
!*   2001-12-20 ceaf => changed order of input arguments; changed types*
!*              from integer to integer(kind=4)                              *
!***********************************************************************

      subroutine bitrev(a,nu,b)

      integer(kind=4) a,b,nu
      integer(kind=4) a1,i,a12

      a1 = a
      b = 0
      do i=1,nu
       a12 = int(a1 / 2)
       b = b * 2 + (a1 - 2 * a12)
       a1 = a12
      enddo ! i

      return
      end

