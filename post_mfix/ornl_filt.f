!***********************************************************************
!* ID => ornl_filt                                                     *
!* Function => filtering routines                                      *
!* Release date => 2001-12-20                                          *
!***********************************************************************
!* Packager => C.E.A. Finney (Oak Ridge National Laboratory)           *
!* Contact information (as of 2001-12-20) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Routine                      Last modified                          *
!* hpf                          2001-12-18 22:00Z                      *
!* lpf                          2001-12-18 20:45Z                      *
!***********************************************************************

!***********************************************************************
!* ID => hpf                                                           *
!* Function => highpass-filters time series                            *
!* Date => 2001-12-14                                                  *
!* Last modified => 2001-12-18 22:00Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-14) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-12-14 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   foo => scratch array, dimensioned at least 1:(iend-ibeg+1)        *
!*   ibeg => first record in TS segment                                *
!*   iend => last record in TS segment                                 *
!*   sr => data sampling rate [samples/time unit]                      *
!*   cf => filter cutoff frequency [/time unit]                        *
!*   fo => filter order (4 is preferred)                               *
!* Output parameters:                                                  *
!*   fts => filtered time series                                       *
!***********************************************************************

      subroutine hpf(TS,foo,ibeg,iend,sr,cf,fo,fts)

      implicit none

      real(kind=8) PI
      parameter (PI=3.14159265359)

      real(kind=8) TS(1:*) !.............................. time series (input)
      real(kind=8) foo(1:*) !........................... scratch array (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      real(kind=8) sr !........ data sampling rate [samples/time unit] (input)
      real(kind=8) cf !.... filter cutoff frequency [cycles/time unit] (input)
      integer(kind=4) fo !............ filter order (number of passes) (input)
      real(kind=8) fts(1:*) !................... filtered time series (output)
      integer(kind=4) i !............................................. counter
      integer(kind=4) m !....................................... order counter
      integer(kind=4) nRec !................. number of records in [ibeg,iend]
      real(kind=8) tau !................................. filter time constant
      real(kind=8) dt !........ data sampling interval (inverse sampling rate)
      real(kind=8) alpha !...................................... filter factor

! --- Sanity and range checks ---
      nRec = iend - ibeg + 1
      if (nRec.le.0) then ! invalid ranges ? abort
       write(*,*) 'ornl.hpf: FATAL: iend < ibeg.'
       return
      endif
      if (cf.le.0) then
       write(*,*) 'ornl.hpf: FATAL: cf <= 0.'
       return
      endif
      if (cf.gt.sr) then
       write(*,*) 'ornl.hpf: FATAL: cf > sr.'
       return
      endif
      if (fo.lt.1) then
       write(*,*) 'ornl.hpf: FATAL: fo < 1.'
       return
      endif
      if (fo.ge.10) then
       write(*,*) 'ornl.hpf: fo >= 10.  WARNING: abnormally high value.'
       return
      endif

! --- Copy data from TS to foo; foo gets overwritten ---
      do i=1,nRec
       foo(i) = TS(ibeg+i-1)
      enddo ! i

! --- Filter constants ---
      tau = 1 / (2 * PI * cf) ! time constant
      dt = 1 / sr ! data inverse sampling rate
      alpha = 1 / (1 + dt / tau) ! filter factor

! --- Filter data ---
      do m=1,fo
       ! apply filter
       do i=1,(nRec-1)
        fts(i+1) = alpha * (foo(i+1) - foo(i) + fts(i))
       enddo ! i
       ! store data for next pass
       do i=1,nRec
        foo(i) = fts(i)
       enddo ! i
      enddo ! m

      return
      end

!***********************************************************************
!* ID => lpf                                                           *
!* Function => lowpass-filters time series                             *
!* Date => 2001-12-05                                                  *
!* Last modified => 2001-12-18 20:45Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-06) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-12-06 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   foo => scratch array, dimensioned at least 1:(iend-ibeg+1)        *
!*   ibeg => first record in TS segment                                *
!*   iend => last record in TS segment                                 *
!*   sr => data sampling rate [samples/time unit]                      *
!*   cf => filter cutoff frequency [/time unit]                        *
!*   fo => filter order (4 is preferred)                               *
!* Output parameters:                                                  *
!*   fts => filtered time series                                       *
!***********************************************************************
!* Changes:                                                            *
!*   2001-12-06 ceaf => added "implicit none"                          *
!*   2001-12-10 ceaf => changed nr to ibeg and iend in argument list   *
!*   2001-12-14 ceaf => added foo scratch array                        *
!*   2001-12-18 ceaf => fixed this comment header                      *
!***********************************************************************

      subroutine lpf(TS,foo,ibeg,iend,sr,cf,fo,fts)

      implicit none

      real(kind=8) PI
      parameter (PI=3.14159265359)

      real(kind=8) TS(1:*) !.............................. time series (input)
      real(kind=8) foo(1:*) !........................... scratch array (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      real(kind=8) sr !........ data sampling rate [samples/time unit] (input)
      real(kind=8) cf !.... filter cutoff frequency [cycles/time unit] (input)
      integer(kind=4) fo !............ filter order (number of passes) (input)
      real(kind=8) fts(1:*) !................... filtered time series (output)
      integer(kind=4) i !............................................. counter
      integer(kind=4) m !....................................... order counter
      integer(kind=4) nRec !................. number of records in [ibeg,iend]
      real(kind=8) tau !................................. filter time constant
      real(kind=8) dt !........ data sampling interval (inverse sampling rate)
      real(kind=8) alpha !...................................... filter factor

! --- Sanity and range checks ---
      nRec = iend - ibeg + 1
      if (nRec.le.0) then ! invalid ranges ? abort
       write(*,*) 'ornl.lpf: FATAL: iend < ibeg.'
       return
      endif
      if (cf.le.0) then
       write(*,*) 'ornl.lpf: FATAL: cf <= 0.'
       return
      endif
      if (cf.gt.(sr/2)) then
       write(*,*) 'ornl.lpf: FATAL: cf > sr/2.'
       return
      endif
      if (fo.lt.1) then
       write(*,*) 'ornl.lpf: FATAL: fo < 1.'
       return
      endif
      if (fo.ge.10) then
       write(*,*) 'ornl.lpf: fo >= 10.  WARNING: abnormally high value.'
       return
      endif

! --- Copy data from TS to foo; foo gets overwritten ---
      do i=1,nRec
       foo(i) = TS(ibeg+i-1)
      enddo ! i

! --- Filter constants ---
      tau = 1 / (2 * PI * cf) ! time constant
      dt = 1 / sr ! data inverse sampling rate
      alpha = exp(-dt / tau) ! filter factor

! --- Filter data ---
      do m=1,fo
       ! apply filter
       do i=1,(nRec-1)
        fts(i+1) = alpha * fts(i) + (1 - alpha) * foo(i)
       enddo ! i
       ! store data for next pass
       do i=1,nRec
        foo(i) = fts(i)
       enddo ! i
      enddo ! m

      return
      end

