!***********************************************************************
!* ID => ornl_sym                                                      *
!* Function => data-symbolization routines                             *
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
!* ts2ss_es                     2001-12-20 22:34Z                      *
!***********************************************************************

!***********************************************************************
!* ID => ts2ss_es                                                      *
!* Function => symbolizes a time series using equispatial partitioning *
!* Date => 2001-02-14                                                  *
!* Last modified => 2001-12-20 22:34Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-01-03) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243 then 865-974-7640                       *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-02-14 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   ibeg => index of first record in TS segment                       *
!*   iend => index of last record in TS segment                        *
!*   nsym => number of symbols (histogram bins)                        *
!* Output parameters:                                                  *
!*   SS => symbolized version of the time series, dimensioned 1:*;     *
!*         symbols taken from set of [0,nsym-1]                        *
!*   part => partition boundaries, including endpoints, dimensioned    *
!*           1:nsym+1                                                  *
!***********************************************************************
!* Changes:                                                            *
!*   2001-02-17 ceaf => added part to argument list                    *
!*   2001-12-20 ceaf => changed TS to real(kind=8)                           *
!***********************************************************************

      subroutine ts2ss_es(TS,ibeg,iend,nsym,SS,part)

      real(kind=8) TS(1:*) !.............................. time series (input)
      integer(kind=4) ibeg !....... index of first record in TS to use (input)
      integer(kind=4) iend !........ index of last record in TS to use (input)
      integer(kind=4) nsym !....... number of symbols (histogram bins) (input)
      integer(kind=4) SS(1:*) !............... symbolized time series (output)
      real(kind=8) part(1:*) !.................. partition boundaries (output)
      integer(kind=4) l,ll !......................................... counters
      integer(kind=4) begidx !...................... first record in SS to use
      integer(kind=4) endidx !....................... last record in SS to use
      real(kind=8) min,max !........................ minimum and maximum of TS
      real(kind=8) dx !.............................................. bin size

! --- Sanity and range checks ---
      if ((iend-ibeg+1).le.0) then ! invalid ranges ? abort
       return
      endif
      if (iend.lt.ibeg) then ! did we make it here from previous if ?
       return
      endif
      if (nsym.lt.2) then ! too few symbols - return all zeros
       ll = 1
       do l=ibeg,iend
        SS(ll) = 0
        ll = ll + 1
       enddo ! l
       return
      endif

! --- Fix indices, mapping first&last to local variables ---
      begidx = ibeg
      endidx = iend

! --- Find time-series extrema and assign bin width ---
      min = TS(begidx)
      max = TS(begidx)
      do l=begidx+1,endidx
       if (TS(l).lt.min) min = TS(l)
       if (TS(l).gt.max) max = TS(l)
      enddo ! l
      dx = (1 + 1.0d-6) * (max - min) / nsym ! 1+1.d-6 prevents overflow

! --- Symbolize the time series ---
      ll = 1 ! initial write index
      do l=begidx,endidx
       SS(ll) = int((TS(l) - min) / dx)
       ll = ll + 1
      enddo ! l

! --- Store partition boundaries ---
      part(1) = min
      do l=2,nsym
       part(l) = part(l-1) + dx
      enddo ! l
      part(nsym+1) = max

      return
      end

