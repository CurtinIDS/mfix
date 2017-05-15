!***********************************************************************
!* ID => ornl_stats                                                    *
!* Function => statistical routines                                    *
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
!* simple_statistics            2001-12-20 22:08Z                      *
!***********************************************************************

!***********************************************************************
!* ID => simple_statistics                                             *
!* Function => computes simple statistics of a time series             *
!* Date => 2001-01-03                                                  *
!* Last modified => 2001-12-20 22:08Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-20) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 1994-2001 by C.E.A. Finney.                           *
!* Placed in the PUBLIC DOMAIN 2001-01-03 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   ibeg => first record in TS for calculating statistics             *
!*   iend => last record in TS for calculating statistics              *
!* Output parameters (here, TS->TS(ibeg:iend):                         *
!*   min => minimal data value in TS                                   *
!*   max => maximal data value in TS                                   *
!*   avg => arithmetic mean (first moment) of TS                       *
!*   var => variance (second moment) of TS                             *
!*   skw => skewness (third moment) of TS                              *
!*   krt => kurtosis (fourth moment) of TS                             *
!*          Following the suggestion in Numerical Recipes, a value of 3*
!*          is subtracted, so that Gaussian data have krt=0.           *
!*   dev => standard deviation of TS (square root of variance)         *
!*   AAD => average absolute deviation about the mean of TS            *
!*          AAD is also known as the mean absolute deviation.          *
!*          AAD = < |TS(i)-avg| >                                      *
!*   tOrb => mean orbital time of TS                                   *
!*           tOrb is defined as the average time (in timesteps) between*
!*           successive upward crossings of the time-series signal     *
!*           through the data mean.                                    *
!*   nOrb => number of orbits in TS                                    *
!*           An orbit is defined as the span between two successive    *
!*           upward crossings of the time-series signal through the    *
!*           data mean.                                                *
!*   tDev => deviation time of TS                                      *
!*           Deviation time is defined as the average time (in         *
!*           timesteps) for the time-series signal to deviate one      *
!*           standard deviation.                                       *
!*           tDev = dev / < |TS(i+1)-TS(i)| >                          *
!***********************************************************************
!* Changes:                                                            *
!*   2001-02-28 ceaf => changed TS to real*4, retained all others as   *
!*              real(kind=8); renamed variable first->ibeg, last->iend       *
!*   2001-03-20 ceaf => changed double to dble                         *
!*   2001-12-20 ceaf => changed TS to real(kind=8)                           *
!***********************************************************************

      subroutine simple_statistics(TS,ibeg,iend,min,max,avg,var,skw,    &
     & krt,dev,AAD,tOrb,nOrb,tDev)

      real(kind=8) TS(1:*) !.............................. time series (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      real(kind=8) min !..................................... minimum (output)
      real(kind=8) max !..................................... maximum (output)
      real(kind=8) avg !........................................ mean (output)
      real(kind=8) var !.................................... variance (output)
      real(kind=8) skw !.................................... skewness (output)
      real(kind=8) krt !.................................... kurtosis (output)
      real(kind=8) dev !.......................... standard deviation (output)
      real(kind=8) AAD !........ average absolute deviation from mean (output)
      real(kind=8) tOrb !....................... average orbital time (output)
      integer(kind=4) nOrb !........................ number of orbits (output)
      real(kind=8) tDev !............................. deviation time (output)
      integer(kind=4) lk !............................................ counter
      integer(kind=4) nCross !......... number of upward crossings through avg
      integer(kind=4) nrec !................. number of records in [ibeg,iend]
      integer(kind=4) TCount !............................... timestep counter
      real(kind=8) old !...... previous point for comparison of mean crossings
      real(kind=8) ratio !............................. interpolation distance
      real(kind=8) sum !.................................................. sum
      real(kind=8) prod !............................................. product
      real(kind=8) tMin !....................... temporal index of first orbit
      real(kind=8) tMax !........................ temporal index of last orbit

! --- Sanity and range checks ---
      nrec = iend - ibeg + 1
      if (nrec.le.0) then ! invalid ranges ? abort
       return
      endif

! --- Determine data mean ---
      sum = 0.
      do lk=ibeg,iend
       sum = sum + TS(lk)
      enddo ! lk
      avg = sum / dble(nrec)

! --- Determine extrema and estimate moments and AAD ---
      min = TS(ibeg)
      max = TS(ibeg)
      var = 0.
      skw = 0.
      krt = 0.
      AAD = 0.
      do lk=ibeg,iend
       if (TS(lk).lt.min) min = TS(lk)
       if (TS(lk).gt.max) max = TS(lk)
       sum = TS(lk) - avg
       AAD = AAD + abs(sum)
       prod = sum * sum
       var = var + prod
       prod = prod * sum
       skw = skw + prod
       prod = prod * sum
       krt = krt + prod
      enddo ! lk
      AAD = AAD / dble(nrec)
      var = var / dble(nrec-1)
      dev = sqrt(var)
      if (var.ne.0.) then
       skw = skw / (dble(nrec) * dev * var)
       krt = krt / (dble(nrec) * var * var) - 3.
      else
       skw = 0.
       krt = 0.
      endif

! --- Determine mean orbital time ---
      nCross = 0
      nOrb = 0
      tCount = 0
      tMin = 0.
      tMax = 0.
      do lk=ibeg,iend
       if ((tCount.gt.0).and.(TS(lk).ge.avg).and.(old.lt.avg)) then
! . . . . . . . . changed TS(lk).gt.avg to TS(lk).ge.avg 2000-01-04
        if (TS(lk).ne.old) then
         ratio = (TS(lk) - avg) / (TS(lk) - old)
        else
         ratio = 0.
        endif
        if (nCross.eq.0) then
         nCross = nCross + 1
         tMin = tCount + ratio
        else
         nCross = nCross + 1
         nOrb = nOrb + 1
         tMax = tCount + ratio
        endif
       endif
       old = TS(lk)
       tCount = tCount + 1
      enddo ! lk
      if (nOrb.gt.0) then
       tOrb = (tMax - tMin) / dble(nOrb)
      else
       tOrb = 0.
      endif

! --- Determine deviation time ---
      sum = 0.
      do lk=ibeg,(iend-1)
       sum = sum + abs(TS(lk+1) - TS(lk))
      enddo ! lk
      tDev = dev / sum

      return
      end

