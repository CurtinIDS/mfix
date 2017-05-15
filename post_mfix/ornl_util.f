!***********************************************************************
!* ID => ornl_util                                                     *
!* Function => utility routines                                        *
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
!* hist_nonculled_i             2001-03-20 14:44Z                      *
!* midpoints_d                  2001-12-20 22:35Z                      *
!* midpoints_s                  2001-02-19 07:55Z                      *
!* sgn_d                        2001-02-07 21:42Z                      *
!***********************************************************************

!***********************************************************************
!* ID => hist_nonculled_i                                              *
!* Function => relative-frequency histogram of an integer vector       *
!* Date => 2001-02-19                                                  *
!* Last modified => 2001-03-20 14:44Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-02-19) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-02-19 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   vec => integer vector, dimensioned 1:*; elements in [0:nsym]      *
!*   ibeg => index of first record in vec segment                      *
!*   iend => index of last record in vec segment                       *
!*   nsym => maximal expected integer in vec                           *
!* Output parameters:                                                  *
!*   hist => frequency histogram, dimensioned 0:nsym                   *
!*   nrec => number of records tallied in hist                         *
!***********************************************************************
!* Changes:                                                            *
!*   2001-02-28 ceaf => trapped out-of-range value in vec (log to      *
!*              stdout)                                                *
!*   2001-03-20 ceaf => changed double to dble                         *
!***********************************************************************

      subroutine hist_nonculled_i(vec,ibeg,iend,nsym,hist,nrec)

      integer(kind=4) vec(1:*) !.......................... data vector (input)
      integer(kind=4) ibeg !...... index of first record in vec to use (input)
      integer(kind=4) iend !....... index of last record in vec to use (input)
      integer(kind=4) nsym !.........  maximal expected integer in vec (input)
      real(kind=8) hist(0:*) !............ frequency histogram of vec (output)
      integer(kind=4) nrec !........ number of records in [ibeg,iend] (output)
      integer(kind=4) i,j !.......................................... counters
      integer(kind=4) begidx !..................... first record in vec to use
      integer(kind=4) endidx !...................... last record in vec to use

! --- Sanity and range checks ---
      nrec = (iend - ibeg + 1)
      if (nrec.le.0) then ! invalid ranges ? abort
       return
      endif

! --- Fix indices, mapping ibeg&iend to local variables ---
      begidx = ibeg
      endidx = iend

! --- Initialize histogram ---
      do i=0,nsym
       hist(i) = 0.
      enddo ! i

! --- Tally observed frequencies of vec ---
      do i=begidx,endidx
       j = vec(i)
       if ((j.lt.0).or.(j.gt.nsym)) then
        nrec = nrec - 1 ! account for out-of-range datum
        write(*,*) 'hist_nonculled_i: out-of-range i,vec(i) ',i,vec(i)
       else
        hist(j) = hist(j) + 1.0d0
       endif
      enddo ! i

! --- Calculate relative frequencies ---
      do i=0,nsym
       hist(i) = hist(i) / dble(nrec)
      enddo ! i

      return
      end

!***********************************************************************
!* ID => midpoints_d                                                   *
!* Function => midpoints of successive vector elements (real(kind=8))        *
!* Date => 2001-12-20                                                  *
!* Last modified => 2001-12-20 22:35Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-20) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-12-20 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   vec => data vector, dimensioned 1:*                               *
!*   lvec => length of vec                                             *
!* Output parameters:                                                  *
!*   mp => midpoints of successive elements of vec                     *
!***********************************************************************

      subroutine midpoints_d(vec,lvec,mp)

      real(kind=8) vec(1:*) !............................. data vector (input)
      integer(kind=4) lvec !............................ length of vec (input)
      real(kind=8) mp(1:*) !. midpoints of successive elements of vec (output)
      integer(kind=4) l !............................................. counter

! --- Find midpoints ---
      do l=1,(lvec-1)
       mp(l) = 0.5 * (vec(l) + vec(l+1))
      enddo ! l

      return
      end

!***********************************************************************
!* ID => midpoints_s                                                   *
!* Function => midpoints of successive vector elements (real(kind=4))        *
!* Date => 2001-02-19                                                  *
!* Last modified => 2001-02-19 07:55Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-02-19) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-02-19 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   vec => data vector, dimensioned 1:*                               *
!*   lvec => length of vec                                             *
!* Output parameters:                                                  *
!*   mp => midpoints of successive elements of vec                     *
!***********************************************************************

      subroutine midpoints_s(vec,lvec,mp)

      real(kind=8) vec(1:*) !............................. data vector (input)
      integer(kind=4) lvec !............................ length of vec (input)
      real(kind=8) mp(1:*) !. midpoints of successive elements of vec (output)
      integer(kind=4) l !............................................. counter

! --- Find midpoints ---
      do l=1,(lvec-1)
       mp(l) = 0.5 * (vec(l) + vec(l+1))
      enddo ! l

      return
      end

!***********************************************************************
!* ID => sgn_d                                                         *
!* Function => returns alegbraic sign (-1 if <0, 0 if 0, +1 if >= 0)   *
!* Date => 2001-01-07                                                  *
!* Last modified => 2001-02-07 21:42Z                                  *
!***********************************************************************
!* Programmer => C.S. Daw (Oak Ridge National Laboratory)              *
!* Contact information (as of 2001-01-23) :                            *
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
!*   a => number                                                       *
!***********************************************************************

      real(kind=8) function sgn_d(a)

      real(kind=8) a

      sgn_d = 1.
      if (a.lt.0.) sgn_d = -1.
      if (a.eq.0.) sgn_d = 0.

      end

