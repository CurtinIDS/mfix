!***********************************************************************
!* ID => ornl_corr                                                     *
!* Function => temporal-correlation routines                           *
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
!* auto_correlation             2001-12-10 01:56Z                      *
!* cross_correlation            2001-12-10 01:59Z                      *
!* mif1                         2001-12-20 21:45Z                      *
!* t3sym                        2001-12-10 01:55Z                      *
!* trevsgn                      2001-12-10 01:50Z                      *
!***********************************************************************
!* Changes:                                                            *
!*   2001-03-09 ceaf => added trevsgn_s (erroneously omitted before)   *
!*   2001-12-07 ceaf => added trevsgn                                  *
!*   2001-12-10 ceaf => removed trevsgn_s                              *
!***********************************************************************

!***********************************************************************
!* ID => auto_correlation                                              *
!* Function => estimates autocorrelation function of a time series     *
!* Date => 2001-01-03                                                  *
!* Last modified => 2001-12-10 01:56Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-10) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 1994-1996 by C.E.A. Finney.                           *
!* Placed in the PUBLIC DOMAIN 2001-01-03 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   ibeg => first record in TS segment                                *
!*   iend => last record in TS segment                                 *
!*   lmin => minimal value of lag for computation of autocorrelation   *
!*             function; should be in [0,lmax]                         *
!*   lmax => maximal value of lag for computation of autocorrelation   *
!*             function; should be in [lmin,length(TS)-lmin]           *
!*   lstp => lag increment                                             *
!* Output parameters:                                                  *
!*   lag => lags between lmin and lmax step lstp, dimensioned 1:*      *
!*   acf => autocorrelation values at each lag contained in lag,       *
!*          dimensioned 1:*                                            *
!***********************************************************************
!* Changes:                                                            *
!*   2001-02-25 ceaf => changed TS to real(kind=4), retained all others as   *
!*              real(kind=8); renamed variable first->ibeg, last->iend,      *
!*              minlag->lmin, maxlag->lmax; added parameter lstp;      *
!*              split acf into lag and acf                             *
!*   2001-02-28 ceaf => fixed ck denominator (added +1)                *
!*   2001-03-20 ceaf => changed double to dble                         *
!*   2001-12-06 sp+ceaf => added "implicit none"; renamed erroneous    *
!*              variable lagstep to lstp; changed TS to real(kind=8);        *
!*              changed denominator in variance calculation from       *
!*              nRec-1 to nRec                                         *
!***********************************************************************

      subroutine auto_correlation(TS,ibeg,iend,lmin,lmax,lstp,lag,acf)

      implicit none

      real(kind=8) TS(1:*) !.............................. time series (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      integer(kind=4) lmin !.......... minimal value of lag to compute (input)
      integer(kind=4) lmax !.......... maximal value of lag to compute (input)
      integer(kind=4) lstp !............................ lag increment (input)
      integer(kind=4) lag(1:*) !........................ list of lags (output)
      real(kind=8) acf(1:*) !................. autocorrelation values (output)
      integer(kind=4) lk !............................................ counter
      integer(kind=4) ll !...................................... write counter
      integer(kind=4) k !......................................... lag counter
      integer(kind=4) nRec !................. number of records in [ibeg,iend]
      real(kind=8) avg !................................................. mean
      real(kind=8) ck !................................. lag-k correlation sum
      real(kind=8) var !............................................. variance
      real(kind=8) sum !.................................................. sum
      integer(kind=4) begidx !...................... first record in TS to use
      integer(kind=4) endidx !....................... last record in TS to use

! --- Sanity and range checks ---
      nRec = iend - ibeg + 1
      if (nRec.le.0) then ! invalid ranges ? abort
       return
      endif
      if ((lmin.lt.0).or.(lmax.lt.0)) then
       return ! no negative lags allowed
      endif
      if ((lmin.gt.lmax).and.(lstp.gt.0)) then
       lstp = -lstp ! decrement through lag loop
      endif

! --- Determine data mean & variance ---
      sum = 0.
      do lk=ibeg,iend
       sum = sum + TS(lk)
      enddo ! lk
      avg = sum / dble(nRec)
      var = 0.
      do lk=ibeg,iend
       var = var + (TS(lk) - avg)**2
      enddo ! lk
      var = var / dble(nRec) ! was (nRec-1) <- 2001-12-06 ceaf
      if (var.eq.0.) then ! no variance -> perfect correlation
       ll = 1 ! initial write index
       do k=lmin,lmax,lstp
        lag(ll) = k
        acf(ll) = 1. ! force correlation coefficient to 1
        ll = ll + 1
       enddo ! k
       return ! no need to calculate below, so return here
      endif

! --- Center time series (avg=0) ---
      do lk=ibeg,iend
       TS(lk) = TS(lk) - avg
      enddo ! lk

! --- Fix indices, mapping ibeg&iend to local variables ---
      begidx = ibeg
      endidx = iend

! --- Estimate autocorrelation function ---
      ll = 1 ! initial write index
      do k=lmin,lmax,lstp ! was lagstep (undefined) <- 2001-12-06 ceaf
       ck = 0.
       do lk=begidx,(endidx-k)
        ck = ck + TS(lk) * TS(lk+k)
       enddo ! lk
       lag(ll) = k
       acf(ll) = (ck / dble(endidx - k - begidx + 1)) / var
       ll = ll + 1
      enddo ! k

      return
      end

!***********************************************************************
!* ID => cross_correlation                                             *
!* Function => estimates crosscorrelation function of two time series  *
!* Date => 2001-01-04                                                  *
!* Last modified => 2001-12-10 01:59Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-10) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 1994-1996 by C.E.A. Finney.                           *
!* Placed in the PUBLIC DOMAIN 2001-01-04 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TSx => array with reference time-series data, dimensioned 1:*     *
!*   TSy => array with test time-series data, dimensioned 1:*          *
!*   ibeg => first record in TS segment                                *
!*   iend => last record in TS segment                                 *
!*   lmin => minimal value of lag for computation of crosscorrelation  *
!*             function; may be negative                               *
!*   lmax => maximal value of lag for computation of crosscorrelation  *
!*             function; may be negative                               *
!*   lstp => lag increment                                             *
!* Output parameters:                                                  *
!*   lag => lags between lmin and lmax step lstp, dimensioned 1:*      *
!*   ccf => crosscorrelation values at each lag contained in lag,      *
!*          dimensioned 1:*                                            *
!***********************************************************************
!* Changes:                                                            *
!*   2001-02-25 ceaf => changed TS to real(kind=4), retained all others as   *
!*              real(kind=8); renamed variable first->ibeg, last->iend,      *
!*              minlag->lmin, maxlag->lmax; added parameter lstp;      *
!*              split ccf into lag and ccf                             *
!*   2001-02-28 ceaf => fixed ck denominator (added +1)                *
!*   2001-03-20 ceaf => changed double to dble                         *
!*   2001-12-06 sp+ceaf => added "implicit none"; added lstp to input  *
!*              argument list; changed TS to real(kind=8); changed           *
!*              denominator in variance calculation from nRec-1 to     *
!*              nRec; changed 1 to ibeg in assignment of begidx;       *
!*              removed **2 from scaling of ck by number of differences*
!***********************************************************************

      subroutine cross_correlation(TSx,TSy,ibeg,iend,lmin,lmax,lstp,lag,ccf)

      implicit none

      real(kind=8) TSx(1:*) !................... reference time series (input)
      real(kind=8) TSy(1:*) !........................ test time series (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      integer(kind=4) lmin !.......... minimal value of lag to compute (input)
      integer(kind=4) lmax !.......... maximal value of lag to compute (input)
      integer(kind=4) lstp !............................ lag increment (input)
      integer(kind=4) lag(1:*) !........................ list of lags (output)
      real(kind=8) ccf(1:*) !................ crosscorrelation values (output)
      integer(kind=4) lk !............................................ counter
      integer(kind=4) ll !...................................... write counter
      integer(kind=4) k !......................................... lag counter
      integer(kind=4) nrec !................. number of records in [ibeg,iend]
      real(kind=8) avgx !......................................... mean of TSx
      real(kind=8) avgy !......................................... mean of TSy
      real(kind=8) ck !................................. lag-k correlation sum
      real(kind=8) varx !..................................... variance of TSx
      real(kind=8) vary !..................................... variance of TSy
      real(kind=8) scale !..................... scale factor (sqrt(varx*vary))
      real(kind=8) sum !.................................................. sum
      integer(kind=4) begidx !...................... first record in TS to use
      integer(kind=4) endidx !....................... last record in TS to use

! --- Sanity and range checks ---
      nrec = iend - ibeg + 1
      if (nrec.le.0) then ! invalid ranges ? abort
       return
      endif
      if ((lmin.gt.lmax).and.(lstp.gt.0)) then
       lstp = -lstp ! decrement through lag loop
      endif

! --- Determine data mean & variance ---
! ... Reference time series ...
      sum = 0.
      do lk=ibeg,iend
       sum = sum + TSx(lk)
      enddo ! lk
      avgx = sum / dble(nrec)
      varx = 0.
      do lk=ibeg,iend
       varx = varx + (TSx(lk) - avgx)**2
      enddo ! lk
      varx = varx / dble(nrec) ! was (nrec-1) <- 2001-12-06 ceaf
! ... Test time series ...
      sum = 0.
      do lk=ibeg,iend
       sum = sum + TSy(lk)
      enddo ! lk
      avgy = sum / dble(nrec) ! was (nrec-1) <- 2001-12-06 ceaf
      vary = 0.
      do lk=ibeg,iend
       vary = vary + (TSy(lk) - avgy)**2
      enddo ! lk
      vary = vary / dble(nrec)
      scale = (varx * vary) ** 0.5
! !!! Add: if (scale .lt. eps) return !!!
!     (avoid overflow for near-zero variance)

! --- Center time series (avg=0) ---
      do lk=ibeg,iend
       TSx(lk) = TSx(lk) - avgx
       TSy(lk) = TSy(lk) - avgy
      enddo ! lk

! --- Fix indices, mapping ibeg&iend to local variables ---
      begidx = max(ibeg,ibeg-lmin) ! changed 1->ibeg <- 2001-12-06 ceaf
      endidx = min(iend,iend-lmax)

! --- Estimate crosscorrelation function ---
      ll = 1 ! initial write index
      do k=lmin,lmax,lstp
       ck = 0.
       do lk=begidx,(endidx-k)
        ck = ck + TSx(lk) * TSy(lk+k)
       enddo ! lk
       lag(ll) = k
       ! The following line was
       !   ccf(ll) = (ck / (dble(endidx - k - begidx + 1)**2)) / scale
       ! and gave wrong scaling. <- 2001-12-06 ceaf
       ccf(ll) = ck / (dble(endidx - k - begidx + 1) * scale)
       ll = ll + 1
      enddo ! k

      return
      end

!***********************************************************************
!* ID => mif1                                                          *
!* Function => estimates mutual information function of a single       *
!*   "symbolized" time series                                          *
!* Date => 2001-01-04                                                  *
!* Last modified => 2001-12-20 21:45Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-10) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney.                                *
!* Placed in the PUBLIC DOMAIN 2001-01-04 for distribution with MFIX.  *
!***********************************************************************
!* Mutual information is typically reported in bits.  This function    *
!* uses a normalization scheme after: Darbellay GA (1998).             *
!* "Predictability: an information-theoretic perspective".  Chapter    *
!* 18 in _Signal Analysis and Prediction_, Prochazka, Uhlir, Rayner,   *
!* Kingsbury, eds.  Boston: Birkhauser, ISBN 0-8176-4042-8.            *
!***********************************************************************
!* Input parameters:                                                   *
!*   SS => array with "symbolized" time-series data, dimensioned 1:*;  *
!*         symbols should be in the set of [0,nbin-1]                  *
!*   ibeg => index of first record in SS segment                       *
!*   iend => index of last record in SS segment                        *
!*   lmin => minimal value of lag for computation of MI function       *
!*   lmax => maximal value of lag for computation of MI function       *
!*   lstp => lag increment                                             *
!*   nbin => number of bins for the observed probabilities (histograms)*
!* Output parameters:                                                  *
!*   lag => lags between lmin and lmax step lstp, dimensioned 1:*      *
!*   mif => mutual-information values at each lag in lag,              *
!*          dimensioned 1:*                                            *
!***********************************************************************
!* Changes:                                                            *
!*   2001-02-14 ceaf => changed mif from two to one-dimensional array  *
!*              and stored lags in separate array; changed real(kind=8) to   *
!*              real(kind=4); changed output units from bits to dimensionless*
!*              after Darbellay; externalized "symbolization"          *
!*   2001-02-15 ceaf => changed probability matrices to be indexed     *
!*              from 0                                                 *
!*   2001-03-01 ceaf => added/changed documentation                    *
!*   2001-12-10 ceaf => added "implicit none"                          *
!*   2001-12-20 ceaf => added range check for nbin>MAXBIN; added write *
!*              statements for range-check errors                      *
!***********************************************************************

      subroutine mif1(SS,ibeg,iend,lmin,lmax,lstp,nbin,lag,mif)

      implicit none

      integer(kind=4) MAXBIN
      parameter (MAXBIN=100)

      integer(kind=4) SS(1:*) !................ symbolized time series (input)
      integer(kind=4) ibeg !....... index of first record in SS to use (input)
      integer(kind=4) iend !........ index of last record in SS to use (input)
      integer(kind=4) lmin !.......... minimal value of lag to compute (input)
      integer(kind=4) lmax !.......... maximal value of lag to com ute (input)
      integer(kind=4) lstp !.......................... lag incrementer (input)
      integer(kind=4) nbin ! number of bins for observed probabilities (input)
      integer(kind=4) lag(1:*) !........................ list of lags (output)
      real(kind=4) mif(1:*) !.............. mutual information values (output)
      real(kind=4) mi !.................................... mutual information
      integer(kind=4) i,j,k,l,ll !................................... counters
      integer(kind=4) begidx !...................... first record in SS to use
      integer(kind=4) endidx !....................... last record in SS to use
      real(kind=4) pi(0:MAXBIN-1) !.. observed probabilities of leading points
      real(kind=4) pj(0:MAXBIN-1) !. observed probabilities of trailing points
      real(kind=4) pij(0:MAXBIN-1,0:MAXBIN-1) !.. joint observed probabilities
      real(kind=4) hi,hj,hij !......................... entropies of pi,pj,pij

! --- Sanity and range checks ---
      if ((iend-ibeg+1).le.0) then ! invalid ranges ? abort
       write(*,*) 'ornl.mif1: FATAL: iend < ibeg.'
       return
      endif
      if ((lmin.lt.0).or.(lmax.lt.0)) then
       write(*,*) 'ornl.mif1: FATAL: lags may not be negative.'
       return ! no negative lags allowed
      endif
      if ((lmin.gt.lmax).and.(lstp.gt.0)) then
       lstp = -lstp ! decrement through lag loop
      endif
      if (nbin.gt.MAXBIN) then
       write(*,*) 'ornl.mif1: FATAL: nbin > MAXBIN.'
       return
      endif

! --- Fix indices, mapping ibeg&iend to local variables ---
      begidx = ibeg
      endidx = iend

! --- Estimate mutual information function ---
      ll = 1 ! initial write index
      do k=lmin,lmax,lstp
! .... Initialize probability matrices
       do i=0,(nbin-1)
        pi(i) = 0.
        pj(i) = 0.
        do j=0,(nbin-1)
         pij(i,j) = 0.
        enddo ! j
       enddo ! i
! .... Leading-point probability
       do l=begidx,endidx-k
        pi(SS(l)) = pi(SS(l)) + 1.
       enddo ! l
! .... Trailing-point probability
       do l=begidx+k,endidx
        pj(SS(l)) = pj(SS(l)) + 1.
       enddo ! l
! .... Joint probability
       do l=begidx,endidx-k
        pij(SS(l),SS(l+k)) = pij(SS(l),SS(l+k)) + 1.
       enddo ! l
! .... Normalize the probabilities
       do i=0,(nbin-1)
        pi(i) = pi(i) / (endidx - k - begidx + 1)
        pj(i) = pj(i) / (endidx - k - begidx + 1)
        do j=0,(nbin-1)
         pij(i,j) = pij(i,j) / (endidx - k - begidx + 1)
        enddo ! j
       enddo ! i
! .... Calculate the entropies
       hi = 0.
       hj = 0.
       hij = 0.
       do i=0,(nbin-1)
        if (pi(i).gt.0.) hi = hi - pi(i) * log(pi(i))
        if (pj(i).gt.0.) hj = hj - pj(i) * log(pj(i))
        do j=0,(nbin-1)
         if (pij(i,j).gt.0.) hij = hij - pij(i,j) * log(pij(i,j))
        enddo ! j
       enddo ! i
       mi = hi + hj - hij ! mutual information in nats
       lag(ll) = k
       mif(ll) = sqrt(1 - exp(-2 * mi)) ! after Darbellay
       ll = ll + 1
      enddo ! k

      return
      end

!***********************************************************************
!* ID => t3sym                                                         *
!* Function => cubed lagged-difference temporal asymmetry              *
!* Date => 2000-12-21                                                  *
!* Last modified => 2001-12-10 01:55Z                                  *
!***********************************************************************
!* Programmer => C.S. Daw (Oak Ridge National Laboratory)              *
!* Maintainer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-10) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2000 by C.S. Daw.                                     *
!* Placed in the PUBLIC DOMAIN 2001-01-07 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   ibeg => first record in TS segment                                *
!*   iend => last record in TS segment                                 *
!*   lmin => minimal value of lag for computation of t3 function       *
!*   lmax => maximal value of lag for computation of t3 function       *
!*   lstp => lag increment                                             *
!* Output parameters:                                                  *
!*   lag => list of lags in lmin,lmax,lstp                             *
!*   tsym => t3 functional value at each lag in lag                    *
!***********************************************************************
!* Changes:                                                            *
!*   2001-12-10 ceaf => added "implicit none"; changed TS to real(kind=8)    *
!***********************************************************************

      subroutine t3sym(TS,ibeg,iend,lmin,lmax,lstp,lag,tsym)

      implicit none

      real(kind=8) TS(1:*) !.............................. time series (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      integer(kind=4) lmin !.......... minimal value of lag to compute (input)
      integer(kind=4) lmax !.......... maximal value of lag to compute (input)
      integer(kind=4) lstp !............................ lag increment (input)
      real(kind=8) lag(1:*) !......... time lags matching tsym values (output)
      real(kind=8) tsym(1:*) !........... asymmetry for all intervals (output)
      real(kind=8) d !...................................... lagged difference
      real(kind=8) s2, s3 !......................... lagged-difference summers
      integer(kind=4) i,j,k !........................................ counters

      if (lstp.eq.0) lstp = 1
      k = 1 ! initialize write index
      do i=lmin,lmax,lstp
       s2 = 0.d0
       s3 = 0.d0
       do j=ibeg,iend-i,1
        d = TS(j+i) - TS(j)
        s2 = s2 + d**2
        s3 = s3 + d**3
       enddo ! j
       lag(k) = i
       tsym(k) = 0.d0
       if (s2.gt.0.) then ! trap division by 0 for perfectly periodic
        tsym(k) = s3 / (s2**1.5)
       endif
       k = k + 1
      enddo ! i

      return
      end

!***********************************************************************
!* ID => trevsgn                                                       *
!* Function => signed lagged-difference temporal asymmetry             *
!* Date => 2001-12-07                                                  *
!* Last modified => 2001-12-10 01:50Z                                  *
!***********************************************************************
!* Programmer => C.E.A. Finney (Oak Ridge National Laboratory)         *
!* Contact information (as of 2001-12-07) :                            *
!*   Email => <mfix@ceafinney.com>                                     *
!*   Telephone => 865-946-1243                                         *
!*   Post => National Transportation Research Center                   *
!*           2360 Cherahala Boulevard                                  *
!*           Knoxville  TN  37932-6472  USA                            *
!***********************************************************************
!* Copyright (C) 2001 by C.E.A. Finney                                 *
!* Placed in the PUBLIC DOMAIN 2001-12-07 for distribution with MFIX.  *
!***********************************************************************
!* Input parameters:                                                   *
!*   TS => array with time-series data, dimensioned 1:*                *
!*   ibeg => index of first record in TS segment                       *
!*   iend => index of last record in TS segment                        *
!*   lmin => minimal value of lag for computation of tsgn function     *
!*   lmax => maximal value of lag for computation of tsgn function     *
!*   lstp => lag increment                                             *
!* Output parameters:                                                  *
!*   lag => list of lags in lmin,lmax,lstp                             *
!*   tsgn => tsgn value at each lag in lag                             *
!***********************************************************************
!*   2001-12-10 ceaf => added "implicit none" and "external sgn_d"     *
!*   2001-12-20 ceaf => added declaration for sgn_d                    *
!***********************************************************************

      subroutine trevsgn(TS,ibeg,iend,lmin,lmax,lstp,lag,tsgn)

      implicit none

      external sgn_d
      real(kind=8) sgn_d !........................... algebraic sign of double

      real(kind=8) TS(1:*) !.............................. time series (input)
      integer(kind=4) ibeg !................ first record in TS to use (input)
      integer(kind=4) iend !................. last record in TS to use (input)
      integer(kind=4) lmin !.......... minimal value of lag to compute (input)
      integer(kind=4) lmax !.......... maximal value of lag to compute (input)
      integer(kind=4) lstp !.......................... lag incrementer (input)
      integer(kind=4) lag(1:*) !...... time lags matching tsgn values (output)
      real(kind=8) tsgn(1:*) !........... asymmetry for all intervals (output)
      real(kind=8) sum !............................. lagged-difference summer
      integer(kind=4) j,k,l !........................................ counters
      integer(kind=4) begidx !...................... first record in SS to use
      integer(kind=4) endidx !....................... last record in SS to use

! --- Sanity and range checks ---
      if ((iend-ibeg+1).le.0) then ! invalid ranges ? abort
       return
      endif
      if ((lmin.lt.0).or.(lmax.lt.0)) then
       return ! no negative lags allowed
      endif
      if (lstp.eq.0) lstp = 1 ! prevent infinite loop
      if ((lmin.gt.lmax).and.(lstp.gt.0)) then
       lstp = -lstp ! decrement through lag loop
      endif

! --- Fix indices, mapping ibeg&iend to local variables ---
      begidx = ibeg
      endidx = iend

! --- Estimate temporal-symmetry function ---
      k = 1 ! initialize write index
      do l=lmin,lmax,lstp
       sum = 0.d0
       do j=begidx,(endidx-l),1
        sum = sum + sgn_d(TS(j+l) - TS(j))
       enddo ! j
       lag(k) = l
       tsgn(k) = sum / dble(endidx - begidx - l + 1)
       k = k + 1
      enddo ! l

      return
      end

