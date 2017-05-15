!***********************************************************************
!* ID => ornl_pca                                                      *
!* Function => principal-components analysis routines                  *
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
!* matdiag                      2001-01-24 06:50Z                      *
!* neigen                       2001-12-20 22:29Z                      *
!* traeigen                     2001-12-20 22:30Z                      *
!***********************************************************************

!***********************************************************************
!* ID => matdiag                                                       *
!* Function => diagnolizes a real symmetric matrix using method of     *
!*   Jacobi as described by J. Greenstadt in "Mathematical Methods for *
!*   Digital Computers", Ralston and Wilf, eds., John Wiley & Sons,    *
!*   1960, chapter 7.  The diagonals of the returned matrix A are the  *
!*   eigenvalues of the original matrix and the matrix eigvec contains *
!*   the corresponding eigenvectors.                                   *
!* Date => 2001-01-09                                                  *
!* Last modified => 2001-01-24 06:50Z ceaf                             *
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
!*   A => real, square symmetric matrix                                *
!*   N => dimension of A                                               *
!* Output parameters:                                                  *
!*   B => diagonalized version of A                                    *
!*   eigvec => eigenvectors                                            *
!***********************************************************************

      subroutine matdiag(A,N,B,eigvec)

      implicit none

      integer(kind=4) MAXDIM !.................... maximal embedding dimension
      parameter (MAXDIM=100)

! ... NOTICE: Keep all reals at least double precision. ...
      real(kind=8) sgn_d !............. external SGN function (algebraic sign)
      real(kind=8) A(MAXDIM,MAXDIM) !........... real symmetric matrix (input)
      integer N !.............................. matrix dimension (input)
      real(kind=8) B(MAXDIM,MAXDIM) !..... transformation of matrix A (output)
      real(kind=8) eigvec(MAXDIM,MAXDIM) !................ eigenvectors (output)
      real(kind=8) RHO !................ small value to represent off-diagonal
      real(kind=8) nu, nufinal !........ dummy values used in norm computation
      integer(kind=4) i,j,q,p !.............................. counting indices
      integer(kind=4) offdiag !...... flag to signal significant off-diagonals
! ... Intermediate quantities in computation ...
      real(kind=8) app,aqq,apq,lambda,mu,temp,sine,cosine,omega
      real(kind=8) co2,si2,sico,sico2

      RHO = 1.0d-16

! --- Generate an identity matrix having off-diagonal elements = 0
!     and diagonal elements = 1.0 ---
      do i=1,N
       do j=1,N
        if (i.ne.j) then
         eigvec(i,j) = 0.
        else
         eigvec(i,j) = 1.
        endif
       enddo ! j
      enddo ! i

! --- Compute initial norm of matrix A
      nu = 0.
      do i = 1,N
       do j = 1,N
        if (i.ne.j) nu = nu + a(i,j) * a(i,j)
       enddo ! j
      enddo ! i
      nu = sqrt(nu)

      nufinal = nu * RHO / N

! --- Begin off-diagonal elimination loop.  Continue loop until
!     nu <= nufinal and no off-diagonals found.
  40  continue ! outer loop
      nu = nu / N
      do q = 2,N
      p = 1
  45  continue ! inner loop until p > q-1
      offdiag = 0 ! initialize off-diagonal flag
      if (abs(a(p,q)).gt.nu) then ! if a
       offdiag = 1
       ! save elements of pivotal set
       app = a(p, p)
       aqq = a(q, q)
       apq = a(p, q)
       lambda = -apq
       mu = (app - aqq) / 2.
       if ((mu.eq.0.).or.(mu.lt.RHO)) then
        omega = sgn_d(lambda)
       else
        omega = lambda / sqrt(lambda * lambda + mu * mu)
        omega = omega * sgn_d(mu)
       endif
       sine = omega / sqrt(2 * (1 + sqrt(1 - omega ** 2)))
       cosine = sqrt(1 - sine ** 2)

       do i=1,N
        temp = a(i,p) * cosine - a(i,q) * sine
        a(i,q) = a(i,p) * sine + a(i,q) * cosine
        a(i,p) = temp
        temp = eigvec(i, p) * cosine - eigvec(i,q) * sine
        eigvec(i,q) = eigvec(i,p) * sine + eigvec(i,q) * cosine
        eigvec(i,p) = temp
       enddo ! i

! ...  now do the diagonal elements ...
       co2 = cosine ** 2
       si2 = sine ** 2
       sico = sine * cosine
       sico2 = 2 * sico * apq
       a(p, p) = app * co2 + aqq * si2 - sico2
       a(q, q) = app * si2 + aqq * co2 + sico2
       a(p, q) = (app - aqq) * sico + (co2 - si2) * apq
       a(q, p) = a(p, q)
       do i=1,N
        a(p,i) = a(i,p)
        a(q,i) = a(i,q)
       enddo ! i
      endif ! if a
      p = p + 1
      if (p.le.(q-1)) goto 45  ! inner loop until p > (q - 1)
      enddo ! q

      if ((nu.le.nufinal).and.offdiag.eq.0) then ! outer loop
       goto 41  ! exit loop
      else
       goto 40  ! continue loop
      endif
  41  continue

! --- Arrange eigenvalues and eigenvectors in decreasing order ---
      do i=2,N
       if (a(i,i).gt.a(i-1,i-1)) then
        temp = a(i,i)
        a(i,i) = a(i-1,i-1)
        a(i-1,i-1) = temp
        do j=1,N
         temp = eigvec(j,i)
         eigvec(j,i) = eigvec(j,i-1)
         eigvec(j,i-1) = temp
        enddo ! j
       endif
      enddo ! i

      return

      end

!***********************************************************************
!* ID => neigen                                                        *
!*   (The name "neigen" is historical, when "normalized eigenvalues"   *
!*   (standardized by the largest magnitude) were returned.)           *
!* Function => principal-components analysis of a time series          *
!* Date => 2001-01-07                                                  *
!* Last modified => 2001-12-20 22:29Z ceaf                             *
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
!*   TS => array with time-series data, dimensioned 1:*                *
!*   ibeg => temporal index for starting analysis                      *
!*   iend => temporal index for stopping analysis                      *
!*   embdim => embedding dimension                                     *
!*   embdel => embedding delay                                         *
!* Output parameters:                                                  *
!*   eigvec => eigenvectors                                            *
!*   eigval => eigenvalues                                             *
!***********************************************************************
!* Changes:                                                            *
!*   2001-03-06 ceaf => changed TS to real*4, retained all others as   *
!*              real(kind=8); renamed variable first->ibeg, last->iend       *
!*   2001-12-20 ceaf => changed TS to real(kind=8)                           *
!***********************************************************************

      subroutine neigen(TS,ibeg,iend,embdim,embdel,eigvec,eigval)

      implicit none
      external matdiag !................. matrix diagonalization routine

      integer(kind=4) MAXDIM !.......... maximal allowable embedding dimension
      parameter (MAXDIM=100) ! used for eigvec,eigval,covin,covout
      integer(kind=4) MAXREC !........... maximal allowable time-series length
      parameter (MAXREC=100000) ! used for TSc

      real(kind=8) TS(*) !................................ time series (input)
      integer(kind=4) embdim !.................... embedding dimension (input)
      integer(kind=4) embdel !........................ embedding delay (input)
      integer(kind=4) ibeg !..... temporal index for starting analysis (input)
      integer(kind=4) iend !..... temporal index for stopping analysis (input)
      real(kind=8) eigvec(MAXDIM,MAXDIM) !.............. eigenvectors (output)
      real(kind=8) eigval(MAXDIM) !........... eigenvalues (sigma**2) (output)
      integer(kind=4) i,j,k !......................................... indices
      integer(kind=4) npact !.......... actual number of data points processed
      real(kind=8) avg !................................................. mean
      real(kind=8) TSc(MAXREC) !......................... centered time series
      integer(kind=4) nwin !............. number of embedded trajectory points
! ... NOTICE: Keep the covariance matrices at least real(kind=8).
      real(kind=8) covin(MAXDIM,MAXDIM) !............ covariance matrix, input
      real(kind=8) covout(MAXDIM,MAXDIM) !.......... covariance matrix, output
      integer(kind=4) ipoint,index1,index2 !......................... pointers

! --- Center time series ---
      npact = iend - ibeg + 1
      avg = 0.
      do i=ibeg,iend
       avg = avg + TS(i)
      enddo ! i
      avg = avg / npact
      do j=ibeg,iend
       TSc(j - ibeg + 1) = TS(j) - avg
      enddo ! j

! --- Normalize covariance matrix ---
      do i=1,embdim
       do j=1,embdim
        covin(i,j) = 0.
       enddo ! j
      enddo ! i

! --- Compute covariance matrix ---
      nwin = npact - ((embdim - 1) * embdel)
      do i=1,nwin
       ipoint = ibeg + i - 1
       do j=1,embdim
        do k=j,embdim
         index1 = ipoint + (j - 1) * embdel
         index2 = ipoint + (k - 1) * embdel
         covin(j,k) = covin(j,k) + TSc(index1) * TSc(index2)
        enddo ! k
        if (k.ne.j) covin(k,j) = covin(j,k) ! force symmetry
       enddo ! j
      enddo ! i

      do j=1,embdim
       do k=j,embdim
        covin(j,k) = covin(j,k) / (nwin - 1) ! normalize covariance
        if (k.ne.j) covin(k,j) = covin(j,k) ! force symmetry
       enddo ! k
      enddo ! j

! --- Diagonalize the covariance matrix ---
      call matdiag(covin,embdim,covout,eigvec)

! --- Save eigenvalues (on diagonal) in vector ---
      do i=1,embdim
       eigval(i) = covout(i,i)
      enddo ! i

      return

      end

!***********************************************************************
!* ID => traeigen                                                      *
!* Function => transforms time series into embedded trajectory and     *
!*             saves selected component                                *
!* Date => 2001-01-07                                                  *
!* Last modified => 2001-12-20 22:30Z ceaf                             *
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
!*   TS => array with time-series data, dimensioned 1:*                *
!*      NOTE: TS must have zero mean                                   *
!*   ibeg => temporal index for starting analysis                      *
!*   iend => temporal index for stopping analysis                      *
!*   embdim => embedding dimension                                     *
!*   embdel => embedding delay                                         *
!*   eigvec => eigenvectors                                            *
!*   comp => index of principal component to output in TSc             *
!* Output parameters:                                                  *
!*   TSc => time series of principal component number comp             *
!***********************************************************************
!* Changes:                                                            *
!*   2001-03-06 ceaf => changed TS to real*4, retained all others as   *
!*              real(kind=8); renamed variable first->ibeg, last->iend       *
!*   2001-03-20 ceaf => incorporated centering of TS for calling ease  *
!*   2001-12-20 ceaf => changed TS and TSc to real(kind=8)                   *
!***********************************************************************

      subroutine traeigen(TS,ibeg,iend,embdim,embdel,eigvec,comp,TSc)

      implicit none

      integer(kind=4) MAXDIM !.......... maximal allowable embedding dimension
      parameter (MAXDIM=100)
      integer(kind=4) MAXREC !........... maximal allowable time-series length
      parameter (MAXREC=100000)

      real(kind=8) TS(*) !................................ time series (input)
      integer(kind=4) ibeg !..... temporal index for starting analysis (input)
      integer(kind=4) iend !..... temporal index for stopping analysis (input)
      integer(kind=4) embdim !.................... embedding dimension (input)
      integer(kind=4) embdel !........................ embedding delay (input)
      real(kind=8) eigvec(MAXDIM,MAXDIM) !............... eigenvectors (input)
      integer(kind=4) comp !..... selected principal component for TSc (input)
      real(kind=8) TSc(MAXREC) !..... output principal component comp (output)
      integer(kind=4) itpts !... number of trajectory component points to save
      integer(kind=4) tidx !. temporal index in TS of current embedding window
      integer(kind=4) icount,j !............................. counting indices
      real(kind=8) avg !..................................... time-series mean

! --- Determine number of embedded data points ---
      itpts = (iend - ibeg + 1) - (embdim - 1) * embdel

! --- Calculate time-series mean ---
      avg = 0.
      do j=ibeg,iend
       avg = avg + TS(j)
      enddo ! j
      avg = avg / dble(iend - ibeg + 1)

! --- Transform time series into PC trajectory ---
      do icount=1,itpts
       tidx = ibeg + icount - 1
       TSc(icount) = 0.
       do j=1,embdim  ! multiply input vector times eigvec()
        TSc(icount) = TSc(icount) + (TS((j-1) * embdel + tidx) - avg)   &
     &                * eigvec(j,comp)
       enddo ! j
      enddo ! icount

      return
      end

