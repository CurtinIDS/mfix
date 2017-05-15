!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: call mass_mobility(s,mi,ni,rho,zeta0,theta,nu,DF)
!
!  author:  C. Hrenya, Feb 2009
!
!  Purpose: find mass mobility coefficient according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine mass_mobility(s,mi,ni,rho,zeta0,theta,nu,DF)
      Implicit NONE

      integer s, indx(s)

      double precision mi(s),ni(s),rho, &
                      zeta0,theta(s),nu(s,s),DF(s,s)

      integer i,j,kk
      double precision kronecker(s,s),Amat(s,s),bmat(s,s), &
                      Amat0(s,s),bmat0(s)
      double precision d

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve

      do i=1,s
         do j=1,s
            if (i.eq.j) then
               kronecker(i,j) = 1.d0
            else
               kronecker(i,j) = 0.d0
            endif
         enddo
      enddo

      do i=1,s
         do j=1,s
             Amat(i,j) = (nu(i,j)+0.5d0*zeta0*kronecker(i,j))    !A matrix for solution of DF (p 8 CMH notes)
             bmat(i,j) = -ni(i)*mi(i)/mi(j)*(kronecker(i,j)-ni(j) &
                *mi(j)/rho)                                      !b matrix for solution of DF (p 8 CMH notes)
         enddo
      enddo

! this extra kk loop and addition of Amat0 and bmat0 is necessary
! since x & b in Ax=b are s by s matrices rather than vectors of length s,
! whereas LUBSKB is specific to x & b vectors of length s
      do kk=1,s
         do i=1,s
            do j=1,s
                Amat0(i,j) = Amat(i,j)
            enddo
         enddo
         do i=1,s
            bmat0(i) = bmat(i,kk)
         enddo

         CALL LUDCMP(Amat0, s, NP, indx, d, 'mass_mobility') ! solve system of s linear equations using
         CALL LUBKSB(Amat0, s, NP, indx, bmat0) ! LU decomposition

         do i=1,s
            DF(i,kk) = bmat0(i)
         enddo

      enddo

      return
      end subroutine mass_mobility
