!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: thermal_diffusivity(s,alpha,ni,mi,rho,v0,mu,sigma,chi,
!                                       zeta0,theta,Ti,p,DT,nu)
!
!  author:  C. Hrenya, Jan 2009
!
!  Purpose: find thermal diffusivity according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine thermal_diffusivity(s,alpha,ni,mi,rho,v0,mu,sigma, &
                                     chi,zeta0,theta,Ti,p,DT,nu)

      Implicit NONE

      integer s, indx(s)

      double precision alpha(s,s),ni(s),mi(s),rho,v0,mu(s,s), &
                       sigma(s,s),chi(s,s),zeta0,theta(s), &
                       Ti(s),p,DT(s)

      integer i,j
      double precision sum1(s),sum2(s),nu(s,s),Amat(s,s),bmat(s)
      double precision d

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve
      double precision pi
      parameter (pi=3.14159265458979323846d0)

      do i=1,s
         sum1(i) = 0.d0     !calculate summation used in nu(i,i) - p 7 CMH notes
         sum2(i) = 0.d0     !calculate summation used in b vector - p 7 CMH notes
      enddo

      do i=1,s
         do j=1,s
            if (i .ne. j) then
               sum1(i) = sum1(i) + ni(j)*sigma(i,j)**2*chi(i,j) &
                    *mu(j,i)*v0*(1.d0+alpha(i,j)) &
                    *dsqrt((theta(i)+theta(j))/(theta(i)*theta(j)))
            endif
            sum2(i) = sum2(i) + ni(j)*mu(i,j)*chi(i,j)* &
                      sigma(i,j)**3*Ti(j)*(1.d0+alpha(i,j))
         enddo
      enddo

      do i=1,s
         do j=1,s
            if (i .eq. j) then
               nu(i,i) = 4.d0*dsqrt(pi)/3.d0*sum1(i)

               Amat(i,j) = nu(i,j)-zeta0                    !A matrix for solution of DT (p 7 CMH notes)
            else
               nu(i,j) = -4.d0*dsqrt(pi)/3.d0*ni(i)*sigma(i,j)**2 &
                  *chi(i,j)*mu(i,j)*v0*(1.d0+alpha(i,j))   &
                  *dsqrt((theta(i)+theta(j))/(theta(i)*theta(j)))

               Amat(i,j) = nu(i,j)                          !A matrix for solution of DT (p 7 CMH notes)
            endif
         enddo
         bmat(i) = -p*ni(i)*mi(i)/rho**2*(1.d0-rho*Ti(i) &  !b matrix for solution of DT (p 7 CMH notes)
                      /(mi(i)*p))+2d0*pi/3d0*ni(i)/rho*sum2(i)
      enddo

      CALL LUDCMP(Amat, s, NP, indx, d, 'thermal_diffusivity')     ! solve system of s linear equations using
      CALL LUBKSB(Amat, s, NP, indx, bmat)  ! LU decomposition

      do i=1,s
         DT(i) = bmat(i)
      enddo

      return
      end subroutine thermal_diffusivity

