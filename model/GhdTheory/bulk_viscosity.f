!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: bulk_viscosity(s,mi,alpha,ni,v0,mu,sigma,chi,theta,kappa)
!
!  author:  C. Hrenya, Jan 2009
!
!  Purpose: find bulk viscosity according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine bulk_viscosity(s,mi,alpha,ni,v0,mu,sigma,chi, &
                                theta,kappa)
      Implicit NONE

      integer s

      double precision mi(s),alpha(s,s),ni(s),v0,mu(s,s),sigma(s,s), &
                       chi(s,s),theta(s),kappa

      integer i,j

      double precision pi
      parameter (pi=3.14159265458979323846d0)

      kappa = 0.d0
      do i=1,s
         do j=1,s
            kappa = kappa + mu(i,j)*mi(j)*ni(i)*ni(j)*v0* &
                      sigma(i,j)**4*chi(i,j)*(1.d0+alpha(i,j)) &
                  *dsqrt((theta(i)+theta(j))/(theta(i)*theta(j)))
         enddo
      enddo
      kappa = 2.d0*dsqrt(pi)/9.d0*kappa        !p 12 CMH notes

      return
      end subroutine bulk_viscosity

