!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: bulk_viscosity(s,mi,sigmai,alpha,ni,v0,mu,sigma,chi,
!                                  beta,zeta0,theta,Ti,kappa,eta)
!
!  author:  C. Hrenya, Jan 2009
!
!  Purpose: find shear viscosity according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine shear_viscosity(s,mi,sigmai,alpha,ni,v0,mu,sigma,chi, &
                                 beta,zeta0,theta,Ti,kappa,eta)
      Implicit NONE

      integer s, indx(s)

      double precision mi(s),sigmai(s),alpha(s,s),ni(s),v0,mu(s,s), &
                       sigma(s,s),chi(s,s),beta(s,s),zeta0,theta(s), &
                       Ti(s),kappa,eta

      integer i,j
      double precision sum1(s),sum2(s),tau(s,s),Amat(s,s),bmat(s), &
                       etajk(s),etakin,etacol
      double precision d

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve
      double precision pi
      parameter (pi=3.14159265458979323846d0)

      do i=1,s
         sum1(i) = 0.d0     !calculate summation used in tau(i,i) - p 10 CMH notes
         sum2(i) = 0.d0     !calculate summation used in b vector - p 10 CMH notes
      enddo

      do i=1,s
         do j=1,s
            if (i .ne. j) then
               sum1(i) = sum1(i) + ni(j)*chi(i,j)*sigma(i,j)**2* &
                 mu(j,i)*(1.d0+alpha(i,j))*theta(i)**1.5d0 &
                 /dsqrt(theta(j)) * (6.d0*beta(i,j)/ &
                 (theta(i)**2*dsqrt(theta(i)+theta(j)))  &
                 + (9.d0-3.d0*alpha(i,j))/2.d0*mu(j,i)/theta(i)**2* &
                 dsqrt(theta(i)+theta(j))   &
                 + 5.d0/(theta(i)*dsqrt(theta(i)+theta(j))))
           endif
           sum2(i) = sum2(i) + 2.d0*pi/15.d0*mi(i)*ni(i)*ni(j)*mu(j,i)* &
             sigma(i,j)**3*chi(i,j)*(1d0+alpha(i,j)) *(mu(j,i)* &
            (3.d0*alpha(i,j)-1.d0)*(Ti(i)/mi(i)+Ti(j)/mi(j)) &
            -4.d0*(Ti(i)-Ti(j))/(mi(i)+mi(j)))
         enddo
      enddo

! calculate tau (p 10-11 CMH notes)
      do i=1,s
         do j=1,s
            if (i .eq. j) then
               tau(i,i) = 4.d0*dsqrt(pi)/15.d0*v0 *(ni(i)*sigmai(i)**2* &
                 chi(i,i)/dsqrt(2.d0*theta(i))*(9.d0-3.d0*alpha(i,i))* &
                 (1.d0+alpha(i,i)) + 2.d0*sum1(i))

               Amat(i,j) = tau(i,j) - 0.5d0*zeta0     !A matrix for solution of etajk (p 10 CMH notes)
            else
               tau(i,j) = 8.d0*dsqrt(pi)/15.d0*v0   &
                 * ni(i)*chi(i,j)*sigma(i,j)**2*mu(i,j)*theta(j)**1.5d0 &
                 /dsqrt(theta(i))*(1.d0+alpha(i,j))  * (6.d0*beta(i,j)/ &
                 (theta(j)**2*dsqrt(theta(i)+theta(j)))   &
                 + (9.d0-3.d0*alpha(i,j))/2.d0*mu(j,i)/theta(j)**2* &
                 dsqrt(theta(i)+theta(j))   &
                 - 5.d0/(theta(j)*dsqrt(theta(i)+theta(j))))

               Amat(i,j) = tau(i,j)                 !A matrix for solution of etajk (p 10 CMH notes)
            endif
!
         enddo
         bmat(i) = ni(i)*Ti(i)+sum2(i)              !b matrix for solution of etajk (p 10 CMH notes)
      enddo

      CALL LUDCMP(Amat, s, NP, indx, d, 'shear_viscosity')     ! solve system of s linear equations using
      CALL LUBKSB(Amat, s, NP, indx, bmat)  ! LU decomposition

      etakin = 0.d0
      do i=1,s
         etajk(i) = bmat(i)
         etakin = etakin + etajk(i)     !kinetic shear viscosity (p 11 CMH notes)
      enddo

      etacol = 0                        !collisional shear viscosity (p 11 CMH notes)
      do i=1,s
         do j=1,s
            etacol = etacol + 4.d0*(pi)/15.d0*ni(j)*sigma(i,j)**3 &
                        *chi(i,j)*mu(j,i)*(1d0+alpha(i,j))*etajk(i)
         enddo
      enddo
      etacol = etacol + 0.6d0*kappa

      eta = etakin + etacol           !total shear viscosity (p 11 CMH notes)

      return
      end subroutine shear_viscosity

