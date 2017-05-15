!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: cooling_rate_tc(s,mi,sigmai,alpha,ni,n,v0,mu,sigma,chi,T,
!                                   zeta0,theta,Ti,zetau)
!
!  author:  C. Hrenya, Feb 2009
!
!  Purpose: find transport coefficient for cooling rate according to
!           GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine cooling_rate_tc(s,mi,sigmai,alpha,ni,n,v0,mu,sigma, &
                                chi,T,zeta0,theta,Ti,zetau)
      Implicit NONE

      integer s, indx(s)

      double precision mi(s),sigmai(s),alpha(s,s),ni(s),n,v0,mu(s,s), &
                      sigma(s,s),chi(s,s),T,zeta0,theta(s), &
                      Ti(s),zetau

      integer i,j
      double precision kronecker(s,s),sum1(s),psi(s,s),eid_bar(s), &
                      Amat(s,s),bmat(s),ejd(s),zeta10,zeta11
      double precision d

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve

      double precision pi
      parameter (pi=3.14159265458979323846d0)

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
         sum1(i) = 0.d0     !calculate summation used in psi(i,i) - p 5 CMH notes
      enddo
      do i=1,s
         do j=1,s
            if (i.ne.j) then
               sum1(i) = sum1(i)+dsqrt(pi)/2.d0*ni(j)*chi(i,j)* &
                sigma(i,j)**2*v0*mu(j,i)*(1.d0+alpha(i,j))/ &
                dsqrt(theta(i)*theta(j))/(theta(i)+theta(j))**2.5d0  &
                *(-2.d0*(90.d0*theta(j)**3+231.d0*theta(i)* &
                theta(j)**2+184.d0*theta(i)**2*theta(j)+40.d0* &
                theta(i)**3)+3.d0*mu(j,i)*(1.d0+alpha(i,j))* &
                (theta(i)+theta(j))*(70.d0*theta(j)**2+117.d0* &
                theta(i)*theta(j)+44.d0*theta(i)**2)-24.d0* &
                mu(j,i)**2*(1.d0+alpha(i,j))**2* &
                (theta(i)+theta(j))**2*(5.d0*theta(j)+4.d0*theta(i))  &
                +30.d0*mu(j,i)**3*(1.d0+alpha(i,j))**3*(theta(i)+ &
                theta(j))**3+5.d0*theta(j)*(theta(i)+theta(j))* &
                (2.d0*(4.d0*theta(i)+3.d0*theta(j))-3.d0*mu(j,i)* &
                (1.d0+alpha(i,j))*(theta(i)+theta(j))))
            endif
         enddo
      enddo

! find psi (p 5 CMH notes) and eid_bar (p 4) values
      do i=1,s
         eid_bar(i) = 0.d0
      enddo
      do i=1,s
         do j=1,s
            if (i.eq.j) then
               psi(i,i) = -2.d0/15.d0*(sum1(i)+dsqrt(2.d0*pi) &
                 /32.d0*ni(i)*chi(i,i)*sigmai(i)**2*(1.d0+ &
                 alpha(i,i))/dsqrt(theta(i))*v0*(30.d0*alpha(i,i)**3 &
                 -126.d0*alpha(i,i)**2+177.d0*alpha(i,i)+48.d0* &
                 (3.d0*alpha(i,i)-7.d0)-137.d0)+3.d0/32.d0* &
                 dsqrt(2.d0*pi)*ni(i)*chi(i,i)*sigmai(i)**2* &
                 (1.d0+alpha(i,i))/dsqrt(theta(i))*v0  &
                 *(10.d0*alpha(i,i)**3+22.d0*alpha(i,i)**2+ &
                 11.d0*alpha(i,i)-3.d0))
             else
                psi(i,j) = -2.d0/15.d0*(dsqrt(pi)/2.d0* &
                  ni(j)*chi(i,j)*sigma(i,j)**2*v0*mu(j,i)*(1.d0+ &
                  alpha(i,j))*(theta(i)/theta(j))**1.5d0/(theta(i)+ &
                  theta(j))**2.5d0*((2.d0*theta(j)+5.d0*theta(i))* &
                  (2.d0*theta(j)+3.d0*mu(j,i)*(1.d0+alpha(i,j))  &
                  *(theta(i)+theta(j)))-24.d0*mu(j,i)**2*(1.d0+ &
                  alpha(i,j))**2*(theta(i)+theta(j))**2+30.d0* &
                  mu(j,i)**3*(1.d0+alpha(i,j))**3*(theta(i)+ &
                  theta(j))**3/theta(j)-5.d0*(theta(i)+theta(j))* &
                  (2.d0*theta(j)+3.d0*mu(j,i)*(1.d0+alpha(i,j))* &
                  (theta(i)+theta(j)))))
             endif
             eid_bar(i) = eid_bar(i) - pi/6.d0*ni(j)*chi(i,j)* &
                  sigma(i,j)**3*mu(j,i)*(1.d0+alpha(i,j))   &
                  *(40.d0*(mu(i,j)-1.d0)+4.d0*(19.d0+9.d0* &
                  alpha(i,j))*mu(j,i)-48.d0*mu(j,i)**2*(theta(i)+ &
                  theta(j))/theta(j)*(1.d0+alpha(i,j))**2   &
                  +15.d0*mu(j,i)**3*(theta(i)+theta(j))**2/ &
                  theta(j)**2*(1.d0+alpha(i,j))**3)
         enddo
      enddo

      do i=1,s
         do j=1,s
            Amat(i,j) = psi(i,j) - 1.5d0*zeta0*kronecker(i,j)     !A matrix for solution of ejd (p 4 CMH notes)
         enddo
         bmat(i) = 2.d0*eid_bar(i)/15.d0                  !b matrix for solution of ejd (p 4 CMH notes)
      enddo

      CALL LUDCMP(Amat, s, NP, indx, d, 'cooling_rate_tc')     ! solve system of s linear equations using
      CALL LUBKSB(Amat, s, NP, indx, bmat)  ! LU decomposition

      do i=1,s
         ejd(i) = bmat(i)
      enddo

! determine summations for cooling rate transport coefficients (p 4 CMH notes)
      zeta10 = 0.d0
      zeta11 = 0.d0
      do i=1,s
         do j=1,s
            zeta10 = zeta10 - 2.d0*pi/(3.d0*n*T)*ni(i)*ni(j)*mu(j,i) &
              *sigma(i,j)**3*chi(i,j)*(1.d0-alpha(i,j)**2)*Ti(i)
            zeta11 = zeta11 + dsqrt(pi)*v0**3/(2.d0*n*T)*ejd(i)* &
              ni(i)*ni(j)*sigma(i,j)**2*chi(i,j)*mi(j)  &
              *mu(i,j)*(1.d0-alpha(i,j)**2)*dsqrt(theta(j)/ &
              (theta(i)**3*(theta(i)+theta(j))))
         enddo
      enddo

      zetau = zeta10 + zeta11         !transport coefficient for cooling rate (p 4 CMH notes)

      return
      end subroutine cooling_rate_tc

