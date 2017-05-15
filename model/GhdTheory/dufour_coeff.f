!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name:  dufour_coeff(s,mi,alpha,T,ni,rho,v0,
!            mu,sigma,chi,beta,zeta0,theta,Ti,Dij,lambdai,gammaij,
!            omega,I_ilj,dTl_dnj,dzeta0_dnj,dchi0il_dnj,Dq)
!
!  author:  C. Hrenya, Feb 2009
!
!  Purpose: find dufour coefficient according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!  Modifications:
!     Sof: removed divisions by ni.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine dufour_coeff(s,mi,alpha,T,ni,rho,v0, &
             mu,sigma,chi,beta,zeta0,theta,Ti,Dij,lambdai,gammaij, &
             omega,I_ilj,dTl_dnj,dzeta0_dnj,dchi0il_dnj,Dq)
      Implicit NONE

      integer s, indx(s)

      double precision mi(s),alpha(s,s),T,ni(s),rho,v0,mu(s,s), &
                      sigma(s,s),chi(s,s),beta(s,s),zeta0,theta(s), &
                      Ti(s),Dij(s,s),lambdai(s),gammaij(s,s), &
                      omega(s,s),I_ilj(s,s,s),Dtl_dnj(s,s), &
                      dzeta0_dnj(s),dchi0il_dnj(s,s,s),Dq(s,s)

      integer i,j,l,p,kk
      double precision kronecker(s,s),integ1(s,s,s),integ2(s,s,s), &
                      integ(s,s,s),sum1(s,s),sum2(s,s), &
                      dq_bar(s,s),dqlj(s,s), &
                      Dqkin(s,s),CipjT(s,s,s),Dqcol(s,s), &
                      Amat(s,s),bmat(s,s), &
                      Amat0(s,s),bmat0(s)
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

!determination of integral - p 18 of CMH notes (eq. 3) - via summation of two smaller integrals
! Note that the new I_ilj(i,l,j) is now defined as old I_ilj(i,l,j)*ni(l)/ni(j)
      do i=1,s
         do l=1,s
            do j=1,s
               !integral 1 - p 18.1 of CMH notes
               integ1(i,l,j) = (kronecker(l,j)*ni(l) + 0.5d0*(ni(l)/chi(i,l) &
                 *dchi0il_dnj(i,l,j)*ni(l)+I_ilj(i,l,j)*ni(j)))  &
                 *pi*mi(i)*ni(i)*chi(i,l)*sigma(i,l)**3*mu(l,i) &
                 *(1.d0+alpha(i,l))  &
                 *((11.d0*mu(i,l)**2+(13.d0-9.d0*alpha(i,l))*mu(i,l) &
                 *mu(l,i)+(5.d0+3.d0*alpha(i,l)**2-3.d0*alpha(i,l)) &
                 *mu(l,i)**2)*Ti(i)**2/mi(i)**2   &
                 +3.d0*mu(l,i)**2*(1.d0+alpha(i,l))**2*Ti(l)**2 &
                 /mi(l)**2+(5.d0*mu(i,l)**2+(1.d0-9.d0*alpha(i,l)) &
                 *mu(i,l)*mu(l,i)+(2.d0+3.d0*alpha(i,l)+6.d0 &
                 *alpha(i,l)**2)*mu(l,i)**2)*Ti(i)*Ti(l)/mi(i)/mi(l)  &
                 -5.d0*(Ti(i)/mi(i)+Ti(l)/mi(l))*Ti(i)/mi(i))
               !integral 2 - p 18.2 of CMH notes
               integ2(i,l,j) = 0.5d0*ni(j)/Ti(l)*dTl_dnj(l,j)  &
                  *2.d0*pi*ni(i)*ni(l)*mu(i,l)*chi(i,l)*sigma(i,l)**3 &
                  *Ti(l)*(1.d0+alpha(i,l))   &
                  *(Ti(i)/mi(i)*(5.d0*(mu(i,l)**2-1.d0)  &
                  +(1.d0-9.d0*alpha(i,l))*mu(i,l)*mu(l,i)  &
                  +(2.d0+3.d0*alpha(i,l)+6.d0*alpha(i,l)**2) &
                  *mu(l,i)**2)  &
                  +6.d0*Ti(l)/mi(l)*mu(l,i)**2*(1.d0+alpha(i,l))**2)
               !summation
               integ(i,l,j) = integ1(i,l,j) + integ2(i,l,j)
            enddo
         enddo
      enddo

      sum1(1:s,1:s) = 0.d0           !calculate 1st summation used in dq_bar - p 18 CMH notes
      sum2(1:s,1:s) = 0.d0           !calculate 2nd summation used in dq_bar - p 18 CMH notes
      do i=1,s
         do j=1,s
            do l=1,s  ! modification to ni(l) > 0 same as in thermal_mobility.f
               if(ni(l) > 0d0) sum1(i,j) = sum1(i,j) + mi(l)*(omega(i,l)-zeta0 &
                          *kronecker(i,l))/ni(l)/Ti(l)*Dij(l,j)
               sum2(i,j) = sum2(i,j) + integ(i,l,j)
            enddo
            dq_bar(i,j) = -2.5d0*ni(j)*Ti(i)**3/mi(i)/T**2   &
                        *(mi(j)/rho/Ti(i)*sum1(i,j)*ni(i)  &
                        -0.4d0*mi(i)*T/Ti(i)**3*dzeta0_dnj(j) &           ! ni(i)/ni(i) cancels
                        *lambdai(i) - dTl_dnj(i,j)/3.d0/Ti(i)**2*ni(i))  &
                        +sum2(i,j)/3.d0/T**2
         enddo
      enddo

!find dq,lj (p 18 CMH notes), which will later be used to find the kinetic contribution to Dq
      do i=1,s
         do l=1,s
            Amat(i,l) = gammaij(i,l)-1.5d0*zeta0*kronecker(i,l)      !A matrix for solution of dq,lj (p 18 CMH notes)
            bmat(i,l) = dq_bar(i,l)                                  !b matrix for solution of dq,lj (p 18 CMH notes)
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

         CALL LUDCMP(Amat0, s, NP, indx, d, 'dufour_coeff') ! solve system of s linear equations using
         CALL LUBKSB(Amat0, s, NP, indx, bmat0) ! LU decomposition

         do i=1,s
            dqlj(i,kk) = bmat0(i)
         enddo
      enddo

!calculate kinetic contribution to the Dufour coefficient - p18.2 CMH notes
      do l=1,s
         do j=1,s
            Dqkin(l,j) = dqlj(l,j)   &
              +2.5d0/T**2*mi(j)*ni(j)*Ti(l)/rho*Dij(l,j)
         enddo
      enddo

!evaluate coefficient CipjT needed to evaluate collisional component - p 18.3 CMH notes
      do i=1,s
         do p=1,s
            do j=1,s
               CipjT(i,p,j) = 8.d0*dsqrt(pi)/3.d0*ni(i)*ni(p)*v0**3   &
                 /dsqrt(theta(i)+theta(p))/(theta(i)*theta(p))**1.5d0  &
                 *(kronecker(j,p)*beta(i,p)*(theta(i)+theta(p))  &
                 -0.5d0*theta(i)*theta(p)*(1.d0+(mu(p,i)* &
                 (theta(i)+theta(p))  &
                 -2.d0*beta(i,p))/theta(p))*ni(j)/Ti(p)*dTl_dnj(p,j))  &
                 +2.d0*dsqrt(pi)/3.d0*ni(i)*ni(p)*v0**3* &
                 (1.d0-alpha(i,p))  &
                 *(mu(p,i)-mu(i,p))*((theta(i)+theta(p))/theta(i) &
                 /theta(p))**1.5d0   &
                 *(kronecker(j,p)+1.5d0*theta(i)/(theta(i)+theta(p))  &
                 *ni(j)/Ti(p)*dTl_dnj(p,j))
            enddo
         enddo
      enddo

!find collisional contribution - p 18.3 CMH notes
      Dqcol(1:s,1:s) = 0.d0
      do i=1,s
         do j=1,s
            do p=1,s
               Dqcol(i,j) = Dqcol(i,j)  &
                 + (1.d0+alpha(i,p))/8.d0*mi(p)*mu(i,p)*sigma(i,p)**3 &
                 *chi(i,p)*(4.d0*pi/5.d0*(1.d0-alpha(i,p))* &
                 (mu(i,p)-mu(p,i))  &
                 *ni(i)*(2.d0/mi(p)*Dqkin(p,j)  &
                 +5.d0*Ti(i)/T**2*mi(j)*ni(j)/rho/mi(i)*Dij(p,j))  &
                 +16.d0*pi/5.d0*ni(i)*(2.d0*mu(p,i)/mi(p)*Dqkin(p,j)  &
                 -5.d0*(2.d0*mu(i,p)-mu(p,i))*Ti(i)/T**2*ni(j)*mi(j) &
                 /rho/mi(i)*Dij(p,j)) &
                 -sigma(i,j)/T**2*CipjT(i,p,j))
            enddo
            !final total Dufour coefficient - p 18 CMH notes
            Dq(i,j) = Dqkin(i,j) + Dqcol(i,j)
         enddo
      enddo

      return
      end subroutine dufour_coeff

