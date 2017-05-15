!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: thermal_conductivity(s,mi,T,sigmai,alpha,ni,rho,v0,mu,
!                     sigma,chi,beta,zeta0,theta,Ti,DT,lambda,
!                     omega,gamma,lambdai)
!
!  author:  C. Hrenya, Jan 2009
!
!  Purpose: find thermal conductivity according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine thermal_conductivity(s,mi,T,sigmai,alpha,ni,rho,v0,mu, &
                                sigma,chi,beta,zeta0,theta,Ti,DT,lambda, &
                               omega,gamma,lambdai)
      Implicit NONE

      integer s, indx(s)

      double precision mi(s),sigmai(s),T,alpha(s,s),ni(s),rho,v0, &
                       mu(s,s),sigma(s,s),chi(s,s),beta(s,s),zeta0, &
                       theta(s),Ti(s),DT(s),lambda

      integer i,j
      double precision E(s,s),A(s,s),C(s,s),F(s,s),sum1(s),sum2(s), &
                       sum3(s),sum4(s),omega(s,s),lambdai_bar(s), &
                       gamma(s,s),Amat(s,s),bmat(s),lambdai(s), &
                       lambdakini(s),lambdakin,lambdacol,CT(s,s)
      double precision d

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve
      double precision pi
      parameter (pi=3.14159265458979323846d0)

!calculate quantities on p.14-15 of CMH notes
      do i=1,s
         do j=1,s
            ! CMH edit (9/18/09):  re-introduced mu(j,i)**2 instead of mu(j,i)
            ! in first line to correct error
            E(i,j) = 2.d0*mu(j,i)**2/theta(i)**2*(theta(i) &
                    +theta(j))**2  &
                    *(2.d0*alpha(i,j)**2-3.d0*alpha(i,j)+4.d0) &
                    *(5.d0*theta(i)+8.d0*theta(j)) &
                    -mu(j,i)*(theta(i)+theta(j))  &
                    *(beta(i,j)/theta(i)**2*(5.d0*theta(i)+8.d0*theta(j)) &
                    *(14.d0*alpha(i,j)-22.d0)-theta(j)/theta(i) &
                    *(20.d0+3.d0*(15.d0-7.d0*alpha(i,j))+9.d0 &
                    *(1.d0-alpha(i,j))-28.d0*alpha(i,j))  &
                    -25.d0*(1d0-alpha(i,j))) +18.d0*(beta(i,j)/theta(i))**2 &
                    *(5.d0*theta(i)+8.d0*theta(j))+2.d0*beta(i,j)  &
                    /theta(i)*(25.d0*theta(i)+66.d0*theta(j))+5.d0*theta(j)  &
                    /theta(i)*(11.d0*theta(i)+6.d0*theta(j))-5.d0*(theta(i) &
                    +theta(j))*theta(j)/theta(i)**2*(5.d0*theta(i) &
                    +6.d0*theta(j))

            ! CMH edit (9/18/09): corrected error in final line,
            ! by making theta(j)/theta(i) instead of its reciprocal
            A(i,j) = 5.d0*(2.d0*beta(i,j)+theta(j))+mu(j,i) &
                    *(theta(i)+theta(j))*(5.d0*(1.d0-alpha(i,j)) &
                    -(14.d0*alpha(i,j)-22.d0)*beta(i,j)/theta(i)) &
                    +18.d0*beta(i,j)**2/theta(i)+2.d0*mu(j,i)**2 &
                    *(2.d0*alpha(i,j)**2-3.d0*alpha(i,j)+4.d0) &
                    *(theta(i)+theta(j))**2/theta(i) &
                    -5.d0*theta(j)/theta(i)*(theta(i)+theta(j))

            C(i,j) = 5.d0*(2.d0*beta(i,j)-theta(i))+mu(j,i) &
                     *(theta(i)+theta(j))*(5.d0*(1.d0-alpha(i,j)) &
                     +(14.d0*alpha(i,j)-22.d0)*beta(i,j)/theta(j)) &
                     -18.d0*beta(i,j)**2/theta(j)-2.d0*mu(j,i)**2 &
                     *(2.d0*alpha(i,j)**2-3.d0*alpha(i,j)+4.d0) &
                     *(theta(i)+theta(j))**2/theta(j) &
                     +5.d0*(theta(i)+theta(j))

            F(i,j) = 2.d0*(mu(j,i)*(theta(i)+theta(j)) &
                    /theta(j))**2*(2.d0*alpha(i,j)**2-3.d0 &
                    *alpha(i,j)+4.d0)*(8.d0*theta(i)+5.d0*theta(j)) &
                    -mu(j,i)*(theta(i)+theta(j))*(beta(i,j) &
                    /theta(j)**2*(8.d0*theta(i)+5.d0*theta(j)) &
                    *(14.d0*alpha(i,j)-22.d0)+theta(i)/theta(j) &
                    *(20.d0+3.d0*(15.d0-7.d0*alpha(i,j))+9.d0 &
                    *(1.d0-alpha(i,j))-28.d0*alpha(i,j))+25.d0 &
                    *(1.d0-alpha(i,j)))+18.d0 &
                    *(beta(i,j)/theta(j))**2*(6.d0*theta(i) &
                    +5.d0*theta(j)) -2.d0*beta(i,j)/theta(j)*(66.d0 &
                    *theta(i)+25.d0*theta(j))+5.d0*theta(i)/theta(j) &
                    *(6.d0*theta(i)+11.d0*theta(j))-5.d0*(theta(i) &
                    +theta(j))/theta(j)*(6.d0*theta(i)+5.d0*theta(j))
         enddo
      enddo

! calculate summations used in gamma(i,i) and omega(i,i) - p. 13 and 16 of CMH notes, respectively
      do i=1,s
         sum1(i) = 0.d0     !calculate summation used in gamma(i,i)
         sum2(i) = 0.d0     !calculate summation used in omega(i,i)
      enddo

      do i=1,s
         do j=1,s
            if (j .ne. i) then
               sum1(i) = sum1(i)+ni(j)*chi(i,j)*sigma(i,j)**2 &
                  *v0*mu(j,i)*(1.d0+alpha(i,j))*(theta(i) &
                  /(theta(j)*(theta(i)+theta(j))))**1.5d0 &
                  *(E(i,j)-5.d0*(theta(i)+theta(j))/theta(i) &
                  *A(i,j))
               sum2(i) = sum2(i) + ni(j)*chi(i,j)*sigma(i,j)**2 &
                  *v0*mu(j,i)*(1.d0+alpha(i,j))*dsqrt(theta(i) &
                  /((theta(i)+theta(j))*theta(j)**3))*A(i,j)
            endif
         enddo
      enddo

! calculate summations used in lambda_bar - p. 15 of CMH notes
      do i=1,s
         sum3(i) = 0.d0     !calculate 1st summation used in b vector (lambda_bar)
         sum4(i) = 0.d0     !calculate 2nd summation used in b vector (lambda_bar)
      enddo

      do i=1,s
         do j=1,s
            !first find omega's (p 16 CMH notes)
            if (j .eq. i) then
               omega(i,i) = 4.d0/3.d0*dsqrt(pi/2.d0)*sigmai(i)**2 &
                 *ni(i)*chi(i,i)*v0*(1.d0-alpha(i,i)**2) &
                 /dsqrt(theta(i)) + 4.d0*dsqrt(pi)/15d0*sum2(i)

               if(ni(j)>0d0) sum3(i) = sum3(i)+  &
                             (omega(i,j)-zeta0)/(ni(j)*Ti(j))*DT(j)
            else
               omega(i,j) = 4.d0*dsqrt(pi)/15.d0*ni(j)*chi(i,j) &
                 *sigma(i,j)**2*v0*mu(j,i)*(1.d0+alpha(i,j))   &
                 *dsqrt(theta(i)/((theta(i)+theta(j)) &
                 *theta(j)**3))*C(i,j)

               if(ni(j)>0d0) sum3(i) = sum3(i)+  &
                             omega(i,j)/(ni(j)*Ti(j))*DT(j)
            endif
            sum4(i) = sum4(i) + 2.d0*pi*ni(i)*ni(j)*mu(i,j) &
                 *chi(i,j)*sigma(i,j)**3*Ti(j)*(1.d0+alpha(i,j))   &
                 *(Ti(i)/mi(i)*(5.d0*(mu(i,j)**2-1.d0)  &
                 +(1.d0-9.d0*alpha(i,j))*mu(i,j)*mu(j,i)  &
                 +(2.d0+3.d0*alpha(i,j)+6.d0*alpha(i,j)**2) &
                 *mu(j,i)**2) +6.d0*Ti(j)/mi(j)*mu(j,i)**2 &
                 *(1.d0+alpha(i,j))**2)
         enddo
      enddo

!calculate gamma (p13 of CMH notes) an lambda_bar (p 15 of CMH notes)
      do i=1,s
         lambdai_bar(i) = -2.5d0*rho*ni(i)*Ti(i)**3 &
            /(Ti(i)*T*mi(i))*sum3(i)+2.5d0*ni(i)*Ti(i)**2 &
            /(mi(i)*T)+sum4(i)/(6.d0*T)
         do j=1,s
            if (i .eq. j) then
               gamma(i,i) = 16.d0/15.d0*dsqrt(pi)*sigmai(i)**2 &
                  *ni(i)*chi(i,i)*v0/dsqrt(2.d0*theta(i)) &
                  *(1.d0+alpha(i,i))*(1.d0+33.d0/16.d0 &
                  *(1.d0-alpha(i,i)))  &
                   +2.d0/15.d0*dsqrt(pi)*sum1(i)

               Amat(i,j) = gamma(i,j) - 2.d0*zeta0 !A matrix for solution of lambdai (p 13 CMH notes)
            else
               gamma(i,j) = -2.d0/15.d0*dsqrt(pi)*ni(i)*chi(i,j) &
                  *sigma(i,j)**2*v0*mu(i,j)*(1.d0+alpha(i,j))   &
                  *(theta(j)/(theta(i)*(theta(i)+theta(j))))**1.5d0   &
                  *(F(i,j)+5.d0*((theta(i)+theta(j))/theta(j)) &
                  *C(i,j))

               Amat(i,j) = gamma(i,j)             !A matrix for solution of lambdai (p 13 CMH notes)
            endif
         enddo
         bmat(i) = lambdai_bar(i)                 !B vector for solution of lambdai (p 13 CMH notes)
      enddo

      CALL LUDCMP(Amat, s, NP, indx, d, 'thermal_conductivity')     ! solve system of s linear equations using
      CALL LUBKSB(Amat, s, NP, indx, bmat)  ! LU decomposition

      do i=1,s
         lambdai(i) = bmat(i)
      enddo

      lambdakin = 0             !kinetic contribution to conductivity (p 16 CMH notes)
      do i = 1,s
         lambdakini(i) = lambdai(i) + 2.5d0/T*rho*Ti(i)*DT(i)/mi(i)
         lambdakin = lambdakin + lambdakini(i)
      enddo

      lambdacol = 0d0            !collisional contribution to conductivity (p 17 CMH notes)
      do i=1,s
         do j=1,s
            CT(i,j) = -4.d0/3.d0*dsqrt(pi)*ni(i)*ni(j)*v0**3 &
               /dsqrt((theta(i)+theta(j))*(theta(i)*theta(j))**3)  &
               *(2.d0*beta(i,j)**2+theta(i)*theta(j)  &
               + (theta(i)+theta(j)) *((theta(i)+theta(j))*mu(i,j) &
               *mu(j,i)+beta(i,j)*(1.d0+mu(j,i))))-dsqrt(pi)*ni(i) &
               *ni(j)*v0**3*(1.d0-alpha(i,j))*(mu(j,i)-mu(i,j))  &
               *((theta(i)+theta(j))/theta(i)/theta(j))**1.5d0 &
               *(mu(j,i)+beta(i,j)/(theta(i)+theta(j)))
             lambdacol = lambdacol + (1.d0+alpha(i,j))/8.d0*mi(j) &
               *mu(i,j)*sigma(i,j)**3*chi(i,j)  *(4.d0*pi/5.d0 &
               *(1.d0-alpha(i,j))*(mu(i,j)-mu(j,i))*ni(i)  &
               *(2.d0/mi(j)*lambdakini(j)+5.d0*Ti(i)/(mi(i)*mi(j) &
               *T)*rho*DT(j)) +48.d0*pi/15.d0*ni(i)*(2.d0*mu(i,j) &
               /mi(j)*lambdakini(j)-5.d0*Ti(i)/(mi(i)*mi(j)*T) &
               *(2.d0*mu(i,j)-mu(j,i))*rho*DT(j))  &
               -sigma(i,j)*CT(i,j)/T)
         enddo
      enddo

      lambda = lambdakin + lambdacol                  !conductivity (p 13 CMH notes)

      return
      end subroutine thermal_conductivity
