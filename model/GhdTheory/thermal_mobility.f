!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: thermal_mobility(s,mi,alpha,ni,mu,sigma,chi,zeta0,
!                                theta,Ti,DF,gammaij,omega,Lij)
!
!  author:  C. Hrenya, Feb 2009
!
!  Purpose: find thermal mobility coefficient according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!  Modifications:
!     Sof: removed /ni in sum1, which yield zero based on Bill's mods.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine thermal_mobility(s,mi,alpha,ni,mu,sigma,chi,zeta0, &
                                 theta,Ti,DF,gammaij,omega,Lij)
      Implicit NONE

      integer s, indx(s)

      double precision mi(s),alpha(s,s),ni(s),mu(s,s),sigma(s,s), &
                      chi(s,s),zeta0,theta(s),Ti(s),DF(s,s), &
                      gammaij(s,s),omega(s,s),Lij(s,s)

      integer i,j,k,p,kk
      double precision kronecker(s,s),sum1(s,s),l_bar(s,s),lkj(s,s), &
                      Lkin(s,s),Lcol(s,s),Amat(s,s),bmat(s,s), &
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

!calculate summation used in l_bar (b vector) - p. 19 of CMH notes
      do i=1,s
         do j=1,s
            sum1(i,j) = 0.d0
         enddo
      enddo
      do i=1,s
         do j=1,s
            do k=1,s
! in case ni(k) == 0d0, the contribution to sum1 can be neglected based on earlier code by Bill Holloway.
               if(ni(k) > 0d0) sum1(i,j) = sum1(i,j) + (omega(i,k)-zeta0* &
                                        kronecker(i,k))/(ni(k)*Ti(k))*DF(k,j)
            enddo
         enddo
      enddo

!calculate l_bar (p 19 of CMH notes)
      do i=1,s
         do j=1,s
            l_bar(i,j) = -2.5d0*ni(i)*Ti(i)**2/mi(i)*sum1(i,j)
            Amat(i,j) = gammaij(i,j) - 0.5d0*zeta0*kronecker(i,j)     !A matrix for solution of lkj (p 19 CMH notes)
            bmat(i,j) = l_bar(i,j)                                    !B matrix for solution of lkj (p 19 CMH notes)
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

         CALL LUDCMP(Amat0, s, NP, indx, d, 'thermal_mobility') ! solve system of s linear equations using
         CALL LUBKSB(Amat0, s, NP, indx, bmat0) ! LU decomposition

         do i=1,s
            lkj(i,kk) = bmat0(i)
         enddo
      enddo

!kinetic contribution to thermal mobility (p 19 CMH notes)
      do i = 1,s
         do j=1,s
            Lkin(i,j) = lkj(i,j) + 2.5d0*Ti(i)/mi(i)*DF(i,j)
         enddo
      enddo

!collisional contribution to thermal mobility (p 20 CMH notes)
      do i=1,s
         do j=1,s
            Lcol(i,j) = 0.d0
         enddo
      enddo
      do i=1,s
         do j=1,s
            do p=1,s
               Lcol(i,j) = Lcol(i,j) + (1.d0+alpha(i,p))/8.d0*mi(p)* &
                 mu(i,p)*sigma(i,p)**3*chi(i,p)*(4.d0*pi/5.d0* &
                 (1.d0-alpha(i,p))*(mu(i,p)-mu(p,i))*ni(i)  &
                 *(2.d0/mi(p)*Lkin(p,j)+5.d0*Ti(i)/mi(i)/mi(p)* &
                 DF(p,j))+48.d0*pi/15.d0*ni(i)*(2.d0*mu(p,i)/mi(p)* &
                 Lkin(p,j)-5.d0*(2.d0*mu(i,p)-mu(p,i))*Ti(i)/mi(i)/ &
                 mi(p)*DF(p,j)))
            enddo
            Lij(i,j) = Lkin(i,j) + Lcol(i,j)          !thermal mobility (p 19 CMH notes)
         enddo
      enddo

      return
      end subroutine thermal_mobility
