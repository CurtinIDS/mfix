!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: cooling_rate(s,mi,ni,n,m,T,chi,sigmai,alpha,rhoi,xvec)
!
!  author:  C. Hrenya, Jan 2009
!
!  Purpose: find theta, Ti, and zeta0 according to GHD polydisperse KT
!
!  Literature/References:
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine cooling_rate(s,mi,ni,n,m,Ti,T,chi,sigmai,alpha,rhoi,xvec)
      Implicit NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer, INTENT(IN) :: s
      double precision :: mi(s)
      double precision :: ni(s)
      double precision :: n, m
      double precision :: Ti(s)
      double precision :: T
      double precision :: chi(s,s)
      double precision :: sigmai(s)
      double precision :: alpha(s,s)
      double precision :: rhoi(s)
      double precision :: xvec(s)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: L
      integer :: ntrial
      double precision :: tolx, tolf
!-----------------------------------------------

      ntrial = 100
      tolx = 1d-14
      tolf = 1d-14

! Initial guess for theta
      DO L = 1, s
          xvec(L) = mi(L) / m
      ENDDO

      CALL MNEWT(ntrial, xvec, s, tolx, tolf, &
                 mi,ni,n,m,chi,sigmai,alpha,rhoi)

      return
      end subroutine cooling_rate


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: MNEWT(ntrial, x, mmax, tolx, tolf, mi,ni,n,m,chi,sigmai,alpha,rhoi)
!  Purpose: use Newton-Raphson method to solve set of non-linear
!           equations to be used for GHD polydisperse KT
!
!  Literature/Document References:
!     Numerical Recipies in Fortran 77, page 374
!
!  Note:  this subroutine has been modified by CMH from Sofiane original
!         to accommodate the inputs/outputs of the stand-along program ghd
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine MNEWT(ntrial, x, s, tolx, tolf, &
                       mi,ni,n,m,chi,sigmai,alpha,rhoi)
      Implicit NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, intent(in) :: s
      INTEGER, intent(in) :: ntrial
      DOUBLE PRECISION, intent(in) :: tolx, tolf
      DOUBLE PRECISION :: mi(s)
      DOUBLE PRECISION :: ni(s)
      DOUBLE PRECISION :: n, m
      DOUBLE PRECISION :: chi(s,s)
      DOUBLE PRECISION :: sigmai(s)
      DOUBLE PRECISION :: alpha(s,s)
      DOUBLE PRECISION :: rhoi(s)
      DOUBLE PRECISION :: X(s)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER :: NP
      PARAMETER (NP=15)  ! solves up to NP variable/equations;
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, K, INDX(s)
      DOUBLE PRECISION :: D, ERRF, ERRX
      DOUBLE PRECISION :: FJAC(s,s), FVEC(s)
      DOUBLE PRECISION :: P(s)
!-----------------------------------------------

      DO K = 1, NTRIAL
        CALL FUNC_JACOBI_EVAL(X, s, FVEC, FJAC, &
                              mi,ni,n,m,chi,sigmai,alpha,rhoi)
        ERRF = 0d0
        DO I = 1, s
             ERRF = ERRF + DABS(FVEC(I))
        ENDDO

        IF(ERRF <= TOLF) RETURN
        DO I = 1, s
            P(I) = -FVEC(I) ! RHS of linear equations.
        ENDDO

        CALL LUDCMP(fjac, s, NP, indx, d, 'MNewt')  ! solve system of s linear equations using
        CALL LUBKSB(fjac, s, NP, indx, p)  ! LU decomposition

        ERRX = 0d0
        DO I = 1, s
            ERRX = ERRX + DABS(P(I))
            X(I) = X(I) + P(I)
        ENDDO
        IF(ERRX <= TOLX) RETURN
      ENDDO  ! for K trials

      RETURN
      END SUBROUTINE MNEWT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: FUNC_JACOBI_EVAL(X, s, FVEC, FJAC, mi,ni,n,m,chi,sigmai,alpha,rhoi)
!  Purpose: computes values of functions and Jacobians of homogenous
!           cooling rates to get values of theta_i (or Ti)
!  Literature/Document References:
!     C. Hrenya handnotes and S. Benyahia derivation, Dec-2008
!
!  Note:  this subroutine has been modified by CMH from Sofiane original
!         to accommodate the inputs/outputs of the stand-along program ghd
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE FUNC_JACOBI_EVAL(X, s, FVEC, FJAC, &
                                  mi,ni,n,m,chi,sigmai,alpha,rhoi)
      Implicit NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, intent(in) :: s
      DOUBLE PRECISION :: X(s)
! vector function and matrix jacobian
      DOUBLE PRECISION :: FVEC(s), FJAC(s,s)

      DOUBLE PRECISION :: mi(s)
      DOUBLE PRECISION :: ni(s)
      DOUBLE PRECISION :: n, m
      DOUBLE PRECISION :: chi(s,s), sigmai(s), alpha(s,s), &
                          rhoi(s)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      double precision :: one
      parameter (one=1.d0)
      double precision :: zero
      parameter (zero=0.d0)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K
! various quantities
      DOUBLE PRECISION :: eta(s,s), lambda(s,s)
      DOUBLE PRECISION :: gamma(s,s), sqrtXiXj(s,s)
      DOUBLE PRECISION :: oneMGama(s,s)
!-----------------------------------------------

      DO i = 1, s
        DO j = 1, s
             ETA(i,j) = (one + alpha(i,j))/2.d0
             lambda(i,j) = chi(i,j)* ni(j)*(mi(j)/(mi(j)+mi(i)))*  &
                         ((sigmai(i)+sigmai(j))/2d0)**2 * (ETA(I,J)*2d0)
! eta = (1+e)/2
             gamma(i,j) = (mi(j)/(mi(j)+mi(i)))* (ETA(I,J))
! The "2" in mu_ji/2 and 2 * eta will cancel
             sqrtXiXj(i,j) = dsqrt((x(i)+x(j))/(x(i)*x(j)))
             oneMGama(i,j) = (one-gamma(i,j)*((x(i)+x(j))/x(j)))
        ENDDO
      ENDDO


! Start computing values of the function FVEC(i)
      DO i = 1, s
        IF(i==1) THEN
          FVEC(1) = ZERO
          DO k = 1, s
            FVEC(1) = FVEC(1) + (rhoi(k)/x(k))
          ENDDO
          FVEC(1) = FVEC(1) - n*m
        ELSE
          FVEC(i) = ZERO
          DO k = 1, s
            FVEC(i) = FVEC(i)  &
                  + lambda(i,k) * sqrtXiXj(i,k) * oneMGama(i,k)  &
                  - lambda(1,k) * sqrtXiXj(1,k) * oneMGama(1,k)
          ENDDO
        ENDIF
      ENDDO
! End of computing values of the function FVEC(i)

! Evaluation of functions FVEC(i) is done above, and their Jacobian is computed below.


! Start computing jacobian (J_ij) of the functions FVEC(i).
      DO i = 1, s
        DO j = 1, s
          IF(i==1) THEN
            FJAC(1,j) = -rhoi(j)/x(j)**2

          ELSEIF(i==j) THEN
            FJAC(j,j) = ZERO

            DO k = 1, s
              IF(k==i) THEN
                FJAC(i,j) = FJAC(i,j)  &
                             - lambda(1,k)*(one/(x(1)*x(j))-(x(j)+x(1))/ &
                          (x(1)*x(j)**2))*(one-gamma(1,k)*(x(j)+x(1))/ &
                             x(j))/dsqrt((x(j)+x(1))/(x(1)*x(j)))/2.d0 &
                              -lambda(1,k)*dsqrt((x(j)+x(1))/  &
                           (x(1)*x(j)))*(gamma(1,k)*(x(j)+x(1))/x(j)**2 &
                 -gamma(1,k)/ x(j))-dsqrt(2.d0)*(one-2.d0*gamma(i,k))* &
                            lambda(i,k)*x(j)**((-3.d0)/2.d0)/2.d0

              ELSE
                FJAC(i,j) = FJAC(i,j)  &
                        + lambda(i,k)*(one-gamma(i,k)*(x(k)+x(j))/x(k))* &
                                (one/(x(j)*x(k))-(x(k)+x(j))/ &
                    (x(j)**2*x(k)))/ dsqrt((x(k)+x(j))/(x(j)*x(k)))/2.d0 &
                            - gamma(i,k)*lambda(i,k)* &
                                dsqrt((x(k)+x(j))/(x(j)*x(k)))/x(k)
              ENDIF
            ENDDO

          ELSEIF(j==1) THEN ! for i /= 1 and i /= j and j == 1
            k = 1  ! k =1 is a special case for computing Jacobioan
            FJAC(i,j)  = -lambda(i,k)*oneMGama(i,1)/ &
                             (x(1)**2*sqrtXiXj(i,1)*2.0d0)  &
                    + lambda(i,k)*sqrtXiXj(i,1)*gamma(i,k)*x(i)/x(1)**2  &
                         + dsqrt(2d0)*(one-2d0*gamma(1,k))* &
                            lambda(1,k)*x(1)**((-3.0d0)/2.0d0)/2.0d0
            DO k = 2, s
                FJAC(i,j) = FJAC(i,j)  &
                            + gamma(1,k)*lambda(1,k)*sqrtXiXj(k,1)/x(k)  &
                            + lambda(1,k)*oneMGama(1,k)/ &
                              (x(1)**2*sqrtXiXj(k,1)*2.0d0)
            ENDDO

          ELSE ! for i /= 1 and i /= j and j /= 1
! no do loop for k here, just one term k = j, the rest will be zero.
            k = j
            FJAC(i,j) = - lambda(i,k)*oneMGama(i,k)/ &
                          (x(j)**2*sqrtXiXj(j,i)*2d0)  &
                + lambda(i,k)*sqrtXiXj(j,i)*gamma(i,k)*x(i)/x(j)**2  &
                + lambda(1,k)*oneMGama(1,k)/(x(j)**2*sqrtXiXj(j,1)*2d0)  &
                - lambda(1,k)*sqrtXiXj(j,1)*gamma(1,k)*x(1)/x(j)**2

          ENDIF ! for checks on i
        ENDDO ! for j
      ENDDO ! for i
! End of computing jacobian (J_ij).

      RETURN
      END SUBROUTINE FUNC_JACOBI_EVAL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: ludcmp(a,n,np,indx,d, calledFrom)
!  Purpose: Replaces matrix a (n,n) by the LU decomposition of a rowwise
!           permutation of itself. Used in combination with lubksb.
!
!  Literature/Document References:
!     Numerical Recipies in Fortran 77, page 38-39
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      subroutine ludcmp(a,n,np,indx,d,calledFrom)
      USE compar
      USE exit, only: mfix_exit

      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer, intent(in) :: n
      double precision, intent(inout) :: a(n,n)
      integer :: np
      integer, intent(out) :: indx(n)
      double precision, intent(out) :: d
      CHARACTER(len=*), intent(in) :: calledFrom
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      integer :: nmax
      double precision :: TINY
      parameter (NMAX=500, TINY=1.0D-20)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: i, j, k, imax
      double precision :: vv(NMAX)
      double precision :: aamax, sum, dum
!-----------------------------------------------

      d = 1.0d0
      do i = 1,n
         aamax=0.0d0
         do j = 1,n
            if (dabs(a(i,j)).gt.aamax) aamax = dabs(a(i,j))
         enddo
         if (aamax.eq.0.0d0) then
           if(myPE==PE_IO) write(*,*) &
              'Singular Matrix in ludcmp called from ', calledFrom
           call mfix_exit(myPE)
         endif
         vv(i) = 1.d0/aamax
      enddo
      do j = 1,n
         do i = 1,j-1
            sum = a(i,j)
            do k = 1,i-1
               sum = sum-a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
         enddo
         aamax = 0.0d0
         do i = j,n
            sum = a(i,j)
            do k = 1,j-1
               sum = sum-a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
            dum = vv(i)*dabs(sum)
            if (dum.ge.aamax) then
               imax = i
               aamax = dum
            endif
         enddo
         if (j.ne.imax) then
            do k = 1,n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
            enddo
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j).eq.0.0d0) a(j,j) = TINY
         if (j.ne.n) then
            dum = 1.0d0/a(j,j)
            do i = j+1,n
               a(i,j) = a(i,j)*dum
            enddo
         endif
      enddo

      return
      end subroutine ludcmp


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  subroutine name: lubksb (a,n,np,indx,b)
!  Purpose: solves the set of n linear equations A(n,n).X(n) = B(n).
!
!
!  Literature/Document References:
!     Numerical Recipies in Fortran 77, page 39.
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine lubksb (a,n,np,indx,b)

      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer, intent(in) :: n
      double precision, intent(in) :: a(n,n)
      integer :: np
      integer, intent(in) :: indx(n)
      double precision, intent(inout) :: b(n)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: i, ii, j, ll
      double precision :: sum
!-----------------------------------------------

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
         elseif (sum.ne.0.d0) then
           ii=i
         endif
         b(i)=sum
      enddo
       do i=n,1,-1
         sum=b(i)
         if (i.lt.n) then
           do j=i+1,n
             sum=sum-a(i,j)*b(j)
           enddo
         endif
         b(i)=sum/a(i,i)
       enddo
       return
       end subroutine lubksb

