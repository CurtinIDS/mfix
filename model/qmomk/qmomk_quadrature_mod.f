!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: qmomk_quadrature                                            C
!  Purpose: Routines to perform the quadrature calculations            C
!  Contains the following subroutines:                                 C
!     quadrature_bounded, moments_twenty_eight_nodes,                  C
!     check_moments_twenty, eight_node_3d, bind_theta                  C
!                                                                      C
!  Author: A. Passalacqua                      Date: 18-Jun-2008       C
!  Reviewer:                                Date: dd-mmm-yy            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

#include "version.inc"

MODULE qmomk_quadrature

  USE qmomk_tools
  USE qmomk_parameters

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: MOMENTS_TWENTY_EIGHT_NODES
  PUBLIC :: EIGHT_NODE_3D
  PUBLIC :: BIND_THETA

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Two-nodes quadrature                                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE QUADRATURE_BOUNDED (m, w, u)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(4) :: m
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2) :: w
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2) :: u

    DOUBLE PRECISION :: eps = 1.D-8
    DOUBLE PRECISION :: smax = 1.D6

    DOUBLE PRECISION :: rho = 0.D0
    DOUBLE PRECISION :: um, sigma, x, q, s

    rho = m(1)
    um = m(2)/rho

    sigma = (m(3)*rho - m(2)*m(2))/(rho**2)

    IF (sigma <= 0.D0) THEN
       PRINT *,'QMOMK: Negative variance in quadrature'
       sigma = 0.D0
       x = 0
    ELSE
       sigma = SQRT(sigma)
       q = m(4)/rho - um**3 - 3.D0*um*(sigma**2)
       s = MIN(smax, ABS(q)/(sigma**3))
       x = 0.D0
       IF ((s**2) < eps) THEN
          x = 0.D0
       ELSE
          x = 0.5 * SIGN(1.D0, q)/SQRT(1.D0 + 4.D0/(s**2))
       END IF
    END IF

    w(1) = rho*(0.5 + x)
    u(1) = um - SQRT((0.5-x)/(0.5+x))*sigma

    w(2) = rho*(0.5-x)
    u(2) = um + SQRT((0.5+x)/(0.5-x))*sigma

    IF (ABS(x) >= 0.5D0) THEN
       PRINT *,'QMOMK: Warning: x > 0.5'
       PRINT *,x
    END IF
  END SUBROUTINE QUADRATURE_BOUNDED

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Calculation of the twenty moments from weights and abscissas        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE MOMENTS_TWENTY_EIGHT_NODES (n, u, v, w, mom)

    USE param1, only: small_number
    IMPLICIT NONE

    INTEGER :: i

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: u
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: v
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: w
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NMOM) :: mom

    DOUBLE PRECISION, dimension(QMOMK_NMOM) :: M

    M=0.D0

    DO i = 1, QMOMK_NN
       M(1) = M(1) + n(i)
       M(2) = M(2) + n(i)*u(i)
       M(3) = M(3) + n(i)*v(i)
       M(4) = M(4) + n(i)*w(i)
       M(5) = M(5) + n(i)*(u(i)**2)
       M(6) = M(6) + n(i)*u(i)*v(i)
       M(7) = M(7) + n(i)*u(i)*w(i)
       M(8) = M(8) + n(i)*(v(i)**2)
       M(9) = M(9) + n(i)*v(i)*w(i)
       M(10)= M(10)+ n(i)*w(i)**2
       M(11)= M(11)+ n(i)*(u(i)**3)
       M(12)= M(12)+ n(i)*(u(i)**2)*v(i)
       M(13)= M(13)+ n(i)*(u(i)**2)*w(i)
       M(14)= M(14)+ n(i)*u(i)*(v(i)**2)
       M(15)= M(15)+ n(i)*u(i)*v(i)*w(i)
       M(16)= M(16)+ n(i)*u(i)*(w(i)**2)
       M(17)= M(17)+ n(i)*(v(i)**3)
       M(18)= M(18)+ n(i)*(v(i)**2)*w(i)
       M(19)= M(19)+ n(i)*v(i)*(w(i)**2)
       M(20)= M(20)+ n(i)*(w(i)**3)
    END DO

    IF (M(1) < SMALL_NUMBER) THEN
      PRINT *,'QMOMK: Zero order moment is very small!'
    END IF

    mom = M
  END SUBROUTINE MOMENTS_TWENTY_EIGHT_NODES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Moments check                                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE CHECK_MOMENTS_TWENTY(mom, n, u, v, w, Check)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NMOM) :: mom
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: u
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: v
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: w
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NMOM) :: Check

    DOUBLE PRECISION, DIMENSION(QMOMK_NMOM) :: F

    INTEGER :: i

    F = 0.D0

    F(1) = - mom(1)
    F(2) = - mom(2)
    F(3) = - mom(3)
    F(4) = - mom(4)
    F(5) = - mom(5)
    F(6) = - mom(6)
    F(7) = - mom(7)
    F(8) = - mom(8)
    F(9) = - mom(9)
    F(10)= - mom(10)
    F(11)= - mom(11)
    F(12)= - mom(12)
    F(13)= - mom(13)
    F(14)= - mom(14)
    F(15)= - mom(15)
    F(16)= - mom(16)
    F(17)= - mom(17)
    F(18)= - mom(18)
    F(19)= - mom(19)
    F(20)= - mom(20)

    DO i = 1, QMOMK_NN
       F(1) = F(1) + n(i)
       F(2) = F(2) + n(i)*u(i)
       F(3) = F(3) + n(i)*v(i)
       F(4) = F(4) + n(i)*w(i)
       F(5) = F(5) + n(i)*(u(i)**2)
       F(6) = F(6) + n(i)*u(i)*v(i)
       F(7) = F(7) + n(i)*u(i)*w(i)
       F(8) = F(8) + n(i)*(v(i)**2)
       F(9) = F(9) + n(i)*v(i)*w(i)
       F(10)= F(10)+ n(i)*(w(i)**2)
       F(11)= F(11)+ n(i)*(u(i)**3)
       F(12)= F(12)+ n(i)*(u(i)**2)*v(i)
       F(13)= F(13)+ n(i)*(u(i)**2)*w(i)
       F(14)= F(14)+ n(i)*u(i)*(v(i)**2)
       F(15)= F(15)+ n(i)*u(i)*v(i)*w(i)
       F(16)= F(16)+ n(i)*u(i)*(w(i)**2)
       F(17)= F(17)+ n(i)*(v(i)**3)
       F(18)= F(18)+ n(i)*(v(i)**2)*w(i)
       F(19)= F(19)+ n(i)*v(i)*(w(i)**2)
       F(20)= F(20)+ n(i)*(w(i)**3)
    END DO

    DO i = 1, QMOMK_NMOM
       F(i) = ABS(F(i))
    END DO

    Check = F
  END SUBROUTINE CHECK_MOMENTS_TWENTY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Calculation of weights and abscissas                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE EIGHT_NODE_3D (mom, n, u, v, w)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: mom
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NN) :: n
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NN) :: u
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NN) :: v
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NN) :: w

    DOUBLE PRECISION :: um, vm, wm, s1, s2, s3, s12, s13, s23
    DOUBLE PRECISION :: q000, q101, q110, q011, q200, q020, q002, q210, q300
    DOUBLE PRECISION :: q201, q120, q111, q102, q030, q021, q012, q003

    DOUBLE PRECISION :: b11, b12, b13, b21, b22, b23, b31, b32, b33, Fmax

    DOUBLE PRECISION :: m000, m100, m010, m001, m200, m020, m002, m300, m030, m003, m111

    DOUBLE PRECISION, DIMENSION (QMOMK_NCOV, QMOMK_NCOV) :: sig
    DOUBLE PRECISION, DIMENSION (QMOMK_NCOV, QMOMK_NCOV) :: chol
    DOUBLE PRECISION, DIMENSION (QMOMK_NCOV, QMOMK_NCOV) :: scal, invScale
    DOUBLE PRECISION, DIMENSION (QMOMK_NCOV, QMOMK_NCOV) :: Umat, Vmat

    DOUBLE PRECISION, DIMENSION (4) :: alpha1, alpha2, alpha3
    DOUBLE PRECISION, DIMENSION (2) :: aq1, aq2, aq3, wq1, wq2, wq3

    DOUBLE PRECISION :: a, b, c, u1, u2, v1, v2, w1, w2, rho123

    DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: nstar, ones

    DOUBLE PRECISION, DIMENSION(32) :: x

    DOUBLE PRECISION, DIMENSION(QMOMK_NMOM) :: F

    q000 = mom(1)
    um   = mom(2)/q000
    vm   = mom(3)/q000
    wm   = mom(4)/q000
    q200 = mom(5)/q000
    q110 = mom(6)/q000
    q101 = mom(7)/q000
    q020 = mom(8)/q000
    q011 = mom(9)/q000
    q002 = mom(10)/q000
    q300 = mom(11)/q000
    q210 = mom(12)/q000
    q201 = mom(13)/q000
    q120 = mom(14)/q000
    q111 = mom(15)/q000
    q102 = mom(16)/q000
    q030 = mom(17)/q000
    q021 = mom(18)/q000
    q012 = mom(19)/q000
    q003 = mom(20)/q000

    sig = 0.D0

    sig(1,1) = q200 - um**2
    sig(2,2) = q020 - vm**2
    sig(3,3) = q002 - wm**2

    IF ( sig(1,1) < 0.D0 ) THEN ! check that u variance is nonnegative
       PRINT *,'QMOMK: Negative u variance',sig
       sig(1,1) = 0.D0
       ERROR_STOP
    END IF

    IF ( sig(2,2) < 0.D0 ) THEN ! check that v variance is nonnegative
       PRINT *,'QMOMK: Negative v variance',sig
       sig(2,2) = 0.D0
       ERROR_STOP
    END IF

    IF ( sig(3,3) < 0.D0 ) THEN ! check that v variance is nonnegative
       PRINT *,'QMOMK: Negative w variance',sig
       sig(3,3) = 0.D0
       ERROR_STOP
    END IF

    s1 = SQRT(sig(1,1))
    s2 = SQRT(sig(2,2))
    s3 = SQRT(sig(3,3))

    sig(1,2) = q110 - um*vm
    sig(2,1) = sig(1,2)
    s12 = sig(1,2)

    sig(1,3) = q101 - um*wm
    sig(3,1) = sig(1,3)
    s13 = sig(1,3)

    sig(2,3) = q011 - vm*wm
    sig(3,2) = sig(2,3)
    s23 = sig(2,3)

    IF ( s1*s2 < ABS(s12) ) THEN ! check that covariance matrix is nonnegative
       PRINT *,'QMOMK: unphysical uv correlation', sig(1,2), s1, s2
       sig(1,2) = SIGN(1.D0, s12)*s1*s2
       sig(2,1) = sig(1,2)

       !q110 = 0.99*q110
       !mom(6) = q000*q110
       !sig(1,2) = q110 - um*vm
       !sig(2,1) = sig(1,2)
       PRINT *,'Sigma = ',sig
       PRINT *,'M = ', mom
       PRINT *,'N = ', n
       PRINT *,'U = ', u
       PRINT *,'V = ', v
       PRINT *,'W = ', w
       ERROR_STOP
    END IF

    IF ( s1*s3 < ABS(s13) ) THEN ! check that covariance matrix is nonnegative
       sig(1,3) = SIGN(1.D0, s13)*s1*s3
       sig(3,1) = sig(1,3)
       PRINT *,'QMOMK: unphysical uw correlation',sig(1,3)
       ERROR_STOP
    END IF

    IF ( s2*s3 < ABS(s23) ) THEN ! check that covariance matrix is nonnegative
       sig(2,3) = SIGN(1.D0, s23)*s2*s3
       sig(3,2) = sig(2,3)
       PRINT *,'QMOMK: unphysical vw correlation'
       ERROR_STOP
    END IF

    CALL CHOLESKY3(sig, chol) ! Bottom cholesky decomposition
    CALL DIAG3(s1, s2, s3, scal)
    CALL INV3(scal, invScale)
    CALL MULTIPLYMATRIX3(chol, invScale, Umat)
    CALL INV3(Umat, Vmat)

    b11 = Vmat(1,1)
    b12 = Vmat(1,2)
    b13 = Vmat(1,3)

    b21 = Vmat(2,1)
    b22 = Vmat(2,2)
    b23 = Vmat(2,3)

    b31 = Vmat(3,1)
    b32 = Vmat(3,2)
    b33 = Vmat(3,3)

    m000 = 1
    m100 = 0
    m010 = 0
    m001 = 0

    m200 = 2.D0*b11*b12*q110 - 2.D0*b11*um*b12*vm + 2.D0*b11*b13*q101 - &
         2.D0*b11*um*b13*wm + 2.D0*b12*b13*q011 - 2.D0*b12*vm*b13*wm - &
         (b11**2)*(um**2) - (b12**2)*(vm**2) - (b13**2)*(wm**2) + &
         (b11**2)*q200 + (b12**2)*q020 + (b13**2)*q002

    m020 = - (b21**2)*(um**2) - (b22**2)*(vm**2) - (b23**2)*(wm**2) + &
         2.D0*b21*b23*q101 - 2.D0*b21*um*b23*wm + 2.D0*b22*b23*q011 - &
         2.D0*b22*vm*b23*wm - 2.D0*b21*um*b22*vm + 2.D0*b21*b22*q110 + &
         (b21**2)*q200 + (b22**2)*q020 + (b23**2)*q002

    m002 = - (b33**2)*(wm**2) - (b31**2)*(um**2) - (b32**2)*(vm**2) - &
         2.D0*b31*um*b32*vm + 2.D0*b31*b32*q110 + (b31**2)*q200 + &
         (b32**2)*q020 + (b33**2)*q002 + 2.D0*b31*b33*q101 - 2.D0*b31*um*b33*wm + &
         2.D0*b32*b33*q011 - 2.D0*b32*vm*b33*wm

    m300 = - 3.D0*(b13**3)*wm*q002 + 3.D0*b11*(b12**2)*q120 + 3.D0*(b12**2)*b13*q021 + &
         3.D0*(b11**2)*b12*q210 + 3.D0*(b11**2)*b13*q201 + 3.D0*b11*(b13**2)*q102 + &
         3.D0*b12*(b13**2)*q012 + 6.D0*b12*vm*(b13**2)*(wm**2) - 3.D0*(b12**3)*vm*q020 - &
         3.D0*(b11**3)*um*q200 - 6.D0*b11*b12*b13*wm*q110 - 6.D0*b11*b12*vm*b13*q101 + &
         (b11**3)*q300 + 2.D0*(b11**3)*(um**3) - 6.D0*b11*um*b12*b13*q011 + 2.D0*(b12**3)*(vm**3) + &
         2.D0*(b13**3)*(wm**3) + (b13**3)*q003 + (b12**3)*q030 - 3.D0*b11*um*(b12**2)*q020 - &
         3.D0*b11*um*(b13**2)*q002 - 3.D0*(b12**2)*b13*wm*q020 + 6.D0*(b11**2)*(um**2)*b12*vm + &
         6.D0*(b12**2)*(vm**2)*b13*wm + 6.D0*(b11**2)*(um**2)*b13*wm + 6.D0*b11*um*(b12**2)*(vm**2) + &
         6.D0*b11*um*(b13**2)*(wm**2) - 6.D0*b11*(b12**2)*vm*q110 - 6.D0*b11*(b13**2)*wm*q101 - &
         6.D0*(b12**2)*vm*b13*q011 + 12.D0*b11*um*b12*vm*b13*wm - 6.D0*b12*(b13**2)*wm*q011 - &
         6.D0*(b11**2)*um*b13*q101 + 6.D0*b11*b12*b13*q111 - 3.D0*b12*vm*(b13**2)*q002 - &
         6.D0*(b11**2)*um*b12*q110 - 3.D0*(b11**2)*b13*wm*q200 - 3.D0*(b11**2)*b12*vm*q200

    m030 = 6.D0*(b21**2)*(um**2)*b22*vm + 6.D0*(b22**2)*(vm**2)*b23*wm - 3.D0*(b21**3)*um*q200 + &
         2.D0*(b22**3)*(vm**3) + 3.D0*(b21**2)*b23*q201 + 2.D0*(b21**3)*(um**3) - 3.D0*(b23**3)*wm*q002 - &
         3.D0*(b21**2)*b22*vm*q200 + 6.D0*b21*um*(b22**2)*(vm**2) + 3.D0*b21*(b22**2)*q120 + &
         3.D0*(b21**2)*b22*q210 + 6.D0*(b21**2)*(um**2)*b23*wm - 3.D0*b21*um*(b22**2)*q020 + &
         3.D0*b21*(b23**2)*q102 + 3.D0*(b22**2)*b23*q021 + 6.D0*b22*vm*(b23**2)*(wm**2) - &
         3.D0*(b21**2)*b23*wm*q200 + (b22**3)*q030 + (b21**3)*q300 - 3.D0*(b22**3)*vm*q020 + &
         6.D0*b21*b22*b23*q111 - 6.D0*(b21**2)*um*b22*q110 + (b23**3)*q003 - 3.D0*b22*vm*(b23**2)*q002 - &
         6.D0*(b21**2)*um*b23*q101 + 3.D0*b22*(b23**2)*q012 - 6.D0*(b22**2)*vm*b23*q011 - &
         3.D0*b21*um*(b23**2)*q002 - 3.D0*(b22**2)*b23*wm*q020 - 6.D0*b22*(b23**2)*wm*q011 - &
         6.D0*b21*b22*b23*wm*q110 - 6.D0*b21*um*b22*b23*q011 + 6.D0*b21*um*(b23**2)*(wm**2) + &
         2.D0*(b23**3)*(wm**3) - 6.D0*b21*(b23**2)*wm*q101 - 6.D0*b21*b22*vm*b23*q101 + &
         12.D0*b21*um*b22*vm*b23*wm - 6.D0*b21*(b22**2)*vm*q110

    m003 = (b31**3)*q300 + 6.D0*b31*um*(b33**2)*(wm**2) + (b32**3)*q030 + (b33**3)*q003 - &
         3.D0*b31*um*(b33**2)*q002  + 6.D0*b32*vm*(b33**2)*(wm**2) + 6.D0*b31*um*(b32**2)*(vm**2) + &
         2.D0*(b31**3)*(um**3) + 2.D0*(b32**3)*(vm**3) + 3.D0*(b31**2)*b33*q201 - &
         3.D0*b32*vm*(b33**2)*q002 - 6.D0*b31*(b32**2)*vm*q110 - 3.D0*(b31**2)*b32*vm*q200 - &
         3.D0*(b31**2)*b33*wm*q200 - 6.D0*b31*(b33**2)*wm*q101 - 3.D0*(b33**3)*wm*q002 + &
         6.D0*(b31**2)*(um**2)*b33*wm + 3.D0*b32*(b33**2)*q012 + 6.D0*b31*b32*b33*q111 - &
         6.D0*(b31**2)*um*b33*q101 - 6.D0*b32*(b33**2)*wm*q011 + 6.D0*(b31**2)*(um**2)*b32*vm - &
         3.D0*b31*um*(b32**2)*q020 - 6.D0*b31*b32*b33*wm*q110 + 3.D0*b31*(b33**2)*q102 - &
         3.D0*(b32**2)*b33*wm*q020 + 12.D0*b31*um*b32*vm*b33*wm + 3.D0*(b31**2)*b32*q210 - &
         6.D0*(b32**2)*vm*b33*q011 + 3.D0*b31*(b32**2)*q120 + 6.D0*(b32**2)*(vm**2)*b33*wm + &
         3.D0*(b32**2)*b33*q021 - 3.D0*(b31**3)*um*q200 - 3.D0*(b32**3)*vm*q020 - &
         6.D0*b31*um*b32*b33*q011 + 2.D0*(b33**3)*(wm**3) - 6.D0*b31*b32*vm*b33*q101 - &
         6.D0*(b31**2)*um*b32*q110

    m111 = - b12*b23*wm*b31*q110 + 2.D0*b11*(um**2)*b22*vm*b31 - b12*b21*b33*wm*q110 - b12*vm*b23*b33*q002 + &
         2.D0*b12*vm*b23*wm*b31*um - b12*vm*b21*b31*q200 - b13*b23*b31*um*q002 - b12*b22*b31*um*q020 + &
         b12*b22*b31*q120 - b11*b21*b33*wm*q200 + 2.D0*b12*(vm**2)*b23*wm*b32 - b12*b23*wm*b32*q020 - &
         b11*b21*b32*vm*q200 - b12*b21*um*b32*q020 - 2.D0*b12*b22*vm*b33*q011 - b11*um*b23*b33*q002 - &
         2.D0*b12*b22*vm*b31*q110 - 3.D0*b11*b21*b31*um*q200 + b11*b21*b32*q210 + b11*b21*b31*q300 - &
         2.D0*b12*b23*b33*wm*q011 + b13*b23*b33*q003 + b12*b23*b33*q012 + b12*b22*b33*q021 - &
         b12*b23*b31*um*q011 - b12*b22*b33*wm*q020 + b13*b23*b31*q102 - 2.D0*b12*b21*b31*um*q110 + &
         b12*b21*b32*q120 + b12*b23*b31*q111 - 2.D0*b12*b21*b32*vm*q110 + 2.D0*b12*(vm**2)*b21*um*b32 - &
         3.D0*b12*b22*b32*vm*q020 + 2.D0*b12*(vm**2)*b22*b31*um + 2.D0*b12*vm*b21*(um**2)*b31 - &
         2.D0*b12*b23*b32*vm*q011 - b12*b21*um*b33*q011 + 2.D0*b12*(vm**2)*b22*b33*wm + b12*b23*b32*q021 + &
         b11*b21*b33*q201 - 3.D0*b13*b23*b33*wm*q002 - 2.D0*b11*b21*um*b32*q110 + 2.D0*b12*vm*b23*(wm**2)*b33 - &
         2.D0*b13*b21*b33*wm*q101 + 2.D0*b13*(wm**2)*b22*vm*b33 - b13*b21*um*b32*q011 - b12*vm*b21*b33*q101 - &
         2.D0*b11*b21*um*b33*q101 - b13*b22*vm*b33*q002 - b13*b21*um*b33*q002 + 2.D0*b11*(um**3)*b21*b31 + &
         2.D0*b11*b21*(um**2)*b32*vm + b13*b21*b31*q201 + b13*b22*b33*q012 + b13*b22*b31*q111 + &
         2.D0*b13*b23*(wm**2)*b31*um - b13*b23*b32*vm*q002 + b11*b23*b31*q201 - b11*b22*vm*b31*q200 + &
         2.D0*b13*wm*b21*(um**2)*b31 - b13*b22*vm*b31*q101 + 2.D0*b11*b21*(um**2)*b33*wm + b13*b23*b32*q012 - &
         2.D0*b11*b22*b31*um*q110 - 2.D0*b13*b22*b32*vm*q011 + b12*b21*b31*q210 - b13*wm*b22*b32*q020 - &
         2.D0*b13*b23*wm*b31*q101 + b11*b22*b33*q111 + 2.D0*b13*(wm**2)*b21*um*b33 - b13*b21*b32*vm*q101 + &
         b12*b22*b32*q030 + 2.D0*b13*b23*(wm**2)*b32*vm - b11*b22*b33*wm*q110 - b12*vm*b23*b31*q101 - &
         b13*wm*b21*b32*q110 + 2.D0*b11*um*b23*wm*b32*vm + 2.D0*b11*um*b22*(vm**2)*b32 + 2*b13*(wm**3)*b23*b33 + &
         2.D0*b12*(vm**3)*b22*b32 - 2.D0*b11*b23*b31*um*q101 - b11*b22*vm*b33*q101 - b11*um*b22*b32*q020 + &
         2.D0*b13*wm*b22*(vm**2)*b32 + b13*b21*b33*q102 + b13*b21*b32*q111 - 2.D0*b13*b21*b31*um*q101 - &
         b13*b22*b31*um*q011 + b11*b22*b32*q120 - b13*wm*b22*b31*q110 + 2.D0*b11*um*b23*(wm**2)*b33 + &
         2.D0*b11*(um**2)*b23*wm*b31 + b11*b23*b32*q111 + b13*b22*b32*q021 + 2.D0*b12*vm*b21*um*b33*wm + &
         2.D0*b13*wm*b22*vm*b31*um - 2.D0*b13*b23*wm*b32*q011 + b12*b21*b33*q111 - 2.D0*b11*b23*b33*wm*q101 + &
         2.D0*b11*um*b22*vm*b33*wm - b11*um*b23*b32*q011 - b11*um*b22*b33*q011 - 2.D0*b11*b22*b32*vm*q110 + &
         b11*b22*b31*q210 - b13*wm*b21*b31*q200 - 2.D0*b13*b22*b33*wm*q011 - b11*b23*wm*b32*q110 + &
         b11*b23*b33*q102 - b11*b23*wm*b31*q200 + 2.D0*b13*wm*b21*um*b32*vm - b11*b23*b32*vm*q101

    alpha1(1) = m000
    alpha1(2) = m100
    alpha1(3) = m200
    alpha1(4) = m300

    CALL QUADRATURE_BOUNDED(alpha1, wq1, aq1)

    alpha2(1) = m000
    alpha2(2) = m010
    alpha2(3) = m020
    alpha2(4) = m030

    call QUADRATURE_BOUNDED(alpha2, wq2, aq2)

    alpha3(1) = m000
    alpha3(2) = m001
    alpha3(3) = m002
    alpha3(4) = m003

    call QUADRATURE_BOUNDED(alpha3, wq3, aq3)

    a = wq1(1)
    b = wq2(1)
    c = wq3(1)

    u1 = aq1(1)
    u2 = aq1(2)

    v1 = aq2(1)
    v2 = aq2(2)

    w1 = aq3(1)
    w2 = aq3(2)

    rho123 = SQRT(a*(1-a)*b*(1-b)*c*(1-c))*m111/SQRT(m200*m020*m002)

    IF (rho123 > 0.D0) THEN
       rho123 = MIN(rho123 , a*b*c , (1-a)*(1-b)*c , (1-a)*b*(1-c) , a*(1-b)*(1-c))
    ELSE
       rho123 = -MIN(-rho123 , (1-a)*b*c , a*(1-b)*c , a*b*(1-c) , (1-a)*(1-b)*(1-c))
    END IF

    IF (a > 1.) THEN
       PRINT *,a
    END IF
    IF (b > 1.) THEN
       PRINT *,b
    END IF
    IF (c > 1.) THEN
       PRINT *,c
    END IF

    nstar(1) = a*b*c                    - rho123
    nstar(2) = (1.D0-a)*b*c             + rho123
    nstar(3) = a*(1.D0-b)*c             + rho123
    nstar(4) = (1.D0-a)*(1.D0-b)*c      - rho123
    nstar(5) = a*b*(1.D0-c)             + rho123
    nstar(6) = (1.D0-a)*b*(1.D0-c)      - rho123
    nstar(7) = a*(1.D0-b)*(1.D0-c)      - rho123
    nstar(8) = (1.D0-a)*(1.D0-b)*(1.D0-c)   + rho123

    IF (MINVAL(nstar) < -1.D-16) THEN
       PRINT *,'QMOMK: Negative weight in quadrature'
       PRINT *,nstar
       nstar = ABS(nstar)
       ERROR_STOP
    END IF

    x(1:8) = nstar
    x(9:16) = (/u1, u2, u1, u2, u1, u2, u1, u2/)
    x(17:24) = (/v1, v1, v2, v2, v1, v1, v2, v2/)
    x(25:32) = (/w1, w1, w1, w1, w2, w2, w2, w2/)

    ones = (/1.,1.,1.,1.,1.,1.,1.,1./)

    n = x(1:8)*q000
    u = Umat(1,1)*x(9:16) + Umat(1,2)*x(17:24) + Umat(1,3)*x(25:32) + um*ones(1:8)
    v = Umat(2,1)*x(9:16) + Umat(2,2)*x(17:24) + Umat(2,3)*x(25:32) + vm*ones(1:8)
    w = Umat(3,1)*x(9:16) + Umat(3,2)*x(17:24) + Umat(3,3)*x(25:32) + wm*ones(1:8)

    CALL CHECK_MOMENTS_TWENTY(mom, n, u, v, w, F)
    Fmax = MAXVAL(F(1:10))

    IF (Fmax > 1.D-10) THEN
       PRINT *,'QMOMK: First 10 moments not satisfied.'
       PRINT *,'F(1:10) = ',F(1:10)
       PRINT *,'M(1:10) = ',mom(1:10)
       ERROR_STOP
    ENDIF
  END SUBROUTINE EIGHT_NODE_3D


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Temperature limiter                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE BIND_THETA (mom, tau_p)

    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: mom
    DOUBLE PRECISION, INTENT(IN) :: tau_p

    DOUBLE PRECISION :: sigma2

    sigma2 = MIN(MINIMUM_THETA/tau_p, MAXIMUM_SIGMA)

    ! Second order moments
    mom(5) = mom(5) + 2.0*sigma2*mom(1)
    mom(8) = mom(8) + 2.0*sigma2*mom(1)
    mom(10) = mom(10) + 2.0*sigma2*mom(1)

    ! Third order moments
    mom(11) = mom(11) + 6.0*sigma2*mom(2)
    mom(12) = mom(12) + 2.0*sigma2*mom(3)
    mom(13) = mom(13) + 2.0*sigma2*mom(4)
    mom(14) = mom(14) + 2.0*sigma2*mom(2)

    mom(16) = mom(16) + 2.0*sigma2*mom(2)
    mom(17) = mom(17) + 6.0*sigma2*mom(3)
    mom(18) = mom(18) + 2.0*sigma2*mom(4)

    mom(19) = mom(19) + 2.0*sigma2*mom(3)
    mom(20) = mom(20) + 6.0*sigma2*mom(4)

  END SUBROUTINE BIND_THETA

END MODULE qmomk_quadrature
