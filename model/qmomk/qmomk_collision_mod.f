!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module:  qmomk_collisions                                           C
!  Purpose: Collision operators for the kinetic equation               C
!  Contains the following subroutines and functions:                   C
!      collisions_istantaneous, collisions_bgk,                        C
!      compute_collision_time, collisions_boltzmann_one_specie,        C
!      collisions_boltzmann_two_specie,                                C
!      solve_boltzmann_collisions_one_specie,                          C
!      solve_boltzmann_collisions_two_specie, radial_g0                C
!                                                                      C
!  Author: A. Passalacqua                           Date:              C
!  Reviewer:                                        Date:              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

#include "version.inc"

MODULE qmomk_collision


  USE qmomk_parameters
  USE qmomk_quadrature

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: COLLISIONS_ISTANTANEOUS
  PUBLIC :: COLLISIONS_BGK
  PUBLIC :: COMPUTE_COLLISION_TIME
  !PUBLIC :: GLOBAL_COLLISION_TIME
  PUBLIC :: SOLVE_BOLTZMANN_COLLISIONS_ONE_SPECIE
  PUBLIC :: SOLVE_BOLTZMANN_COLLISIONS_TWO_SPECIES

  PUBLIC :: RADIAL_G0

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Instantaneous collisions                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE COLLISIONS_ISTANTANEOUS(M, dp)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: M
    DOUBLE PRECISION, INTENT(IN) :: dp

    DOUBLE PRECISION :: m0, m1, m2, m3, m11, m22, m33

    DOUBLE PRECISION :: up1, up2, up3
    DOUBLE PRECISION :: sig1, sig2, sig3, seq

    DOUBLE PRECISION :: d11, d12, d13, d22, d23, d33
    DOUBLE PRECISION :: d111, d112, d113, d122, d123, d133, d222, d223, d233, d333

    m0 = M(1)
    m1 = M(2)
    m2 = M(3)
    m3 = M(4)
    m11 = M(5)
    m22 = M(8)
    m33 = M(10)

    up1 = m1/m0
    up2 = m2/m0
    up3 = m3/m0
    sig1 = max(0.D0, m11/m0 - up1**2)
    sig2 = max(0.D0, m22/m0 - up2**2)
    sig3 = max(0.D0, m33/m0 - up3**2)
    seq = ( sig1 + sig2 + sig3 )/3.D0

    d11 = m0*(seq + up1**2)
    d12 = m0*up1*up2
    d13 = m0*up1*up3
    d22 = m0*(seq + up2**2)
    d23 = m0*up2*up3
    d33 = m0*(seq + up3**2)
    d111 = m1*(3.D0*seq + up1**2)
    d112 = m2*(seq + up1**2)
    d113 = m3*(seq + up1**2)
    d122 = m1*(seq + up2**2)
    d123 = m0*up1*up2*up3
    d133 = m1*(seq + up3**2)
    d222 = m2*(3.D0*seq + up2**2)
    d223 = m3*(seq + up2**2)
    d233 = m2*(seq + up3**2)
    d333 = m3*(3.D0*seq + up3**2)

    M(5)  =  d11
    M(6)  =  d12
    M(7)  =  d13
    M(8)  =  d22
    M(9)  =  d23
    M(10) =  d33
    M(11) =  d111
    M(12) =  d112
    M(13) =  d113
    M(14) =  d122
    M(15) =  d123
    M(16) =  d133
    M(17) =  d222
    M(18) =  d223
    M(19) =  d233
    M(20) =  d333

  END SUBROUTINE COLLISIONS_ISTANTANEOUS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Bathnagar-Gross-Krook (BGK) collision operator                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE COLLISIONS_BGK(M, Dt, tcol, dp, e)

    USE param1, only: small_number
    USE constant, only: Pi, EP_STAR

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: M
    DOUBLE PRECISION, INTENT(IN) :: Dt, tcol, dp, e

    DOUBLE PRECISION :: m0, m1, m2, m3, m11, m12, m13, m22, m23, m33
    DOUBLE PRECISION :: m111, m112, m113, m122, m123, m133, m222
    DOUBLE PRECISION :: m223, m233, m333

    DOUBLE PRECISION :: up1, up2, up3, b, a1, b1, gamma, om
    DOUBLE PRECISION :: sig1, sig2, sig3, sig11, sig12, sig13, sig22, sig23
    DOUBLE PRECISION :: sig33, seq, T

    DOUBLE PRECISION :: d11, d12, d13, d22, d23, d33
    DOUBLE PRECISION :: d111, d112, d113, d122, d123, d133, d222, d223, d233, d333

    DOUBLE PRECISION :: kap, ALPHA_MAX

    m0   = M(1)
    m1   = M(2)
    m2   = M(3)
    m3   = M(4)
    m11  = M(5)
    m12  = M(6)
    m13  = M(7)
    m22  = M(8)
    m23  = M(9)
    m33  = M(10)
    m111 = M(11)
    m112 = M(12)
    m113 = M(13)
    m122 = M(14)
    m123 = M(15)
    m133 = M(16)
    m222 = M(17)
    m223 = M(18)
    m233 = M(19)
    m333 = M(20)

    up1 = m1/m0
    up2 = m2/m0
    up3 = m3/m0
    sig1 = MAX(0.D0, m11/m0 - up1**2)
    sig2 = MAX(0.D0, m22/m0 - up2**2)
    sig3 = MAX(0.D0, m33/m0 - up3**2)
    T = ( sig1 + sig2 + sig3 )/3.D0

    !A.P. This parameter is hardcoded:
    !     b = 0 for BGK, b = -1/2 for ES-BGK
    b = 0.D0
    gamma = 1.D0 - b

    om = (1 + e)/2.D0
    a1 = gamma*(om**2)
    b1 = a1 - 2.D0*gamma*om + 1.D0

    sig11 = a1*T + b1*sig1
    sig22 = a1*T + b1*sig2
    sig33 = a1*T + b1*sig3

    sig12 = b1*(m12/m0 - up1*up2)
    sig13 = b1*(m13/m0 - up1*up3)
    sig23 = b1*(m23/m0 - up2*up3)

    seq = ( sig11 + sig22 + sig33 )/3.D0

    d11 = m0*(sig11 + up1**2)
    d12 = m0*(sig12 + up1*up2)
    d13 = m0*(sig13 + up1*up3)
    d22 = m0*(sig22 + up2**2)
    d23 = m0*(sig23 + up2*up3)
    d33 = m0*(sig33 + up3**2)

    d111 = -(-3.D0*d11 + 2.D0*m0*up1**2)*up1
    d112 = -(-d11*up2 - (-2.D0*d12 + 2.D0*up1*up2*m0)*up1)

    d113 = -(-d11*up3+(-2.D0*d13+2.D0*up1*up3*m0)*up1)
    d122 = -(-2.D0*d12*up2+(2.D0*(up2**2)*m0-d22)*up1)
    d123 = -(-d12*up3-d13*up2+(2.D0*up2*up3*m0-d23)*up1)
    d133 = -(-2.D0*d12*up3+(2.D0*(up3**2)*m0-d33)*up1)
    d222 = -(-3.D0*d22+2.D0*(up2**2)*m0)*up2
    d223 = -(-d22*up2+(-2.D0*d23+2.D0*up2*up3*m0)*up2)
    d233 = -(-2.D0*d23*up3+(2.D0*(up3**2)*m0-d33)*up2)
    d333 = -(-3.D0*d33+2.D0*(up3**2)*m0)*up3

! JEC: Probable BUG: alpha_max was not defined. include definition as found in
! other places
    ALPHA_MAX = 1. - EP_STAR

    IF (m0 > ALPHA_MAX - SMALL_NUMBER) THEN
      kap = 0.
    ELSE
      kap = EXP(-12.0*RADIAL_G0(m0)*m0*Dt*sqrt(T/pi)/(dp*gamma))
    END IF

    M(5)  =  d11  - kap*( d11 - m11 )
    M(6)  =  d12  - kap*( d12 - m12 )
    M(7)  =  d13  - kap*( d13 - m13 )
    M(8)  =  d22  - kap*( d22 - m22 )
    M(9)  =  d23  - kap*( d23 - m23 )
    M(10) =  d33  - kap*( d33 - m33 )
    M(11) =  d111 - kap*( d111 - m111 )
    M(12) =  d112 - kap*( d112 - m112 )
    M(13) =  d113 - kap*( d113 - m113 )
    M(14) =  d122 - kap*( d122 - m122 )
    M(15) =  d123 - kap*( d123 - m123 )
    M(16) =  d133 - kap*( d133 - m133 )
    M(17) =  d222 - kap*( d222 - m222 )
    M(18) =  d223 - kap*( d223 - m223 )
    M(19) =  d233 - kap*( d233 - m233 )
    M(20) =  d333 - kap*( d333 - m333 )
  END SUBROUTINE COLLISIONS_BGK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Collision time                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE COMPUTE_COLLISION_TIME (M, dp, theta, tcol, dt)

    USE param1, only: small_number
    USE constant, only: Pi, EP_STAR
! JEC: should ep_star be replaced by ep_star_array?

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NMOM) :: M
    DOUBLE PRECISION, INTENT(IN) :: theta, dp, dt
    DOUBLE PRECISION, INTENT(OUT) :: tcol

    DOUBLE PRECISION :: ALPHA_MAX

    !A.P.: Assuming m0 = M(1) = solids volume fraction due to the moments scaling
    ALPHA_MAX = 1. - EP_STAR

    IF (M(1) > ALPHA_MAX - SMALL_NUMBER) THEN
      tcol = dt ! If this condition is satisfied, moments are set to
                ! equilibrium. No need to integrate on more than one time step
    ELSE
      tcol = dp*SQRT(pi/(theta+SMALL_NUMBER))/(12.0*M(1)*RADIAL_G0(M(1))+SMALL_NUMBER)
    END IF
  END SUBROUTINE COMPUTE_COLLISION_TIME

  !     Global collision time based on the mean volume fraction (1 specie)
  !SUBROUTINE GLOBAL_COLLISION_TIME (alpha, dp, theta, tcol, dt)
  !  DOUBLE PRECISION, INTENT(IN) :: alpha
  !  DOUBLE PRECISION, INTENT(IN) :: theta, dp, dt
  !  DOUBLE PRECISION, INTENT(OUT) :: tcol

  !  DOUBLE PRECISION :: ALPHA_MAX

  !  ALPHA_MAX = 1. - EP_STAR

  ! IF (alpha > ALPHA_MAX - SMALL_NUMBER) THEN
  !    tcol = dt ! If this condition is satisfied, moments are set to
                ! equilibrium. No need to integrate on more than one time step
  !  ELSE
  !    tcol = dp*SQRT(pi/(theta+SMALL_NUMBER))/(12.0*alpha*RADIAL_G0(M(1))+SMALL_NUMBER)
  !  END IF
  !END SUBROUTINE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Boltzmann collision operator for the monodispersed case             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE COLLISIONS_BOLTZMANN_ONE_SPECIE (N, U, V, W, dp, e, Coll)

    USE constant, only: Pi

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: N, U, V, W
    DOUBLE PRECISION, INTENT(IN) :: dp, e
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NMOM) :: Coll

    DOUBLE PRECISION, DIMENSION(3) :: a1, a, vr, v1r, v2r, g1, g2

    DOUBLE PRECISION :: C200, C110, C101, C020, C011, C002, C300, C210, C201
    DOUBLE PRECISION :: C120, C111, C102, C030, C021, C012, C003
    DOUBLE PRECISION :: d, g, gsqr, o, ni, ui, vi, wi, nj, uj, vj, wj, cross
    DOUBLE PRECISION :: denom, L11, L12, L13, L21, L22, L23, L31, L32, L33

    INTEGER :: i, j

    a1 = (/0.0462,0.0971,0.8235/) ! chosen at random (cannot be colinear with vr)
    a = 1/((a1(1))**2 + (a1(2))**2 + (a1(3))**2) * a1;

    C200 = 0.
    C110 = 0.
    C101 = 0.
    C020 = 0.
    C011 = 0.
    C002 = 0.
    C300 = 0.
    C210 = 0.
    C201 = 0.
    C120 = 0.
    C111 = 0.
    C102 = 0.
    C030 = 0.
    C021 = 0.
    C012 = 0.
    C003 = 0.

    o = (1+e)/2
    DO i= 1, QMOMK_NN
       ni = N(i)
       ui = U(i)
       vi = V(i)
       wi = W(i)
       v1r(1) = ui
       v1r(2) = vi
       v1r(3) = wi
       DO j = 1, 8
          IF (i == j) THEN
             CYCLE ! vi - vj = 0 always => nothing to add to source terms
          END IF
          nj = N(j)
          uj = U(j)
          vj = V(j)
          wj = W(j)

          v2r(1) = uj
          v2r(2) = vj
          v2r(3) = wj

          vr(1) = v1r(1) - v2r(1)
          vr(2) = v1r(2) - v2r(2)
          vr(3) = v1r(3) - v2r(3)

          g = SQRT( vr(1)**2 + vr(2)**2 + vr(3)**2 )

          g1(1) = vr(1)/g
          g1(2) = vr(2)/g
          g1(3) = vr(3)/g

          cross = a(1)*g1(1) + a(2)*g1(2) + a(3)*g1(3)

          g2(1) = a(1) - cross*g1(1)
          g2(2) = a(2) - cross*g1(2)
          g2(3) = a(3) - cross*g1(3)

          denom = SQRT( g2(1)**2 + g2(2)**2 + g2(3)**2 )

          L11 = g2(1)/denom
          L12 = g2(2)/denom
          L13 = g2(3)/denom
          L21 = ( g1(3)*a(2)-g1(2)*a(3) )/denom
          L22 = ( g1(1)*a(3)-g1(3)*a(1) )/denom
          L23 = ( g1(2)*a(1)-g1(1)*a(2) )/denom
          L31 = g1(1)
          L32 = g1(2)
          L33 = g1(3)

          gsqr = g**2

          C200 = C200 + ni*nj*(pi*gsqr*o*(4*g*o*(L31**2)+6*uj*L31-6*ui*L31+g*o*(L21**2)+g*o*(L11**2)))/6

          C110 = C110 + ni*nj*(pi*gsqr*o*(4*g*o*L31*L32+3*uj*L32-3*ui*L32+3*vj*L31-3*vi*L31+g*o*L21*L22 + &
               g*o*L11*L12))/6

          C101 = C101 + ni*nj*(pi*gsqr*o*(4*g*o*L31*L33+3*uj*L33-3*ui*L33+3*wj*L31-3*wi*L31+g*o*L21*L23 + &
               g*o*L11*L13))/6

          C020 = C020 + ni*nj*(pi*gsqr*o*(4*g*o*(L32**2)+6*vj*L32-6*vi*L32+g*o*(L22**2)+g*o*(L12**2)))/6

          C011 = C011 + ni*nj*(pi*gsqr*o*(4*g*o*L32*L33+3*vj*L33-3*vi*L33+3*wj*L32-3*wi*L32+g*o*L22*L23 + &
               g*o*L12*L13))/6

          C002 = C002 + ni*nj*(pi*gsqr*o*(4*g*o*(L33**2)+6*wj*L33-6*wi*L33+g*o*(L23**2)+g*o*(L13**2)))/6

          C300 = C300 + ni*nj*(pi*gsqr*o*(uj+ui)*(4*g*o*(L31**2)+6*uj*L31-6*ui*L31+g*o*(L21**2)+g*o*(L11**2)))/4

          C210 = C210 + ni*nj*(pi*gsqr*o*(8*g*o*uj*L31*L32+8*g*o*ui*L31*L32+6*(uj**2)*L32-6*(ui**2)*L32 + &
               4*g*o*vj*(L31**2)+4*g*o*vi*(L31**2)+12*uj*vj*L31-12*ui*vi*L31+2*g*o*uj*L21*L22 + &
               2*g*o*ui*L21*L22+g*o*vj*(L21**2)+g*o*vi*(L21**2)+2*g*o*uj*L11*L12+2*g*o*ui*L11*L12 + &
               g*o*vj*(L11**2)+g*o*vi*(L11**2)))/12

          C201 = C201 + ni*nj*(((8*pi*(g**3)*(o**2)*uj+8*pi*(g**3)*(o**2)*ui)*L31+6*pi*gsqr*o*(uj**2)-6*pi*gsqr*o*(ui**2))*L33 + &
               (4*pi*(g**3)*(o**2)*wj+4*pi*(g**3)*(o**2)*wi)*(L31**2)+(12*pi*(g**2)*o*uj*wj-12*pi*(g**2)*o*ui*wi)*L31 + &
               (2*pi*(g**3)*(o**2)*uj+2*pi*(g**3)*(o**2)*ui)*L21*L23+(pi*(g**3)*(o**2)*wj+pi*(g**3)*(o**2)*wi)*(L21**2) + &
               (2*pi*(g**3)*(o**2)*uj+2*pi*(g**3)*(o**2)*ui)*L11*L13+(pi*(g**3)*(o**2)*wj+pi*(g**3)*(o**2)*wi)*(L11**2))/12

          C120 = C120 + ni*nj*(pi*gsqr*o*(4*g*o*uj*(L32**2)+4*g*o*ui*(L32**2)+8*g*o*vj*L31*L32+8*g*o*vi*L31*L32 + &
               12*uj*vj*L32-12*ui*vi*L32+6*(vj**2)*L31-6*(vi**2)*L31+g*o*uj*(L22**2)+g*o*ui*(L22**2) + &
               2*g*o*vj*L21*L22+2*g*o*vi*L21*L22+g*o*uj*(L12**2)+g*o*ui*(L12**2)+2*g*o*vj*L11*L12 + &
               2*g*o*vi*L11*L12))/12

          C111 = C111 + ni*nj*(pi*gsqr*o*(4*g*o*uj*L32*L33+4*g*o*ui*L32*L33+4*g*o*vj*L31*L33+4*g*o*vi*L31*L33 + &
               6*uj*vj*L33-6*ui*vi*L33+4*g*o*wj*L31*L32+4*g*o*wi*L31*L32+6*uj*wj*L32-6*ui*wi*L32 + &
               6*vj*wj*L31-6*vi*wi*L31+g*o*uj*L22*L23+g*o*ui*L22*L23+g*o*vj*L21*L23+g*o*vi*L21*L23 + &
               g*o*wj*L21*L22+g*o*wi*L21*L22+g*o*uj*L12*L13+g*o*ui*L12*L13+g*o*vj*L11*L13+g*o*vi*L11*L13 + &
               g*o*wj*L11*L12+g*o*wi*L11*L12))/12

          C102 = C102 + ni*nj*(gsqr*o*pi*(4*g*o*vj*(L33**2)+4*g*o*vi*(L33**2)+8*g*o*wj*L32*L33+8*g*o*wi*L32*L33 + &
               12*vj*wj*L33-12*vi*wi*L33+6*(wj**2)*L32-6*(wi**2)*L32+g*o*vj*(L23**2)+g*o*vi*(L23**2)+2*g*o*wj*L22*L23 + &
               2*g*o*wi*L22*L23+g*o*vj*(L13**2)+g*o*vi*(L13**2)+2*g*o*wj*L12*L13+2*g*o*wi*L12*L13))/12

          C030 = C030 + ni*nj*(pi*gsqr*o*(vj+vi)*(4*g*o*(L32**2)+6*vj*L32-6*vi*L32+g*o*(L22**2)+g*o*(L12**2)))/4

          C021 = C021 + ni*nj*(pi*gsqr*o*(8*g*o*vj*L32*L33+8*g*o*vi*L32*L33+6*(vj**2)*L33-6*(vi**2)*L33 + &
               4*g*o*wj*(L32**2)+4*g*o*wi*(L32**2)+12*vj*wj*L32-12*vi*wi*L32+2*g*o*vj*L22*L23+2*g*o*vi*L22*L23 + &
               g*o*wj*(L22**2)+g*o*wi*(L22**2)+2*g*o*vj*L12*L13+2*g*o*vi*L12*L13+g*o*wj*(L12**2)+g*o*wi*(L12**2)))/12

          C012 = C012 + ni*nj*(gsqr*o*pi*(4*g*o*vj*(L33**2)+4*g*o*vi*(L33**2)+8*g*o*wj*L32*L33+8*g*o*wi*L32*L33 + &
               12*vj*wj*L33-12*vi*wi*L33+6*(wj**2)*L32-6*(wi**2)*L32+g*o*vj*(L23**2)+g*o*vi*(L23**2)+2*g*o*wj*L22*L23 + &
               2*g*o*wi*L22*L23+g*o*vj*(L13**2)+g*o*vi*(L13**2)+2*g*o*wj*L12*L13+2*g*o*wi*L12*L13))/12

          C003 = C003 + ni*nj*(pi*gsqr*o*(wj+wi)*(4*g*o*(L33**2)+6*wj*L33-6*wi*L33+g*o*(L23**2)+g*o*(L13**2)))/4
       END DO
    END DO

    d = dp*dp/2 ! Scaling
    Coll(1) = 0
    Coll(2) = 0
    Coll(3) = 0
    Coll(4) = 0
    Coll(5) = d*C200
    Coll(6) = d*C110
    Coll(7) = d*C101
    Coll(8) = d*C020
    Coll(9) = d*C011
    Coll(10) = d*C002
    Coll(11) = d*C300
    Coll(12) = d*C210
    Coll(13) = d*C201
    Coll(14) = d*C120
    Coll(15) = d*C111
    Coll(16) = d*C102
    Coll(17) = d*C030
    Coll(18) = d*C021
    Coll(19) = d*C012
    Coll(20) = d*C003

  END SUBROUTINE COLLISIONS_BOLTZMANN_ONE_SPECIE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Boltzmann collision operator for bidisperse case                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE COLLISIONS_BOLTZMANN_TWO_SPECIES (N1, U1, V1, W1, N2, U2, V2, W2, m1, m2, dp1, dp2, e, Coll)

    USE constant, only: Pi

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: N1, U1, V1, W1, N2, U2, V2, W2
    DOUBLE PRECISION, INTENT(IN) :: m1, m2, dp1, dp2, e
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(QMOMK_NMOM) :: Coll

    DOUBLE PRECISION, DIMENSION(3) :: a1, a, vr, v1r, v2r, g1, g2

    DOUBLE PRECISION :: C100, C010, C001,  C200, C110, C101, C020, C011
    DOUBLE PRECISION :: C002, C300, C210, C201, C120, C111, C102, C030, C021, C012, C003

    DOUBLE PRECISION :: dp, d, g, gsqr, o, mu,  ni, ui, vi, wi, nj, uj, vj, wj
    DOUBLE PRECISION :: cross, denom, L11, L12, L13, L21, L22, L23, L31, L32, L33

    INTEGER :: i, j

    a1 = (/0.0462,0.0971,0.8235/) ! chosen at random (cannot be colinear with vr)
    a = 1/((a1(1))**2 + (a1(2))**2 + (a1(3))**2) * a1;

    C100 = 0.D0
    C010 = 0.D0
    C001 = 0.D0
    C200 = 0.D0
    C110 = 0.D0
    C101 = 0.D0
    C020 = 0.D0
    C011 = 0.D0
    C002 = 0.D0
    C300 = 0.D0
    C210 = 0.D0
    C201 = 0.D0
    C120 = 0.D0
    C111 = 0.D0
    C102 = 0.D0
    C030 = 0.D0
    C021 = 0.D0
    C012 = 0.D0
    C003 = 0.D0

    ! Reduced mass
    mu = m1*m2/(m1+m2)
    o = -(1+e)*mu/m1
    ! Collision diameter
    dp = (dp1 + dp2)/2
    ! Scaling factor
    d = (dp**2)/m2

    DO i = 1, QMOMK_NN
       ui = U1(i)
       vi = V1(i)
       wi = W1(i)
       ni = N1(i)

       v1r(1) = ui
       v1r(2) = vi
       v1r(3) = wi

       DO j = 1, QMOMK_NN
          IF (i == j) THEN
             CYCLE ! vi - vj = 0 always => nothing to add to source term
          END IF

          nj = N2(j)
          uj = U2(j)
          vj = V2(j)
          wj = W2(j)

          v2r(1) = uj
          v2r(2) = vj
          v2r(3) = wj

          vr(1) = v1r(1) - v2r(1)
          vr(2) = v1r(2) - v2r(2)
          vr(3) = v1r(3) - v2r(3)

          g = SQRT(vr(1)*vr(1) + vr(2)*vr(2) + vr(3)*vr(3))

          IF (g == 0) THEN ! Avoiding division by zero in the transformation matrix
             CYCLE
          END IF

          g1(1) = vr(1)/g
          g1(2) = vr(2)/g
          g1(3) = vr(3)/g

          cross = a(1)*g1(1) + a(2)*g1(2) + a(3)*g1(3)

          g2(1) = a(1) - cross*g1(1)
          g2(2) = a(2) - cross*g1(2)
          g2(3) = a(3) - cross*g1(3)

          denom = SQRT( g2(1)*g2(1) + g2(2)*g2(2) + g2(3)*g2(3))

          L11 = g2(1)/denom
          L12 = g2(2)/denom
          L13 = g2(3)/denom
          L21 = ( g1(3)*a(2)-g1(2)*a(3) )/denom
          L22 = ( g1(1)*a(3)-g1(3)*a(1) )/denom
          L23 = ( g1(2)*a(1)-g1(1)*a(2) )/denom
          L31 = g1(1)
          L32 = g1(2)
          L33 = g1(3)

          gsqr = g**2

          C100 = C100 + ni*nj*(pi*gsqr*o*L31)/2
          C010 = C010 + ni*nj*(pi*gsqr*o*L32)/2
          C001 = C001 + ni*nj*(pi*gsqr*o*L33)/2

          C200 = C200 + ni*nj*(4*pi*(g**3)*(o**2)*(L31**2)+12*pi*gsqr*o*ui*L31 &
               + pi*(g**3)*(o**2)*(L21**2)+pi*(g**3)*(o**2)*(L11**2))/12
          C110 = C110 + ni*nj*((4*pi*(g**3)*(o**2)*L31+6*pi*gsqr*o*ui)*L32+6*pi*gsqr*o*vi*L31+pi*(g**3)*(o**2)*L21*L22 + &
               pi*(g**3)*(o**2)*L11*L12)/12

          C101 = C101 + ni*nj*((4*pi*(g**3)*(o**2)*L31+6*pi*gsqr*o*ui)*L33+6*pi*gsqr*o*wi*L31+pi*(g**3)*(o**2)*L21*L23 + &
               pi*(g**3)*(o**2)*L11*L13)/12

          C020 = C020 + ni*nj*(4*pi*(g**3)*(o**2)*(L32**2)+12*pi*gsqr*o*vi*L32 &
               + pi*(g**3)*(o**2)*(L22**2)+pi*(g**3)*(o**2)*(L12**2))/12

          C011 = C011 + ni*nj*((4*pi*(g**3)*(o**2)*L32+6*pi*gsqr*o*vi)*L33+6*pi*gsqr*o*wi*L32+pi*(g**3)*(o**2)*L22*L23 + &
               pi*(g**3)*(o**2)*L12*L13)/12

          C002 = C002 + ni*nj*(4*pi*(g**3)*(o**2)*(L33**2)+12*pi*gsqr*o*wi*L33 &
               + pi*(g**3)*(o**2)*(L23**2)+pi*(g**3)*(o**2)*(L13**2))/12

          C300 = C300 + ni*nj*(2*pi*(g**4)*(o**3)*(L31**3)+8*pi*(g**3)*(o**2)*ui*(L31**2) &
               + (pi*(g**4)*(o**3)*(L21**2)+pi*(g**4)*(o**3)*(L11**2) + 12*pi*gsqr*o*(ui**2))*L31 &
               + 2*pi*(g**3)*(o**2)*ui*(L21**2)+2*pi*(g**3)*(o**2)*ui*(L11**2))/8

          C210 = C210 + ni*nj*((6*pi*(g**4)*(o**3)*(L31**2)+16*pi*(g**3)*(o**2)*ui*L31 &
               + pi*(g**4)*(o**3)*(L21**2)+pi*(g**4)*(o**3)*(L11**2) + 12*pi*gsqr*o*(ui**2))*L32 &
               + 8*pi*(g**3)*(o**2)*vi*(L31**2)+(2*pi*(g**4)*(o**3)*L21*L22+2*pi*(g**4)*(o**3)*L11*L12 &
               + 24*pi*(g**2)*o*ui*vi)*L31+4*pi*(g**3)*(o**2)*ui*L21*L22+2*pi*(g**3)*(o**2)*vi*(L21**2) &
               + 4*pi*(g**3)*(o**2)*ui*L11*L12+2*pi*(g**3)*(o**2)*vi*(L11**2))/(24)

          C201 = C201 + ni*nj*((6*pi*(g**4)*(o**3)*(L31**2)+16*pi*(g**3)*(o**2)*ui*L31 &
               + pi*(g**4)*(o**3)*(L21**2)+pi*(g**4)*(o**3)*(L11**2) + 12*pi*gsqr*o*(ui**2))*L33 &
               + 8*pi*(g**3)*(o**2)*wi*(L31**2)+(2*pi*(g**4)*(o**3)*L21*L23+2*pi*(g**4)*(o**3)*L11*L13 &
               + 24*pi*gsqr*o*ui*wi)*L31+4*pi*(g**3)*(o**2)*ui*L21*L23+2*pi*(g**3)*(o**2)*wi*(L21**2) &
               + 4*pi*(g**3)*(o**2)*ui*L11*L13+2*pi*(g**3)*(o**2)*wi*(L11**2))/(24)

          C120 = C120 + ni*nj*((6*pi*(g**4)*(o**3)*L31+8*pi*(g**3)*(o**2)*ui)*(L32**2)+(16*pi*(g**3)*(o**2)*vi*L31 + &
               2*pi*(g**4)*(o**3)*L21*L22+2*pi*(g**4)*(o**3)*L11*L12+24*pi*gsqr*o*ui*vi)*L32+(pi*(g**4)*(o**3)*(L22**2) + &
               pi*(g**4)*(o**3)*(L12**2)+12*pi*gsqr*o*(vi**2))*L31+2*pi*(g**3)*(o**2)*ui*(L22**2)+4*pi*(g**3)*(o**2)*vi*L21*L22 + &
               2*pi*(g**3)*(o**2)*ui*(L12**2)+4*pi*(g**3)*(o**2)*vi*L11*L12)/(24)

          C111 = C111 + ni*nj*(((6*pi*(g**4)*(o**3)*L31+8*pi*(g**3)*(o**2)*ui)*L32 &
               + 8*pi*(g**3)*(o**2)*vi*L31+pi*(g**4)*(o**3)*L21*L22 + pi*(g**4)*(o**3)*L11*L12 &
               + 12*pi*gsqr*o*ui*vi)*L33+(8*pi*(g**3)*(o**2)*wi*L31+pi*(g**4)*(o**3)*L21*L23 + &
               pi*(g**4)*(o**3)*L11*L13+12*pi*gsqr*o*ui*wi)*L32+(pi*(g**4)*(o**3)*L22*L23+pi*(g**4)*(o**3)*L12*L13+ &
               12*pi*gsqr*o*vi*wi)*L31+(2*pi*(g**3)*(o**2)*ui*L22+2*pi*(g**3)*(o**2)*vi*L21)*L23 + &
               2*pi*(g**3)*(o**2)*wi*L21*L22+(2*pi*(g**3)*(o**2)*ui*L12+2*pi*(g**3)*(o**2)*vi*L11)*L13 + &
               2*pi*(g**3)*(o**2)*wi*L11*L12)/(24)

          C102 = C102 + ni*nj*((6*pi*(g**4)*(o**3)*L31+8*pi*(g**3)*(o**2)*ui)*(L33**2)+(16*pi*(g**3)*(o**2)*wi*L31 + &
               2*pi*(g**4)*(o**3)*L21*L23+2*pi*(g**4)*(o**3)*L11*L13+24*pi*(g**2)*o*ui*wi)*L33 + &
               (pi*(g**4)*(o**3)*(L23**2)+pi*(g**4)*(o**3)*(L13**2)+12*pi*(g**2)*o*(wi**2))*L31+2*pi*(g**3)*(o**2)*ui*(L23**2)+ &
               4*pi*(g**3)*(o**2)*wi*L21*L23+2*pi*(g**3)*(o**2)*ui*(L13**2)+4*pi*(g**3)*(o**2)*wi*L11*L13)/(24)

          C030 = C030 + ni*nj*(2*pi*(g**4)*(o**3)*(L32**3)+8*pi*(g**3)*(o**2)*vi*(L32**2) &
               + (pi*(g**4)*(o**3)*(L22**2)+pi*(g**4)*(o**3)*(L12**2) + 12*pi*(g**2)*o*(vi**2))*L32 &
               + 2*pi*(g**3)*(o**2)*vi*(L22**2)+2*pi*(g**3)*(o**2)*vi*(L12**2))/8

          C021 = C021 + ni*nj*((6*pi*(g**4)*(o**3)*(L32**2)+16*pi*(g**3)*(o**2)*vi*L32 &
               + pi*(g**4)*(o**3)*(L22**2)+pi*(g**4)*(o**3)*(L12**2) + 12*pi*gsqr*o*(vi**2))*L33 &
               + 8*pi*(g**3)*(o**2)*wi*(L32**2)+(2*pi*(g**4)*(o**3)*L22*L23+2*pi*(g**4)*(o**3)*L12*L13 &
               + 24*pi*(g**2)*o*vi*wi)*L32+4*pi*(g**3)*(o**2)*vi*L22*L23+2*pi*(g**3)*(o**2)*wi*(L22**2) &
               + 4*pi*(g**3)*(o**2)*vi*L12*L13+2*pi*(g**3)*(o**2)*wi*(L12**2))/(24)

          C012 = C012 + ni*nj*((6*pi*(g**4)*(o**3)*L32+8*pi*(g**3)*(o**2)*vi)*(L33**2)+(16*pi*(g**3)*(o**2)*wi*L32 + &
               2*pi*(g**4)*(o**3)*L22*L23+2*pi*(g**4)*(o**3)*L12*L13+24*pi*gsqr*o*vi*wi)*L33+(pi*(g**4)*(o**3)*(L23**2) + &
               pi*(g**4)*(o**3)*(L13**2)+12*pi*gsqr*o*(wi**2))*L32+2*pi*(g**3)*(o**2)*vi*(L23**2)+4*pi*(g**3)*(o**2)*wi*L22*L23 + &
               2*pi*(g**3)*(o**2)*vi*(L13**2)+4*pi*(g**3)*(o**2)*wi*L12*L13)/(24)

          C003 = C003 + ni*nj*(2*pi*(g**4)*(o**3)*(L33**3)+8*pi*(g**3)*(o**2)*wi*(L33**2) &
               + (pi*(g**4)*(o**3)*(L23**2)+pi*(g**4)*(o**3)*(L13**2) + 12*pi*gsqr*o*(wi**2))*L33 &
               + 2*pi*(g**3)*(o**2)*wi*(L23**2)+2*pi*(g**3)*(o**2)*wi*(L13**2))/8;

       END DO
    END DO

    Coll(1) = 0;
    Coll(2) = d*C100;
    Coll(3) = d*C010;
    Coll(4) = d*C001;
    Coll(5) = d*C200;
    Coll(6) = d*C110;
    Coll(7) = d*C101;
    Coll(8) = d*C020;
    Coll(9) = d*C011;
    Coll(10) = d*C002;
    Coll(11) = d*C300;
    Coll(12) = d*C210;
    Coll(13) = d*C201;
    Coll(14) = d*C120;
    Coll(15) = d*C111;
    Coll(16) = d*C102;
    Coll(17) = d*C030;
    Coll(18) = d*C021;
    Coll(19) = d*C012;
    Coll(20) = d*C003;

  END SUBROUTINE COLLISIONS_BOLTZMANN_TWO_SPECIES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Calculates the rate of change in the moments due to collisions     C
!   (1 specie)                                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE SOLVE_BOLTZMANN_COLLISIONS_ONE_SPECIE(M, N, U, V, W, dt, &
     e, dp, order)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: M
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: N, U, V, W
    DOUBLE PRECISION, INTENT(IN) ::  dt, e, dp
    INTEGER, INTENT(IN) :: order

    DOUBLE PRECISION, DIMENSION(QMOMK_NMOM) :: Coll, Mtmp
    DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: Ntmp, Utmp, Vtmp, Wtmp
    DOUBLE PRECISION :: h


    IF (order == 0) THEN      ! Euler method - First order
       h = dt
       CALL COLLISIONS_BOLTZMANN_ONE_SPECIE (N, U, V, W, dp, e, Coll)
       Mtmp = M + h*Coll
    ELSE IF (order == 1) THEN ! Runge - Kutta second order
       h = dt ;
       CALL COLLISIONS_BOLTZMANN_ONE_SPECIE (N, U, V, W, dp, e, Coll)
       Mtmp = M + 0.5*h*Coll
       CALL EIGHT_NODE_3D (Mtmp, Ntmp, Utmp, Vtmp, Wtmp) ! Projection step
       CALL COLLISIONS_BOLTZMANN_ONE_SPECIE (Ntmp, Utmp, Vtmp, Wtmp, dp, e, Coll)
       Mtmp = M + h*Coll
    ELSE
       PRINT *,'QMOMK: Wrong order in Boltzmann collisions'
       ERROR_STOP
    END IF
    M= Mtmp
  END SUBROUTINE SOLVE_BOLTZMANN_COLLISIONS_ONE_SPECIE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Calculates the rate of change in the moments due to collisions     C
!   (2 specie)                                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE SOLVE_BOLTZMANN_COLLISIONS_TWO_SPECIES(M1, N1, U1, V1, &
     W1, M2, N2, U2, V2, W2, dt, m_1, m_2, dp1, dp2, e11, e12, order)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: M1
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NMOM) ::  M2
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: N1, U1, V1, W1, N2, U2, V2, W2
    DOUBLE PRECISION, INTENT(IN) ::  dt, dp1, dp2, m_1, m_2, e11, e12
    INTEGER, INTENT(IN) :: order

    DOUBLE PRECISION, DIMENSION(QMOMK_NMOM) :: Coll, Colltmp, M1tmp
    DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: N1tmp, U1tmp, V1tmp, W1tmp
    DOUBLE PRECISION :: h

    IF (order == 0) THEN      ! Euler method - First order
       h = dt ;
       !     Collisions inside the same specie
       CALL COLLISIONS_BOLTZMANN_TWO_SPECIES (N1, U1, V1, W1, N1, U1, V1, W1, m_1, m_1, dp1, dp1, e11, Colltmp)
       Coll = Colltmp
       !     Collisions between the two species
       CALL COLLISIONS_BOLTZMANN_TWO_SPECIES (N1, U1, V1, W1, N2, U2, V2, W2, m_1, m_2, dp1, dp2, e12, Colltmp)
       Coll = Coll + Colltmp
       M1tmp = M1 + h*Coll
    ELSE IF (order == 1) THEN     ! Runge-Kutta second order method
       h = dt ;
       !     Collisions inside the same specie
       CALL COLLISIONS_BOLTZMANN_TWO_SPECIES (N1, U1, V1, W1, N1, U1, V1, W1, m_1, m_1, dp1, dp1, e11, Colltmp)
       Coll = Colltmp
       !     Collisions between the two species
       CALL COLLISIONS_BOLTZMANN_TWO_SPECIES (N1, U1, V1, W1, N2, U2, V2, W2, m_1, m_2, dp1, dp2, e12, Colltmp)
       Coll = Coll + Colltmp
       M1tmp = M1 + 0.5*h*Coll ;
       CALL EIGHT_NODE_3D (M1tmp, N1tmp, U1tmp, V1tmp, W1tmp) ! Projection step
       !     Collisions inside the same specie

       CALL COLLISIONS_BOLTZMANN_TWO_SPECIES (N1tmp, U1tmp, V1tmp, W1tmp, &
            N1tmp, U1tmp, V1tmp, W1tmp, m_1, m_1, dp1, dp1, e11, Colltmp)

       Coll = Colltmp
       !     Collisions between the two species
       CALL COLLISIONS_BOLTZMANN_TWO_SPECIES (N1tmp, U1tmp, V1tmp, W1tmp, N2, U2, V2, W2, m_1, m_2, dp1, dp2, e12, Colltmp)
       Coll = Coll + Colltmp
       M1tmp = M1 + h*Coll;
    ELSE
       PRINT *,'QMOMK: Wrong order in Boltzmann collisions'
       ERROR_STOP
    END IF
    M1 = M1tmp;
  END SUBROUTINE SOLVE_BOLTZMANN_COLLISIONS_TWO_SPECIES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  DOUBLE PRECISION FUNCTION RADIAL_G0(ALPHA)

    USE param1, only: small_number, large_number
    USE constant, only: EP_STAR

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: ALPHA

    DOUBLE PRECISION ALPHA_MAX

    ALPHA_MAX = 1. - EP_STAR

    RADIAL_G0 = 1.
    IF (ALPHA < ALPHA_MAX - SMALL_NUMBER) THEN
      RADIAL_G0 = 1./(1.0-ALPHA) + 3.0*ALPHA/(2.*(1.-ALPHA)**2) + (ALPHA**2)/(2*(1.0-ALPHA)**3)
    ELSE
      RADIAL_G0 = LARGE_NUMBER
    END IF

  END FUNCTION RADIAL_G0

END MODULE qmomk_collision
