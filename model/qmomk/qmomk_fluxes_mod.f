!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!       Module name:  qmomk_fluxes                                     C
!       Purpose: Calculation of the spatial fluxes for QMOM            C
!                                                                      C
!       Author: A. Passalacqua                      Date: 18-Jun-2008  C
!       Reviewer:                                   Date: dd-mmm-yy    C
!                                                                      C
!       Revision Number:                                               C
!       Purpose:                                                       C
!       Author:                                     Date: dd-mmm-yy    C
!       Reviewer:                                   Date: dd-mmm-yy    C
!                                                                      C
!       Literature/Document References:                                C
!                                                                      C
!       Variables referenced:                                          C
!       Variables modified:                                            C
!                                                                      C
!       Local variables:                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE qmomk_fluxes

  USE qmomk_parameters
  USE qmomk_collision

  IMPLICIT NONE

  PRIVATE

  ! Subroutines for kinetic fluxes
  PUBLIC :: KINETIC_FLUX_X_TWENTY_EIGHT_NODES
  PUBLIC :: KINETIC_FLUX_Y_TWENTY_EIGHT_NODES
  PUBLIC :: KINETIC_FLUX_Z_TWENTY_EIGHT_NODES

CONTAINS

  !     Flux calculation in the x direction
  SUBROUTINE KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nl, Ul, Vl, Wl, Nr, Ur, Vr, Wr, f)

    implicit none

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Nl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Ul
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Vl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Wl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Nr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Ur
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Vr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Wr
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: f

    INTEGER :: i

    f = 0.D0

    DO i = 1, QMOMK_NN
       f(1) = f(1) + Nl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*MIN(Ur(i),0.D0)

       f(2) = f(2) + Nl(i)*Ul(i) *MAX(Ul(i),0.D0) &
            + Nr(i)*Ur(i)*MIN(Ur(i),0.D0)

       f(3) = f(3) + Nl(i)*Vl(i) *MAX(Ul(i),0.D0) &
            + Nr(i)*Vr(i)*MIN(Ur(i),0.D0)

       f(4) = f(4) + Nl(i)*Wl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*Wr(i)*MIN(Ur(i),0.D0)

       f(5) = f(5) + Nl(i)*(Ul(i)**2)*MAX(Ul(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*MIN(Ur(i),0.D0)

       f(6) = f(6) + Nl(i)*Ul(i)*Vl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*Ur(i)*Vr(i)*MIN(Ur(i),0.D0)

       f(7) = f(7) + Nl(i)*Ul(i)*Wl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*Ur(i)*Wr(i)*MIN(Ur(i),0.D0)

       f(8) = f(8) + Nl(i)*(Vl(i)**2)*MAX(Ul(i),0.D0) &
            + Nr(i)*(Vr(i)**2)*MIN(Ur(i),0.D0)

       f(9) = f(9) + Nl(i)*Vl(i)*Wl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*Vr(i)*Wr(i)*MIN(Ur(i),0.D0)

       f(10)= f(10)+ Nl(i)*(Wl(i)**2) *MAX(Ul(i),0.D0) &
            + Nr(i)*(Wr(i)**2)*MIN(Ur(i),0.D0)

       f(11)= f(11)+ Nl(i)*(Ul(i)**3)*MAX(Ul(i),0.D0) &
            + Nr(i)*(Ur(i)**3)*MIN(Ur(i),0.D0)

       f(12)= f(12)+ Nl(i)*(Ul(i)**2)*Vl(i) *MAX(Ul(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*Vr(i) *MIN(Ur(i),0.D0)

       f(13)= f(13)+ Nl(i)*(Ul(i)**2)*Wl(i) *MAX(Ul(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*Wr(i)*MIN(Ur(i),0.D0)

       f(14)= f(14)+ Nl(i)*Ul(i)*(Vl(i)**2) *MAX(Ul(i),0.D0) &
            + Nr(i)*Ur(i)*(Vr(i)**2)*MIN(Ur(i),0.D0)

       f(15)= f(15)+ Nl(i)*Ul(i)*Vl(i)*Wl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*Ur(i)*Vr(i)*Wr(i)*MIN(Ur(i),0.D0)

       f(16)= f(14)+ Nl(i)*Ul(i)*(Wl(i)**2)*MAX(Ul(i),0.D0) &
            + Nr(i)*Ur(i)*(Wr(i)**2)*MIN(Ur(i),0.D0)

       f(17)= f(17)+ Nl(i)*(Vl(i)**3)*MAX(Ul(i),0.D0) &
            + Nr(i)*(Vr(i)**3)*MIN(Ur(i),0.D0)

       f(18)= f(18)+ Nl(i)*(Vl(i)**2)*Wl(i)*MAX(Ul(i),0.D0) &
            + Nr(i)*(Vr(i)**2)*Wr(i)*MIN(Ur(i),0.D0)

       f(19)= f(19)+ Nl(i)*Vl(i)*(Wl(i)**2)*MAX(Ul(i),0.D0) &
            + Nr(i)*Vr(i)*(Wr(i)**2)*MIN(Ur(i),0.D0)

       f(20)= f(20)+ Nl(i)*(Wl(i)**3)*MAX(Ul(i),0.D0) &
            + Nr(i)*(Wr(i)**3)*MIN(Ur(i),0.D0)
    END DO
  END SUBROUTINE KINETIC_FLUX_X_TWENTY_EIGHT_NODES

  !     Flux calculation in the y direction
  SUBROUTINE KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nl, Ul, Vl, Wl, Nr, Ur, Vr, Wr, f)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Nl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Ul
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Vl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Wl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Nr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Ur
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Vr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Wr
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: f

    INTEGER :: i

    f = 0.D0

    DO i = 1, QMOMK_NN
       f(1) = f(1) + Nl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*MIN(Vr(i),0.D0)

       f(2) = f(2) + Nl(i)*Ul(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Ur(i)*MIN(Vr(i),0.D0)

       f(3) = f(3) + Nl(i)*Vl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Vr(i)*MIN(Vr(i),0.D0)

       f(4) = f(4) + Nl(i)*Wl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Wr(i)*MIN(Vr(i),0.D0)

       f(5) = f(5) + Nl(i)*(Ul(i)**2)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*MIN(Vr(i),0.D0)

       f(6) = f(6) + Nl(i)*Ul(i)*Vl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Ur(i)*Vr(i)*MIN(Vr(i),0.D0)

       f(7) = f(7) + Nl(i)*Ul(i)*Wl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Ur(i)*Wr(i)*MIN(Vr(i),0.D0)

       f(8) = f(8) + Nl(i)*(Vl(i)**2)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Vr(i)**2)*MIN(Vr(i),0.D0)

       f(9) = f(9) + Nl(i)*Vl(i)*Wl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Vr(i)*Wr(i)*MIN(Vr(i),0.D0)

       f(10)= f(10)+ Nl(i)*(Wl(i)**2)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Wr(i)**2)*MIN(Vr(i),0.D0)

       f(11)= f(11)+ Nl(i)*(Ul(i)**3)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Ur(i)**3)*MIN(Vr(i),0.D0)

       f(12)= f(12)+ Nl(i)*(Ul(i)**2)*Vl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*Vr(i)*MIN(Vr(i),0.D0)

       f(13)= f(13)+ Nl(i)*(Ul(i)**2)*Wl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*Wr(i)*MIN(Vr(i),0.D0)

       f(14)= f(14)+ Nl(i)*Ul(i)*(Vl(i)**2)*MAX(Vl(i),0.D0) &
            + Nr(i)*Ur(i)*(Vr(i)**2)*MIN(Vr(i),0.D0)

       f(15)= f(15)+ Nl(i)*Ul(i)*Vl(i)*Wl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*Ur(i)*Vr(i)*Wr(i)*MIN(Vr(i),0.D0)

       f(16)= f(14)+ Nl(i)*Ul(i)*(Wl(i)**2)*MAX(Vl(i),0.D0) &
            + Nr(i)*Ur(i)*(Wr(i)**2)*MIN(Vr(i),0.D0)

       f(17)= f(17)+ Nl(i)*(Vl(i)**3)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Vr(i)**3)*MIN(Vr(i),0.D0)

       f(18)= f(18)+ Nl(i)*(Vl(i)**2)*Wl(i)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Vr(i)**2)*Wr(i)*MIN(Vr(i),0.D0)

       f(19)= f(19)+ Nl(i)*Vl(i)*(Wl(i)**2)*MAX(Vl(i),0.D0) &
            + Nr(i)*Vr(i)*(Wr(i)**2)*MIN(Vr(i),0.D0)

       f(20)= f(20)+ Nl(i)*(Wl(i)**3)*MAX(Vl(i),0.D0) &
            + Nr(i)*(Wr(i)**3)*MIN(Vr(i),0.D0)

    END DO
  END SUBROUTINE KINETIC_FLUX_Y_TWENTY_EIGHT_NODES

  !     FLUX CALCULATION IN the z direction
  SUBROUTINE KINETIC_FLUX_Z_TWENTY_EIGHT_NODES (Nl, Ul, Vl, Wl, Nr, Ur, Vr, Wr, f)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Nl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Ul
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Vl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Wl
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Nr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Ur
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Vr
    DOUBLE PRECISION, INTENT(IN), DIMENSION(QMOMK_NN) :: Wr
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(QMOMK_NMOM) :: f

    INTEGER :: i

    f = 0.D0

    DO i = 1, QMOMK_NN
       f(1) = f(1) + Nl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*MIN(Wr(i),0.D0)

       f(2) = f(2) + Nl(i)*Ul(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Ur(i)*MIN(Wr(i),0.D0)

       f(3) = f(3) + Nl(i)*Vl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Vr(i)*MIN(Wr(i),0.D0)

       f(4) = f(4) + Nl(i)*Wl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Wr(i)*MIN(Wr(i),0.D0)

       f(5) = f(5) + Nl(i)*(Ul(i)**2)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*MIN(Wr(i),0.D0)

       f(6) = f(6) + Nl(i)*Ul(i)*Vl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Ur(i)*Vr(i)*MIN(Wr(i),0.D0)

       f(7) = f(7) + Nl(i)*Ul(i)*Wl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Ur(i)*Wr(i)*MIN(Wr(i),0.D0)

       f(8) = f(8) + Nl(i)*(Vl(i)**2)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Vr(i)**2)*MIN(Wr(i),0.D0)

       f(9) = f(9) + Nl(i)*Vl(i)*Wl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Vr(i)*Wr(i)*MIN(Wr(i),0.D0)

       f(10)= f(10)+ Nl(i)*(Wl(i)**2)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Wr(i)**2)*MIN(Wr(i),0.D0)

       f(11)= f(11)+ Nl(i)*(Ul(i)**3)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Ur(i)**3)*MIN(Wr(i),0.D0)

       f(12)= f(12)+ Nl(i)*(Ul(i)**2)*Vl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*Vr(i)*MIN(Wr(i),0.D0)

       f(13)= f(13)+ Nl(i)*(Ul(i)**2)*Wl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Ur(i)**2)*Wr(i)*MIN(Wr(i),0.D0)

       f(14)= f(14)+ Nl(i)*Ul(i)*(Vl(i)**2)*MAX(Wl(i),0.D0) &
            + Nr(i)*Ur(i)*(Vr(i)**2)*MIN(Wr(i),0.D0)

       f(15)= f(15)+ Nl(i)*Ul(i)*Vl(i)*Wl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*Ur(i)*Vr(i)*Wr(i)*MIN(Wr(i),0.D0)

       f(16)= f(14)+ Nl(i)*Ul(i)*(Wl(i)**2)*MAX(Wl(i),0.D0) &
            + Nr(i)*Ur(i)*(Wr(i)**2)*MIN(Wr(i),0.D0)

       f(17)= f(17)+ Nl(i)*(Vl(i)**3)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Vr(i)**3)*MIN(Wr(i),0.D0)

       f(18)= f(18)+ Nl(i)*(Vl(i)**2)*Wl(i)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Vr(i)**2)*Wr(i)*MIN(Wr(i),0.D0)

       f(19)= f(19)+ Nl(i)*Vl(i)*(Wl(i)**2)*MAX(Wl(i),0.D0) &
            + Nr(i)*Vr(i)*(Wr(i)**2)*MIN(Wr(i),0.D0)

       f(20)= f(20)+ Nl(i)*(Wl(i)**3)*MAX(Wl(i),0.D0) &
            + Nr(i)*(Wr(i)**3)*MIN(Wr(i),0.D0)

    END DO
  END SUBROUTINE KINETIC_FLUX_Z_TWENTY_EIGHT_NODES

END MODULE qmomk_fluxes
