!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: scales_mod.f                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      MODULE scales

! reference pressure
      DOUBLE PRECISION :: P_ref

! pressure scale
      DOUBLE PRECISION :: P_scale

      CONTAINS

      DOUBLE PRECISION FUNCTION SCALE_PRESSURE(XXX)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXX
      SCALE_PRESSURE   = (XXX - P_ref) / P_scale
      END FUNCTION SCALE_PRESSURE

      DOUBLE PRECISION FUNCTION UNSCALE_PRESSURE(XXX)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXX
      UNSCALE_PRESSURE = (XXX * P_scale + P_ref)
      END FUNCTION UNSCALE_PRESSURE

      END MODULE scales
