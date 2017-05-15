!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: tolerance.inc                                          C
!  Purpose: Specify all tolerance parameters                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUL-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE toleranc

!      PUBLIC :: COMPARE

      Use param
      Use param1

!                      Minimum value of solids volume fraction tracked
      DOUBLE PRECISION, PARAMETER          ::  ZERO_EP_s = 1.0D-10
!
!                      Small value for species mass fraction for disregarding
!                      residual calculation
      DOUBLE PRECISION, PARAMETER          ::  ZERO_X_gs = 1.0D-7

!                      Dilute flow threshold.  When the volume fraction of a
!                      certain phase in a cell is smaller than this value the
!                      momentum equation for that phase is not solved in the cell.
      DOUBLE PRECISION, PARAMETER          ::  DIL_EP_s = 1.0D-8
!
!                      Tolerance used for comparing two numbers for equality
!                      in function compare(a, b)
      DOUBLE PRECISION, PARAMETER          ::  TOL_COM = 1.0D-4
!
!                      Upper bound for temperatures
      DOUBLE PRECISION, PARAMETER          ::  TMAX = 4000.D0
!
!                      Lower bound for temperatures
      DOUBLE PRECISION, PARAMETER          ::  TMIN = 250.D0
!
!                      Reciprocal of a maximum molecular weight
      DOUBLE PRECISION, PARAMETER          ::  oMW_MAX = (ONE/500.D0)
!
!                      Maximum value of velocities set to avoid divergence problems.
      DOUBLE PRECISION MAX_INLET_VEL
!
!                      User definable factor used to scale MAX_INLET_VEL.  Default value is 1.
      DOUBLE PRECISION MAX_INLET_VEL_FAC
!
!                      Maximum allowed velocity of gas or solids in case no inlet velocities
!                      (or zero velocities) are defined at inlet (see function check_vel_bound)
      DOUBLE PRECISION, PARAMETER          ::  MAX_ALLOWED_VEL = 500.0D+2
!
!  The following quantities can be specified through the input data file, with
!  namelist inputs of the same name.
!
!                      Tolerance in residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID
!
!                      Tolerance in energy eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_T
!
!                      Tolerance in species eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_X
!
!                      Tolerance in scalr eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_Scalar
!
!                      Tolerance in K & Epsilon eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_K_Epsilon
!
!                      Tolerance in Granular Temperature eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_Th
!
!                      Minimum residual for declaring divergence
      DOUBLE PRECISION TOL_DIVERGE
!
!                      Factor for normalizing the residual of gas cont. eq.
      DOUBLE PRECISION NORM_g
!
!                      Factor for normalizing the residual of solids cont. eq.
      DOUBLE PRECISION NORM_s
!
!                      Detect negative Rho_g in physical_prop to reduce DT in iterate
      LOGICAL  ::      Neg_RHO_G = .FALSE.

      CONTAINS

      LOGICAL FUNCTION COMPARE (V1, V2)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Values to be compared
      DOUBLE PRECISION V1, V2
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
      IF (ABS(V1) <= SMALL_NUMBER) THEN
         IF (ABS(V2) <= SMALL_NUMBER) THEN
            COMPARE = .TRUE.
         ELSE
            COMPARE = .FALSE.
         ENDIF
      ELSE
         IF (ABS(V2/V1 - ONE) <= TOL_COM) THEN
            COMPARE = .TRUE.
         ELSE
            COMPARE = .FALSE.
         ENDIF
      ENDIF
      RETURN
      END FUNCTION COMPARE

      LOGICAL FUNCTION IS_SMALL (V, TOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Tolerance value for small
      DOUBLE PRECISION TOL
!
!                      Field vriable array
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: V
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJK
!-----------------------------------------------

      IS_SMALL = .FALSE.
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IF (ABS(V(IJK)) > TOL) RETURN
         ENDIF
      END DO
      IS_SMALL = .TRUE.
!
      RETURN
      END FUNCTION IS_SMALL

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3


      END MODULE toleranc
