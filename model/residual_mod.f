! -*- f90 -*-
MODULE residual

      Use param, only: dim_n, dim_m

!     residual.inc

      INTEGER, PARAMETER :: MAX_RESID_INDEX = 8    !for printing; don't change this

      INTEGER, PARAMETER :: RESID_p  = 1     !pressure
      INTEGER, PARAMETER :: RESID_ro = 2     !density, volume fraction
      INTEGER, PARAMETER :: RESID_u  = 3     !u-velocity
      INTEGER, PARAMETER :: RESID_v  = 4     !v-velocity
      INTEGER, PARAMETER :: RESID_w  = 5     !w-velocity
      INTEGER, PARAMETER :: RESID_t  = 6     !temperature
      INTEGER, PARAMETER :: RESID_th = 7     !granular temperature
      INTEGER, PARAMETER :: RESID_sc = 8     !user-defined scalar
      INTEGER, PARAMETER :: NRESID   = 8 + DIM_N
      INTEGER, PARAMETER :: RESID_ke = 9     !k-epsilon equations
      INTEGER, PARAMETER :: RESID_x  = 10    !mass fraction (keep this the last)
      INTEGER, PARAMETER :: NPREFIX  = 10
!
!    Group Resisuals by equation
      INTEGER, PARAMETER :: HYDRO_GRP   = 1     !hydrodynamics
      INTEGER, PARAMETER :: THETA_GRP   = 2     !Granular Energy
      INTEGER, PARAMETER :: ENERGY_GRP  = 3     !Energy
      INTEGER, PARAMETER :: SPECIES_GRP = 4     !Species
      INTEGER, PARAMETER :: SCALAR_GRP  = 5     !Scalars
      INTEGER, PARAMETER :: KE_GRP      = 6     !K-Epsilon

!                      prefix of Residuals string
      CHARACTER, PARAMETER, DIMENSION(NPREFIX) :: RESID_PREFIX = &
        (/ 'P', 'R', 'U', 'V', 'W', 'T', 'G', 'S', 'K', 'X' /)

!
!                      Average residual
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RESID
!
!                      Maximum residual
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  MAX_RESID
!
!                      sum of residuals every 5 iterations
      DOUBLE PRECISION SUM5_RESID
!
!                      IJK location of maximum residual
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IJK_RESID

!                      Residual Numerator
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: NUM_RESID
!
!                      Residual Denominator
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DEN_RESID
!
!                      Residual Packing for Global Operations
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RESID_PACK
!
!                      Residual sum within a group of equations
      LOGICAL          :: GROUP_RESID
      DOUBLE PRECISION :: RESID_GRP(6)
!
!                      Residuals to be printed out
      CHARACTER(LEN=4)      RESID_STRING(MAX_RESID_INDEX)
      CHARACTER(LEN=8)      RESID_GRP_STRING(6)
!
!                      Indices of residuals to be printed out
      INTEGER          RESID_INDEX(MAX_RESID_INDEX, 2)
!

!                        fluid and solids accumulation, for checking the over-all fluid mass balance
      DOUBLE PRECISION accum_resid_g, accum_resid_s(DIM_M)

   CONTAINS

      FUNCTION GET_RESID_STRING(INDEX)
         IMPLICIT NONE
         CHARACTER(LEN=4) :: GET_RESID_STRING
         INTEGER, INTENT(IN) :: INDEX

         GET_RESID_STRING = RESID_STRING(INDEX)

      END FUNCTION GET_RESID_STRING

      FUNCTION GET_RESID_GRP_STRING(INDEX)
         IMPLICIT NONE
         CHARACTER(LEN=8) :: GET_RESID_GRP_STRING
         INTEGER, INTENT(IN) :: INDEX

         GET_RESID_GRP_STRING = RESID_GRP_STRING(INDEX)

      END FUNCTION GET_RESID_GRP_STRING

      FUNCTION GET_RESID(INDEX)
         IMPLICIT NONE
         DOUBLE PRECISION :: GET_RESID
         INTEGER, INTENT(IN) :: INDEX
         INTEGER :: RI, RI2

         IF (INDEX > SIZE(RESID_INDEX,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR INDEX ",INDEX
            GET_RESID = 0.0
            RETURN
         ENDIF
         RI = RESID_INDEX(INDEX,1)
         RI2 = RESID_INDEX(INDEX,2)
         IF (RI > SIZE(RESID,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR RESID_INDEX 1 ",RI
            GET_RESID = 0.0
            RETURN
         ENDIF
         IF (RI2 > SIZE(RESID,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR RESID_INDEX 2 ",RI2
            GET_RESID = 0.0
            RETURN
         ENDIF
         GET_RESID = RESID(RI,RI2)

      END FUNCTION GET_RESID

      FUNCTION GET_RESID_GRP(INDEX)
         IMPLICIT NONE
         DOUBLE PRECISION :: GET_RESID_GRP
         INTEGER, INTENT(IN) :: INDEX

          GET_RESID_GRP = RESID_GRP(INDEX)

      END FUNCTION GET_RESID_GRP

      END MODULE residual
