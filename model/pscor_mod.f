      MODULE pscor


      Use param
      Use param1


      INTEGER          P_star_bdry
      PARAMETER        (P_star_bdry = 100)

!
!  Variables for solids pressure correction equation
!
!
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  e_e
!
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  e_n
!
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  e_t
!
!                      dPodEP_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_cp
!
!                      Solids volume fraction correction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  EPp
!
!                      Indicates the phase used for solids pressure correction
      INTEGER, DIMENSION(:), ALLOCATABLE ::           PHASE_4_P_s
!
!                      Indicates whether different phases were used for solids
!                      pressure correction
      LOGICAL          SWITCH_4_P_s
!
!                      Indicates whether solids pressure correction is needed.
      LOGICAL          DO_P_s
!
!                      Index for close-packed solids phase
      INTEGER          Mcp
!


!!!HPF$ align e_e(:) with TT(:)
!!!HPF$ align e_n(:) with TT(:)
!!!HPF$ align e_t(:) with TT(:)
!!!HPF$ align K_cp(:) with TT(:)
!!!HPF$ align EPp(:) with TT(:)
!!!HPF$ align PHASE_4_P_s(:) with TT(:)


      END MODULE pscor
