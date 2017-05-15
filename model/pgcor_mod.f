      MODULE pgcor


      Use param
      Use param1


!
!  Variables for gas pressure correction equation
!
!
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  d_e
!
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  d_n
!
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  d_t
!
!                      gas pressure correction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Pp_g
!
!                      Indicates the phase used for gas pressure correction
      INTEGER, DIMENSION(:), ALLOCATABLE ::           PHASE_4_P_g
!
!                      Indicates whether different phases were used for gas
!                      pressure correction
      LOGICAL          SWITCH_4_P_g
!


!!!HPF$ align d_e(:, *) with TT(:)
!!!HPF$ align d_n(:, *) with TT(:)
!!!HPF$ align d_t(:, *) with TT(:)
!!!HPF$ align Pp_g(:) with TT(:)
!!!HPF$ align PHASE_4_P_g(:) with TT(:)

      END MODULE pgcor
