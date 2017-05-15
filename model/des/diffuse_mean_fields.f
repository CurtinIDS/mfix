!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose: Given the field variable PHI (e.g., volume fraction),      !
!  diffuse it across the Eulerian grid such that the Full Width at     !
!  Half Maximum (FWHM) equals the user specified diffusion width.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIFFUSE_MEAN_FIELD(PHI, VNAME)

! Global Variables:
!---------------------------------------------------------------------//
! Max bound for array sizes.
      use geometry, only: IJKMAX2
! Coefficient matrix and force vector.
      use ambm, only: A_M, B_M
! Method to solve linear system and max iterations
      use leqsol, only: LEQ_METHOD, LEQ_IT
! Preconditioner, sweep method, convergence tolerance
      use leqsol, only: LEQ_PC, LEQ_SWEEP, LEQ_TOL
! Size of fluid variable arrays.
      use param, only: DIMENSION_3

! Module procedures
!---------------------------------------------------------------------//
! Routines to manage messages to user.
      use error_manager
      use machine, only: wall_time

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Variable to diffuse
      DOUBLE PRECISION, INTENT(INOUT) :: PHI(DIMENSION_3)
! Name of variable to diffuse
      CHARACTER(LEN=*), INTENT(IN) :: VNAME

! Local Variables:
!---------------------------------------------------------------------//
! Integer error flag
      INTEGER :: IER
! Linear equation solver method and iterations
      INTEGER :: LEQM, LEQI
! Start, stop and step size of diffusion time march
      DOUBLE PRECISION :: DIF_TIME, DIF_STOP, DIF_DT
! wall time at start
      DOUBLE PRECISION :: WALL_START
! Local flag to print debug messages
      LOGICAL, PARAMETER :: setDBG = .TRUE.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIF
!......................................................................!

      IF(setDBG) THEN
         WRITE(ERR_MSG, "(/3x,'Diffusing Variable: ',A)") VNAME
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         WALL_START = WALL_TIME()
      ENDIF

! Populate the diffusion coefficients
      CALL CALC_DIF_DES(DIF, setDBG, IER)

      DIF_STOP = 1.0d0
      DIF_TIME = 0.0d0
      DIF_DT = DIF_STOP/5.0

! Integrate the diffusion equation (time, space)
      DO WHILE(DIF_TIME < DIF_STOP)
! Initialize the coefficient matrix and force vector
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)
! Calculate the coefficients
         CALL DIF_PHI_DES(0, DIF, A_M, B_M)
! Apply zero-flux BC at all walls
         CALL DIF_PHI_BC_DES(PHI, 0, A_M, B_M)
! Collect the center coefficient and force vector
         CALL DIF_PHI_SOURCE_DES(PHI, 0, A_M, B_M, DIF_DT)
! Set the local method and iterations.
         CALL ADJUST_LEQ(0.0d0, LEQ_IT(10), LEQ_METHOD(10), LEQI, LEQM)
! Solve the linear system.
         CALL SOLVE_LIN_EQ (VNAME, 10, PHI, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(10), LEQ_TOL(10), LEQ_PC(10), IER)
! Advance time.
         DIF_TIME = DIF_TIME + DIF_DT
      ENDDO

! Debugging information
      IF(setDBG) THEN
         WRITE(ERR_MSG, 9001) WALL_TIME() - WALL_START
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 9001 FORMAT(5x,'Wall Time: ',g11.4)

      RETURN
      END SUBROUTINE DIFFUSE_MEAN_FIELD

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DIF_DES                                            !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose: Calculate the diffusion coefficient for diffusing mean     !
!  fields. Presently the diffusion coefficient is constant.            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DIF_DES(DIF, lDBG, IER)

      use param
      use param1

      use compar, only: IJKStart3, IJKEnd3
      use particle_filter, only: DES_DIFFUSE_WIDTH
      use functions, only: FLUID_AT

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: DIF(DIMENSION_3)
      LOGICAL, INTENT(IN) :: lDBG
      INTEGER, INTENT(INOUT) :: IER

! Fluid Cell indices
      INTEGER :: IJK

      DOUBLE PRECISION :: lDIF

      IER = 0

! The diffusion coefficient is set so that over one second, the
! quantity diffuses such that the Full Width at Half Maximum (FWHM)
! equals what the user specified as the "filter wideth."
      lDIF = ((0.5*DES_DIFFUSE_WIDTH)**2) / &
         (2.0*sqrt(2.0*log(2.0)))

! Store the diffusion coefficient in all fluid cells.
      DO IJK = IJKStart3, IJKEnd3
         DIF(IJK) = ZERO
         IF(FLUID_AT(IJK)) DIF(IJK) = lDIF
      ENDDO

! Information included for debugging.
      IF(lDBG) THEN
         WRITE(ERR_MSG, 9100) iVal(lDIF)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 9100 FORMAT(/3x,'Diffusion Coefficient: ',A)

      RETURN
      END SUBROUTINE CALC_DIF_DES
