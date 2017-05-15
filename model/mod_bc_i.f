!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MOD_BC_I                                                C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!                                                                      C
!  Purpose: modify the "I" values for the b.c. plane                   C
!     This is a yz plane                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE MOD_BC_I(BCV)

      use bc, only: BC_I_W, BC_I_E
      use bc, only: BC_J_S, BC_J_N
      use bc, only: BC_K_B, BC_K_T
      use bc, only: BC_PLANE

      USE geometry, only: ICBC_FLAG

      USE compar
      USE mpi_utility

      use error_manager
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! boundary condition index
      INTEGER, INTENT(IN) :: BCV

! i cell indices defining location of yz plane
      INTEGER :: I_w, I_e
! south/bottom j,k cell indices of yz plane
      INTEGER :: J_s, K_b

!-----------------------------------------------
! Local variables
!-----------------------------------------------
!     'IJK' indices
      INTEGER :: OWNER
      INTEGER :: J, K
      INTEGER :: IJK , IPJK

      INTEGER :: IER
      LOGICAL :: ERROR
      INTEGER :: I_FLUID, IJK_FLUID
      INTEGER :: I_WALL,  IJK_WALL

      CALL INIT_ERR_MSG("MOD_BC_I")

      I_W = BC_I_W(BCV)
      I_E = BC_I_E(BCV)

      J_S = BC_J_S(BCV)
      K_B = BC_K_B(BCV)

! Establish the OWNER of the BC
      OWNER = merge(myPE, 0, IS_ON_myPE_owns(I_W, J_S, K_B))
      CALL GLOBAL_ALL_SUM(OWNER)

      IF(myPE == OWNER) THEN

         IJK  = FUNIJK(I_W,   J_S, K_B)
         IPJK = FUNIJK(I_W+1, J_S, K_B)

! Flow on west boundary (fluid cell on east).
         IF(WALL_ICBC_FLAG(IJK) .AND. ICBC_FLAG(IPJK)(1:1)=='.') THEN
            I_W = I_W
            I_E = I_E
            BC_PLANE(BCV) = 'E'

! Flow on east boundary (fluid cell on west).
         ELSEIF(WALL_ICBC_FLAG(IPJK) .AND. ICBC_FLAG(IJK)(1:1)=='.') THEN
            I_W = I_W + 1
            I_E = I_E + 1
            BC_PLANE(BCV) = 'W'

! Set the plane of a value we know to be wrong so we can detect the error.
         ELSE
            BC_PLANE(BCV) = '.'
         ENDIF
      ENDIF

! The owner distributes the new Iw/Ie coordinates to the other ranks.
      CALL BCAST(I_W,   OWNER)
      CALL BCAST(I_E,   OWNER)
      CALL BCAST(BC_PLANE(BCV), OWNER)


! If there is an error, send IJK/IPJK to all ranks. Report and exit.
      IF(BC_PLANE(BCV) == '.') THEN
         CALL BCAST(IPJK,OWNER)
         CALL BCAST(IJK, OWNER)

         WRITE(ERR_MSG, 1100) BCV, I_W, I_E, J_S, K_B,                 &
            IJK, ICBC_FLAG(IJK),  IPJK, ICBC_FLAG(IPJK)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'I West   =  ',I6,' I East   = ',I6,/&
         3x,'J South  =  ',I6,' K Bottom = ',I6,2/' The following ',   &
         'should conttain a wall cell and fluid cell:',/3x,'IJK  ',I9, &
         ' :: ',A3,/3x,'IPJK ',I9,' :: ',A3,2/' Maybe no IC was ',     &
         'specified for the fluid cell.')

! Store the new values in the global data array.
      BC_I_W(BCV) = I_W
      BC_I_E(BCV) = I_E

! Set up the I-indices for checking the entire BC region.
      I_WALL = BC_I_W(BCV)
      I_FLUID = merge(I_WALL-1, I_WALL+1, BC_PLANE(BCV)=='W')


! First pass through all of the BC region and verify that you have
! inflow/outflow specified against a wall. Flag any errors.
      ERROR = .FALSE.
      DO K = BC_K_B(BCV), BC_K_T(BCV)
      DO J = BC_J_S(BCV), BC_J_N(BCV)

         IF(.NOT.IS_ON_myPE_plus2layers(I_FLUID,J,K)) CYCLE
         IF(.NOT.IS_ON_myPE_plus2layers(I_WALL, J,K)) CYCLE
         IF(DEAD_CELL_AT(I_FLUID,J,K)) CYCLE
         IF(DEAD_CELL_AT(I_WALL, J,K)) CYCLE

         IJK_WALL = FUNIJK(I_WALL,J,K)
         IJK_FLUID = FUNIJK(I_FLUID,J,K)

! Verify that the the fluid and wall cells match the ICBC_FLAG.
         IF(.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND.                       &
            ICBC_FLAG(IJK_FLUID)(1:1) == '.')) ERROR = .TRUE.

      ENDDO
      ENDDO

! Sync up the error flag across all processes.
      CALL GLOBAL_ALL_OR(ERROR)

! If an error is detected, have each rank open a log file and write
! it's own message. Otherwise, we need to send all the data back to
! PE_IO and that's too much work!
      IF(ERROR) THEN

         CALL OPEN_PE_LOG(IER)

         WRITE(ERR_MSG, 1200) BCV
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)

            IF(.NOT.IS_ON_myPE_plus2layers(I_FLUID,J,K)) CYCLE
            IF(.NOT.IS_ON_myPE_plus2layers(I_WALL, J,K)) CYCLE
            IF(DEAD_CELL_AT(I_FLUID,J,K)) CYCLE
            IF(DEAD_CELL_AT(I_WALL, J,K)) CYCLE

            IJK_WALL = FUNIJK(I_WALL,J,K)
            IJK_FLUID = FUNIJK(I_FLUID,J,K)

            IF(.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND.                    &
               ICBC_FLAG(IJK_FLUID)(1:1) == '.')) THEN

               WRITE(ERR_MSG, 1201) &
                  I_WALL,  J, K, IJK_WALL, ICBC_FLAG(IJK_WALL),        &
                  I_FLUID, J, K, IJK_FLUID, ICBC_FLAG(IJK_FLUID)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF
         ENDDO
         ENDDO

         WRITE(ERR_MSG,"('Please correct the mfix.dat file.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ENDIF


 1200 FORMAT('Error 1200: Illegal geometry for boundary condition:',I3)

 1201 FORMAT(' ',/14X,'I',7X,'J',7X,'K',7X,'IJK',4x,'FLAG',/3x,        &
         'WALL ',3(2x,I6),2x,I9,3x,A,/3x,'FLUID',3(2x,I6),2x,I9,3x,A)


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MOD_BC_I

