!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MOD_BC_J(BC, I_w, J_s, J_n, K_b, PLANE)                !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: modify the "J" values for the b.c. plane                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MOD_BC_J(BCV)

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

! boundary condition index
      INTEGER, INTENT(in) :: BCV

! Calculated cell indices in I,J,K directions
      INTEGER :: J_s, J_n
      INTEGER :: I_w, K_b

      INTEGER :: OWNER
      INTEGER :: I, K
      INTEGER :: IJK , IJPK

      INTEGER :: IER
      LOGICAL :: ERROR
      INTEGER :: J_FLUID, IJK_FLUID
      INTEGER :: J_WALL,  IJK_WALL


!-----------------------------------------------

      CALL INIT_ERR_MSG("MOD_BC_J")

      J_S = BC_J_S(BCV)
      J_N = BC_J_N(BCV)

      I_W = BC_I_W(BCV)
      K_B = BC_K_B(BCV)

! Establish the OWNER of the BC
      OWNER = merge(myPE, 0, IS_ON_myPE_owns(I_W,J_S,K_B))
      CALL GLOBAL_ALL_SUM(OWNER)

      IF(myPE == OWNER)THEN

         IJK  = FUNIJK(I_W, J_S,   K_B)
         IJPK = FUNIJK(I_W, J_S+1, K_B)

         IF (WALL_ICBC_FLAG(IJK) .AND. ICBC_FLAG(IJPK)(1:1)=='.') THEN
            J_S = J_S
            J_N = J_N
            BC_PLANE(BCV) = 'N'

         ELSE IF (WALL_ICBC_FLAG(IJPK) .AND. ICBC_FLAG(IJK)(1:1)=='.') THEN
            J_S = J_S + 1
            J_N = J_N + 1
            BC_PLANE(BCV) = 'S'

         ELSE
            BC_PLANE(BCV) = '.'
         ENDIF
      ENDIF

      CALL BCAST(J_S,OWNER)
      CALL BCAST(J_N,OWNER)
      CALL BCAST(BC_PLANE(BCV),OWNER)

! If there is an error, send IJK/IPJK to all ranks. Report and exit.
      IF(BC_PLANE(BCV) == '.') THEN
         CALL BCAST(IJPK,OWNER)
         CALL BCAST(IJK, OWNER)

         WRITE(ERR_MSG, 1100) BCV, J_S, J_N, I_W, K_B,                 &
            IJK, ICBC_FLAG(IJK),  IJPK, ICBC_FLAG(IJPK)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'J South  =  ',I6,' J North  = ',I6,/&
         3x,'I West   =  ',I6,' K Bottom = ',I6,2/' The following ',   &
         'should conttain a wall cell and fluid cell:',/3x,'IJK  ',I9, &
         ' :: ',A3,/3x,'IJPK ',I9,' :: ',A3,2/' Maybe no IC was ',     &
         'specified for the fluid cell.')

! Store the new values in the global data array.
      BC_J_S(BCV) = J_S
      BC_J_N(BCV) = J_N


      J_WALL = BC_J_S(BCV)
      J_FLUID = merge(J_WALL-1, J_WALL+1, BC_PLANE(BCV)=='S')


! First pass through all of the BC region and verify that you have
! inflow/outflow specified against a wall. Flag any errors.
      ERROR = .FALSE.
      DO K = BC_K_B(BCV), BC_K_T(BCV)
      DO I = BC_I_W(BCV), BC_I_E(BCV)
         IF(.NOT.IS_ON_myPE_plus2layers(I,J_FLUID,K)) CYCLE
         IF(.NOT.IS_ON_myPE_plus2layers(I,J_WALL, K)) CYCLE
         IF(DEAD_CELL_AT(I,J_FLUID,K)) CYCLE
         IF(DEAD_CELL_AT(I,J_WALL, K)) CYCLE

         IJK_WALL  = FUNIJK(I,J_WALL, K)
         IJK_FLUID = FUNIJK(I,J_FLUID,K)

          IF(.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND. &
             ICBC_FLAG(IJK_FLUID)(1:1)=='.')) ERROR = .TRUE.

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

 1200 FORMAT('Error 1200: Illegal geometry for boundary condition:',I3)

         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)

            IF(.NOT.IS_ON_myPE_plus2layers(I,J_FLUID,K)) CYCLE
            IF(.NOT.IS_ON_myPE_plus2layers(I,J_WALL, K)) CYCLE
            IF(DEAD_CELL_AT(I, J_FLUID,K)) CYCLE
            IF(DEAD_CELL_AT(I, J_WALL, K)) CYCLE

            IJK_WALL  = FUNIJK(I,J_WALL ,K)
            IJK_FLUID = FUNIJK(I,J_FLUID,K)

            IF(.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND.                    &
               ICBC_FLAG(IJK_FLUID)(1:1)=='.')) THEN

               WRITE(ERR_MSG, 1201) &
                  I, J_WALL,  K, IJK_WALL, ICBC_FLAG(IJK_WALL),        &
                  I, J_FLUID, K, IJK_FLUID, ICBC_FLAG(IJK_FLUID)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

 1201 FORMAT(' ',/14X,'I',7X,'J',7X,'K',7X,'IJK',4x,'FLAG',/3x,        &
         'WALL ',3(2x,I6),2x,I9,3x,A,/3x,'FLUID',3(2x,I6),2x,I9,3x,A)

         ENDDO
         ENDDO

         WRITE(ERR_MSG,"('Please correct the mfix.dat file.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MOD_BC_J
