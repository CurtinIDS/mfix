!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MOD_BC_K (BC, I_w, J_s, K_b, K_t, PLANE)               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: modify the "K" values for the b.c. plane                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MOD_BC_K(BCV)

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

! calculated cell indices in I,J,K directions
      INTEGER :: K_b, K_t
      INTEGER :: I_w, J_s

      INTEGER :: OWNER

      INTEGER :: I, J
      INTEGER :: IJK, IJKP

      INTEGER :: IER
      LOGICAL :: ERROR
      INTEGER :: K_FLUID, IJK_FLUID
      INTEGER :: K_WALL,  IJK_WALL


!-----------------------------------------------

      CALL INIT_ERR_MSG("MOD_BC_K")

      K_B = BC_K_B(BCV)
      K_T = BC_K_T(BCV)

      I_W = BC_I_W(BCV)
      J_S = BC_J_S(BCV)


! Establish the OWNER of the BC
      OWNER = merge(myPE, 0, IS_ON_myPE_owns(I_W, J_S, K_B))
      CALL GLOBAL_ALL_SUM(OWNER)

      IF(myPE == OWNER) THEN

         IJK  = FUNIJK(I_W, J_S, K_B)
         IJKP = FUNIJK(I_W, J_S, K_B+1)

         IF(WALL_ICBC_FLAG(IJK) .AND. ICBC_FLAG(IJKP)(1:1)=='.')THEN
            K_B = K_B
            K_T = K_T
            BC_PLANE(BCV) = 'T'
         ELSEIF(WALL_ICBC_FLAG(IJKP) .AND. ICBC_FLAG(IJK)(1:1)=='.')THEN
            K_B = K_B + 1
            K_T = K_T + 1
            BC_PLANE(BCV) = 'B'
         ELSE
            BC_PLANE(BCV) = '.'
         ENDIF
      ENDIF

! The owner distributes the new Iw/Ie coordinates to the other ranks.
      CALL BCAST(K_B,OWNER)
      CALL BCAST(K_T,OWNER)
      CALL BCAST(BC_PLANE(BCV),OWNER)

! If there is an error, send IJK/IPJK to all ranks. Report and exit.
      IF(BC_PLANE(BCV) == '.') THEN
         CALL BCAST(IJKP,OWNER)
         CALL BCAST(IJK, OWNER)

         WRITE(ERR_MSG, 1100) BCV, K_B, K_T, I_W, J_S,                 &
            IJK, ICBC_FLAG(IJK),  IJKP, ICBC_FLAG(IJKP)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Cannot locate flow plane for boundary ',     &
         'condition ',I3,'.',2/3x,'K Bottom =  ',I6,' K Top    = ',I6,/&
         3x,'I West   =  ',I6,' J South  = ',I6,2/' The following ',   &
         'should conttain a wall cell and fluid cell:',/3x,'IJK  ',I9, &
         ' :: ',A3,/3x,'IJKP ',I9,' :: ',A3,2/' Maybe no IC was ',     &
         'specified for the fluid cell.')

! Store the new values in the global data array.
      BC_K_B(BCV) = K_B
      BC_K_T(BCV) = K_T

! Set up the I-indices for checking the entire BC region.
      K_WALL = BC_K_B(BCV)
      K_FLUID = merge(K_WALL-1, K_WALL+1, BC_PLANE(BCV)=='B')


      ERROR = .FALSE.
      DO J = BC_J_S(BCV), BC_J_N(BCV)
      DO I = BC_I_W(BCV), BC_I_E(BCV)

! Only check cells that you own and contain fluid.
         IF(.NOT.IS_ON_myPE_plus2layers(I,J,K_FLUID)) CYCLE
         IF(.NOT.IS_ON_myPE_plus2layers(I,J,K_WALL )) CYCLE
         IF(DEAD_CELL_AT(I,J,K_FLUID)) CYCLE
         IF(DEAD_CELL_AT(I,J,K_WALL )) CYCLE

         IJK_WALL = FUNIJK(I,J,K_WALL)
         IJK_FLUID = FUNIJK(I,J,K_FLUID)

         IF(.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND.                       &
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

         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)

! Only check cells that you own and contain fluid.
            IF(.NOT.IS_ON_myPE_plus2layers(I,J,K_FLUID)) CYCLE
            IF(.NOT.IS_ON_myPE_plus2layers(I,J,K_WALL )) CYCLE
            IF(DEAD_CELL_AT(I,J,K_FLUID)) CYCLE
            IF(DEAD_CELL_AT(I,J,K_WALL )) CYCLE

            IJK_WALL = FUNIJK(I,J,K_WALL)
            IJK_FLUID = FUNIJK(I,J,K_FLUID)

            IF(.NOT.(WALL_ICBC_FLAG(IJK_WALL) .AND.                    &
               ICBC_FLAG(IJK_FLUID)(1:1)=='.')) THEN

               WRITE(ERR_MSG, 1201) &
                  I, J, K_WALL,  IJK_WALL,  ICBC_FLAG(IJK_WALL),       &
                  I, J, K_FLUID, IJK_FLUID, ICBC_FLAG(IJK_FLUID)
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
      END SUBROUTINE MOD_BC_K
