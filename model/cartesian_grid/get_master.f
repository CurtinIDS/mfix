!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_U_MASTER_CELLS                                     C
!  Purpose: Identify master cells for wall U-Momentum cells            C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_U_MASTER_CELLS

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE bc
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,IJKC,D,DIR(18),DMAX
      LOGICAL :: U_NODE,V_NODE,W_NODE,VEL_NODE,MASTER_FOUND
      INTEGER :: NC,L,BCV

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)'FINDING MASTER CELLS FOR U-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

!======================================================================
! Loop though all BC's and create a default NSW BC in case it is needed
! when the search for a master cell fails
!======================================================================
      DO L = 1, DIMENSION_BC
         IF (.NOT. (BC_DEFINED(L).OR.IS_CG(BC_TYPE_ENUM(L)))) THEN
            BC_TYPE_ENUM(L)=CG_NSW
            NSW_GHOST_BC_ID = L
            EXIT
         ENDIF
      ENDDO

!======================================================================
!  For each wall_u_cell: Probe neighboors and identify the master cell
!  as the first one having all velocity components available
!  A velocity component is available is the cell is not blocked nor a wall cell
!======================================================================

      NC = 0

      DO IJK = IJKSTART3, IJKEND3

        IF(WALL_U_AT(IJK)) THEN

            MASTER_FOUND = .FALSE.
            NC = NC + 1

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            DIR(1) = EAST_OF(IJK)
            DIR(2) = WEST_OF(IJK)

            DIR(3) = NORTH_OF(IJK)
            DIR(4) = SOUTH_OF(IJK)


            DIR(5) = EAST_OF(DIR(3))   ! NORTH-EAST
            DIR(6) = EAST_OF(DIR(4))   ! SOUTH-EAST

            DIR(7) = WEST_OF(DIR(3))   ! NORTH-WEST
            DIR(8) = WEST_OF(DIR(4))   ! SOUTH-WEST

            DIR(9) = TOP_OF(IJK)
            DIR(10) = BOTTOM_OF(IJK)

            DIR(11) = NORTH_OF(DIR(9)) ! NORTH-TOP
            DIR(12) = SOUTH_OF(DIR(9)) ! SOUTH-TOP

            DIR(13) = NORTH_OF(DIR(10)) ! NORTH-BOTTOM
            DIR(14) = SOUTH_OF(DIR(10)) ! SOUTH-BOTTOM

            DIR(15) = EAST_OF(DIR(9)) ! EAST-TOP
            DIR(16) = WEST_OF(DIR(9)) ! WEST-TOP

            DIR(17) = EAST_OF(DIR(10)) ! EAST-BOTTOM
            DIR(18) = WEST_OF(DIR(10)) ! WEST-BOTTOM

            IF(NO_K)   THEN
               DMAX = 4 !8         ! In 2D, probe E,W,N,S,NE,SE,NW,SW
            ELSE
               DMAX = 18        ! In 3D, probe T,B,NT,ST,NB,SB,ET,WT,EB,WB as well
            ENDIF

            DO D = 1,DMAX

               IJKC = DIR(D)

               U_NODE = ((.NOT.BLOCKED_U_CELL_AT(IJKC)).AND.(.NOT.WALL_U_AT(IJKC)))
               V_NODE = ((.NOT.BLOCKED_V_CELL_AT(IJKC)).AND.(.NOT.WALL_V_AT(IJKC)))

               IF(NO_K)   THEN
                  VEL_NODE = ((U_NODE).AND.(V_NODE))
               ELSE
                  W_NODE = ((.NOT.BLOCKED_W_CELL_AT(IJKC)).AND.(.NOT.WALL_W_AT(IJKC)))
                  VEL_NODE = ((U_NODE).AND.(V_NODE).AND.(W_NODE))
               ENDIF

               IF(U_NODE) THEN
                  U_MASTER_OF(IJK) = IJKC
                  MASTER_FOUND = .TRUE.
                  EXIT
               ENDIF
            ENDDO

            IF(.NOT.MASTER_FOUND) THEN
               BCV = BC_U_ID(IJK)
               IF(BCV>0) THEN
                  IF(BC_TYPE_ENUM(BCV) == CG_FSW) THEN
                     WRITE(*,*) ' WARNING IN SUBROUTINE: GET_U_MASTER_CELLS:'
                     WRITE(*,*) ' NO MASTER CELL FOUND FOR U_MOMENTUM WALL CELL:', IJK,I,J,K
                     WRITE(*,*) ' REVERTING TO NO SLIP WALL BOUNDARY CONDITION IN THIS CELL'
                     BC_U_ID(IJK) = NSW_GHOST_BC_ID
                     WRITE(*,*) ' BC_U_ID(IJK) = ', BC_U_ID(IJK)
                  ENDIF
!               WRITE(*,*) ' ERROR IN SUBROUTINE: GET_U_MASTER_CELLS:'
!               WRITE(*,*) ' NO MASTER CELL FOUND FOR U_MOMENTUM WALL CELL:', IJK,I,J,K
!               WRITE(*,*) ' MFIX WILL EXIT NOW.'
!               CALL MFIX_EXIT(myPE)

               ENDIF
            ENDIF
         ENDIF

      END DO



      RETURN

      END SUBROUTINE GET_U_MASTER_CELLS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_V_MASTER_CELLS                                     C
!  Purpose: Identify master cells for wall V-Momentum cells            C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_V_MASTER_CELLS

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE bc
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,IJKC,D,DIR(18),DMAX
      LOGICAL :: U_NODE,V_NODE,W_NODE,VEL_NODE,MASTER_FOUND
      INTEGER :: NC,BCV

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'FINDING MASTER CELLS FOR V-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

!======================================================================
!  For each wall_v_cell: Probe neighboors and identify the master cell
!  as the first one having all velocity components available
!  A velocity component is available is the cell is not blocked nor a wall cell
!======================================================================

      NC = 0

      DO IJK = IJKSTART3, IJKEND3

        IF(WALL_V_AT(IJK)) THEN

            MASTER_FOUND = .FALSE.
            NC = NC + 1

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)


            DIR(1) = EAST_OF(IJK)
            DIR(2) = WEST_OF(IJK)

            DIR(3) = NORTH_OF(IJK)
            DIR(4) = SOUTH_OF(IJK)


            DIR(5) = EAST_OF(DIR(3))   ! NORTH-EAST
            DIR(6) = EAST_OF(DIR(4))   ! SOUTH-EAST

            DIR(7) = WEST_OF(DIR(3))   ! NORTH-WEST
            DIR(8) = WEST_OF(DIR(4))   ! SOUTH-WEST

            DIR(9) = TOP_OF(IJK)
            DIR(10) = BOTTOM_OF(IJK)

            DIR(11) = NORTH_OF(DIR(9)) ! NORTH-TOP
            DIR(12) = SOUTH_OF(DIR(9)) ! SOUTH-TOP

            DIR(13) = NORTH_OF(DIR(10)) ! NORTH-BOTTOM
            DIR(14) = SOUTH_OF(DIR(10)) ! SOUTH-BOTTOM

            DIR(15) = EAST_OF(DIR(9)) ! EAST-TOP
            DIR(16) = WEST_OF(DIR(9)) ! WEST-TOP

            DIR(17) = EAST_OF(DIR(10)) ! EAST-BOTTOM
            DIR(18) = WEST_OF(DIR(10)) ! WEST-BOTTOM

            IF(NO_K)   THEN
               DMAX = 4 !8         ! In 2D, probe E,W,N,S,NE,SE,NW,SW
            ELSE
               DMAX = 18        ! In 3D, probe T,B,NT,ST,NB,SB,ET,WT,EB,WB as well
            ENDIF

            DO D = 1,DMAX

               IJKC = DIR(D)

               U_NODE = ((.NOT.BLOCKED_U_CELL_AT(IJKC)).AND.(.NOT.WALL_U_AT(IJKC)))
               V_NODE = ((.NOT.BLOCKED_V_CELL_AT(IJKC)).AND.(.NOT.WALL_V_AT(IJKC)))


               IF(NO_K)   THEN
                  VEL_NODE = ((U_NODE).AND.(V_NODE))
               ELSE
                  W_NODE = ((.NOT.BLOCKED_W_CELL_AT(IJKC)).AND.(.NOT.WALL_W_AT(IJKC)))
                  VEL_NODE = ((U_NODE).AND.(V_NODE).AND.(W_NODE))
               ENDIF

               IF(V_NODE) THEN
                  V_MASTER_OF(IJK) = IJKC
                  MASTER_FOUND = .TRUE.
                  EXIT
               ENDIF
            ENDDO

            IF(.NOT.MASTER_FOUND) THEN
               BCV = BC_V_ID(IJK)
               IF(BCV>0) THEN
                  IF(BC_TYPE_ENUM(BCV) == CG_FSW) THEN
                     WRITE(*,*) ' WARNING IN SUBROUTINE: GET_V_MASTER_CELLS:'
                     WRITE(*,*) ' NO MASTER CELL FOUND FOR V_MOMENTUM WALL CELL:', IJK,I,J,K
                     WRITE(*,*) ' REVERTING TO NO SLIP WALL BOUNDARY CONDITION IN THIS CELL'
                     BC_V_ID(IJK) = NSW_GHOST_BC_ID
                     WRITE(*,*) ' BC_V_ID(IJK) = ', BC_V_ID(IJK)
                  ENDIF
!               WRITE(*,*) ' ERROR IN SUBROUTINE: GET_V_MASTER_CELLS:'
!               WRITE(*,*) ' NO MASTER CELL FOUND FOR V_MOMENTUM WALL CELL:', IJK,I,J,K
!               WRITE(*,*) ' MFIX WILL EXIT NOW.'
!               CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDIF

      END DO



      RETURN

      END SUBROUTINE GET_V_MASTER_CELLS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_W_MASTER_CELLS                                     C
!  Purpose: Identify master cells for wall W-Momentum cells            C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_W_MASTER_CELLS

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE bc
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,IJKC,D,DIR(18),DMAX
      LOGICAL :: U_NODE,V_NODE,W_NODE,VEL_NODE,MASTER_FOUND
      INTEGER :: NC,BCV

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'FINDING MASTER CELLS FOR W-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

!======================================================================
!  For each wall_w_cell: Probe neighboors and identify the master cell
!  as the first one having all velocity components available
!  A velocity component is available is the cell is not blocked nor a wall cell
!======================================================================

      NC = 0

      DO IJK = IJKSTART3, IJKEND3

         IF(WALL_W_AT(IJK)) THEN

            MASTER_FOUND = .FALSE.
            NC = NC + 1

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)


            DIR(1) = EAST_OF(IJK)
            DIR(2) = WEST_OF(IJK)

            DIR(3) = NORTH_OF(IJK)
            DIR(4) = SOUTH_OF(IJK)


            DIR(5) = EAST_OF(DIR(3))   ! NORTH-EAST
            DIR(6) = EAST_OF(DIR(4))   ! SOUTH-EAST

            DIR(7) = WEST_OF(DIR(3))   ! NORTH-WEST
            DIR(8) = WEST_OF(DIR(4))   ! SOUTH-WEST

            DIR(9) = TOP_OF(IJK)
            DIR(10) = BOTTOM_OF(IJK)

            DIR(11) = NORTH_OF(DIR(9)) ! NORTH-TOP
            DIR(12) = SOUTH_OF(DIR(9)) ! SOUTH-TOP

            DIR(13) = NORTH_OF(DIR(10)) ! NORTH-BOTTOM
            DIR(14) = SOUTH_OF(DIR(10)) ! SOUTH-BOTTOM

            DIR(15) = EAST_OF(DIR(9)) ! EAST-TOP
            DIR(16) = WEST_OF(DIR(9)) ! WEST-TOP

            DIR(17) = EAST_OF(DIR(10)) ! EAST-BOTTOM
            DIR(18) = WEST_OF(DIR(10)) ! WEST-BOTTOM

            IF(NO_K)   THEN
               DMAX = 8         ! In 2D, probe E,W,N,S,NE,SE,NW,SW
            ELSE
               DMAX = 18        ! In 3D, probe T,B,NT,ST,NB,SB,ET,WT,EB,WB as well
            ENDIF

            DO D = 1,DMAX

               IJKC = DIR(D)

               U_NODE = ((.NOT.BLOCKED_U_CELL_AT(IJKC)).AND.(.NOT.WALL_U_AT(IJKC)))
               V_NODE = ((.NOT.BLOCKED_V_CELL_AT(IJKC)).AND.(.NOT.WALL_V_AT(IJKC)))


               IF(NO_K)   THEN
                  VEL_NODE = ((U_NODE).AND.(V_NODE))
               ELSE
                  W_NODE = ((.NOT.BLOCKED_W_CELL_AT(IJKC)).AND.(.NOT.WALL_W_AT(IJKC)))
                  VEL_NODE = ((U_NODE).AND.(V_NODE).AND.(W_NODE))
               ENDIF

               IF(W_NODE) THEN
                  W_MASTER_OF(IJK) = IJKC
                  MASTER_FOUND = .TRUE.
                  EXIT
               ENDIF
            ENDDO

            IF(.NOT.MASTER_FOUND) THEN
               BCV = BC_W_ID(IJK)
               IF(BCV>0) THEN
                  IF(BC_TYPE_ENUM(BCV) == CG_FSW) THEN
                     WRITE(*,*) ' WARNING IN SUBROUTINE: GET_W_MASTER_CELLS:'
                     WRITE(*,*) ' NO MASTER CELL FOUND FOR W_MOMENTUM WALL CELL:', IJK,I,J,K
                     WRITE(*,*) ' REVERTING TO NO SLIP WALL BOUNDARY CONDITION IN THIS CELL'
                     BC_W_ID(IJK) = NSW_GHOST_BC_ID
                     WRITE(*,*) ' BC_W_ID(IJK) = ', BC_W_ID(IJK)
                  ENDIF
!               WRITE(*,*) ' ERROR IN SUBROUTINE: GET_W_MASTER_CELLS:'
!               WRITE(*,*) ' NO MASTER CELL FOUND FOR W_MOMENTUM WALL CELL:', IJK,I,J,K
!               WRITE(*,*) ' MFIX WILL EXIT NOW.'
!               CALL MFIX_EXIT(myPE)
               ENDIF

            ENDIF
         ENDIF
      END DO


      RETURN

      END SUBROUTINE GET_W_MASTER_CELLS
