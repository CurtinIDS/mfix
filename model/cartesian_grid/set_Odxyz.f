!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_ODXYZ_U_CUT_CELL                                   C
!  Purpose: Set 1/dx, 1/dy, and 1/dz for U-Momentum cell               C
!           (only when cartesian grid is used)                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_ODXYZ_U_CUT_CELL

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE vtk
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: IP,JP,KP
      INTEGER :: IJK1,IJK2

      IF(MyPE == PE_IO) THEN
         IF(NO_K) THEN
            WRITE(*,10)'COMPUTING 1/DX, 1/DY FOR U-MOMENTUM CELLS...'
         ELSE
            WRITE(*,10)'COMPUTING 1/DX, 1/DY, 1/DZ FOR U-MOMENTUM CELLS...'
         ENDIF
      ENDIF
10    FORMAT(1X,A)

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

!======================================================================
!  1/dx at East face of U-Momentume cell
!======================================================================

         IF(I == ISTART3) THEN
            IJK2 = FUNIJK(ISTART1,J,K)
            IJK1 = FUNIJK(ISTART1-1,J,K)
         ELSEIF(I == IEND3) THEN
            IJK2 = FUNIJK(IEND3,J,K)
            IJK1 = FUNIJK(IEND3-1,J,K)
         ELSE
            IP = I + 1
            IJK2 = FUNIJK(IP,J,K)
            IJK1 = IJK
         ENDIF

         IF(X_U(IJK2)/=X_U(IJK1)) THEN
            ONEoDX_E_U(IJK) = ONE / (X_U(IJK2)-X_U(IJK1))
         ELSE
            ONEoDX_E_U(IJK) = ZERO
         ENDIF

!======================================================================
!  1/dy at North face of U_Momentum cell
!======================================================================

         IF(J == JSTART3) THEN
            IJK2 = FUNIJK(I,JSTART1,K)
            IJK1 = FUNIJK(I,JSTART1-1,K)
         ELSEIF(J == JEND3) THEN
            IJK2 = FUNIJK(I,JEND3,K)
            IJK1 = FUNIJK(I,JEND3-1,K)
         ELSE
            JP = J + 1
            IJK2 = FUNIJK(I,JP,K)
            IJK1 = IJK
         ENDIF

         IF(Y_U(IJK2)/=Y_U(IJK1)) THEN
            ONEoDY_N_U(IJK) = ONE / (Y_U(IJK2)-Y_U(IJK1))
         ELSE
            ONEoDY_N_U(IJK) = ZERO
         ENDIF

!======================================================================
!  1/dz at Top face of U_Momentum cell
!======================================================================

         IF(DO_K) THEN
            IF(K == KSTART3) THEN
               IJK2 = FUNIJK(I,J,KSTART1)
               IJK1 = FUNIJK(I,J,KSTART1-1)
            ELSEIF(K == KEND3) THEN
               IJK2 = FUNIJK(I,J,KEND3)
               IJK1 = FUNIJK(I,J,KEND3-1)
            ELSE
               KP = K + 1
               IJK2 = FUNIJK(I,J,KP)
               IJK1 = IJK
            ENDIF

            IF(Z_U(IJK2)/=Z_U(IJK1)) THEN
               ONEoDZ_T_U(IJK) = ONE / (Z_U(IJK2)-Z_U(IJK1))
            ELSE
               ONEoDZ_T_U(IJK) = ZERO
            ENDIF

         ELSE
            ONEoDZ_T_U = ONE / ZLENGTH
         ENDIF

      END DO


      RETURN


      END SUBROUTINE SET_ODXYZ_U_CUT_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_ODXYZ_V_CUT_CELL                                   C
!  Purpose: Set 1/dx, 1/dy, and 1/dz for V-Momentum cell               C
!           (only when cartesian grid is used)                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_ODXYZ_V_CUT_CELL

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE vtk
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: IP,JP,KP
      INTEGER :: IJK1,IJK2

      IF(MyPE == PE_IO) THEN
         IF(NO_K) THEN
            WRITE(*,10)'COMPUTING 1/DX, 1/DY FOR V-MOMENTUM CELLS...'
         ELSE
            WRITE(*,10)'COMPUTING 1/DX, 1/DY, 1/DZ FOR V-MOMENTUM CELLS...'
         ENDIF
      ENDIF
10    FORMAT(1X,A)

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

!======================================================================
!  1/dx at East face of V-Momentume cell
!======================================================================

         IF(I == ISTART3) THEN
            IJK2 = FUNIJK(ISTART1,J,K)
            IJK1 = FUNIJK(ISTART1-1,J,K)
         ELSEIF(I == IEND3) THEN
            IJK2 = FUNIJK(IEND3,J,K)
            IJK1 = FUNIJK(IEND3-1,J,K)
         ELSE
            IP = I + 1
            IJK2 = FUNIJK(IP,J,K)
            IJK1 = IJK
         ENDIF

         IF(X_V(IJK2)/=X_V(IJK1)) THEN
            ONEoDX_E_V(IJK) = ONE / (X_V(IJK2)-X_V(IJK1))
         ELSE
            ONEoDX_E_V(IJK) = ZERO
         ENDIF

!======================================================================
!  1/dy at North face of V_Momentum cell
!======================================================================

         IF(J == JSTART3) THEN
            IJK2 = FUNIJK(I,JSTART1,K)
            IJK1 = FUNIJK(I,JSTART1-1,K)
         ELSEIF(J == JEND3) THEN
            IJK2 = FUNIJK(I,JEND3,K)
            IJK1 = FUNIJK(I,JEND3-1,K)
         ELSE
            JP = J + 1
            IJK2 = FUNIJK(I,JP,K)
            IJK1 = IJK
         ENDIF

         IF(Y_V(IJK2)/=Y_V(IJK1)) THEN
            ONEoDY_N_V(IJK) = ONE / (Y_V(IJK2)-Y_V(IJK1))
         ELSE
            ONEoDY_N_V(IJK) = ZERO
         ENDIF

!======================================================================
!  1/dz at Top face of V_Momentum cell
!======================================================================

         IF(DO_K) THEN
            IF(K == KSTART3) THEN
               IJK2 = FUNIJK(I,J,KSTART1)
               IJK1 = FUNIJK(I,J,KSTART1-1)
            ELSEIF(K == KEND3) THEN
               IJK2 = FUNIJK(I,J,KEND3)
               IJK1 = FUNIJK(I,J,KEND3-1)
            ELSE
               KP = K + 1
               IJK2 = FUNIJK(I,J,KP)
               IJK1 = IJK
            ENDIF

            IF(Z_V(IJK2)/=Z_V(IJK1)) THEN
               ONEoDZ_T_V(IJK) = ONE / (Z_V(IJK2)-Z_V(IJK1))
            ELSE
               ONEoDZ_T_V(IJK) = ZERO
            ENDIF

         ELSE

            ONEoDZ_T_V = ONE / ZLENGTH

         ENDIF

      END DO



      RETURN


      END SUBROUTINE SET_ODXYZ_V_CUT_CELL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_ODXYZ_W_CUT_CELL                                   C
!  Purpose: Set 1/dx, 1/dy, and 1/dz for W-Momentum cell               C
!           (only when cartesian grid is used)                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_ODXYZ_W_CUT_CELL

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE vtk
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: IP,JP,KP
      INTEGER :: IJK1,IJK2

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'COMPUTING 1/DX, 1/DY, 1/DZ FOR W-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

!======================================================================
!  1/dx at East face of V-Momentume cell
!======================================================================

         IF(I == ISTART3) THEN
            IJK2 = FUNIJK(ISTART1,J,K)
            IJK1 = FUNIJK(ISTART1-1,J,K)
         ELSEIF(I == IEND3) THEN
            IJK2 = FUNIJK(IEND3,J,K)
            IJK1 = FUNIJK(IEND3-1,J,K)
         ELSE
            IP = I + 1
            IJK2 = FUNIJK(IP,J,K)
            IJK1 =IJK
         ENDIF

         IF(X_W(IJK2)/=X_W(IJK1)) THEN
            ONEoDX_E_W(IJK) = ONE / (X_W(IJK2)-X_W(IJK1))
         ELSE
            ONEoDX_E_W(IJK) = ZERO
         ENDIF
!======================================================================
!  1/dy at North face of V_Momentum cell
!======================================================================

         IF(J == JSTART3) THEN
            IJK2 = FUNIJK(I,JSTART1,K)
            IJK1 = FUNIJK(I,JSTART1-1,K)
         ELSEIF(J == JEND3) THEN
            IJK2 = FUNIJK(I,JEND3,K)
            IJK1 = FUNIJK(I,JEND3-1,K)
         ELSE
            JP = J + 1
            IJK2 = FUNIJK(I,JP,K)
            IJK1 = IJK
         ENDIF

         IF(Y_W(IJK2)/=Y_W(IJK1)) THEN
            ONEoDY_N_W(IJK) = ONE / (Y_W(IJK2)-Y_W(IJK1))
         ELSE
            ONEoDY_N_W(IJK) = ZERO
         ENDIF
!         print*,'myPE,IJK,ONEoDY_N_W(IJK)=',myPE,IJK,ONEoDY_N_W(IJK)

!======================================================================
!  1/dz at Top face of V_Momentum cell
!======================================================================

         IF(K == KSTART3) THEN
            IJK2 = FUNIJK(I,J,KSTART1)
            IJK1 = FUNIJK(I,J,KSTART1-1)
         ELSEIF(K == KEND3) THEN
            IJK2 = FUNIJK(I,J,KEND3)
            IJK1 = FUNIJK(I,J,KEND3-1)
         ELSE
            KP = K + 1
            IJK2 = FUNIJK(I,J,KP)
            IJK1 = IJK
         ENDIF

         IF(Z_W(IJK2)/=Z_W(IJK1)) THEN
            ONEoDZ_T_W(IJK) = ONE / (Z_W(IJK2)-Z_W(IJK1))
         ELSE
            ONEoDZ_T_W(IJK) = ZERO
         ENDIF

      END DO



      RETURN


      END SUBROUTINE SET_ODXYZ_W_CUT_CELL
