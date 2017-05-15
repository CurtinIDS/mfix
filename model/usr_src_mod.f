! or put the actual routine calc_usr_source here? move logic variable
!  usr_source into run or something...
      module usr_src
      use param, only: dim_eqs

! If .TRUE. call user-defined source subroutines
      LOGICAL :: CALL_USR_SOURCE(DIM_EQS)

      ENUM, BIND(C)
         ENUMERATOR :: Pressure_Correction, Solids_Correction
         ENUMERATOR :: Gas_Continuity, Solids_Continuity
         ENUMERATOR :: Gas_U_Mom, Solids_U_Mom, Gas_V_Mom, Solids_V_Mom
         ENUMERATOR :: Gas_W_Mom, Solids_W_Mom
         ENUMERATOR :: Gas_Energy, Solids_Energy
         ENUMERATOR :: Gas_Species, Solids_Species
         ENUMERATOR :: Gran_Energy
         ENUMERATOR :: Usr_Scalar, K_Epsilon_K, K_Epsilon_E
         ENUMERATOR :: BLANK
         ENUMERATOR :: DUMMY
      END ENUM

      CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_USR_SOURCE                                         C
!  Purpose: Driver routine to calculate user defined source terms      C
!  for each of the various equations                                   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_USR_SOURCE(lEQ_NO, A_M, B_M, lB_MMAX, lM, lN)

! Modules
!-----------------------------------------------
      use compar, only: ijkstart3, ijkend3
      use fldvar, only: ep_s, ep_g
      use fun_avg, only: avg_x, avg_y, avg_z
      use functions, only: fluid_at, wall_at
      use functions, only: ip_at_e, ip_at_n, ip_at_t
      use functions, only: sip_at_e, sip_at_n, sip_at_t
      use functions, only: east_of, north_of, top_of
      use indices, only: i_of, j_of, k_of
      use param, only: dimension_3, dimension_m
      use param1, only: undefined_i, zero
      use pgcor, only: phase_4_p_g
      use physprop, only: smax, mmax
      use pscor, only: phase_4_p_s
      use run, only: kt_type_enum, ghd_2007
      use run, only: momentum_x_eq, momentum_y_eq, momentum_z_eq
      use toleranc, only: dil_ep_s
      use error_manager
      IMPLICIT NONE
! Dummy arguments
!-----------------------------------------------
! Equation number
! Table relating the eq_no set in the mfix.dat to the corresponding
! pde/variable of interest
!     1, PP_g (pressure correction eq)
!     2, EPP (solids correction eq)
!     2, rop_g (gas continuity eq)
!     2, rop_s (solids continuity eq)
!     3, u_g (gas u-momentum eq)
!     3, u_s (solids u-momentum eq)
!     4, v_g (gas v-momentum eq)
!     4, v_s (solids v-momentum eq)
!     5, w_g (gas w-momentum eq)
!     5, w_s (solids w-momentum eq)
!     6, T_g (gas energy eq)
!     6, T_s (solids energy eq)
!     7, x_g (gas species eq)
!     7, x_s (solids species eq)
!     8, theta_m (granular energy eq)
!     9, scalar_eq (scalar eq)
!     9, k_epsilon k (k_epsilon eqs)
!     9, k_epsilon e (k_epsilon_eqs)
      INTEGER, INTENT(IN) :: lEQ_NO

! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Vector b_mmax
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT) :: lB_mmax(DIMENSION_3, 0:DIMENSION_M)
! Phase index
      INTEGER, OPTIONAL, INTENT(IN) :: lM
! Species index OR scalar equation number
      INTEGER, OPTIONAL, INTENT(IN) :: lN

! Local variables
!-----------------------------------------------
! source terms which appear appear in the
! center coefficient (lhs) - part of a_m matrix
! source vector (rhs) - part of b_m vector
      DOUBLE PRECISION :: sourcelhs, sourcerhs
! indices
      INTEGER :: IJK, I, J, K, IJKE, IJKN, IJKT
! tmp indices
      INTEGER :: L, ll, M, N
! volume fractions
      DOUBLE PRECISION :: EPGA, EPSA, EPtmp
!-----------------------------------------------

! set local values for phase index and species or scalar index
      IF (.NOT.present(lM)) THEN
         M = UNDEFINED_I
      ELSE
         M = lM
      ENDIF
      IF (.NOT.present(lN)) THEN
         N = UNDEFINED_I
      ELSE
         N = lN
      ENDIF

      SELECT CASE(lEQ_NO)

! gas pressure correction equation
      CASE(Pressure_Correction)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
               CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
               A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
               B_M(IJK,M) = B_M(IJK,M) - sourcerhs
               lB_MMAX(IJK,M) = max(abs(lB_MMAX(IJK,M)), abs(B_M(IJK,M)))
            ENDIF
         ENDDO

! solids correction currently only called when smax==1 and
! mcp/=undefined_i
      CASE(Solids_correction)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
               CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
               A_M(IJK,0,0) = A_M(IJK,0,0) - sourcelhs
               B_M(IJK,0) = B_M(IJK,0) - sourcerhs
               lB_MMAX(IJK,0) = max(abs(lB_MMAX(IJK,0)), abs(B_M(IJK,0)))
            ENDIF
         ENDDO

! gas continuity
      CASE(Gas_continuity)
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!                            .AND. PHASE_4_P_G(IJK)/=M) THEN
! have not included phase_4_p_g logic...which would require unique
! structure here
               CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
               A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
               B_M(IJK,M) = B_M(IJK,M) - sourcerhs
            ENDIF
         ENDDO

! soldis continuity
      CASE(Solids_continuity)
! have not included phase_4_p_g/p_s logic...which would require unique
! structure here
         DO IJK=IJKSTART3,IJKEND3
            IF (FLUID_AT(IJK)) THEN
!               .AND. PHASE_4_P_G(IJK)/=M .AND. PHASE_4_P_S(IJK)/=M) THEN
               CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
               A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
               B_M(IJK,M) = B_M(IJK,M) - sourcerhs
            ENDIF
         ENDDO

! u-momentum equations
      CASE(Gas_U_Mom)
         M=0
         IF(.NOT.MOMENTUM_X_EQ(M)) RETURN
         DO IJK=IJKSTART3,IJKEND3
            IF (.NOT.FLUID_AT(IJK)) CYCLE
            I = I_OF(IJK)
            IJKE = EAST_OF(IJK)
            EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)
            IF (IP_AT_E(IJK) .OR. EPGA <= DIL_EP_S) CYCLE
            CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
            A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
            B_M(IJK,M) = B_M(IJK,M) - sourcerhs
         ENDDO

! u-momentum equations
      CASE(Solids_U_Mom)
         DO L=1,MMAX
            IF (KT_TYPE_ENUM /= GHD_2007 .OR. &
               (KT_TYPE_ENUM == GHD_2007 .AND. L==MMAX)) THEN
               IF(.NOT.MOMENTUM_X_EQ(L)) RETURN
               DO IJK=IJKSTART3,IJKEND3
                  IF (.NOT.FLUID_AT(IJK)) CYCLE
                  I = I_OF(IJK)
                  IJKE = EAST_OF(IJK)
                  EPtmp = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                     DO lL = 1, SMAX
                        EPtmp = EPtmp + AVG_X(EP_S(IJK,lL),EP_S(IJKE,lL),I)
                     ENDDO
                     EPSA = EPtmp
                  ELSE
                     EPSA = AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I)
                  ENDIF
                  IF (IP_AT_E(IJK) .OR. SIP_AT_E(IJK) .OR. &
                      EPSA <= DIL_EP_S) CYCLE
                  CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, L, N)
                  A_M(IJK,0,L) = A_M(IJK,0,L) - sourcelhs
                  B_M(IJK,L) = B_M(IJK,L) - sourcerhs
               ENDDO    ! enddo ijk
            ENDIF
         ENDDO   ! enddo mmax

! v-momentum equations
      CASE(Gas_V_Mom)
         M=0
         IF(.NOT.MOMENTUM_Y_EQ(M)) RETURN
         DO IJK=IJKSTART3,IJKEND3
            IF (.NOT.FLUID_AT(IJK)) CYCLE
            J = J_OF(IJK)
            IJKN = NORTH_OF(IJK)
            EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)
            IF (IP_AT_N(IJK) .OR. EPGA <= DIL_EP_S) CYCLE
            CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
            A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
            B_M(IJK,M) = B_M(IJK,M) - sourcerhs
         ENDDO

! v-momentum equations
      CASE(Solids_V_Mom)
         DO L=1,MMAX
            IF (KT_TYPE_ENUM /= GHD_2007 .OR. &
                (KT_TYPE_ENUM == GHD_2007 .AND. L==MMAX)) THEN
               IF(.NOT.MOMENTUM_Y_EQ(L)) RETURN
               DO IJK=IJKSTART3,IJKEND3
                  IF (.NOT.FLUID_AT(IJK)) CYCLE
                  J = J_OF(IJK)
                  IJKN = NORTH_OF(IJK)
                  EPtmp = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                     DO lL = 1, SMAX
                        EPtmp = EPtmp + AVG_Y(EP_S(IJK,lL),EP_S(IJKN,lL),J)
                     ENDDO
                     EPSA = EPtmp
                  ELSE
                     EPSA = AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J)
                  ENDIF
                  IF (IP_AT_N(IJK) .OR. SIP_AT_N(IJK) .OR. &
                      EPSA <= DIL_EP_S) CYCLE
                  CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, L, N)
                  A_M(IJK,0,L) = A_M(IJK,0,L) - sourcelhs
                  B_M(IJK,L) = B_M(IJK,L) - sourcerhs
               ENDDO   ! enddo ijk
            ENDIF
         ENDDO   ! enddo mmax

! w-momentum equations
      CASE (Gas_W_Mom)
         M=0
         IF(.NOT.MOMENTUM_Z_EQ(M)) RETURN
         DO IJK=IJKSTART3,IJKEND3
            IF (.NOT.FLUID_AT(IJK)) CYCLE
            K = K_OF(IJK)
            IJKT = TOP_OF(IJK)
            EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
            IF (IP_AT_T(IJK) .OR. EPGA <= DIL_EP_S) CYCLE
            CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
            A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
            B_M(IJK,M) = B_M(IJK,M) - sourcerhs
         ENDDO

! w-momentum equations
      CASE (Solids_W_Mom)
         DO L=1,MMAX
            IF (KT_TYPE_ENUM /= GHD_2007 .OR. &
                (KT_TYPE_ENUM == GHD_2007 .AND. L==MMAX)) THEN
               IF(.NOT.MOMENTUM_Z_EQ(L)) RETURN
               DO IJK=IJKSTART3,IJKEND3
                  IF (.NOT.FLUID_AT(IJK)) CYCLE
                  K = K_OF(IJK)
                  IJKT = TOP_OF(IJK)
                  EPtmp = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                     DO lL = 1, SMAX
                        EPtmp = EPtmp + AVG_Z(EP_S(IJK,lL),EP_S(IJKT,lL),K)
                     ENDDO
                     EPSA = EPtmp
                  ELSE
                     EPSA = AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)
                  ENDIF
                  IF (IP_AT_T(IJK) .OR. SIP_AT_T(IJK) .OR. &
                      EPSA <= DIL_EP_S) CYCLE
                  CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, L, N)
                  A_M(IJK,0,L) = A_M(IJK,0,L) - sourcelhs
                  B_M(IJK,L) = B_M(IJK,L) - sourcerhs
               ENDDO   ! enddo ijk
            ENDIF
         ENDDO   ! enddo mmax


! gas and solids energy equations (6)
! gas and solids species mass fractions (7)
! scalars, k_epsilon k and e (9)
      CASE (Gas_Energy, Solids_Energy, Gas_Species, Solids_Species,&
            Usr_Scalar, K_Epsilon_K, K_Epsilon_E )
          DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               IF (M==0) THEN
                  eptmp=ep_g(ijk)
               ELSE
                  eptmp=ep_s(ijk,M)
               ENDIF
               IF (eptmp <= DIL_EP_S) CYCLE
               CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
               A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
               B_M(IJK,M) = B_M(IJK,M) - sourcerhs
            ENDIF
         ENDDO


      CASE (Gran_Energy)   ! unique due to ghd
! granular temperature
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               IF (KT_TYPE_ENUM == GHD_2007) THEN
                  eptmp = zero
                  DO lL=1,SMAX
                     eptmp=eptmp+EP_S(IJK,lL)
                  ENDDO
               ELSE
                  eptmp = ep_s(IJK,M)
               ENDIF
            ENDIF
            IF (eptmp<=dil_ep_s) cycle
            CALL USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)
            A_M(IJK,0,M) = A_M(IJK,0,M) - sourcelhs
            B_M(IJK,M) = B_M(IJK,M) - sourcerhs
         ENDDO


! error out
      CASE DEFAULT
! should never hit this
! Initialize the error manager.
         CALL INIT_ERR_MSG("CALC_USR_SOURCE")
         WRITE(ERR_MSG, 1001) ival(leq_no)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1001 FORMAT('Error 1101: Unknown Equation= ', A)
         CALL FINL_ERR_MSG
      END SELECT   ! end selection of usr_source equation

      RETURN
      END SUBROUTINE CALC_USR_SOURCE

      end module usr_src

