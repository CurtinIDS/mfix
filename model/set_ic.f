!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: This module sets all the initial conditions.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE SET_IC

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE constant
      USE physprop
      USE ic
      USE fldvar
      USE visc_g
      USE indices
      USE scales
      USE energy
      USE scalars
      USE compar
      USE run
      USE sendrecv
      USE solids_pressure
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, IJK, M
! Local index for initial condition
      INTEGER :: L
! Temporary variable for storing IC_EP_g
      DOUBLE PRECISION :: EPGX
! Temporary variable for storing IC_P_g
      DOUBLE PRECISION :: PGX
! Temporary variable for storing P_star
      DOUBLE PRECISION :: PSX
! Temporary variable for storing IC_T_g
      DOUBLE PRECISION :: TGX
! Temporary variable for storing IC_U_g
      DOUBLE PRECISION :: UGX
! Temporary variable for storing IC_V_g
      DOUBLE PRECISION :: VGX
! Temporary variable for storing IC_W_g
      DOUBLE PRECISION :: WGX
! Temporary variable for storing IC_ROP_s
      DOUBLE PRECISION :: ROPSX (DIMENSION_M)
! Temporary variable for storing IC_T_s
      DOUBLE PRECISION :: TSX (DIMENSION_M)
! Temporary variable for storing IC_U_s
      DOUBLE PRECISION :: USX (DIMENSION_M)
! Temporary variable for storing IC_V_s
      DOUBLE PRECISION :: VSX (DIMENSION_M)
! Temporary variable for storing IC_W_s
      DOUBLE PRECISION :: WSX (DIMENSION_M)
! number density for GHD theory
      DOUBLE PRECISION :: nM, nTOT
!-----------------------------------------------

!  Set the initial conditions.
      DO L = 1, DIMENSION_IC
         IF (IC_DEFINED(L)) THEN

            EPGX = IC_EP_G(L)
            PGX = IC_P_G(L)
            PSX = IC_P_STAR(L)
            IF (PSX==UNDEFINED .AND. IC_TYPE(L)/='PATCH') PSX = ZERO
            TGX = IC_T_G(L)
            UGX = IC_U_G(L)
            VGX = IC_V_G(L)
            WGX = IC_W_G(L)

            M = 1
            IF (MMAX > 0) THEN
              ROPSX(:MMAX) = IC_ROP_S(L,:MMAX)
              TSX(:MMAX) = IC_T_S(L,:MMAX)
              USX(:MMAX) = IC_U_S(L,:MMAX)
              VSX(:MMAX) = IC_V_S(L,:MMAX)
              WSX(:MMAX) = IC_W_S(L,:MMAX)
              M = MMAX + 1
            ENDIF

            DO K = IC_K_B(L), IC_K_T(L)
            DO J = IC_J_S(L), IC_J_N(L)
            DO I = IC_I_W(L), IC_I_E(L)
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)

               IF (.NOT.WALL_AT(IJK)) THEN
                  IF (EPGX /= UNDEFINED) EP_G(IJK) = EPGX

                  IF (IC_TYPE(L) == 'PATCH') THEN
                      IF (PGX /= UNDEFINED) P_G(IJK) = SCALE_PRESSURE(PGX)
                  ELSE
                     P_G(IJK) = merge(SCALE_PRESSURE(PGX), UNDEFINED,           &
                        PGX /= UNDEFINED)
                  ENDIF

                  IF (PSX /= UNDEFINED) P_STAR(IJK) = PSX
                  IF (TGX /= UNDEFINED) T_G(IJK) = TGX
                  IF (IC_L_SCALE(L) /= UNDEFINED) L_SCALE(IJK) =       &
                      IC_L_SCALE(L)

                  IF (NMAX(0) > 0) THEN
                     WHERE (IC_X_G(L,:NMAX(0)) /= UNDEFINED)           &
                        X_G(IJK,:NMAX(0)) = IC_X_G(L,:NMAX(0))
                  ENDIF

                  IF (NScalar > 0) THEN
                     WHERE (IC_Scalar(L,:NScalar) /= UNDEFINED)        &
                        Scalar(IJK,:NScalar) = IC_Scalar(L,:NScalar)
                  ENDIF

                  IF (K_Epsilon) THEN
                      IF (IC_K_Turb_G(L) /= UNDEFINED)                 &
                          K_Turb_G(IJK) = IC_K_Turb_G(L)
                      IF (IC_E_Turb_G(L) /= UNDEFINED)                 &
                          E_Turb_G(IJK) = IC_E_Turb_G(L)
                  ENDIF

                  IF (UGX /= UNDEFINED) U_G(IJK) = UGX
                  IF (VGX /= UNDEFINED) V_G(IJK) = VGX
                  IF (WGX /= UNDEFINED) W_G(IJK) = WGX

                  GAMA_RG(IJK) = IC_GAMA_RG(L)
                  T_RG(IJK) = merge(IC_T_RG(L), ZERO,                  &
                     IC_T_RG(L) /= UNDEFINED)

                  DO M = 1, MMAX
                    IF (ROPSX(M) /= UNDEFINED) ROP_S(IJK,M) = ROPSX(M)
                    IF (TSX(M) /= UNDEFINED) T_S(IJK,M) = TSX(M)
                    IF (IC_THETA_M(L,M) /= UNDEFINED)                  &
                       THETA_M(IJK,M) = IC_THETA_M(L,M)
                    IF (USX(M) /= UNDEFINED) U_S(IJK,M) = USX(M)
                    IF (VSX(M) /= UNDEFINED) V_S(IJK,M) = VSX(M)
                    IF (WSX(M) /= UNDEFINED) W_S(IJK,M) = WSX(M)

                    GAMA_RS(IJK,M) = IC_GAMA_RS(L,M)
                    T_RS(IJK,M) = merge(IC_T_RS(L,M),ZERO,             &
                       IC_T_RS(L,M) /= UNDEFINED)

                    IF (NMAX(M) > 0) THEN
                        WHERE (IC_X_S(L,M,:NMAX(M)) /= UNDEFINED)      &
                           X_S(IJK,M,:NMAX(M)) = IC_X_S(L,M,:NMAX(M))
                    ENDIF
                  ENDDO

! for GHD theory to compute mixture IC of velocity and density
                  IF(KT_TYPE_ENUM == GHD_2007) THEN
                     ROP_S(IJK,MMAX) = ZERO
                     U_S(IJK,MMAX) = ZERO
                     V_S(IJK,MMAX) = ZERO
                     W_S(IJK,MMAX) = ZERO
                     THETA_M(IJK,MMAX) = ZERO
                     nTOT = ZERO
                     nM = ZERO
                     DO M = 1, SMAX
                        IF (ROPSX(M) /= UNDEFINED) THEN
                           ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + ROPSX(M)
                           nM = ROPSX(M)*6d0 /                         &
                              (PI*D_p(IJK,M)**3*RO_S(IJK,M))
                           nTOT = nTOT + nM
                        ENDIF
                        IF (IC_THETA_M(L,M) /= UNDEFINED)              &
                           THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) +     &

                           nM*IC_THETA_M(L,M)
                        IF(USX(M) /= UNDEFINED .AND.                   &
                           ROPSX(M) /= UNDEFINED)  U_S(IJK,MMAX) =     &
                           U_S(IJK,MMAX) + ROPSX(M)*USX(M)

                        IF(VSX(M) /= UNDEFINED .AND.                   &
                           ROPSX(M) /= UNDEFINED) V_S(IJK,MMAX) =      &
                           V_S(IJK,MMAX) +  ROPSX(M)*VSX(M)

                        IF(WSX(M) /= UNDEFINED .AND.                   &
                           ROPSX(M) /= UNDEFINED) W_S(IJK,MMAX) =      &
                           W_S(IJK,MMAX) +  ROPSX(M)*WSX(M)
                     ENDDO

! If ropsTotal > 0 then RoN_T > 0
                     IF(ROP_S(IJK,MMAX) > ZERO) THEN
                        U_S(IJK,MMAX) = U_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                        V_S(IJK,MMAX) = V_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                        W_S(IJK,MMAX) = W_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                        THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) / nTOT
! For initially empty bed:
                     ELSE
                        U_S(IJK,MMAX) = U_S(IJK,MMAX)
                        V_S(IJK,MMAX) = V_S(IJK,MMAX)
                        W_S(IJK,MMAX) = W_S(IJK,MMAX)
! Set T > 0 in case Ti > 0
                        DO M = 1, SMAX
                           THETA_M(IJK,MMAX) = THETA_M(IJK,M)
                        ENDDO
                        IF(THETA_M(IJK,MMAX)==ZERO)                     &
                           THETA_M(IJK,MMAX) = small_number
                    ENDIF
                  ENDIF
! end of modifications for GHD theory
               ENDIF     ! Fluid at
            ENDDO   ! over i
            ENDDO   ! over j
            ENDDO   ! over k
         ENDIF   ! if (ic_defined)
      ENDDO   ! over dimension_ic

      CALL SEND_RECV(L_SCALE,2)

      RETURN
      END SUBROUTINE SET_IC
