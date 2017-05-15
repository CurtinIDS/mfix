!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: K_Epsilon_PROP                                          C
!  Purpose: Calculate diffusion coefficeint and sources for K &        C
!           Epsilon equations                                          C
!                                                                      C
!  Author:                                       Date:                 C
!  Modified: S. Benyahia                         Date:May-13-04        C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Wilcox, D.C., Turbulence Modeling for CFD. DCW Industries, Inc.     C
!     La Canada, Ca. 1994.                                             C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE K_Epsilon_PROP()

! Modules
!---------------------------------------------------------------------//
      USE param
      USE param1
      USE parallel
      USE physprop
      USE drag
      USE run
      USE output
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE indices
      USE constant
      Use vshear
      USE turb
      USE toleranc
      USE compar
      USE TAU_G
      USE sendrecv

      USE cutcell
      USE fun_avg
      USE functions

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Loop indices
      INTEGER :: I1, J1, K1
! Solids phase index
      INTEGER :: M
! U_g at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION :: UgE, UgW
! U_g at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION :: UgN, UgS
! U_g at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION :: UgT, UgB
! U_g at the center of a scalar cell (i, j, k)
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION :: UgcC
! Cell center value of U_g
      DOUBLE PRECISION UGC
! V_g at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION :: VgE, VgW
! V_g at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION :: VgN, VgS
! V_g at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION :: VgT, VgB
! Cell center value of V_g
      DOUBLE PRECISION VGC
! W_g at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION :: WgE, WgW
! W_g at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION :: WgN, WgS
! W_g at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION :: WgT, WgB
! W_g at the center of a scalar cell (i, j, k).
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION :: WgcC
! Cell center value of W_g
      DOUBLE PRECISION WGC
! trace_g and eddy viscosity
      DOUBLE PRECISION :: Trace_G, Mu_gas_t
      DOUBLE PRECISION :: C_Eps_Pope, Xsi_Pope
! gradient of velocity
      DOUBLE PRECISION :: DelV_G(3,3)
! rate of strain tensor
      DOUBLE PRECISION :: D_g(3,3)

!  Production of Turb. Due to shear, Turb Visc, and Constants.
!  See Wilcox PP. 89
      DOUBLE PRECISION Tauij_gDUi_gODxj, C_MU, Sigma_k, Sigma_E, Kappa
      DOUBLE PRECISION Pos_Tauij_gDUi_gODxj, Neg_Tauij_gDUi_gODxj
      DOUBLE PRECISION Ceps_1, Ceps_2, C_Eps_3, Check_Log
      DOUBLE PRECISION Pos_PI_kq_2, Neg_PI_kq_2
! Modif. for Sof Local Var.
      INTEGER :: P,Q
!---------------------------------------------------------------------//

      IF( .NOT. K_Epsilon) RETURN

! M should be forced to be equal to one to get some info.
! from solids phase.
      M = 1

! Add constants. Most of these constants have the same names and values
! as the ones defined in Wilcox book (turbulence modeling for CFD).
! Some are necessary only for Simonin turbulence model

      C_MU = 9.0D-02
      Kappa = 0.42D+0
      Sigma_k = 1.0D0
      Sigma_E = 1.3D0
      Ceps_1 = 1.44D0 !should be 1.6 for axisymmetric cases with no Pope Correction.
      Ceps_2 = 1.92D0
      C_Eps_3 = 1.2D0 ! for Simonin model
      C_Eps_Pope = 0.79D0 ! for Pope Correction.

!!!$omp  parallel do private(ijk, L)
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

! Velocity components at the scalar faces
            CALL GET_FACE_VEL_GAS(IJK, uge, ugw, ugn, ugs, ugt, ugb, ugcc, &
                                       vge, vgw, vgn, vgs, vgt, vgb, &
                                       wge, wgw, wgn, wgs, wgt, wgb, wgcc)

! Velocity derivatives (gradient and rate of strain tensor)
            CALL CALC_DERIV_VEL_GAS(IJK, DelV_g, D_g)

! Value of trace is also available via trd_g but the latter is
! calculated only once per time step
            Trace_g = D_G(1,1) + D_G(2,2) + D_G(3,3)

! Pope's correction in 2-D to the Epsilon Equation for a round-Jet
! from Wilcox book, Page 103.
            Xsi_Pope = ZERO
            DO I1 = 1,3
               DO J1 = 1,3
                  DO K1 = 1,3
                     Xsi_Pope = Xsi_Pope + (DelV_g(I1,J1) - DelV_g(J1,I1))* &
                                (DelV_g(J1,K1) - DelV_g(K1,J1))*&
                                (DelV_g(K1,I1) + DelV_g(I1,K1))
                  ENDDO
               ENDDO
            ENDDO
            Xsi_Pope = Xsi_Pope/6. * K_Turb_G(IJK)**2/(E_Turb_G(IJK)+Small_number)

! This IF statment is to ensure that we are using the turbulent viscosity
! and NOT the effective viscosity.
            IF (MU_GT(IJK) .GE. MU_g(IJK)) THEN
               Mu_gas_t = MU_GT(IJK) - MU_g(IJK)
            ELSE
               Mu_gas_t = ZERO
            ENDIF

! Calculate Tau(i,j)*dUi/dXj (production term in the K Equation
            IF(.NOT.CUT_CELL_AT(IJK)) THEN
               Tauij_gDUi_gODxj = 2D0*Mu_gas_t*(             &
                  D_G(1,1) * D_G(1,1) + &
                  D_G(1,2) * (UgN-UgS)*ODY(J) + &
                  D_G(1,3) * ((UgT-UgB)*(OX(I)*ODZ(K))-WgcC*OX(I)) + &
                  D_G(2,1) * (VgE-VgW)*ODX(I) + &
                  D_G(2,2) * D_G(2,2) + &
                  D_G(2,3) * (VgT-VgB)*(OX(I)*ODZ(K)) + &
                  D_G(3,1) * (WgE-WgW)*ODX(I) + &
                  D_G(3,2) * (WgN-WgS)*ODY(J) + &
                  D_G(3,3) * D_G(3,3)) - &
                  F2O3 * RO_G(IJK) * K_Turb_G(IJK)*Trace_g - &
                  F2O3 * Mu_gas_t * Trace_g**2

            ELSE  ! CUT_CELL
! This is actually not used because of wall functions in cut cells
!               Tauij_gDUi_gODxj = 2D0*Mu_gas_t*(                             &
!                  D_G(1,1) * UG(1,1)  +                      &
!                  D_G(1,2) * UG(1,2)  +                      &
!                  D_G(1,3) * UG(1,3)  +                      &
!                  D_G(2,1) * UG(2,1)  +                      &
!                  D_G(2,2) * UG(2,2)  +                      &
!                  D_G(2,3) * UG(2,3)  +                      &
!                  D_G(3,1) * UG(3,1)  +                      &
!                  D_G(3,2) * UG(3,2)  +                      &
!                  D_G(3,3) * UG(3,3)) -                      &
!                  F2O3 * RO_G(IJK) * K_Turb_G(IJK)*Trace_g-  &
!                  F2O3 * Mu_gas_t * Trace_g**2
            ENDIF


! To avoid very small negative numbers
            IF (Tauij_gDUi_gODxj .GE. ZERO) THEN
               Pos_Tauij_gDUi_gODxj = Tauij_gDUi_gODxj
               Neg_Tauij_gDUi_gODxj = ZERO
            ELSE
               Pos_Tauij_gDUi_gODxj = ZERO
               Neg_Tauij_gDUi_gODxj = Tauij_gDUi_gODxj
            ENDIF


! Interaction terms in the K-Epsilon equations FOR USE Simonin and Ahmadi models
            IF(SIMONIN) THEN
               Pos_PI_kq_2 = F_GS(IJK,1)*K_12(IJK)
               Neg_PI_kq_2 = F_GS(IJK,1)*2.0D0* K_Turb_G(IJK)
            ELSE IF(AHMADI) THEN
               Pos_PI_kq_2 = F_GS(IJK,1)*3.0D0*Theta_m(IJK,M)
               Neg_PI_kq_2 = F_GS(IJK,1)*2.0D0* K_Turb_G(IJK)
               C_Eps_3 = zero ! no extra terms in epsilon equation for Ahmadi model !
            ELSE
               Pos_PI_kq_2 = ZERO
               Neg_PI_kq_2 = ZERO
            ENDIF


            IF(K_Turb_G(IJK) > Small_number .AND. &
               E_Turb_G(IJK) > Small_number) THEN

! Start Adding source terms to the K equation
               K_Turb_G_c (IJK) = (EP_g(IJK)*Pos_Tauij_gDUi_gODxj +  &
                                   Pos_PI_kq_2 )
               K_Turb_G_p (IJK) =(-EP_g(IJK)*Neg_Tauij_gDUi_gODxj +  &
                  EP_g(IJK)*RO_G(IJK)*E_Turb_G(IJK)+Neg_PI_kq_2)/ &
                  K_Turb_G(IJK)


! Calculate velocity components at i, j, k to be used in wall functions
            IMJK = IM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK)
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

! Implementing wall functions targeted to fluid cells next to walls...
! Setting Source and sink terms in the K equation since the production
! in the K eq. due to shear should include the LOG law of the wall.
! North/South wall
               IF(WALL_AT(IJPK).OR.WALL_AT(IJMK)) THEN
                  Check_Log = LOG(9.81*C_mu**0.25*                     &
                     RO_G(IJK)*SQRT(K_Turb_G(IJK))*DY(J)/2.0D+0/       &
                     Mu_g(IJK))
                  IF(Check_Log .LE. ZERO) THEN
                     K_Turb_G_c (IJK) = zero
                     K_Turb_G_p (IJK) = zero
                  ELSE
                     K_Turb_G_c (IJK) = SQRT(C_mu) * 2.D+0/DY(J)*      &
                        MAX(ABS(UGC),ABS(WGC))*EP_g(IJK)*RO_G(IJK)*    &
                        K_Turb_G(IJK)/Check_Log
                     K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)*         &
                        C_mu**0.75*K_Turb_G(IJK)**1.5/DY(J)*           &
                        2.0D+0/Kappa))/K_Turb_G(IJK)
                  ENDIF !for check_log less than zero
! Top/Bottom wall
               ELSEIF(WALL_AT(IJKP).OR.WALL_AT(IJKM)) THEN
                  Check_Log = LOG(9.81*C_mu**0.25*                     &
                     RO_G(IJK)*SQRT(K_Turb_G(IJK))/                    &
                     (ODZ(K)*OX(I)*2.0D+0)/Mu_g(IJK))
                  IF(Check_Log .LE. ZERO) THEN
                     K_Turb_G_c (IJK) = zero
                     K_Turb_G_p (IJK) = zero
                  ELSE
                     K_Turb_G_c (IJK) = SQRT(C_mu)*(ODZ(K)*OX(I)*      &
                        2.0D+0)* MAX(ABS(UGC),ABS(VGC))*EP_g(IJK)*     &
                        RO_G(IJK)*K_Turb_G(IJK)/Check_Log
                     K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)*         &
                        C_mu**0.75*K_Turb_G(IJK)**1.5/Kappa*           &
                        (ODZ(K)*OX(I)*2.0D+0)))/K_Turb_G(IJK)
                  ENDIF !for check_log less than zero
               ENDIF !For walls

! For Cylindrical cases, wall_at (IP) is a wall cell, but wall_at (IM) is
! the axis of symmetry where wall functions obviously don't apply.
               IF(CYLINDRICAL) THEN
                  IF (WALL_AT(IPJK))  THEN
                     Check_Log = LOG(9.81*C_mu**0.25* RO_G(IJK)*       &
                     SQRT(K_Turb_G(IJK))*DX(I)/2.0D+0/Mu_g(IJK))
                     IF(Check_Log .LE. ZERO) THEN
                        K_Turb_G_c (IJK) = zero
                        K_Turb_G_p (IJK) = zero
                     ELSE
                        K_Turb_G_c (IJK) = SQRT(C_mu)*2.D+0/DX(I)*     &
                           MAX(ABS(VGC),ABS(WGC)) *EP_g(IJK)*RO_G(IJK)*&
                           K_Turb_G(IJK)/Check_Log
                        K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)*      &
                           C_mu**0.75*K_Turb_G(IJK)**1.5/DX(I)*        &
                           2.0D+0/Kappa))/ K_Turb_G(IJK)
                     ENDIF !for check_log less than zero
                  ENDIF  ! for wall cells in I direction
! East/West wall
               ELSEIF (WALL_AT(IPJK).OR.WALL_AT(IMJK)) THEN
                  Check_Log = LOG(9.81*C_mu**0.25*RO_G(IJK)*           &
                     SQRT(K_Turb_G(IJK))*DX(I)/2.0D+0/Mu_g(IJK))
                  IF(Check_Log .LE. ZERO) THEN
                     K_Turb_G_c (IJK) = zero
                     K_Turb_G_p (IJK) = zero
                  ELSE
                     K_Turb_G_c (IJK) = SQRT(C_mu)*2.D+0/DX(I)*        &
                        MAX(ABS(VGC),ABS(WGC))*EP_g(IJK)*RO_G(IJK)*    &
                        K_Turb_G(IJK)/Check_Log
                     K_Turb_G_p (IJK) = ((EP_g(IJK)*RO_G(IJK)*         &
                        C_mu**0.75*K_Turb_G(IJK)**1.5/DX(I)*           &
                        2.0D+0/Kappa))/K_Turb_G(IJK)
                  ENDIF !for check_log less than zero
               ENDIF ! for cylindrical

               IF(CUT_CELL_AT(IJK)) THEN
                  Check_Log = LOG(9.81*C_mu**0.25* RO_G(IJK)*          &
                  SQRT(K_Turb_G(IJK))*DELH_Scalar(IJK)/Mu_g(IJK))
                  IF(Check_Log .LE. ZERO) THEN
                     K_Turb_G_c (IJK) = zero
                     K_Turb_G_p (IJK) = zero
                  ELSE
!                     K_Turb_G_c (IJK) = SQRT(C_mu)/DELH_Scalar(IJK)*   &
!                        MAX(ABS(UGC),ABS(VGC),ABS(WGC))*               &
!                        EP_g(IJK)*RO_G(IJK)*K_Turb_G(IJK) /Check_Log
                     K_Turb_G_c (IJK) = SQRT(C_mu)/DELH_Scalar(IJK)*   &
                        DSQRT(UGC**2 + VGC**2 + WGC**2) *              &
                        EP_g(IJK)*RO_G(IJK)*K_Turb_G(IJK)/Check_Log

                     K_Turb_G_p (IJK) = (EP_g(IJK)*RO_G(IJK)*          &
                        C_mu**0.75*K_Turb_G(IJK)**1.5)/                &
                       (DELH_Scalar(IJK)*Kappa)/K_Turb_G(IJK)
                  ENDIF !for check_log less than zero
               ENDIF


! Diffusion coefficient for turbulent kinetic energy (K)
               Dif_K_Turb_G(IJK) = EP_g(IJK)* &
                  (MU_G(IJK) + Mu_gas_t /Sigma_k)

! Add here Dissipation of Turbulence source terms
               E_Turb_G_c (IJK) = (Ceps_1 *&
                  EP_g(IJK)*Pos_Tauij_gDUi_gODxj*                      &
                  E_Turb_G(IJK)/K_Turb_G(IJK) +                        &
                  C_Eps_3*Pos_PI_kq_2*E_Turb_G(IJK)/K_Turb_G(IJK))

! Pope Correction in Xsi_Pope, Add it to E_Turb_G_c to use this option.
                  ! + C_Eps_Pope*RO_G(IJK)*EP_g(IJK)*&
                   !  ZMAX(Xsi_Pope)

               E_Turb_G_p (IJK) = -Ceps_1 * &
                  EP_g(IJK)*Neg_Tauij_gDUi_gODxj /K_Turb_G(IJK) +      &
                  Ceps_2 * EP_g(IJK) *RO_G(IJK)*                       &
                  E_Turb_G(IJK)/K_Turb_G(IJK) +                        &
                  C_Eps_3*(Neg_PI_kq_2)/K_Turb_G(IJK)

! Diffusion coefficient for Dissipation of turbulent energy (Epsilon)
               Dif_E_Turb_G(IJK) =EP_g(IJK)* &
                  (MU_G(IJK) + Mu_gas_t /Sigma_E)

            ELSE
               K_Turb_G_c (IJK) = zero
               K_Turb_G_p (IJK) = zero
               E_Turb_G_c (IJK) = zero
               E_Turb_G_p (IJK) = zero
               Dif_K_Turb_G(IJK) = zero
               Dif_E_Turb_G(IJK) = zero
            ENDIF !for K_turb_G and E_Turb_G having very small numbers
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE K_Epsilon_PROP
