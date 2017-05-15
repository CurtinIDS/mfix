!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: QMOMK_TIME_MARCH                                       C
!     Purpose: Called in time_march.f to do QMOMK calculations            C
!                                                                         C
!     Author: Alberto Passalacqua                        Date:            C
!     Reviewer:                                          Date:            C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
SUBROUTINE QMOMK_TIME_MARCH

  USE param
  USE param1
  USE constant
  USE run
  USE output
  USE physprop
  USE fldvar
  USE geometry
  USE cont
  USE tau_g
  USE tau_s
  USE visc_g
  USE visc_s
  USE funits
  USE vshear
  USE scalars
  USE drag
  USE rxns
  USE compar
  USE time_cpu
  USE is
  USE indices
  USE sendrecv
  USE qmom_kinetic_equation
  USE qmomk_fluxes
  USE qmomk_quadrature
  USE qmomk_collision
  USE qmomk_parameters
  USE drag
  USE ur_facs
  USE fun_avg
  USE functions

  IMPLICIT NONE

  INTEGER :: I,J,K,M,M2,IN
  INTEGER :: IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
  DOUBLE PRECISION :: Vmax, min_space_delta
  DOUBLE PRECISION :: vrel, Rep, CD, beta_drag, min_tau_drag
  DOUBLE PRECISION :: UGC, VGC, WGC, mom0, QMOMK_TCOL, drag_exp, QMOMK_OMEGA
  DOUBLE PRECISION QMOMK_TIME, QMOMK_DT_TMP, TMP_DTS
  INTEGER :: TIME_FACTOR, TIME_I
  LOGICAL, SAVE ::  FIRST_PASS = .TRUE.

  DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: Nlminus, Nlplus, Nrminus, Nrplus
  DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: Ulminus, Ulplus, Urminus, Urplus
  DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: Vlminus, Vlplus, Vrminus, Vrplus
  DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: Wlminus, Wlplus, Wrminus, Wrplus
  DOUBLE PRECISION, DIMENSION(QMOMK_NMOM) :: F_x_left, F_x_right, F_y_left
  DOUBLE PRECISION, DIMENSION(QMOMK_NMOM) :: F_y_right, F_z_left, F_z_right

  IF (FIRST_PASS) THEN
     FIRST_PASS = .FALSE.
     ! Setting initial guess for drag time
     QMOMK_TIME = 0.D0
     QMOMK_TAU_DRAG = 1.0D-5

     ! A.P. For restart runs, IC's are already set and would be overwritten
     IF (RUN_TYPE == 'NEW') THEN
       CALL QMOMK_INITIAL_CONDITIONS
       CALL QMOMK_INIT_BC
     ELSE
       CALL QMOMK_INIT_BC
     ENDIF

     QMOMK_M1 = QMOMK_M0
     QMOMK_N1 = QMOMK_N0
     QMOMK_U1 = QMOMK_U0
     QMOMK_V1 = QMOMK_V0
     QMOMK_W1 = QMOMK_W0
  END IF

  ! Finding initial QMOMK time step
  ! CFL condition
  Vmax = 0.d0
  DO M = 1, MMAX
    DO IJK = ijkstart3, ijkend3
        DO I = 1, QMOMK_NN
          IF (Vmax < MAX(ABS(QMOMK_U0(I,IJK,M)), ABS(QMOMK_V0(I,IJK,M)), ABS(QMOMK_W0(I,IJK,M)))) THEN
              Vmax = MAX(ABS(QMOMK_U0(I,IJK,M)), ABS(QMOMK_V0(I,IJK,M)), ABS(QMOMK_W0(I,IJK,M)))
          END IF
        END DO
        CALL COMPUTE_COLLISION_TIME(QMOMK_M0(:,IJK,M), D_p0(M), THETA_M(IJK,M), QMOMK_COLLISION_TIME(IJK,M), DT)
    END DO
  END DO

  min_space_delta = MIN (MINVAL(DX), MINVAL(DY), MINVAL(DZ))

  QMOMK_DT = QMOMK_CFL*min_space_delta/Vmax

  ! Drag time and collision time
  min_tau_drag = 1.D0
  QMOMK_TCOL = 1.D10
  DO IJK = ijkstart3, ijkend3
    IF (FLUID_AT(IJK)) THEN
        IF (min_tau_drag >  MINVAL(QMOMK_TAU_DRAG(:,IJK,:))) THEN
          min_tau_drag = MINVAL(QMOMK_TAU_DRAG(:,IJK,:))
        END IF
        IF (QMOMK_TCOL > MINVAL(QMOMK_COLLISION_TIME(IJK,:))) THEN
           QMOMK_TCOL = MINVAL(QMOMK_COLLISION_TIME(IJK,:))
        END IF
    END IF
  END DO

  IF ( QMOMK_COLLISIONS /= 'NONE') THEN
     QMOMK_DT = MIN(QMOMK_TCOL,QMOMK_DT)
  END IF

  QMOMK_DT = MIN(QMOMK_DT,min_tau_drag/10.D0)

  ! Preparing for QMOMK internal time stepping, if required
  QMOMK_TIME = TIME
  IF (DT >= QMOMK_DT) THEN
     TIME_FACTOR = CEILING(real(DT/QMOMK_DT)) + 1
     QMOMK_DT_TMP = ZERO
  ELSE
     TIME_FACTOR = 1
     QMOMK_DT_TMP = QMOMK_DT
     QMOMK_DT = DT
  END IF

  TMP_DTS = ZERO

  DO TIME_I = 1, TIME_FACTOR  ! Starting QMOMK internal time stepping

     ! Exit if time is greater than flow time
     IF (QMOMK_TIME .GE. (TIME + DT)) EXIT

     IF ((QMOMK_TIME + QMOMK_DT) .GT. (TIME + DT)) THEN
      TMP_DTS = QMOMK_DT
      QMOMK_DT = TIME + DT - QMOMK_TIME
     END IF

     PRINT *,'Flow time step: ',DT
     PRINT *,'QMOM DT = ',QMOMK_DT
     PRINT *,'QMOMK internal time step: ',TIME_I
     PRINT *,'Time factor', TIME_FACTOR

     QMOMK_OMEGA = (1.D0 + C_e)/2.D0

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! First step of Runge-Kutta time integration !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Calculating kinetic based fluxes (Half time step)
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK)) THEN
              I = I_OF(IJK)
              J = J_OF(IJK)
              K = K_OF(IJK)

              IMJK = IM_OF(IJK)
              IJMK = JM_OF(IJK)
              IJKM = KM_OF(IJK)

              IPJK = IP_OF(IJK)
              IJPK = JP_OF(IJK)
              IJKP = KP_OF(IJK)

              ! 2D Case
              IF (NO_K) THEN
                 ! x - direction
                 Nlminus = QMOMK_N0(:,IMJK,M)
                 Ulminus = QMOMK_U0(:,IMJK,M)
                 Vlminus = QMOMK_V0(:,IMJK,M)
                 Wlminus = QMOMK_W0(:,IMJK,M)

                 Nlplus  = QMOMK_N0(:,IJK,M)
                 Ulplus  = QMOMK_U0(:,IJK,M)
                 Vlplus  = QMOMK_V0(:,IJK,M)
                 Wlplus  = QMOMK_W0(:,IJK,M)

                 Nrminus = QMOMK_N0(:,IJK,M)
                 Urminus = QMOMK_U0(:,IJK,M)
                 Vrminus = QMOMK_V0(:,IJK,M)
                 Wrminus = QMOMK_W0(:,IJK,M)

                 Nrplus  = QMOMK_N0(:,IPJK,M)
                 Urplus  = QMOMK_U0(:,IPJK,M)
                 Vrplus  = QMOMK_V0(:,IPJK,M)
                 Wrplus  = QMOMK_W0(:,IPJK,M)

                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_x_left)
                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_x_right)

                 ! y - direction
                 Nlminus = QMOMK_N0(:,IJMK,M)
                 Ulminus = QMOMK_U0(:,IJMK,M)
                 Vlminus = QMOMK_V0(:,IJMK,M)
                 Wlminus = QMOMK_W0(:,IJMK,M)

                 Nlplus  = QMOMK_N0(:,IJK,M)
                 Ulplus  = QMOMK_U0(:,IJK,M)
                 Vlplus  = QMOMK_V0(:,IJK,M)
                 Wlplus  = QMOMK_W0(:,IJK,M)

                 Nrminus = QMOMK_N0(:,IJK,M)
                 Urminus = QMOMK_U0(:,IJK,M)
                 Vrminus = QMOMK_V0(:,IJK,M)
                 Wrminus = QMOMK_W0(:,IJK,M)

                 Nrplus  = QMOMK_N0(:,IJPK,M)
                 Urplus  = QMOMK_U0(:,IJPK,M)
                 Vrplus  = QMOMK_V0(:,IJPK,M)
                 Wrplus  = QMOMK_W0(:,IJPK,M)

                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_y_left)
                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_y_right)

                 QMOMK_M1(:,IJK,M) = QMOMK_M0(:,IJK,M) &
                      - 0.5*(QMOMK_DT/MINVAL(DX))*(F_x_right - F_x_left) - 0.5*(QMOMK_DT/MINVAL(DY))*(F_y_right - F_y_left)


              ELSE ! 3D case
                 ! x - direction
                 Nlminus = QMOMK_N0(:,IMJK,M)
                 Ulminus = QMOMK_U0(:,IMJK,M)
                 Vlminus = QMOMK_V0(:,IMJK,M)
                 Wlminus = QMOMK_W0(:,IMJK,M)

                 Nlplus  = QMOMK_N0(:,IJK,M)
                 Ulplus  = QMOMK_U0(:,IJK,M)
                 Vlplus  = QMOMK_V0(:,IJK,M)
                 Wlplus  = QMOMK_W0(:,IJK,M)

                 Nrminus = QMOMK_N0(:,IJK,M)
                 Urminus = QMOMK_U0(:,IJK,M)
                 Vrminus = QMOMK_V0(:,IJK,M)
                 Wrminus = QMOMK_W0(:,IJK,M)

                 Nrplus  = QMOMK_N0(:,IPJK,M)
                 Urplus  = QMOMK_U0(:,IPJK,M)
                 Vrplus  = QMOMK_V0(:,IPJK,M)
                 Wrplus  = QMOMK_W0(:,IPJK,M)

                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_x_left)
                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_x_right)

                 ! y - direction
                 Nlminus = QMOMK_N0(:,IJMK,M)
                 Ulminus = QMOMK_U0(:,IJMK,M)
                 Vlminus = QMOMK_V0(:,IJMK,M)
                 Wlminus = QMOMK_W0(:,IJMK,M)

                 Nlplus  = QMOMK_N0(:,IJK,M)
                 Ulplus  = QMOMK_U0(:,IJK,M)
                 Vlplus  = QMOMK_V0(:,IJK,M)
                 Wlplus  = QMOMK_W0(:,IJK,M)

                 Nrminus = QMOMK_N0(:,IJK,M)
                 Urminus = QMOMK_U0(:,IJK,M)
                 Vrminus = QMOMK_V0(:,IJK,M)
                 Wrminus = QMOMK_W0(:,IJK,M)

                 Nrplus  = QMOMK_N0(:,IJPK,M)
                 Urplus  = QMOMK_U0(:,IJPK,M)
                 Vrplus  = QMOMK_V0(:,IJPK,M)
                 Wrplus  = QMOMK_W0(:,IJPK,M)

                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_y_left)
                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_y_right)

                 ! z - direction
                 Nlminus = QMOMK_N0(:,IJKM,M)
                 Ulminus = QMOMK_U0(:,IJKM,M)
                 Vlminus = QMOMK_V0(:,IJKM,M)
                 Wlminus = QMOMK_W0(:,IJKM,M)

                 Nlplus  = QMOMK_N0(:,IJK,M)
                 Ulplus  = QMOMK_U0(:,IJK,M)
                 Vlplus  = QMOMK_V0(:,IJK,M)
                 Wlplus  = QMOMK_W0(:,IJK,M)

                 Nrminus = QMOMK_N0(:,IJK,M)
                 Urminus = QMOMK_U0(:,IJK,M)
                 Vrminus = QMOMK_V0(:,IJK,M)
                 Wrminus = QMOMK_W0(:,IJK,M)

                 Nrplus  = QMOMK_N0(:,IJKP,M)
                 Urplus  = QMOMK_U0(:,IJKP,M)
                 Vrplus  = QMOMK_V0(:,IJKP,M)
                 Wrplus  = QMOMK_W0(:,IJKP,M)

                 CALL KINETIC_FLUX_Z_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_z_left)
                 CALL KINETIC_FLUX_Z_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_z_right)

                 QMOMK_M1(:,IJK,M) = QMOMK_M0(:,IJK,M) - 0.5*(QMOMK_DT/MINVAL(DX))*(F_x_right - F_x_left) - &
                      0.5*(QMOMK_DT/MINVAL(DY))*(F_y_right - F_y_left) - 0.5*(QMOMK_DT/MINVAL(DZ))*(F_z_right - F_z_left)
              END IF
           END IF
        END DO
     END DO
     ! End of fluxes calcuation

     ! Update weights and abscissas
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK))  THEN
              CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                   QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
              CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                   QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
           END IF
        END DO
     END DO

     ! RHS of moments equations
     ! Gravity and drag contribution (Half time step)
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
           I = I_OF(IJK)
           J = J_OF(IJK)
           K = K_OF(IJK)
           IMJK = IM_OF(IJK)
           IJMK = JM_OF(IJK)
           IJKM = KM_OF(IJK)

           UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
           VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
           WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

           DO IN = 1, QMOMK_NN
              vrel = SQRT(((QMOMK_U1(IN,IJK,M) - UGC)**2) + ((QMOMK_V1(IN,IJK,M) - VGC)**2) + ((QMOMK_W1(IN,IJK,M) - WGC)**2))

              ! Particle Reynolds number
              Rep = RO_g0*D_p0(M)*vrel/MU_g0
              Cd = 24.D0*(EP_G(IJK)**(-2.65))*(1.D0 + 0.15D0*((EP_G(IJK)*Rep)**0.687))/(Rep + SMALL_NUMBER)
              beta_drag = 3.D0*RO_g0*Cd*vrel/(4.D0*D_p0(M)*RO_s(IJK,M))

              drag_exp = EXP(-0.5*QMOMK_DT*beta_drag)

              IF (beta_drag > SMALL_NUMBER) THEN   !CHECK HERE, it was /= 0
                QMOMK_TAU_DRAG(IN,IJK,M) = 1.D0/beta_drag
              ELSE
                QMOMK_TAU_DRAG(IN,IJK,M) = LARGE_NUMBER
              END IF

              ! Drag calculation for QMOMK
              QMOMK_U1 (IN,IJK,M) = drag_exp*QMOMK_U1 (IN,IJK,M) + (1.D0 - drag_exp)*UGC
              QMOMK_V1 (IN,IJK,M) = drag_exp*QMOMK_V1 (IN,IJK,M) + (1.D0 - drag_exp)*(VGC - GRAVITY*QMOMK_TAU_DRAG(IN,IJK,M))
              QMOMK_W1 (IN,IJK,M) = drag_exp*QMOMK_W1 (IN,IJK,M) + (1.D0 - drag_exp)*WGC
           END DO
           CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
           CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
           CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
          END IF
        END DO
     END DO

     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
          IF (FLUID_AT(IJK)) THEN
            CALL BIND_THETA(QMOMK_M1(:,IJK,M), MINVAL(QMOMK_TAU_DRAG(:,IJK,M)))
            CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                 QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
           CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
          END IF
        END DO
     END DO

     ! End of gravity and drag contribution

     PRINT *,'Before collisions'

     ! Collisions contribution (Half time step)
     IF (QMOMK_COLLISIONS == 'BGK' .AND. MMAX == 1) THEN
        DO M = 1, MMAX
           DO IJK = ijkstart3, ijkend3
              IF (FLUID_AT(IJK))  THEN
                 mom0 = QMOMK_M1(1,IJK,M)
                 IF (mom0 > epsn) THEN
                    CALL COLLISIONS_BGK(QMOMK_M1(:,IJK,M), 0.5*QMOMK_DT, QMOMK_TCOL, D_p0 (M), C_e)
                    CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                         QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
                    CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                         QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
                 END IF
              END IF
           END DO
        END DO
     ELSE IF (QMOMK_COLLISIONS == 'BOLTZMANN' .AND. MMAX == 1) THEN
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK)) THEN
              CALL SOLVE_BOLTZMANN_COLLISIONS_ONE_SPECIE(QMOMK_M1(:,IJK,1), QMOMK_N1(:,IJK,1), &
                   QMOMK_U1(:,IJK,1), QMOMK_V1(:,IJK,1), QMOMK_W1(:,IJK,1), &
                   0.5*QMOMK_DT, C_e, D_p0 (1), QMOMK_COLLISIONS_ORDER)
              CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,1), QMOMK_N1(:,IJK,1), QMOMK_U1(:,IJK,1), QMOMK_V1(:,IJK,1), QMOMK_W1(:,IJK,1))
              CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,1), &
                   QMOMK_U1(:,IJK,1), QMOMK_V1(:,IJK,1), QMOMK_W1(:,IJK,1), QMOMK_M1(:,IJK,1))
           END IF
        END DO
     ELSE
        DO M = 1, MMAX
           DO IJK = ijkstart3, ijkend3
              IF (FLUID_AT(IJK)) THEN
                 DO M2 = 1, MMAX
                    CALL SOLVE_BOLTZMANN_COLLISIONS_TWO_SPECIES(QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                         QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), &
                         QMOMK_M1(:,IJK,M2), QMOMK_N1(:,IJK,M2), QMOMK_U1(:,IJK,M2), QMOMK_V1(:,IJK,M2), QMOMK_W1(:,IJK,M2), &
                         0.5*QMOMK_DT, 1.D0/6.D0*Pi*(D_p0(M)**3)*RO_s(IJK,M), 1.D0/6.D0*Pi*(D_p0(M2)**3)*RO_s(IJK,M2), &
                         D_p0(M), D_p0(M2), C_e, C_e, QMOMK_COLLISIONS_ORDER)
                 END DO
                 CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
                 CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                      QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
              END IF
           END DO
        END DO
     END IF

     PRINT *,'After collisions'

     ! End of collisions contribution

     ! Boundary conditions update
      CALL QMOMK_SET_BC

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Second step of Runge-Kutta time integration !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! M0 -> old values
     ! M1 -> first step RK values

     ! Calculating kinetic based fluxes (Full time step)
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK)) THEN
              I = I_OF(IJK)
              J = J_OF(IJK)
              K = K_OF(IJK)

              IMJK = IM_OF(IJK)
              IJMK = JM_OF(IJK)
              IJKM = KM_OF(IJK)

              IPJK = IP_OF(IJK)
              IJPK = JP_OF(IJK)
              IJKP = KP_OF(IJK)

              ! 2D Case

              IF (NO_K) THEN
                 ! x - direction
                 Nlminus = QMOMK_N1(:,IMJK,M)
                 Ulminus = QMOMK_U1(:,IMJK,M)
                 Vlminus = QMOMK_V1(:,IMJK,M)
                 Wlminus = QMOMK_W1(:,IMJK,M)

                 Nlplus  = QMOMK_N1(:,IJK,M)
                 Ulplus  = QMOMK_U1(:,IJK,M)
                 Vlplus  = QMOMK_V1(:,IJK,M)
                 Wlplus  = QMOMK_W1(:,IJK,M)

                 Nrminus = QMOMK_N1(:,IJK,M)
                 Urminus = QMOMK_U1(:,IJK,M)
                 Vrminus = QMOMK_V1(:,IJK,M)
                 Wrminus = QMOMK_W1(:,IJK,M)

                 Nrplus  = QMOMK_N1(:,IPJK,M)
                 Urplus  = QMOMK_U1(:,IPJK,M)
                 Vrplus  = QMOMK_V1(:,IPJK,M)
                 Wrplus  = QMOMK_W1(:,IPJK,M)

                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_x_left)
                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_x_right)

                 ! y - direction
                 Nlminus = QMOMK_N1(:,IJMK,M)
                 Ulminus = QMOMK_U1(:,IJMK,M)
                 Vlminus = QMOMK_V1(:,IJMK,M)
                 Wlminus = QMOMK_W1(:,IJMK,M)

                 Nlplus  = QMOMK_N1(:,IJK,M)
                 Ulplus  = QMOMK_U1(:,IJK,M)
                 Vlplus  = QMOMK_V1(:,IJK,M)
                 Wlplus  = QMOMK_W1(:,IJK,M)

                 Nrminus = QMOMK_N1(:,IJK,M)
                 Urminus = QMOMK_U1(:,IJK,M)
                 Vrminus = QMOMK_V1(:,IJK,M)
                 Wrminus = QMOMK_W1(:,IJK,M)

                 Nrplus  = QMOMK_N1(:,IJPK,M)
                 Urplus  = QMOMK_U1(:,IJPK,M)
                 Vrplus  = QMOMK_V1(:,IJPK,M)
                 Wrplus  = QMOMK_W1(:,IJPK,M)

                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_y_left)
                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_y_right)

                 QMOMK_M1(:,IJK,M) = QMOMK_M0(:,IJK,M) - (QMOMK_DT/MINVAL(DX))*(F_x_right - F_x_left) &
                      - (QMOMK_DT/MINVAL(DY))*(F_y_right - F_y_left)

                 ! 3D case
              ELSE
                 ! x - direction
                 Nlminus = QMOMK_N1(:,IMJK,M)
                 Ulminus = QMOMK_U1(:,IMJK,M)
                 Vlminus = QMOMK_V1(:,IMJK,M)
                 Wlminus = QMOMK_W1(:,IMJK,M)

                 Nlplus  = QMOMK_N1(:,IJK,M)
                 Ulplus  = QMOMK_U1(:,IJK,M)
                 Vlplus  = QMOMK_V1(:,IJK,M)
                 Wlplus  = QMOMK_W1(:,IJK,M)

                 Nrminus = QMOMK_N1(:,IJK,M)
                 Urminus = QMOMK_U1(:,IJK,M)
                 Vrminus = QMOMK_V1(:,IJK,M)
                 Wrminus = QMOMK_W1(:,IJK,M)

                 Nrplus  = QMOMK_N1(:,IPJK,M)
                 Urplus  = QMOMK_U1(:,IPJK,M)
                 Vrplus  = QMOMK_V1(:,IPJK,M)
                 Wrplus  = QMOMK_W1(:,IPJK,M)

                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_x_left)
                 CALL KINETIC_FLUX_X_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_x_right)

                 ! y - direction
                 Nlminus = QMOMK_N1(:,IJMK,M)
                 Ulminus = QMOMK_U1(:,IJMK,M)
                 Vlminus = QMOMK_V1(:,IJMK,M)
                 Wlminus = QMOMK_W1(:,IJMK,M)

                 Nlplus  = QMOMK_N1(:,IJK,M)
                 Ulplus  = QMOMK_U1(:,IJK,M)
                 Vlplus  = QMOMK_V1(:,IJK,M)
                 Wlplus  = QMOMK_W1(:,IJK,M)

                 Nrminus = QMOMK_N1(:,IJK,M)
                 Urminus = QMOMK_U1(:,IJK,M)
                 Vrminus = QMOMK_V1(:,IJK,M)
                 Wrminus = QMOMK_W1(:,IJK,M)

                 Nrplus  = QMOMK_N1(:,IJPK,M)
                 Urplus  = QMOMK_U1(:,IJPK,M)
                 Vrplus  = QMOMK_V1(:,IJPK,M)
                 Wrplus  = QMOMK_W1(:,IJPK,M)

                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_y_left)
                 CALL KINETIC_FLUX_Y_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_y_right)

                 ! z - direction
                 Nlminus = QMOMK_N1(:,IJKM,M)
                 Ulminus = QMOMK_U1(:,IJKM,M)
                 Vlminus = QMOMK_V1(:,IJKM,M)
                 Wlminus = QMOMK_W1(:,IJKM,M)

                 Nlplus  = QMOMK_N1(:,IJK,M)
                 Ulplus  = QMOMK_U1(:,IJK,M)
                 Vlplus  = QMOMK_V1(:,IJK,M)
                 Wlplus  = QMOMK_W1(:,IJK,M)

                 Nrminus = QMOMK_N1(:,IJK,M)
                 Urminus = QMOMK_U1(:,IJK,M)
                 Vrminus = QMOMK_V1(:,IJK,M)
                 Wrminus = QMOMK_W1(:,IJK,M)

                 Nrplus  = QMOMK_N1(:,IJKP,M)
                 Urplus  = QMOMK_U1(:,IJKP,M)
                 Vrplus  = QMOMK_V1(:,IJKP,M)
                 Wrplus  = QMOMK_W1(:,IJKP,M)

                 CALL KINETIC_FLUX_Z_TWENTY_EIGHT_NODES (Nlminus, Ulminus, Vlminus, Wlminus, &
                      Nlplus, Ulplus, Vlplus, Wlplus, F_z_left)
                 CALL KINETIC_FLUX_Z_TWENTY_EIGHT_NODES (Nrminus, Urminus, Vrminus, Wrminus, &
                      Nrplus, Urplus, Vrplus, Wrplus, F_z_right)

                 QMOMK_M1(:,IJK,M) = QMOMK_M0(:,IJK,M) - (QMOMK_DT/MINVAL(DX))*(F_x_right - F_x_left) - &
                      (QMOMK_DT/MINVAL(DY))*(F_y_right - F_y_left) - (QMOMK_DT/MINVAL(DZ))*(F_z_right - F_z_left)
              END IF
           END IF
        END DO
     END DO

     ! Update weights and abscissas
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK))  THEN
              CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
              CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                   QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
           END IF
        END DO
     END DO

     ! RHS of moments equations
     ! Collisions contribution (Full time step)
     IF (QMOMK_COLLISIONS == 'BGK' .AND. MMAX == 1) THEN
        DO M = 1, MMAX
           DO IJK = ijkstart3, ijkend3
              IF (FLUID_AT(IJK))  THEN
                 mom0 = QMOMK_M1(1,IJK,M)
                 IF (mom0 > epsn) THEN
                    CALL COLLISIONS_BGK(QMOMK_M1(:,IJK,M), QMOMK_DT, QMOMK_TCOL, D_p0 (M), C_e)
                    CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                         QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
                    CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                         QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
                 END IF
              END IF
           END DO
        END DO
     ELSE IF (QMOMK_COLLISIONS == 'BOLTZMANN' .AND. MMAX == 1) THEN
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK)) THEN
              CALL SOLVE_BOLTZMANN_COLLISIONS_ONE_SPECIE(QMOMK_M1(:,IJK,1), QMOMK_N1(:,IJK,1), &
                   QMOMK_U1(:,IJK,1), QMOMK_V1(:,IJK,1), QMOMK_W1(:,IJK,1), &
                   QMOMK_DT, C_e, D_p0 (1), QMOMK_COLLISIONS_ORDER)
              CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,1), QMOMK_N1(:,IJK,1), QMOMK_U1(:,IJK,1), QMOMK_V1(:,IJK,1), QMOMK_W1(:,IJK,1))
              CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,1), &
                   QMOMK_U1(:,IJK,1), QMOMK_V1(:,IJK,1), QMOMK_W1(:,IJK,1), QMOMK_M1(:,IJK,1))
           END IF
        END DO
     ELSE
        DO M = 1, MMAX
           DO IJK = ijkstart3, ijkend3
              IF (FLUID_AT(IJK)) THEN
                 DO M2 = 1, MMAX
                    CALL SOLVE_BOLTZMANN_COLLISIONS_TWO_SPECIES(QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
                         QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), &
                         QMOMK_M1(:,IJK,M2), QMOMK_N1(:,IJK,M2), &
                         QMOMK_U1(:,IJK,M2), QMOMK_V1(:,IJK,M2), QMOMK_W1(:,IJK,M2), &
                         QMOMK_DT, 1.D0/6.D0*Pi*(D_p0(M)**3)*RO_s(IJK,M), &
                         1.D0/6.D0*Pi*(D_p0(M2)**3)*RO_s(IJK,M2),  D_p0(M), D_p0(M2), &
                         C_e, C_e, QMOMK_COLLISIONS_ORDER)
                 END DO
                 CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
                 CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                      QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
              END IF
           END DO
        END DO
     END IF
     ! End of collisions contribution

     ! Gravity and drag contribution (Full time step)
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK)) THEN
              I = I_OF(IJK)

              IMJK = IM_OF(IJK)
              IJMK = JM_OF(IJK)
              IJKM = KM_OF(IJK)

              UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
              VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
              WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

              DO IN = 1, QMOMK_NN
                 vrel = SQRT(((QMOMK_U1(IN,IJK,M) - UGC)**2) + ((QMOMK_V1(IN,IJK,M) - VGC)**2) + ((QMOMK_W1(IN,IJK,M) - WGC)**2))

                 ! Particle Reynolds number
                 Rep = RO_g0*D_p0(M)*vrel/MU_g0
                 Cd = 24.D0*(EP_G(IJK)**(-2.65))*(1.D0 + 0.15*((EP_G(IJK)*Rep)**0.687))/(Rep+SMALL_NUMBER)
                 beta_drag = 3.D0*RO_g0*Cd*vrel/(4.D0*D_p0(M)*RO_s(IJK,M))

                 IF (beta_drag > SMALL_NUMBER) THEN ! CHECK HERE, it was /= 0.
                    QMOMK_TAU_DRAG(IN,IJK,M) = 1.D0/beta_drag
                 ELSE
                    QMOMK_TAU_DRAG(IN,IJK,M) = LARGE_NUMBER
                 END IF

                 drag_exp = EXP(-QMOMK_DT*beta_drag)

                 ! Drag calculation for QMOMK
                 QMOMK_U1 (IN,IJK,M) = drag_exp*QMOMK_U1 (IN,IJK,M) + (1.D0 - drag_exp)*UGC
                 QMOMK_V1 (IN,IJK,M) = drag_exp*QMOMK_V1 (IN,IJK,M) + (1.D0 - drag_exp)*(VGC - GRAVITY*QMOMK_TAU_DRAG(IN,IJK,M))
                 QMOMK_W1 (IN,IJK,M) = drag_exp*QMOMK_W1 (IN,IJK,M) + (1.D0 - drag_exp)*WGC
              END DO
              CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                   QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
              CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
              CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                   QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
           END IF
        END DO
     END DO
     ! End of gravity and drag contribution

     DO M = 1, MMAX
       DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
           CALL BIND_THETA(QMOMK_M1(:,IJK,M), MINVAL(QMOMK_TAU_DRAG(:,IJK,M)))
           CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
           CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
                QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
         END IF
       END DO
     END DO

     ! Boundary conditions update
     CALL QMOMK_SET_BC

     ! Preparing data for next time step
     QMOMK_M0 = QMOMK_M1
     QMOMK_N0 = QMOMK_N1
     QMOMK_U0 = QMOMK_U1
     QMOMK_V0 = QMOMK_V1
     QMOMK_W0 = QMOMK_W1

     PRINT *,'QMOM DT 5 = ',QMOMK_DT

     QMOMK_TIME = QMOMK_TIME + QMOMK_DT

     PRINT *,'QMOM DT 6 = ',QMOMK_DT
  END DO  ! QMOMK internal time stepping

  IF (DT .LT. QMOMK_DT_TMP) THEN
    QMOMK_DT = QMOMK_DT_TMP
  ENDIF

  IF (TMP_DTS .NE. ZERO) THEN
    QMOMK_DT = TMP_DTS
    TMP_DTS = ZERO
  ENDIF

  PRINT *,'Time = ',TIME
  PRINT *,'Time QMOMK = ',QMOMK_TIME

  ! Resetting gas volume fraction to one, in order to update it
  DO IJK = ijkstart3, ijkend3
     IF (FLUID_AT(IJK)) THEN
        EP_G(IJK) = ONE
     END IF
  END DO

  ! Passing values to the fluid solver and for storage
  DO M = 1, MMAX
     DO IJK = ijkstart3, ijkend3
        IF (FLUID_AT(IJK)) THEN
           ! Mean particle velocities
           U_S(IJK, M) = QMOMK_M1(2,IJK,M)/QMOMK_M1(1,IJK,M)
           V_S(IJK, M) = QMOMK_M1(3,IJK,M)/QMOMK_M1(1,IJK,M)
           W_S(IJK, M) = QMOMK_M1(4,IJK,M)/QMOMK_M1(1,IJK,M)

           ! Gas-phase volume fraction
           EP_G(IJK) = EP_G(IJK) - QMOMK_M1(1,IJK,M)

           ! Particle phases "densities"
           ROP_S(IJK, M) = QMOMK_M1(1,IJK,M)*RO_s(IJK,M)

           ! Particle phases granular temperatures
           THETA_M(IJK,M) = ((QMOMK_M1(5,IJK,M)/QMOMK_M1(1,IJK,M) - &
                            (QMOMK_M1(2,IJK,M)/QMOMK_M1(1,IJK,M))*(QMOMK_M1(2,IJK,M)/QMOMK_M1(1,IJK,M))) + &
                            (QMOMK_M1(8,IJK,M)/QMOMK_M1(1,IJK,M) - &
                               (QMOMK_M1(3,IJK,M)/QMOMK_M1(1,IJK,M))*(QMOMK_M1(3,IJK,M)/QMOMK_M1(1,IJK,M))) + &
                            (QMOMK_M1(10,IJK,M)/QMOMK_M1(1,IJK,M) - &
                               (QMOMK_M1(4,IJK,M)/QMOMK_M1(1,IJK,M))*(QMOMK_M1(4,IJK,M)/QMOMK_M1(1,IJK,M))))/3.D0

        END IF
     END DO
  END DO

  ! Calculating momentum exchange terms
  IF (QMOMK_COUPLED) THEN
     DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
           IF (FLUID_AT(IJK)) THEN
              I = I_OF(IJK)

              IMJK = IM_OF(IJK)
              IJMK = JM_OF(IJK)
              IJKM = KM_OF(IJK)

              UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
              VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
              WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

              DO IN = 1, QMOMK_NN
                 vrel = SQRT(((QMOMK_U1(IN,IJK,M) - UGC)**2) + ((QMOMK_V1(IN,IJK,M) - VGC)**2) + ((QMOMK_W1(IN,IJK,M) - WGC)**2))

                 Rep = RO_g0*D_p0(M)*vrel/MU_g0
                 Cd = 24.D0*(EP_G(IJK)**(-2.65))*(1.D0 + 0.15*((EP_G(IJK)*Rep)**0.687))/(Rep+SMALL_NUMBER)
                 beta_drag = 3.D0*QMOMK_N1(IN,IJK,M)*RO_g0*Cd*vrel/(4.D0*D_p0(M))

                 QMOMK_F_GS(IN,IJK,M) = beta_drag
              END DO
           END IF
        END DO
     END DO
  END IF
END SUBROUTINE QMOMK_TIME_MARCH
