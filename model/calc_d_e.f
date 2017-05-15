!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_e                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Note that the D_E coefficients for phases M>0 are generally not     !
!  used unless the solids phase has close_packed=F, in which case the  !
!  D_E coefficients for that phase are employed in a mixture pressure  !
!  correction equation and for correcting velocities.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D_E(A_M, VXF_GS, VXF_SS, D_E, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Number of solids phases.
      use physprop, only: MMAX
! Flag: Solve the X-Momentum Equations
      use run, only: MOMENTUM_X_EQ
! Flag: Coupled DEM, PIC, or Hybrid simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag: TMF-DEM solids hybrid model
      use discretelement, only: DES_CONTINUUM_HYBRID
! Volume x average at momentum cell center drag for DEM/PIC
      use discretelement, only: VXF_GDS, VXF_SDS
! Number of solids phases.
      use physprop, only: MMAX

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3, DIMENSION_M
! Size of MxM upper triangular matrix.
      use param1, only: DIMENSION_LM

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Septadiagonal matrix A_m.  The center coefficient is negative.
      DOUBLE PRECISION, INTENT(IN):: A_m(DIMENSION_3,-3:3,0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
      DOUBLE PRECISION, INTENT(INOUT) :: d_e(DIMENSION_3, 0:DIMENSION_M)


! Local variables:
!---------------------------------------------------------------------//
! Numerator and demoninator needed for evaluation of pressure correction
! coefficient for GAS phase. A temporary variable is used for local
! manipulations.
      DOUBLE PRECISION, dimension(:,:), allocatable :: AM0
! Flag: One or more solids momentum equations are solved.
      LOGICAL :: ANY_SOLIDS_X_MOMENTUM

!......................................................................!

      allocate(AM0(DIMENSION_3, 0:DIMENSION_M))

! Initialize the error flag.
      IER = 0

! Copy the gas phase A_M coefficients into temp array.
      AM0(:,:) = A_M(:,0,:)

! Add DEM temp A_M so they are included in the pressure correciton eq.
      IF(DES_CONTINUUM_COUPLED) THEN
         AM0(:,0) = AM0(:,0) - VXF_GDS(:)
         IF (DES_CONTINUUM_HYBRID) &
            AM0(:,1:MMAX) = AM0(:,1:MMAX) - VXF_SDS(:,1:MMAX)
      ENDIF

      ANY_SOLIDS_X_MOMENTUM = any(MOMENTUM_X_EQ(1:MMAX))

      IF(MOMENTUM_X_EQ(0)) THEN
         IF(ANY_SOLIDS_X_MOMENTUM) THEN
            CALL CALC_D_E_GAS_AND_SOLIDS(AM0, VXF_GS, VXF_SS, D_E)
         ELSE
            CALL CALC_D_E_GAS_ONLY(AM0, VXF_GS, D_E)
         ENDIF

      ELSEIF(ANY_SOLIDS_X_MOMENTUM) THEN
         CALL CALC_D_E_SOLIDS_ONLY(AM0, VXF_GS, VXF_SS, D_E)
      ENDIF

      deallocate(AM0)

      RETURN
      END SUBROUTINE CALC_D_E




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_D_E_GAS_AND_SOLIDS                                 !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  CASE: The gas phase X-momentum equation is solved and at least one  !
!    solids phase (m>0) X-momentum equation is solved.                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D_E_GAS_AND_SOLIDS(AM0, VXF_GS, VXF_SS, D_E)

! Global Variables:
!---------------------------------------------------------------------//
! Volume fractions of gas and solids phases.
      use fldvar, only: EP_G, EP_s
! Number of solids phases.
      use physprop, only: MMAX
! Flag: Solve the X-Momentum Equations
      use run, only: MOMENTUM_X_EQ
! Flag: Use MODEL_B
      use run, only: MODEL_B
! Pressure scale factor
      use scales, only: P_SCALE
! Flag: Cartesian grid simulation
      use cutcell, only: CARTESIAN_GRID
! Area of cell's XZ face.
      use geometry, only: AYZ
! Function to average across X face.
      use fun_avg, only: AVG_X
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to convert to phases to single array, IJK of cell to east
      use functions, only: FUNLM, EAST_OF
! Flags: Impermeable surface and mass flow at east face
      use functions, only: IP_AT_E, MFLOW_AT_E
! Indices: I index of cell
      use indices, only: I_OF

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3, DIMENSION_M
! Size of MxM upper triangular matrix.
      use param1, only: DIMENSION_LM
! Double precision parameters.
      use param1, only: ZERO, SMALL_NUMBER, ONE

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION, INTENT(IN) :: AM0(DIMENSION_3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
      DOUBLE PRECISION, INTENT(INOUT) :: D_E(DIMENSION_3, 0:DIMENSION_M)

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: I, IJK, IJKE
! Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPSA(DIMENSION_M)
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION :: SUM_VXF_GS
! sum of Solid M - All other Solid drag
      DOUBLE PRECISION :: SUM_VXF_SS(DIMENSION_M)
! sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION :: SUM_VXF_SS_wt_M
! numerator and demoninator needed for evaluation of pressure correction
! coefficient for GAS phase
      DOUBLE PRECISION :: DEN_MGas, NUM_MGas
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M. specifically, for when the
! sum L is over solids phase only
      DOUBLE PRECISION :: NUM_MSol_LSol(DIMENSION_M)
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over solids phases only
      DOUBLE PRECISION :: DEN_MSol_LSol(DIMENSION_M)
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
      DOUBLE PRECISION :: NUM_MSol_LGas
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over gas phase only
      DOUBLE PRECISION :: DEN_MSol_LGas
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp

!......................................................................!

!$omp parallel default(none) &
!$omp          private(ijk,area_face,i,ijke,epga,sum_vxf_gs,epsa,lm,sum_vxf_ss,den_mgas,num_mgas,tmpdp, &
!$omp                  num_msol_lgas,den_msol_lgas,den_msol_lsol,num_msol_lsol,sum_vxf_ss_wt_m,lpl) &
!$omp          shared(ijkstart3,ijkend3,d_e,mmax,cartesian_grid,ayz,i_of,ep_g,vxf_gs,vxf_ss,momentum_x_eq,am0,model_b,p_scale)
!$omp do
      DO IJK = ijkstart3, ijkend3

         IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN
            DO M= 0, MMAX
               D_E(IJK,M) = ZERO
            ENDDO
            CYCLE
         ENDIF

         AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

         SUM_VXF_GS = ZERO
         DO M= 1, MMAX
            EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
! Gas - All Solids VolxDrag summation
            SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
            SUM_VXF_SS(M) = ZERO
            DO L = 1, MMAX
               IF (L .NE. M) THEN
                  LM = FUNLM(L,M)
! Solids M - All other Solids VolxDrag summation
                  SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)
               ENDIF
            ENDDO
         ENDDO

! calculating DEN_MGas and NUM_MGas
         DEN_MGas  = ZERO
         NUM_MGas = ZERO
         DO M= 1, MMAX
            IF (MOMENTUM_X_EQ(M)) THEN
               NUM_MGas = NUM_MGas + (EPSA(M)*VXF_GS(IJK,M)/           &
                 ((-AM0(IJK,M)) + VXF_GS(IJK,M) + SUM_VXF_SS(M) +      &
                 SMALL_NUMBER))
               DEN_MGas = DEN_MGas + (VXF_GS(IJK,M)*                   &
                 ((-AM0(IJK,M))+SUM_VXF_SS(M))/                        &
                 ((-AM0(IJK,M))+VXF_GS(IJK,M)+SUM_VXF_SS(M)+           &
                 SMALL_NUMBER))
            ELSE
               DEN_MGas = DEN_MGas + VXF_GS(IJK,M)
            ENDIF
         ENDDO

! Model B
! -------------------------------------------------------------------
         IF (MODEL_B) THEN
! Linking velocity correction coefficient to pressure - GAS Phase
            TMPdp = -AM0(IJK,0) + DEN_MGas
            IF(abs(TMPdp) > SMALL_NUMBER) THEN
               D_E(IJK,0) = P_SCALE*AREA_FACE/TMPdp
            ELSE
               D_E(IJK,0) = ZERO
            ENDIF
! Linking velocity correction coefficient to pressure - SOLIDs Phase
            DO M = 1, MMAX
               IF(MOMENTUM_X_EQ(M)) THEN
                  TMPdp = -AM0(IJK,M) + VXF_GS(IJK,M)
                  IF(abs(TMPdp) > SMALL_NUMBER) THEN
                     D_E(IJK,M) = D_E(IJK,0)*VXF_GS(IJK,M)/TMPdp
                  ELSE
                     D_E(IJK,M) = ZERO
                  ENDIF
               ELSE
                  D_E(IJK,M) = ZERO
               ENDIF
            ENDDO

! Model A
! -------------------------------------------------------------------
         ELSE

! Linking velocity correction coefficient to pressure - GAS Phase
            TMPdp = -AM0(IJK,0) + DEN_MGas
            IF(abs(TMPdp) >SMALL_NUMBER) THEN
               D_E(IJK,0) = P_SCALE*AREA_FACE*(EPGA+NUM_MGas)/TMPdp
            ELSE
               D_E(IJK,0) = ZERO
            ENDIF

            DO M= 1, MMAX
! calculating NUM_MSol_LGas and DEN_MSol_LGas
               NUM_MSol_LGas = VXF_GS(IJK,M)*EPGA/                     &
                  ((-AM0(IJK,0))+SUM_VXF_GS+SMALL_NUMBER)
               DEN_MSol_LGas = VXF_GS(IJK,M)*(                         &
                  ((-AM0(IJK,0)) + (SUM_VXF_GS - VXF_GS(IJK,M)))/      &
                  ((-AM0(IJK,0)) + SUM_VXF_GS + SMALL_NUMBER) )

! calculating NUM_MSol_LSol and DEN_MSol_LSol
               NUM_MSol_LSol(M) = ZERO
               DEN_MSol_LSol(M)  = ZERO
               DO L = 1, MMAX
                  IF (L /= M) THEN
                     LM = FUNLM(L,M)
                     SUM_VXF_SS_wt_M = ZERO
                     DO Lp = 1, MMAX
                        IF((Lp /= L) .AND. (Lp /= M) ) THEN
                           LpL = FUNLM(Lp,L)
! Solids L - All other Solids VolxDrag but M summation
                           SUM_VXF_SS_wt_M =                           &
                              SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)
                        ENDIF
                     ENDDO

                     IF(MOMENTUM_X_EQ(L)) THEN
                        NUM_MSol_LSol(M) = NUM_MSol_LSol(M) +          &
                           ( VXF_SS(IJK,LM)*EPSA(L)/                   &
                           ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+ &
                           SMALL_NUMBER) )

                        DEN_MSol_LSol(M) = DEN_MSol_LSol(M) +          &
                           VXF_SS(IJK,LM)*(((-AM0(IJK,L))+             &
                           VXF_GS(IJK,L)+SUM_VXF_SS_wt_M)/             &
                           ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+ &
                           SMALL_NUMBER) )
                     ELSE
                        DEN_MSol_LSol(M) =                             &
                           DEN_MSol_LSol(M) + VXF_SS(IJK,LM)
                     ENDIF
                  ENDIF  ! end if (l.ne.m)
               ENDDO  ! end do (l=1,mmax)

! Linking velocity correction coefficient to pressure - SOLIDs Phase
               IF(MOMENTUM_X_EQ(M)) THEN
                  TMPdp = -AM0(IJK,M) + DEN_MSol_LGas + DEN_MSol_LSol(M)
                  IF(abs(TMPdp) > SMALL_NUMBER) THEN
                     D_E(IJK,M) = P_SCALE*AREA_FACE * (EPSA(M) +       &
                        NUM_MSol_LSol(M) + NUM_MSol_LGas)/TMPdp
                  ELSE
                     D_E(IJK,M) = ZERO
                  ENDIF
               ELSE
                  D_E(IJK,M) = ZERO
               ENDIF
            ENDDO   ! end do (m=1,mmax)

         ENDIF    !end if/else branch Model_B/Model_A

      ENDDO  ! end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel

      RETURN
      END SUBROUTINE CALC_D_E_GAS_AND_SOLIDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_D_E_GAS_ONLY                                       !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Notes: The gas phase x-momentum equation is solved; NO solids phase !
!  x-momentum equations are solved.                                    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_D_E_GAS_ONLY(AM0, VXF_GS, D_E)

! Global Variables:
!---------------------------------------------------------------------//
! Volume fraction of gas phase.
      use fldvar, only: EP_G
! Number of solids phases.
      use physprop, only: MMAX
! Flag: Use MODEL_B
      use run, only: MODEL_B
! Pressure scale factor
      use scales, only: P_SCALE
! Flag: Cartesian grid simulation
      use cutcell, only: CARTESIAN_GRID
! Area of cell's YZ face.
      use geometry, only: AYZ
! Volume of V-momentum cell.
      use geometry, only: VOL_U
! Function to average across Y face.
      use fun_avg, only: AVG_X
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flags: Impermeable surface and mass flow at east face, IJK of cell to east
      use functions, only: IP_AT_E, MFLOW_AT_E, EAST_OF
! Indices: J index of cell
      use indices, only: I_OF
! Flag and variables for QMOM implementation.
      USE qmom_kinetic_equation, only: QMOMK, QMOMK_NN
      USE qmom_kinetic_equation, only: QMOMK_F_GS, QMOMK_F_GS

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3, DIMENSION_M
! Double precision parameters.
      use param1, only: ZERO, SMALL_NUMBER, ONE

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION, INTENT(IN) :: AM0(DIMENSION_3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
      DOUBLE PRECISION, INTENT(INOUT) :: D_E(DIMENSION_3, 0:DIMENSION_M)

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: M
      INTEGER :: INN
      INTEGER :: I, IJK, IJKE
! Average solid volume fraction at momentum cell centers
      !DOUBLE PRECISION :: EPSA(DIMENSION_M)
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION :: SUM_VXF_GS
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp

!......................................................................!

!$omp parallel default(none) &
!$omp          private(ijk,area_face,i,ijke,epga,sum_vxf_gs,tmpdp) &
!$omp          shared(ijkstart3,ijkend3,d_e,mmax,cartesian_grid,ayz,i_of,ep_g,vxf_gs,am0,model_b,p_scale,qmomk,vol_u,qmomk_f_gs)
!$omp do
      DO IJK = IJKSTART3, IJKEND3

         IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
            D_E(IJK,0) = ZERO
            CYCLE
         ENDIF

         AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)

         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

         SUM_VXF_GS = ZERO
         IF (.NOT. QMOMK) THEN
            DO M= 1, MMAX
! Gas - All Solids VolxDrag summation
               SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)
            ENDDO
         ELSE
            DO INN = 1, QMOMK_NN
               DO M = 1, MMAX
                  SUM_VXF_GS = SUM_VXF_GS + VOL_U(IJK)* &
                     AVG_X(QMOMK_F_GS(INN,IJK,M),QMOMK_F_GS(INN,IJKE,M),I)
               ENDDO
            ENDDO
         ENDIF

         TMPdp = -AM0(IJK,0) + SUM_VXF_GS
         IF(abs(TMPdp) > SMALL_NUMBER) THEN
            IF (MODEL_B) THEN
               D_E(IJK,0) = P_SCALE*AREA_FACE/TMPdp
            ELSE
               D_E(IJK,0) = P_SCALE*AREA_FACE*EPGA/TMPdp
            ENDIF   !end if/else branch Model_B/Model_A
         ELSE
            D_E(IJK,0) = ZERO
         ENDIF

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!$omp end parallel

      RETURN
      END SUBROUTINE CALC_D_E_GAS_ONLY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_D_E_SOLIDS_ONLY                                    !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: Calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!  Notes: The gas phase momentum equations are NOT solved but at least !
!  one solids phase momentum equation is solved.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D_E_SOLIDS_ONLY(AM0, VXF_GS, VXF_SS, D_E)

! Global Variables:
!---------------------------------------------------------------------//
! Volume fractions of gas and solids phases.
     use fldvar, only: EP_s
! Number of solids phases.
     use physprop, only: MMAX
! Flag: Solve the X-Momentum Equations
     use run, only: MOMENTUM_X_EQ
! Flag: Use MODEL_B
     use run, only: MODEL_B
! Pressure scale factor
     use scales, only: P_SCALE
! Flag: Cartesian grid simulation
      use cutcell, only: CARTESIAN_GRID
! Area of cell's YZ face.
      use geometry, only: AYZ
! Function to average across X face.
      use fun_avg, only: AVG_X
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to convert to phases to single array, IJK of cell to east
      use functions, only: FUNLM, EAST_OF
! Flags: Impermeable surface and mass flow at east face
      use functions, only: IP_AT_E, MFLOW_AT_E
! Indices: I index of cell
      use indices, only: I_OF

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3, DIMENSION_M
! Size of MxM upper triangular matrix.
      use param1, only: DIMENSION_LM
! Double precision parameters.
      use param1, only: ZERO, SMALL_NUMBER, ONE

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Error index
!      INTEGER, INTENT(INOUT) :: IER
! Temporary variable to store A matrix for local manipulation
      DOUBLE PRECISION, INTENT(IN) :: AM0(DIMENSION_3, 0:DIMENSION_M)
! Volume x average at momentum cell centers
      DOUBLE PRECISION, INTENT(IN) :: VxF_gs(DIMENSION_3, DIMENSION_M)
! Volume x average at momentum cell centers Solid-Solid Drag
      DOUBLE PRECISION, INTENT(IN) :: VxF_ss(DIMENSION_3, DIMENSION_LM)
! Coefficients for pressure correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
      DOUBLE PRECISION, INTENT(INOUT) :: D_E(DIMENSION_3, 0:DIMENSION_M)


! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: LM, M, L, LPL, LP
      INTEGER :: I, IJK, IJKE
! Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPSA(DIMENSION_M)
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! sum of Solid M - All other Solid drag
      DOUBLE PRECISION :: SUM_VXF_SS(DIMENSION_M)
! sum of Solid L - All other Solid VolxDrag but M
      DOUBLE PRECISION :: SUM_VXF_SS_wt_M
! component of the total numerator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M. specifically, for when the
! sum L is over solids phase only
      DOUBLE PRECISION :: NUM_MSol_LSol(DIMENSION_M)
! component of the total denominator needed for evaluation of pressure
! correction coefficient for SOLIDS phase M.  specifically, for when the
! sum L is over solids phases only
      DOUBLE PRECISION :: DEN_MSol_LSol(DIMENSION_M)
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp

!......................................................................!

!$omp parallel default(none) &
!$omp          private(ijk,area_face,i,ijke,epsa,lm,sum_vxf_ss,tmpdp,den_msol_lsol,num_msol_lsol,sum_vxf_ss_wt_m,lpl) &
!$omp          shared(ijkstart3,ijkend3,d_e,mmax,cartesian_grid,ayz,i_of,vxf_gs,vxf_ss,momentum_x_eq,am0,model_b,p_scale)
!$omp do
      DO IJK = IJKSTART3, IJKEND3
         IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK) .OR. MODEL_B) THEN
            DO M= 1, MMAX
               D_E(IJK,M) = ZERO
            ENDDO
            CYCLE
         ENDIF

         AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)

         DO M= 1, MMAX
            EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
            SUM_VXF_SS(M) = ZERO
            DO L = 1, MMAX
              IF (L .NE. M) THEN
                 LM = FUNLM(L,M)
! Solids M - All other Solids VolxDrag summation
                 SUM_VXF_SS(M) = SUM_VXF_SS(M) + VXF_SS(IJK,LM)
              ENDIF
           ENDDO
         ENDDO

         DO M= 1, MMAX
! calculating NUM_MSol_LSol and DEN_MSol_LSol
            NUM_MSol_LSol(M) = ZERO
            DEN_MSol_LSol(M)  = ZERO
            DO L = 1, MMAX
               IF (L .NE. M) THEN
                  LM = FUNLM(L,M)
                  SUM_VXF_SS_wt_M = ZERO
                  DO Lp = 1, MMAX
                     IF ( (Lp .NE. L) .AND. (Lp .NE. M) ) THEN
                        LpL = FUNLM(Lp,L)
! Solids L - All other Solids VolxDrag but M summation
                        SUM_VXF_SS_wt_M =                              &
                           SUM_VXF_SS_wt_M + VXF_SS(IJK,LpL)
                     ENDIF
                  ENDDO
                  IF (MOMENTUM_X_EQ(L)) THEN
                     NUM_MSol_LSol(M) = NUM_MSol_LSol(M) +             &
                        (VXF_SS(IJK,LM)*EPSA(L) /                      &
                        ((-AM0(IJK,L))+VXF_GS(IJK,L)+SUM_VXF_SS(L)+    &
                        SMALL_NUMBER))

                     DEN_MSol_LSol(M) = DEN_MSol_LSol(M) +             &
                        VXF_SS(IJK,LM)*(((-AM0(IJK,L))+VXF_GS(IJK,L) + &
                        SUM_VXF_SS_wt_M)/((-AM0(IJK,L))+VXF_GS(IJK,L)+ &
                        SUM_VXF_SS(L)+ SMALL_NUMBER ))
                  ELSE
                     DEN_MSol_LSol(M) = DEN_MSol_LSol(M)+VXF_SS(IJK,LM)
                  ENDIF
               ENDIF   ! end if (L.ne.M)
            ENDDO  ! end do (l=1,mmax)

! Linking velocity correction coefficient to pressure - SOLIDs Phase (Model_A only)
            IF (MOMENTUM_X_EQ(M)) THEN
               TMPdp = -AM0(IJK,M)+VXF_GS(IJK,M)+DEN_MSol_LSol(M)
               IF(abs(TMPdp) > SMALL_NUMBER) THEN
                  D_E(IJK,M) = P_SCALE*AREA_FACE*                      &
                     (EPSA(M) + NUM_MSol_LSol(M))/TMPdp
               ELSE
                  D_E(IJK,M) = ZERO
               ENDIF
            ELSE
               D_E(IJK,M) = ZERO
            ENDIF
         ENDDO  ! end do (m=1,mmax)
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!$omp end parallel

      RETURN
      END SUBROUTINE CALC_D_E_SOLIDS_ONLY
