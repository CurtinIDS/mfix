!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Driver routine to call checks for WALL BCs.                 !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS(M_TOT, SKIP, BCV)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Identifies solids model (TFM,DEM,PIC)
      use run, only: SOLIDS_MODEL
! User-input: solids kinetic-theory model.
      use run, only: KT_TYPE_ENUM, GHD_2007

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC being checked.
      INTEGER, INTENT(in) :: BCV
! Total number of solids phases.
      INTEGER, INTENT(in) :: M_TOT
! Flag. Solids not present at this BC (used for flow BCs).
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local Variables:
!---------------------------------------------------------------------//
! Loop/counter variable.
      INTEGER :: M
! Local total number of solids phases
      INTEGER :: MTOT_L
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS")

! Input checks for gas phase.
      CALL CHECK_BC_WALLS_GAS(BCV)

      MTOT_L = merge( M_TOT+1, M_TOT, KT_TYPE_ENUM == GHD_2007)

! Input checks for solid phases.
      DO M=1, MTOT_L
         SELECT CASE(SOLIDS_MODEL(M))
         CASE ('TFM'); CALL CHECK_BC_WALLS_TFM(BCV, M)
         CASE ('DEM'); CALL CHECK_BC_WALLS_DISCRETE(BCV, M)
         CASE ('PIC'); CALL CHECK_BC_WALLS_DISCRETE(BCV, M)
         END SELECT
      ENDDO

! Input checks for user-defined scalar equations.
      CALL CHECK_BC_WALLS_SCALAR_EQ(BCV)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_BC_WALLS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_GAS                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for gas phase WALL BC parameters.          !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_GAS(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! User-input: type of BC
      use bc
! User-Input: gas velocity at wall BCs.
      use bc, only: BC_UW_G, BC_VW_G, BC_WW_G
! User-Input: gas energy eq BCs.
      use bc, only: BC_HW_T_G, BC_TW_G, BC_C_T_G
! User-Input: gas species eq BCs.
      use bc, only: BC_HW_X_G, BC_XW_G, BC_C_X_G
! Total number of speices in each phase.
      use physprop, only: NMAX
! Flag: Solve energy equations.
      use run, only: ENERGY_EQ
! Flag: Solve species equations.
      use run, only: SPECIES_EQ
! Flag: Solve K-th direction (3D)
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants.
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV

      INTEGER :: N
!......................................................................!


! Initialize the error manger.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_GAS")

! The wall velocities are not needed for no-slip or free-slip
      IF(BC_TYPE_ENUM(BCV) == PAR_SLIP_WALL) THEN
         IF(BC_UW_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Uw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_VW_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Vw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_WW_G(BCV) == UNDEFINED) THEN
            IF(DO_K)THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Ww_g',BCV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSE
               BC_WW_G(BCV) = ZERO
            ENDIF
         ENDIF
      ENDIF

! Check energy equation input.
      IF(ENERGY_EQ) THEN
         IF(BC_HW_T_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('BC_HW_T_g',BCV)),           &
               trim(iVal(BC_HW_T_G(BCV)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_T_G(BCV)/=ZERO .AND.                                 &
            BC_TW_G(BCV)==UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Tw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_T_G(BCV)/=UNDEFINED .AND.                            &
            BC_C_T_G(BCV)==UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_C_T_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF


! Check species equation input.
      IF(SPECIES_EQ(0)) THEN
         DO N=1, NMAX(0)
            IF(BC_HW_X_G(BCV,N) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('BC_HW_X_g',BCV,N)),      &
                  trim(iVal(BC_HW_X_G(BCV,N)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_X_G(BCV,N)/=ZERO .AND.                            &
               BC_XW_G(BCV,N)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Xw_g',BCV,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_X_G(BCV,N)/=UNDEFINED .AND.                       &
               BC_C_X_G(BCV,N)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_C_X_g',BCV,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF


! Clear the error manager.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_GAS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_TFM                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for TFM solids WALL BC parameters.         !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_TFM(BCV,M)


! Global Variables:
!---------------------------------------------------------------------//
! User-input: type of BC
      use bc
! User-Input: solids velocity at wall BCs.
      use bc, only: BC_UW_s, BC_VW_s, BC_WW_s
! User-Input: solids energy eq BCs.
      use bc, only: BC_HW_T_s, BC_TW_s, BC_C_T_s
! User-Input: solids species eq BCs.
      use bc, only: BC_HW_X_s, BC_XW_s, BC_C_X_s
! User-Input: granular energy eq BCs.
      use bc, only: BC_HW_THETA_M, BC_ThetaW_M, BC_C_Theta_M
! Total number of solids phases
      use physprop, only: MMAX
! Total number of speices in each phase.
      use physprop, only: NMAX
! Flag: Solve energy equations.
      use run, only: ENERGY_EQ
! Flag: Solve species equations.
      use run, only: SPECIES_EQ
! Flag: Solve Granular energy PDE
      use run, only: GRANULAR_ENERGY
! User-input: solids kinetic-theory model.
      use run, only: KT_TYPE_ENUM, GHD_2007
! Flag: Use johnson and jackson bc
      use bc, only: BC_JJ_PS
! Flag: use revised phip for JJ BC.
      use run, only: bc_jj_m, phip_out_jj
! Flag: use jenkins small friction bc
      use run, only: jenkins
! User input: particle wall restitution coefficient and
! angle of internal friction at walls (degrees)
      use constant, only: e_w, phi_w
! Used by jenkins or revised phip bcs
      use constant, only: tan_phi_w
! Used by revised phip for jj bc
      use constant, only: k4phi, phip0
! User-Input: number of reactionrates
!      use rxns, only: nrr
! Flag: Solve K-th direction (3D)
      use geometry, only: DO_K
! Flag: use cartesian grid
      use cutcell, only: cartesian_grid

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants.
      use param1, only: ONE, ZERO, UNDEFINED, UNDEFINED_I
! parameter values
      use constant, only: pi

      use rxns
! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager
      use funits, only: unit_log
      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC being checked.
      INTEGER, INTENT(in) :: BCV
! Index of solids phase.
      INTEGER, INTENT(in) :: M

! Local Variables:
!---------------------------------------------------------------------//
! Loop/variable counter.
      INTEGER :: NN
! Flag to check momentum eq input.
      LOGICAL :: CHECK_MOMENTUM
! Flag to check scalar eq input.
      LOGICAL :: CHECK_SCALARS
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_TFM")

! Toggle the momentum and scalar input variable checks.
      SELECT CASE(KT_TYPE_ENUM)
      CASE (GHD_2007)
         CHECK_MOMENTUM = (M == MMAX)
         CHECK_SCALARS  = (M /= MMAX)
      CASE DEFAULT
         CHECK_MOMENTUM = .TRUE.
         CHECK_SCALARS  = .TRUE.
      END SELECT

! Set the default specification of Johnson-Jackson BC
      IF(BC_JJ_PS(BCV) == UNDEFINED_I)                                 &
         BC_JJ_PS(BCV) = merge(1,0,GRANULAR_ENERGY)

! specifying bc_jj_ps=1 without granular_energy would cause problem in
! the momentum bc routines
      IF(.NOT.GRANULAR_ENERGY .AND. BC_JJ_PS(BCV) == 1) THEN
         WRITE(ERR_MSG, 1101)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1101 FORMAT('Error 1101: Invoking BC_JJ_PS requires GRANULAR_ENERGY', &
         '=.TRUE.',/ 'Please correct the mfix.dat file.')

! The wall velocities are not needed for no-slip or free-slip
! Wall velocities are needed if johnson-jackson bc model is used
      IF(CHECK_MOMENTUM) THEN
         IF(BC_TYPE_ENUM(BCV) == PAR_SLIP_WALL .OR. &
            BC_JJ_PS(BCV) /= ZERO) THEN
            IF(BC_UW_S(BCV,M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Uw_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(BC_VW_S(BCV,M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Vw_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(BC_WW_S(BCV,M) == UNDEFINED) THEN
               IF(DO_K)THEN
                  WRITE(ERR_MSG,1000) trim(iVar('BC_Ww_s',BCV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ELSE
                  BC_WW_S(BCV,M) = ZERO
               ENDIF
            ENDIF
         ENDIF

         IF(GRANULAR_ENERGY .AND. BC_JJ_PS(BCV)==0) THEN
            IF(BC_HW_THETA_M(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('BC_HW_Theta_M',BCV,M)),  &
                  trim(iVal(BC_HW_Theta_M(BCV,M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_THETA_M(BCV,M)/=ZERO .AND.                        &
               BC_THETAW_M(BCV,M)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_ThetaW_M',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_THETA_M(BCV,M)/=UNDEFINED .AND.                   &
               BC_C_THETA_M(BCV,M)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_C_THETA_M',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ELSE
      ENDIF

      IF(CHECK_SCALARS)THEN
         IF(ENERGY_EQ) THEN
            IF(BC_HW_T_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('BC_HW_T_s',BCV,M)),      &
                  trim(iVal(BC_HW_T_S(BCV,M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_T_S(BCV,M)/=ZERO .AND.                            &
               BC_TW_S(BCV,M)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Tw_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_T_S(BCV,M)/=UNDEFINED .AND.                       &
               BC_C_T_S(BCV,M)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_C_T_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(SPECIES_EQ(M)) THEN
            DO NN=1, NMAX(M)
               IF(BC_HW_X_S(BCV,M,NN) < ZERO) THEN
                  WRITE(ERR_MSG,1001) trim(iVar('BC_HW_X_s',BCV,M,NN)), &
                     trim(iVal(BC_HW_X_S(BCV,M,NN)))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               IF(BC_HW_X_S(BCV,M,NN)/=ZERO .AND.                       &
                  BC_XW_S(BCV,M,NN)==UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('BC_Xw_s',BCV,M,NN))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               IF(BC_HW_X_S(BCV,M,NN)/=UNDEFINED .AND.                  &
                  BC_C_X_S(BCV,M,NN)==UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('BC_C_X_s',BCV,M,NN))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO
         ENDIF ! Species Equation
      ELSE
      ENDIF ! Check Scalars


! might make sense to move these checks out of this routine
! but placed here for now
      IF(GRANULAR_ENERGY .AND. BC_JJ_PS(BCV) == 1) THEN

         IF(KT_TYPE_ENUM == GHD_2007) THEN
            WRITE(ERR_MSG, 1201)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(CARTESIAN_GRID) THEN
! the user should really be warned this is not implemented as
! opposed to running with it
            WRITE(ERR_MSG, 1202)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

! small frictional boundary condition model
         IF(JENKINS) THEN
            IF (BC_JJ_M) THEN
               WRITE(ERR_MSG, 1203)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF (PHI_W == UNDEFINED) THEN
               WRITE(ERR_MSG, 1204)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF (E_W > ONE .OR. E_W < ZERO) THEN
               WRITE(ERR_MSG, 1001) 'E_W', E_W
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
! PHI_W is given in degrees but calculated in radian within
! the fortran codes
            TAN_PHI_W = TAN(PHI_W*PI/180.D0)
         ENDIF

! k4phi, phip0 for variable specularity coefficient
         k4phi = undefined
         IF(BC_JJ_M) THEN
            IF (JENKINS) THEN
               WRITE(ERR_MSG, 1203)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF (PHI_W == UNDEFINED) THEN
               WRITE(ERR_MSG, 1204)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF (E_W > ONE .OR. E_W < ZERO) THEN
               WRITE(ERR_MSG, 1001) 'E_W', E_W
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
! PHI_W is given in degrees but calculated in radian within
! the fortran codes
            TAN_PHI_W = TAN(PHI_W*PI/180.D0)

            k4phi = 7.d0/2.d0*tan_phi_w*(1.d0+e_w)
            IF (phip0 .eq. UNDEFINED) THEN
               phip0 = -0.0012596340709032689 + &
                        0.10645510095633175*k4phi - &
                        0.04281476447854031*k4phi**2 + &
                        0.009759402181229842*k4phi**3 - &
                        0.0012508257938705263*k4phi**4 + &
                        0.00008369829630479206*k4phi**5 - &
                        0.000002269550565981776*k4phi**6
! if k4phi is less than 0.2, the analytical expression for phi is used
! to estimate the phi at r->0
               IF (k4phi .le. 0.2d0) THEN
                  phip0=0.09094568176225006*k4phi
               ENDIF
               WRITE (UNIT_LOG, 1207) phip0
            ENDIF
            IF (phip0 < 0) THEN
               WRITE(ERR_MSG, 1208)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

            IF (PHIP_OUT_JJ) THEN
               IF(nRR < 1) THEN
                  WRITE(ERR_MSG, 1209)
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               WRITE (UNIT_LOG, 1210) phip0
            ENDIF
         ENDIF
      ENDIF   ! if granular_energy and bc_jj_ps = 1


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1201 FORMAT('Error 1201: KT_TYPE = "GHD" cannot be used with ',&
         ' BC_JJ_PS',/'Please correct the mfix.dat file.')

 1202 FORMAT('Error 1202: CARTESIAN_GRID cannot be used with ',&
         ' BC_JJ_PS',/'Please correct the mfix.dat file.')

 1203 FORMAT('Error 1203: JENKINS and BC_JJ_M cannot be used at the',&
         ' same time',/'Please correct the mfix.dat file.')
 1204 FORMAT('Error 1204: Angle of particle-wall friction (PHI_W) not',&
         ' specified.',/'Please correct the mfix.dat file.')

 1208 FORMAT('Error 1208: phip0 less than zero.')
 1209 FORMAT('Error 1209: nRR should be at least 1 for storing ',&
             'specularity.')

 1207 FORMAT(/1X,70('*')//' From: CHECK_BC_WALLS_TFM',/' Message: ',&
         'No input for phip0 available, working expression is used.',/ &
         'phip0=',G12.5,/1X,70('*')/)
 1210 FORMAT(/1X,70('*')//' From: CHECK_BC_WALLS_TFM',/' Message: ',&
         'Specularity will be stored as the first element of ', &
         'ReactionRates',/1X,'array in WALL CELLS. Please avoid ', &
         'overwriting it when reacting flow',/1X,' is simulated.', &
         /1X,70('*')/)

      END SUBROUTINE CHECK_BC_WALLS_TFM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_DISCRETE                                  !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for DEM/PIC solids WALL BC parameters.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_DISCRETE(BCV,M)


! Global Variables:
!---------------------------------------------------------------------//
! User-Input: solids velocity at wall BCs.
      use bc, only: BC_UW_s, BC_VW_s, BC_WW_s
! User-Input: solids energy eq BCs.
      use bc, only: BC_HW_T_s, BC_TW_s, BC_C_T_s
! User-Input: solids species eq BCs.
      use bc, only: BC_HW_X_s, BC_XW_s, BC_C_X_s

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of possible species.
      use param, only: DIM_N_S
! Parameter constants.
      use param1, only: UNDEFINED, ZERO

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC getting checked.
      INTEGER, INTENT(in) :: BCV
! Index of solid phase getting checked.
      INTEGER, INTENT(in) :: M

! Local Variables:
!---------------------------------------------------------------------//
! Loop/variable counter.
      INTEGER :: NN
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_DISCRETE")

! DEM and PIC are restricted to adiabatic walls.
      IF((BC_HW_T_S(BCV,M) /= UNDEFINED) .and. &
      &  (BC_HW_T_S(BCV,M) /= ZERO)) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_HW_T_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF((BC_C_T_S(BCV,M) /= UNDEFINED) .and. &
         &   (BC_C_T_S(BCV,M) /= ZERO))THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_C_T_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: ',A,' should not specified for DEM/PIC',/    &
         'to be non-zero as they are currently limited to constant',/  &
         'constant temperature BCs.',&
         /'Please correct the mfix.dat file.')


! The following checks verify that TFM solids parameters are not
! specified for discrete solids.


! The wall velocities are not needed DEM/PIC solids
      IF(BC_UW_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Uw_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(BC_VW_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Vw_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(BC_WW_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Ww_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! DEM cannot have a species flux at the walls.
      DO NN=1, DIM_N_s
         IF(BC_HW_X_S(BCV,M,nn) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_HW_X_s',BCV,M,nn))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_XW_S(BCV,M,nn) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Xw_s',BCV,M,nn))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_C_X_S(BCV,M,nn) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_C_X_s',BCV,M,nn))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

 1101 FORMAT('Error 1101: Illegal input for boundary condition: ',I3,/ &
         A,' should not be specified for DEM/PIC simulations.',/       &
         'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BC_WALLS_DISCRETE




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_SCALAR_EQ                                 !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for generic scalar eq WALL BC parameters.  !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_SCALAR_EQ(BCV)


! Global Variables:
!---------------------------------------------------------------------//
! User-input: number of generic scalar equations to solve.
      use scalars, only: NSCALAR
! User-Input: generic scalar eq at wall BCs.
      use bc, only: BC_HW_SCALAR, BC_SCALARw, BC_C_SCALAR

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC getting checked.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Loop/counter variable.
      INTEGER :: NN
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_SCALAR_EQ")

      DO NN=1, NSCALAR
         IF(BC_HW_Scalar(BCV,nn) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('BC_HW_SCALAR',BCV,nn)),      &
               trim(iVal(BC_HW_Scalar(BCV,nn)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_SCALAR(BCV,nn) /= ZERO .AND.                          &
            BC_SCALARw(BCV,nn) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_SCALARw',BCV,nn))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_SCALAR(BCV,nn) /= UNDEFINED .AND.                     &
            BC_C_Scalar(BCV,nn) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_C_SCALAR',BCV,nn))
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_SCALAR_EQ
