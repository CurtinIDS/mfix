!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_INFLOW                                          !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on common inflow BC       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_x_g, bc_x_s
      use bc, only: bc_t_g, bc_t_s, bc_theta_m
      use bc, only: bc_k_turb_g, bc_e_turb_g
      use bc, only: bc_scalar
      use param, only: dim_m
      use param1, only: undefined, one, zero
      use physprop, only: inert_species
      use physprop, only: mu_g0
      use physprop, only: mw_avg, nmax
      use physprop, only: ro_g0
      use run, only: energy_eq, granular_energy, k_epsilon
      use run, only: solids_model, solve_ros, species_eq
      use scalars, only: nscalar
      use toleranc, only: compare

      use error_manager
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local variables
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: M, N
      DOUBLE PRECISION SUM_X
! Index of inert species
      INTEGER :: INERT
!---------------------------------------------------------------------//

      CALL INIT_ERR_MSG("CHECK_BC_INFLOW")

! Check temperature dependency.
      IF((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR.MU_G0==UNDEFINED) .AND. &
         BC_T_G(BCV)==UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('BC_T_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Sum together defined gas phase species mass fractions.
      SUM_X = ZERO
      DO N = 1, NMAX(0)
         IF(BC_X_G(BCV,N) /= UNDEFINED) THEN
            SUM_X = SUM_X + BC_X_G(BCV,N)
         ELSE
            BC_X_G(BCV,N) = ZERO
         ENDIF
      ENDDO

! Enforce that the species mass fractions must sum to one.
      IF(.NOT.COMPARE(ONE,SUM_X)) THEN

         IF(SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1110) BCV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1110 FORMAT('Error 1110: BC_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'species equations are solved. Please correct ', &
         'the mfix.dat file.')

         ELSEIF(RO_G0 == UNDEFINED .AND. MW_AVG == UNDEFINED) THEN
            WRITE(ERR_MSG, 1111) BCV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1111 FORMAT('Error 1111: BC_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'is compressible and MW_AVG is UNDEFINED.',/     &
         'Please correct the mfix.dat the mfix.dat file.')

         ELSEIF(.NOT.COMPARE(SUM_X,ZERO)) THEN
            WRITE(ERR_MSG, 1112) BCV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1112 FORMAT('Error 1112: BC_X_g(',I3,',:) do not sum to ONE or ZERO ',&
         'and they',/'are not needed. Please correct the mfix.dat ',   &
         'the mfix.dat file.')

         ELSE
            BC_X_G(BCV,:) = ZERO
            BC_X_G(BCV,1) = ONE
         ENDIF
      ENDIF


! Verify that species mass fractions are defined for mass flow BCs when
! using variable solids density. Needed to calculation RO_s
      DO M = 1, M_TOT

! If this phase is not present, clear out x_s for the BC and
! cycle the solids loop. No need to continue checks.
         IF(SKIP(M)) THEN
           IF(SPECIES_EQ(M))THEN
               BC_X_S(BCV,M,:) = ZERO
               BC_X_S(BCV,M,1) = ONE
            ENDIF
            CYCLE
         ENDIF

! Sum together defined species mass fractions.
         SUM_X = ZERO
         DO N = 1, NMAX(M)
            IF(BC_X_S(BCV,M,N) /= UNDEFINED) THEN
               SUM_X = SUM_X + BC_X_S(BCV,M,N)
            ELSE
               BC_X_S(BCV,M,N) = ZERO
            ENDIF
         ENDDO

! Enforce that the species mass fractions must sum to one.
         IF(.NOT.COMPARE(ONE,SUM_X)) THEN
            IF(SPECIES_EQ(M)) THEN
               WRITE(ERR_MSG, 1210) BCV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1210 FORMAT('Error 1210: BC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'and the solids phase',/'species equations are solved. ',     &
         'Please correct the mfix.dat file.')

            ELSEIF(SOLVE_ROS(M)) THEN
               WRITE(ERR_MSG, 1211) BCV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1211 FORMAT('Error 1211: BC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'and the solids phase',/'density is calculated. Please ',     &
         'correct the mfix.dat file.')

            ELSEIF(.NOT.COMPARE(SUM_X,ZERO)) THEN
               WRITE(ERR_MSG, 1212) BCV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1212 FORMAT('Error 1212: BC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'or ZERO and',/'they are not needed. Please correct the ',    &
         'mfix.dat file.')

            ELSE
               BC_X_S(BCV,M,:) = ZERO
               BC_X_S(BCV,M,1) = ONE
            ENDIF
         ENDIF

! Set the solids density for the BC region.
         IF(SOLVE_ROs(M)) THEN
! Verify that the species mass fraction for the inert material is not
! zero in the IC region when the solids is present.
            INERT = INERT_SPECIES(M)
            IF(BC_X_S(BCV,M,INERT) == ZERO) THEN
               WRITE(ERR_MSG,1213) M, BCV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1213 FORMAT('Error 1213: No inert species for phase ',I2,' in BC ',   &
         'region',I3,'.',/'Unable to calculate solids phase density. ',&
         'Please refer to the Readme',/' file for required variable ', &
         'solids density model input parameters and',/' make the ',   &
         'necessary corrections to the mfix.dat file.')

            ENDIF
         ENDIF
      ENDDO


      DO M = 1, M_TOT
! Check solids phase temperature dependency.
         IF(ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN
            IF(SKIP(M)) THEN
               BC_T_S(BCV,M) = BC_T_G(BCV)
            ELSE
               WRITE(ERR_MSG, 1000) trim(iVar('BC_T_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

! Check granular energy dependency
         IF(GRANULAR_ENERGY) THEN
            IF(BC_THETA_M(BCV,M) == UNDEFINED) THEN
               IF(SKIP(M) .OR. SOLIDS_MODEL(M) /= 'TFM') THEN
                  BC_THETA_M(BCV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG,1000) trim(iVar('BC_Theta_m',BCV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

! Check K-Epsilon BCs.
      IF(K_Epsilon) THEN
         IF(BC_K_Turb_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('BC_K_Turb_G',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(BC_E_Turb_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('BC_E_Turb_G',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Check scalar equation BCs.
      DO N = 1, NScalar
         IF(BC_Scalar(BCV,N) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1001) trim(iVar('BC_Scalar',BCV,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_INFLOW                                     !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on BC                     !
!                                                                      !
! Comments:                                                            !
!     The velocities at the inflow face are fixed and the momentum     !
!     equations are not solved in the inflow cells. Since the flow is  !
!     into the domain all other scalars that are used need to be       !
!     specified (e.g., mass fractions, void fraction, etc.,)           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_ep_g, bc_p_g
      use bc, only: bc_rop_s, bc_ep_s, bc_x_s
      use bc, only: bc_massflow_g
      use eos, only: eoss
      use param, only: dim_m
      use param1, only: undefined, one, zero
      use physprop, only: ro_g0, ro_s0
      use physprop, only: inert_species, x_s0
      use run, only: solve_ros
      use toleranc, only: compare
      use error_manager
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local variables
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: M
      DOUBLE PRECISION :: SUM_EP
! Solids phase density in BC region.
      DOUBLE PRECISION :: BC_ROs(DIM_M)
! Index of inert species
      INTEGER :: INERT

!---------------------------------------------------------------------//

      CALL INIT_ERR_MSG("CHECK_BC_MASS_INFLOW")

! Check gas phase volume fraction.
      IF(BC_EP_G(BCV) == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('BC_EP_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Verify compressible boundary condition variables.
      IF(RO_G0 == UNDEFINED) THEN
         IF(BC_P_G(BCV) == UNDEFINED) THEN
            IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .AND.                   &
               BC_MASSFLOW_G(BCV) /= ZERO) THEN
               WRITE(ERR_MSG, 1100) trim(iVar('BC_P_g',BCV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
 1100 FORMAT('Error 1100: ',A,' must be specified for compressible ',  &
         'flows',/'when specifying BC_MASSFLOW_g to make the ',        &
         'conversion to velocity.',/'Please correct the mfix.dat file.')

         ELSEIF(BC_P_G(BCV) <= ZERO) THEN
            WRITE(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1101 FORMAT('Error 1101: Pressure must be greater than zero for ',    &
         'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
         'correct the mfix.dat file.')
      ENDIF

! Calculate the solids volume fraction from the gas phase if there is
! only one solids phase.
      IF(M_TOT == 1 .AND. BC_EP_S(BCV,1) == UNDEFINED) THEN
         BC_EP_S(BCV,1) = ONE - BC_EP_g(BCV)
      ENDIF

! Bulk density or solids volume fraction must be explicitly defined
! if there are more than one solids phase.
      IF(M_TOT > 1 .AND. .NOT.COMPARE(BC_EP_g(BCV),ONE)) THEN
         DO M = 1, M_TOT
            IF(BC_ROP_S(BCV,M) == UNDEFINED .AND. &
               BC_EP_S(BCV,M) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1200) M, BCV, 'BC_ROP_s and BC_EP_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF
 1200 FORMAT('Error 1200: Insufficient solids phase ',I2,' data ',     &
         'for BC',I3,'. ',/A,' not specified.',/'Please correct the ', &
         'mfix.dat file.')

! Initialize the sum of the total volume fraction.
      SUM_EP = BC_EP_G(BCV)
      DO M = 1, M_TOT

! If this phase is not present, clear out EPs and ROPs for the BC and
! cycle the solids loop. No need to continue checks.
         IF(SKIP(M)) THEN
            BC_EP_S(BCV,M)  = ZERO
            BC_ROP_S(BCV,M) = ZERO
            CYCLE
         ENDIF

! Set the solids density for the BC region.
         IF(.NOT.SOLVE_ROs(M)) THEN
            BC_ROs(M) = RO_s0(M)
         ELSE
! presence of non-zero inert species is checked by bc_inflow
            INERT = INERT_SPECIES(M)
            BC_ROs(M) = EOSS(RO_s0(M), X_s0(M,INERT),&
               BC_X_S(BCV,M,INERT))
         ENDIF

! If both input parameters are defined. Make sure they are equivalent.
         IF(BC_ROP_S(BCV,M) /= UNDEFINED .AND.                         &
            BC_EP_S(BCV,M) /= UNDEFINED) THEN

            IF(.NOT.COMPARE(BC_EP_S(BCV,M)*BC_ROs(M),                  &
               BC_ROP_S(BCV,M))) THEN
               WRITE(ERR_MSG,1214) BCV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
 1214 FORMAT('Error 1214: Illegal initial condition region : ',I3,/    &
         'BC_EP_s and BC_ROP_s are inconsistent. Please correct the ',/&
         'mfix.dat file.')

! Compute BC_EP_s from BC_ROP_s
         ELSEIF(BC_EP_S(BCV,M) == UNDEFINED) THEN
            BC_EP_S(BCV,M) = BC_ROP_S(BCV,M) / BC_ROs(M)

! Compute BC_ROP_s from BC_EP_s and BC_ROs
         ELSEIF(BC_ROP_S(BCV,M) == UNDEFINED) THEN
            BC_ROP_S(BCV,M) = BC_EP_S(BCV,M) * BC_ROs(M)

         ENDIF
! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + BC_EP_S(BCV,M)
      ENDDO

! Verify that the volume fractions sum to one.
      IF(.NOT.COMPARE(SUM_EP,ONE)) THEN
         WRITE(ERR_MSG,1215) BCV, trim(iVal(SUM_EP))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1215 FORMAT('Error 1215: Illegal boundary condition region: ',I3,'. ',&
         'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
         ')',/'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_MASS_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_P_INFLOW                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided detailed error message on bc                       !
!                                                                      !
! Comments:                                                            !
!     Unlike the MI boundary, for the PI boundary the velocities at    !
!     the inflow face are calculated by solving the momentum eqns      !
!     and are not fixed. In this way, the PI is similar to the PO      !
!     except that the flow is into the domain and hence all other      !
!     scalars (e.g., mass fractions, void fraction, temperature,       !
!     etc.,) at the inflow cells need to be specified. To satisfy      !
!     the error routines at the start of the simulation, both the      !
!     tangential and normal components at the inflow also need to      !
!     be specified. The velocities values essentially serve as IC.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_P_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_p_g, bc_rop_s
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use geometry, only: no_i, no_j, no_k
      use param, only: dim_m
      use param1, only: undefined, zero
      use physprop, only: ro_g0
      use error_manager
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: M

!---------------------------------------------------------------------//

      CALL INIT_ERR_MSG("CHECK_BC_P_INFLOW")

! Remove checks on bc_ep_g/bc_rop_s; using routine check_bc_outflow

      IF (BC_P_G(BCV) == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'BC_P_g', BCV
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN
         WRITE(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV)))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1101 FORMAT('Error 1101: Pressure must be greater than zero for ',    &
         'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
         'correct the mfix.dat file.')
      ENDIF

! Check that velocities are also specified. These are essentially used
! as initial conditions for the boundary region. If they are not
! specified then a default value is set here otherwise check_data_20
! will complain and cause MFIX to exit.
      IF(BC_U_G(BCV) == UNDEFINED) THEN
         BC_U_G(BCV) = ZERO
         IF(.NOT.NO_I) WRITE(ERR_MSG, 1300) trim(iVar('BC_U_g',BCV))
      ENDIF

      IF(BC_V_G(BCV) == UNDEFINED) THEN
         BC_V_G(BCV) = ZERO
         IF(.NOT.NO_J) WRITE(ERR_MSG, 1300) trim(iVar('BC_V_g',BCV))
      ENDIF

      IF(BC_W_G(BCV) == UNDEFINED) THEN
         BC_W_G(BCV) = ZERO
         IF(.NOT.NO_K) WRITE(ERR_MSG, 1300) trim(iVar('BC_W_g',BCV))
      ENDIF

      DO M = 1, M_TOT
         IF (SKIP(M)) THEN
            BC_U_S(BCV,M) = ZERO
            BC_V_S(BCV,M) = ZERO
            BC_W_S(BCV,M) = ZERO
         ELSE
            IF(BC_U_S(BCV,M) == UNDEFINED) THEN
               BC_U_S(BCV,M) = ZERO
               IF(BC_ROP_S(BCV,M) /= ZERO .AND. .NOT.NO_I) &
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M))
            ENDIF

            IF(BC_V_S(BCV,M) == UNDEFINED) THEN
               BC_V_S(BCV,M) = ZERO
               IF(BC_ROP_S(BCV,M) /= ZERO .AND. .NOT.NO_J) &
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M))
            ENDIF

            IF(BC_W_S(BCV,M) == UNDEFINED) THEN
               BC_W_S(BCV,M) = ZERO
               IF(BC_ROP_S(BCV,M) /= ZERO .AND. .NOT.NO_K) &
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M))
            ENDIF
         ENDIF
      ENDDO

 1300 FORMAT('Warning 1300: ',A,' was undefined. This variable was ', &
         'set ',/ 'to zero to be used as the initial value in the BC ',&
         'region.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_BC_P_INFLOW
