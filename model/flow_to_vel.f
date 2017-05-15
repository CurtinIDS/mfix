!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: FLOW_TO_VEL_NEW                                         !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert volumetric and mass flow rates to velocities       !
!     A specified mass flow rate is first converted to volumetric      !
!     flow rate. The volumetric flow rate is then converted to a       !
!     velocity.                                                        !
!                                                                      !
!    When both flow rates and velocities are specified, a consistency  !
!    check is done. The first time flow_to_vel is called in by setting !
!    the logical DO_VEL_CHECK to .TRUE.. If cut-cells are not used,    !
!    flow_to_vel is only called once.  When cut-cells are used,        !
!    flow_to_vel is called another time after the cut-cell pre-        !
!    processing stage. During, the second call, the velocity check     !
!    should not be performed, because the velocity assigned suring the !
!    first call will not match the flow rate. Therfore, when called    !
!    from cut_cell_preprocessing.f DO_VEL_CHECK is set to .FALSE..     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE FLOW_TO_VEL_NEW(DO_VEL_CHECK, M_TOT, SKIP, BCV)

      use param, only: DIM_M
      use param1, only: UNDEFINED
      use geometry, only: NO_I, NO_J, NO_K
      use bc, only: BC_MASSFLOW_G
      use bc, only: BC_VOLFLOW_G
      use bc, only: BC_MASSFLOW_S
      use bc, only: BC_VOLFLOW_S
      use run, only: REINITIALIZING

      use error_manager
      use toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      LOGICAL, intent(in) :: DO_VEL_CHECK

! loop/variable indices
      INTEGER, intent(in) :: M_TOT
      LOGICAL, intent(in) :: SKIP(DIM_M)

! loop/variable indices
      INTEGER, intent(in) :: BCV

! Whether any volumetric flow conversion was done
      LOGICAL :: CONVERTED = .FALSE.

! Loop index
      INTEGER :: M


      CALL INIT_ERR_MSG("FLOW_TO_VEL_NEW")

! Mass flows rates are converted to volumetric flow rates.
      IF(BC_MASSFLOW_G(BCV) /= UNDEFINED) &
         CALL GAS_MASSFLOW_TO_VOLFLOW(BCV)

      DO M=1,M_TOT
         IF(BC_MASSFLOW_S(BCV,M) /= UNDEFINED) &
            CALL SOLIDS_MASSFLOW_TO_VOLFLOW(BCV,M,SKIP(M))
      ENDDO

! Volumetric flow rates are converted to velocities.
      IF(BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
         CALL GAS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK, BCV)
! Set the conversion flag.
         CONVERTED = .TRUE.
      ENDIF

      DO M=1,M_TOT
         IF(BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
            CALL SOLIDS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK,BCV,M,SKIP(M))
! Set the conversion flag.
            CONVERTED = .TRUE.
         ENDIF
      ENDDO

      IF(CONVERTED .AND. .NOT.REINITIALIZING .AND. &
         (NO_I.OR.NO_J.OR.NO_K)) THEN
         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT('Warning 1100: Some volumetric or mass flow rates have ', &
         'been converted',/'velocity. Ensure that the third (unused) ',&
         'dimension in 2D simulations',/'is correctly specified (e.g.',&
         ', in axisymmetric cylindrical coordinates',/'ZLENGTH = 2*Pi)')

      END SUBROUTINE FLOW_TO_VEL_NEW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_MASSFLOW_TO_VOLFLOW                                 !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert a gas phase BC input from a mass flow rate to      !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_MASSFLOW_TO_VOLFLOW(BCV)

      use bc, only: BC_MASSFLOW_g
      use bc, only: BC_P_g
      use bc, only: BC_T_g
      use bc, only: BC_VOLFLOW_g
      use bc, only: BC_X_g
      use eos, only: EOSG
      use error_manager
      use param, only: DIMENSION_BC
      use param1, only: UNDEFINED
      use param1, only: ZERO
      use physprop, only: CALC_MW
      use physprop, only: MW_AVG, MW_g
      use physprop, only: NMAX
      use physprop, only: RO_g0
      use scales, only: P_REF
      use toleranc

      IMPLICIT NONE

      INTEGER, INTENT(in) :: BCV

! Volumetric flow rate computed from mass flow rate
      DOUBLE PRECISION :: VOLFLOW
! Average molecular weight
      DOUBLE PRECISION :: MW

      CALL INIT_ERR_MSG("GAS_MASSFLOW_TO_VOLFLOW")

! No need to convert if the mass flow is zero.
      IF(COMPARE(BC_MASSFLOW_G(BCV),ZERO)) THEN
         VOLFLOW = ZERO

! Incompressible gas BC.
      ELSEIF(RO_G0 /= UNDEFINED) THEN
         VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0

! Well-defined compresible gas BC.
      ELSEIF(BC_P_G(BCV)/=UNDEFINED .AND. BC_T_G(BCV)/=UNDEFINED) THEN
         IF(MW_AVG == UNDEFINED) THEN
            MW = CALC_MW(BC_X_G,DIMENSION_BC,BCV,NMAX(0),MW_G)
         ELSE
            MW = MW_AVG
         ENDIF
         VOLFLOW = BC_MASSFLOW_G(BCV) / &
            EOSG(MW,(BC_P_G(BCV)-P_REF),BC_T_G(BCV))

! Fails. This shouldn't happen as previous checks should catch any
! errors leading to this routine.
      ELSE
         WRITE(ERR_MSG, 1100) BCV
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Boundary condition ',I3,' failed sanity ',   &
      'check.',/'Please report this to the MFIX mailing list.')

      ENDIF

! Check that a specified volumetric flow matches the calculated value.
      IF(BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
         IF(.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN
            WRITE(ERR_MSG,1101) trim(iVar('BC_MASSFLOW_g',BCV)), BCV,  &
               VOLFLOW, BC_VOLFLOW_g(BCV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ELSE

! Store the calculated volumetric flow rate.
         BC_VOLFLOW_G(BCV) = VOLFLOW
      ENDIF

 1101 FORMAT('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the mfix.dat file.')

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE GAS_MASSFLOW_TO_VOLFLOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_MASSFLOW_TO_VOLFLOW                              !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert solids phase BC input from a mass flow rate to     !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SOLIDS_MASSFLOW_TO_VOLFLOW(BCV,M, SKIP_M)

      USE bc, only: BC_MASSFLOW_s
      USE bc, only: BC_VOLFLOW_s
      USE bc, only: BC_X_s
      USE param1, only: UNDEFINED, ZERO
      USE physprop, only: INERT_SPECIES
      USE physprop, only: RO_s0
      USE physprop, only: X_s0
      use eos, only: EOSS
      use error_manager
      use toleranc

      IMPLICIT NONE

! loop/variable indices
      INTEGER, INTENT(in) :: BCV, M
      LOGICAL, INTENT(in) :: SKIP_M

! Volumetric flow rate computed from mass flow rate
      DOUBLE PRECISION :: VOLFLOW
! Index of inert species
      INTEGER :: INERT

      CALL INIT_ERR_MSG("SOLIDS_MASSFLOW_TO_VOLFLOW")


      IF(SKIP_M) THEN
         WRITE(ERR_MSG,1100) M, BCV, trim(iVar("BC_MASSFLOW_S",BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Solids phase ',I2,' has a specified mass ',  &
         'flow rate',/'at BC ',I3,', ',A,'. But, both BC_ROP_s and ',&
         'BC_EP_s are zero or undefined.',/'Please correct the ',&
         'mfix.dat file.')

      IF(COMPARE(BC_MASSFLOW_S(BCV,M),ZERO)) THEN
         VOLFLOW = ZERO
! Constant solids density.
      ELSEIF(RO_S0(M) /= UNDEFINED) THEN
         VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S0(M)
      ELSE
! Set an alias for the inert species.
         INERT = INERT_SPECIES(M)
! Variable solids density.
         VOLFLOW = BC_MASSFLOW_S(BCV,M)/EOSS(RO_s0(M),              &
            X_s0(M,INERT), BC_X_S(BCV,M,INERT))
      ENDIF

! If volumetric flow is also specified compare both
      IF(BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
         IF(.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN
            WRITE(ERR_MSG,1101) trim(iVar('BC_MASSFLOW_S',BCV,M)), BCV, &
               VOLFLOW, BC_VOLFLOW_S(BCV,M)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ELSE
         BC_VOLFLOW_S(BCV,M) = VOLFLOW
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1101 FORMAT('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the mfix.dat file.')


      END SUBROUTINE SOLIDS_MASSFLOW_TO_VOLFLOW




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_VOLFLOW_TO_VELOCITY                                 !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert gas phase volumetric rate to a velocity.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK, BCV)

      USE bc
      USE compar
      USE discretelement
      USE error_manager
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE mfix_pic
      USE param
      USE param1
      USE physprop
      USE run
      USE scales
      USE toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER, INTENT(in) :: BCV

! Whether any volumetric flow conversion was done
      LOGICAL, INTENT(in) :: DO_VEL_CHECK


      DOUBLE PRECISION :: SGN, OFF

! Velocity computed from volumetric flow rate
      DOUBLE PRECISION :: VEL
!-----------------------------------------------
      CALL INIT_ERR_MSG("GAS_VOLFLOW_TO_VELOCITY")

      SELECT CASE (BC_TYPE_ENUM(BCV))
      CASE (MASS_INFLOW);  SGN =  ONE; OFF = ZERO
      CASE (MASS_OUTFLOW); SGN = -ONE; OFF = ONE
      CASE DEFAULT
        write(*,*) 'error in GAS_VOLFLOW_TO_VELOCITY'
        call mfix_exit(myPE)
      END SELECT

      SELECT CASE (BC_PLANE(BCV))
      CASE ('W'); SGN = -SGN
      CASE ('S'); SGN = -SGN
      CASE ('B'); SGN = -SGN
      END SELECT

! Calculate the velocity based on the volumetric flow rate,
! BC area and BC volume fraction.
      VEL = SGN*BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV))

! if the user also defined the boundary velocity through the plane, then
! check that the calculated value agrees with the specified value. if
! the user did not define the boundary velocity through the plane, then
! if mass_inflow set the value of the boundary velocity to the
! calculated value. otherwise do nothing.
      IF(BC_PLANE(BCV) == 'W' .OR. BC_PLANE(BCV)== 'E') THEN

         IF(BC_U_G(BCV) /= UNDEFINED .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_U_G(BCV))) THEN
               WRITE(ERR_MSG,1100) BCV, VEL, 'BC_U_g', BC_U_G(BCV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_U_G(BCV) = VEL
            BC_V_G(BCV) = OFF * BC_V_G(BCV)
            BC_W_G(BCV) = OFF * BC_W_G(BCV)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'S' .OR. BC_PLANE(BCV)== 'N') THEN
         IF(BC_V_G(BCV) /= UNDEFINED .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_V_G(BCV))) THEN
               WRITE(ERR_MSG, 1100) BCV, VEL, 'BC_V_g', BC_V_G(BCV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_V_G(BCV) = VEL
            BC_U_G(BCV) = OFF * BC_U_G(BCV)
            BC_W_G(BCV) = OFF * BC_W_G(BCV)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'B' .OR. BC_PLANE(BCV)== 'T') THEN
         IF(BC_W_G(BCV) /= UNDEFINED .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL, BC_W_G(BCV))) THEN
               WRITE(ERR_MSG, 1100) BCV, VEL, 'BC_W_g', BC_W_G(BCV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_W_G(BCV) = VEL
            BC_U_G(BCV) = OFF * BC_U_G(BCV)
            BC_V_G(BCV) = OFF * BC_V_G(BCV)
         ENDIF

      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,') = ',G14.7,/1X,70('*')/)

      END SUBROUTINE GAS_VOLFLOW_TO_VELOCITY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_VOLFLOW_TO_VELOCITY                              !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert volumetric and mass flow rates to velocities       !
!     A specified mass flow rate is first converted to volumetric      !
!     flow rate. The volumetric flow rate is then converted to a       !
!     velocity.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SOLIDS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK, BCV, M, SKIP_M)

      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE bc
      USE scales
      USE indices
      USE funits
      USE compar
      USE discretelement
      USE mfix_pic

      use error_manager
      use toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER, INTENT(in) :: BCV, M
! Logical to preform velocity check.
      LOGICAL, INTENT(in) :: DO_VEL_CHECK, SKIP_M

! Velocity computed from volumetric flow rate
      DOUBLE PRECISION :: VEL

      DOUBLE PRECISION :: SGN, OFF

!-----------------------------------------------

      CALL INIT_ERR_MSG("SOLIDS_VOLFLOW_TO_VELOCITY")

      IF(SKIP_M) THEN
         WRITE(ERR_MSG,1100) M, BCV, trim(iVar("BC_VOLFLOW_S",BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Solids phase ',I2,' has a specified ',       &
         'volumetric flow rate',/'at BC ',I3,', ',A,'. But, both ',&
         'BC_ROP_s and BC_EP_s are zero or undefined.',/'Please ',&
         'the mfix.dat file.')

      SELECT CASE (BC_TYPE_ENUM(BCV))
      CASE (MASS_INFLOW);  SGN =  ONE; OFF = ZERO
      CASE (MASS_OUTFLOW); SGN = -ONE; OFF = ONE
      CASE DEFAULT
        write(*,*) 'error in SOLIDS_VOLFLOW_TO_VELOCITY'
        call mfix_exit(myPE)
      END SELECT

      SELECT CASE (BC_PLANE(BCV))
      CASE ('W'); SGN = -SGN
      CASE ('S'); SGN = -SGN
      CASE ('B'); SGN = -SGN
      END SELECT

      IF(BC_EP_S(BCV,M) /= ZERO) THEN
         VEL = SGN * BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*BC_EP_S(BCV,M))
      ELSE
         IF(BC_VOLFLOW_S(BCV,M) == ZERO) THEN
            VEL = ZERO
         ELSE
            IF(DMP_LOG)WRITE (UNIT_LOG, 1101) BCV, M
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1101 FORMAT('Error 1101: BC No:',I2,' Non-zero vol. or mass flow ',&
         'specified with BC_ROP_s', I1,' = 0.')

      IF(BC_PLANE(BCV) == 'W' .OR. BC_PLANE(BCV)== 'E') THEN
         IF(BC_U_S(BCV,M) /= UNDEFINED .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL, BC_U_S(BCV,M))) THEN
              WRITE(ERR_MSG, 1300) BCV, (-VEL), 'BC_U_s', M, BC_U_S(BCV,M)
              CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_U_S(BCV,M) = VEL
            BC_V_S(BCV,M) = OFF * BC_V_S(BCV,M)
            BC_W_S(BCV,M) = OFF * BC_W_S(BCV,M)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'S' .OR. BC_PLANE(BCV)== 'N') THEN
         IF(BC_V_S(BCV,M) /= UNDEFINED .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_V_S(BCV,M))) THEN
               WRITE(ERR_MSG,1300) BCV, VEL, 'BC_V_s', M, BC_V_S(BCV,M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_V_S(BCV,M) = VEL
            BC_U_S(BCV,M) = OFF * BC_U_S(BCV,M)
            BC_W_S(BCV,M) = OFF * BC_W_S(BCV,M)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'B' .OR. BC_PLANE(BCV)== 'T') THEN
         IF(BC_W_S(BCV,M) /= UNDEFINED .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_W_S(BCV,M))) THEN
               WRITE(ERR_MSG, 1300) BCV, VEL, 'BC_W_s', M, BC_W_S(BCV,M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_W_S(BCV,M) = VEL
            BC_U_S(BCV,M) = OFF * BC_U_S(BCV,M)
            BC_V_S(BCV,M) = OFF * BC_V_S(BCV,M)
         ENDIF
      ENDIF


      CALL FINL_ERR_MSG

      RETURN

 1300 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,I1,') = ',G14.7,/1X,70('*')/)

      END SUBROUTINE SOLIDS_VOLFLOW_TO_VELOCITY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: FLOW_TO_VEL                                             C
!  Purpose: Convert volumetric and mass flow rates to velocities       C
!     A specified mass flow rate is first converted to volumetric      C
!     flow rate. The volumetric flow rate is then converted to a       C
!     velocity.                                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE FLOW_TO_VEL(DO_VEL_CHECK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE discretelement
      USE eos, ONLY: EOSG, EOSS
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE mfix_pic
      USE param
      USE param1
      USE physprop
      USE run
      USE scales
      USE toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER :: BCV, M
! Volumetric flow rate computed from mass flow rate
      DOUBLE PRECISION :: VOLFLOW
! Velocity computed from volumetric flow rate
      DOUBLE PRECISION :: VEL
! Solids phase volume fraction
      DOUBLE PRECISION :: EPS
! Whether any volumetric flow conversion was done
      LOGICAL :: CONVERTED,DO_VEL_CHECK
! Average molecular weight
      DOUBLE PRECISION :: MW
! Index of inert species
      INTEGER :: INERT
! Solids density at BC plane
      DOUBLE PRECISION :: BC_ROs
!-----------------------------------------------

! When both flow rates and velocities are specified, a consistency check is done
! The first time flow_to_vel is called in
! by setting the logical DO_VEL_CHECK to .TRUE.
! If cut-cells are not used, flow_to_vel is only called once.
! When cut-cells are used, flow_to_vel is called another time after
! the cut-cell preprocessing stage. During, the second call, the velocity check
! should not be performed, because the velocity assigned suring the first call
! will not match the flow rate. Therfore, when called from cut_cell_preprocessing.f
! DO_VEL_CHECK is set to .FALSE.

! initialize
      VOLFLOW = UNDEFINED

      CONVERTED = .FALSE.
      DO BCV = 1, DIMENSION_BC
         IF (BC_DEFINED(BCV)) THEN
            IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .OR.&
                BC_TYPE_ENUM(BCV)==MASS_OUTFLOW) THEN

! If gas mass flow is defined convert it to volumetric flow
! ---------------------------------------------------------------->>>
               IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN
                  IF (RO_G0 /= UNDEFINED) THEN
                     VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0
                  ELSE
                     IF (BC_P_G(BCV)/=UNDEFINED .AND. &
                         BC_T_G(BCV)/=UNDEFINED) THEN
                        IF (MW_AVG == UNDEFINED) THEN
                           MW = CALC_MW(BC_X_G,DIMENSION_BC,BCV,NMAX(0),MW_G)
                        ELSE
                           MW = MW_AVG
                        ENDIF
                        VOLFLOW = BC_MASSFLOW_G(BCV)/&
                           EOSG(MW,(BC_P_G(BCV)-P_REF),BC_T_G(BCV))

                     ELSE
! for mass_inflow, check_data_07 has already required that either ro_g0
! be defined or bc_p_g and bc_t_g be defined. So this branch will never
! be entered when mass_inflow. if mass_outflow, ro_g0 and either bc_p_g,
! or bc_t_g, or both, must be undefined.
! if the code comes through this branch without exiting and massflow is
! non-zero, then volflow will remain undefined.

                        IF (BC_TYPE_ENUM(BCV) == MASS_OUTFLOW) THEN ! this check seems unnecessary

! if no mass flow through the boundary, the volume flow is zero.
! otherwise check that the value of velocity component through the
! boundary plane is defined, and is non-zero (otherwise would be caught
! by bc_massflow_g == zero branch)
                           IF (BC_MASSFLOW_G(BCV) == ZERO) THEN
                              VOLFLOW = ZERO
                           ELSEIF (BC_PLANE(BCV)=='W' .OR. &
                                   BC_PLANE(BCV)=='E') THEN
                              IF (BC_U_G(BCV)==UNDEFINED .OR. &
                                  BC_U_G(BCV)/=ZERO) THEN
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1010)&
                                    BCV, 'BC_U_g'
                                 call mfix_exit(myPE)
                              ENDIF
                           ELSEIF (BC_PLANE(BCV)=='N' .OR. &
                                   BC_PLANE(BCV)=='S') THEN
                              IF (BC_V_G(BCV)==UNDEFINED .OR. &
                                  BC_V_G(BCV)/=ZERO) THEN
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1010) &
                                    BCV, 'BC_V_g'
                                 call mfix_exit(myPE)
                              ENDIF
                           ELSEIF (BC_PLANE(BCV)=='T' .OR. &
                                   BC_PLANE(BCV)=='B') THEN
                              IF (BC_W_G(BCV)==UNDEFINED .OR. &
                                  BC_W_G(BCV)/=ZERO) THEN
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1010) &
                                    BCV, 'BC_W_g'
                                 call mfix_exit(myPE)
                              ENDIF
                           ENDIF
                        ELSE   ! else branch if(bc_type=mass_outflow)
! not sure how this branch will be reached by mass_inflow
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1020) BCV
                           call mfix_exit(myPE)
                        ENDIF   ! end if (bc_type_enum(bcv)==mass_outflow)
                     ENDIF
                  ENDIF   ! end if/else (ro_g0 /=undefined)

! If volumetric flow is also specified compare both
                  IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
! volflow may be undefined for mass_outflow boundaries wherein ro_g0 and
! either bc_p_g, or bc_t_g, or both, were undefined.
                     IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, &
                           VOLFLOW, BC_VOLFLOW_G(BCV)
                        call mfix_exit(myPE)
                     ENDIF
                  ELSE
                     BC_VOLFLOW_G(BCV) = VOLFLOW
                  ENDIF
               ENDIF   ! end if (bc_massflow_g(bcv) /= undefined)
! end gas mass flow conversion to volumetric flow
! ----------------------------------------------------------------<<<

! If gas volumetric flow is defined convert it to velocity
! ---------------------------------------------------------------->>>
               IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN

                  IF (BC_EP_G(BCV) /= UNDEFINED) THEN
! volumetric flow rate and void fraction at the boundary are specified
! (known) so that the corresponding gas velocity through the boundary
! plane may be calculated.
                     VEL = BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV))
                  ELSE
! for mass_inflow, check_data_07 has already required that bc_ep_g be
! defined (i.e., this section will not happen for MI). For mass_outflow
! the routine will exit here if bc_ep_g is not defined.  However, for
! this MO the boundary velocities must also be set or mfix will exit due
! to later checks in check_data_07.
                     RETURN
                  ENDIF

! if the user also defined the boundary velocity through the plane, then
! check that the calculated value agrees with the specified value. if
! the user did not define the boundary velocity through the plane, then
! if mass_inflow set the value of the boundary velocity to the
! calculated value. otherwise do nothing.
                  CONVERTED = .TRUE.
                  SELECT CASE (TRIM(BC_PLANE(BCV)))
                  CASE ('W')
                     IF (BC_U_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                        IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                            .NOT.COMPARE((-VEL),BC_U_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_U_g', BC_U_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                        IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                            .NOT.COMPARE(VEL,BC_U_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, &
                              VEL, 'BC_U_g', BC_U_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                           BC_U_G(BCV) = -VEL
                           BC_V_G(BCV) = ZERO
                           BC_W_G(BCV) = ZERO
                        ELSE
                           BC_U_G(BCV) = VEL
                        ENDIF
                     ENDIF
                  CASE ('E')
                     IF (BC_U_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                        IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                            .NOT.COMPARE(VEL,BC_U_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_U_g', BC_U_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                        IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                            .NOT.COMPARE((-VEL),BC_U_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, &
                              (-VEL), 'BC_U_g', BC_U_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                           BC_U_G(BCV) = VEL
                           BC_V_G(BCV) = ZERO
                           BC_W_G(BCV) = ZERO
                        ELSE
                           BC_U_G(BCV) = -VEL
                        ENDIF
                     ENDIF
                  CASE ('S')
                     IF (BC_V_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                        IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                            .NOT.COMPARE((-VEL),BC_V_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_V_g', BC_V_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                        IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                            .NOT.COMPARE(VEL,BC_V_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_V_g', BC_V_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                           BC_U_G(BCV) = ZERO
                           BC_V_G(BCV) = -VEL
                           BC_W_G(BCV) = ZERO
                        ELSE
                           BC_V_G(BCV) = VEL
                        ENDIF
                     ENDIF
                  CASE ('N')
                     IF (BC_V_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                        IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                            .NOT.COMPARE(VEL,BC_V_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_V_g', BC_V_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                        IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                            .NOT.COMPARE((-VEL),BC_V_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, &
                              (-VEL), 'BC_V_g', BC_V_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                           BC_U_G(BCV) = ZERO
                           BC_V_G(BCV) = VEL
                           BC_W_G(BCV) = ZERO
                        ELSE
                           BC_V_G(BCV) = -VEL
                        ENDIF
                     ENDIF
                  CASE ('B')
                     IF (BC_W_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                        IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND.  &
                            .NOT.COMPARE((-VEL),BC_W_G(BCV))) THEN
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_W_g', BC_W_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                        IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                            .NOT.COMPARE(VEL,BC_W_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_W_g', BC_W_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                           BC_U_G(BCV) = ZERO
                           BC_V_G(BCV) = ZERO
                           BC_W_G(BCV) = -VEL
                        ELSE
                           BC_W_G(BCV) = VEL
                        ENDIF
                     ENDIF
                  CASE ('T')
                     IF (BC_W_G(BCV) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                        IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                            .NOT.COMPARE(VEL,BC_W_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV, VEL,&
                              'BC_W_g', BC_W_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                        IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                            .NOT.COMPARE((-VEL),BC_W_G(BCV))) THEN
                           IF(DMP_LOG) WRITE (UNIT_LOG, 1100) BCV,&
                              (-VEL), 'BC_W_g', BC_W_G(BCV)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                           BC_U_G(BCV) = ZERO
                           BC_V_G(BCV) = ZERO
                           BC_W_G(BCV) = VEL
                        ELSE
                           BC_W_G(BCV) = -VEL
                        ENDIF
                     ENDIF
                  END SELECT    ! end select (trim(bc_plane(bcv))
               ENDIF   ! end if (bc_volflow_g(bcv)/=undefined)
! end gas volumetric flow conversion to velocity
! ----------------------------------------------------------------<<<

               IF (.NOT.DISCRETE_ELEMENT .OR. (DISCRETE_ELEMENT &
                   .AND. DES_CONTINUUM_HYBRID).OR. (DISCRETE_ELEMENT &
                   .AND. MPPIC)) THEN
! The following quantities should not be required for DEM simulations
! To ensure this is the case leave them undefined in mfix.dat
! MPPIC BC's are based on Two Fluid Model based specification.
! So call the below for MPPIC.

! Do flow conversions for solids phases
               DO M = 1, SMAX

! initialize
                  VOLFLOW = UNDEFINED
                  BC_ROs = UNDEFINED

                  IF((BC_MASSFLOW_S(BCV,M) /= UNDEFINED) .OR.          &
                     (BC_VOLFLOW_S(BCV,M)  /= UNDEFINED)) THEN

! Calculate the solid density.
                     BC_ROs = UNDEFINED
                     IF(SOLVE_ROs(M))THEN
                        INERT = INERT_SPECIES(M)
! Verify that the species mass fraction for the inert material is not
! zero in the IC region when the solids is present.
                        IF(BC_X_S(BCV,M,INERT) == ZERO) THEN
                           IF(BC_ROP_S(BCV,M) /= ZERO) THEN
                              IF(DMP_LOG) THEN
                                 WRITE(*,1401) M, BCV
                                 WRITE(UNIT_LOG,1401) M, BCV
                              ENDIF
                              CALL MFIX_EXIT(myPE)
                           ELSE
! If the solids isn't present, give it the baseline density.
                              BC_ROs = RO_s0(M)
                           ENDIF
                        ELSE
! Calculate the solids density.
                           BC_ROs = EOSS(RO_s0(M),X_s0(M,INERT),    &
                              BC_X_S(BCV,M,INERT))
                        ENDIF
                     ELSE
                        BC_ROs = RO_S0(M)
                     ENDIF
                  ELSE
! Set a generic value to pass through sanity check.
                     BC_ROs = RO_S0(M)
                  ENDIF


! If solids mass flow is defined convert it to volumetric flow
! ---------------------------------------------------------------->>>
                  IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED) THEN

! Sanity check on solids phase density.
                     IF(BC_ROs <= ZERO .OR. BC_ROs==UNDEFINED) THEN
                        IF(DMP_LOG)THEN
                           WRITE(*,1401) M, BCV
                           WRITE(UNIT_LOG,1401) M, BCV
                        ENDIF
                        CALL MFIX_EXIT(myPE)
                     ENDIF

                     VOLFLOW = BC_MASSFLOW_S(BCV,M)/BC_ROs

! If volumetric flow is also specified compare both
                     IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
                        IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN
                           IF(DMP_LOG) WRITE(UNIT_LOG,1200) &
                              BCV,VOLFLOW,M,BC_VOLFLOW_S(BCV,M)
                           call mfix_exit(myPE)
                        ENDIF
                     ELSE
                        BC_VOLFLOW_S(BCV,M) = VOLFLOW
                     ENDIF
                  ENDIF   ! end if (bc_massflow_s(bcv,m)/=undefined)
! end solids mass flow conversion to volumetric flow
! ----------------------------------------------------------------<<<


! if possible, define bulk density based on ep_g. note if mass_outflow,
! bc_ep_g may still be undefined at this point wherein bc_rop_s would
! become set as 1-undefined. to avoid this issue, a check for bc_ep_g
! defined was added. note if mass_inflow, check_data_07 later performs
! the same calculations as below but with additional checks.
                  IF (BC_ROP_S(BCV,M) == UNDEFINED .AND. &
                     BC_EP_G(BCV) /= UNDEFINED ) THEN
                     IF (BC_EP_G(BCV) == ONE) THEN
                         BC_ROP_S(BCV,M) = ZERO
                     ELSEIF (SMAX == 1 .AND. &
                              .NOT.DES_CONTINUUM_HYBRID) THEN

! Sanity check on solids phase density.
                        IF(BC_ROs <= ZERO .OR. BC_ROs==UNDEFINED) THEN
                           IF(DMP_LOG)THEN
                              WRITE(*,1401) M, BCV
                              WRITE(UNIT_LOG,1401) M, BCV
                           ENDIF
                           CALL MFIX_EXIT(myPE)
                        ENDIF


! bulk density must be explicitly defined for hybrid model and cannot be
! defined from 1-bc_ep_g
                         BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*BC_ROs
                     ENDIF
                  ENDIF
! note bc_rop_s may still be undefined at this point


! If solids volumetric flow is defined convert it to velocity
! ---------------------------------------------------------------->>>
                  IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN

! Sanity check on solids phase density.
                     IF(BC_ROs <= ZERO .OR. BC_ROs==UNDEFINED) THEN
                        IF(DMP_LOG)THEN
                           WRITE(*,1401) M, BCV
                           WRITE(UNIT_LOG,1401) M, BCV
                        ENDIF
                        CALL MFIX_EXIT(myPE)
                     ENDIF

                     IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN
                        EPS = BC_ROP_S(BCV,M)/BC_ROs
! volumetric flow rate and solids volume fraction at the boundary are
! specified (known) so that the corresponding solids velocity through
! the boundary plane may be calculated.
                        IF (EPS /= ZERO) THEN
                           VEL = BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*EPS)
                        ELSE
                           IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN
                              VEL = ZERO
                           ELSE
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1250) BCV, M
                              call mfix_exit(myPE)
                           ENDIF
                        ENDIF
                     ELSE   ! bc_rop_s is undefined
                        IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN
                           VEL = ZERO
                        ELSE
! if bc_rop_s is undefined and bc_volflow_s is not zero exit MFIX with
! error
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1260) BCV, M
                           call mfix_exit(myPE)
                        ENDIF
                     ENDIF

! if the user also defined the boundary velocity through the plane, then
! check that the calculated value agrees with the specified value. if
! the user did not define the boundary velocity through the plane, then
! if mass_inflow set the value of the boundary velocity to the
! calculated value. otherwise do nothing.
                     CONVERTED = .TRUE.
                     SELECT CASE (TRIM(BC_PLANE(BCV)))
                     CASE ('W')
                        IF (BC_U_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                           IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                               .NOT.COMPARE((-VEL),BC_U_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_U_s', M, BC_U_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                           IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                               .NOT.COMPARE(VEL,BC_U_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_U_s', M, BC_U_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                        ELSE
                           IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                              BC_U_S(BCV,M) = -VEL
                              BC_V_S(BCV,M) = ZERO
                              BC_W_S(BCV,M) = ZERO
                           ELSE
                              BC_U_S(BCV,M) = VEL
                           ENDIF
                        ENDIF
                     CASE ('E')
                        IF (BC_U_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                           IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                               .NOT.COMPARE(VEL,BC_U_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_U_s', M, BC_U_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                           IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                               .NOT.COMPARE((-VEL),BC_U_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_U_s', M, BC_U_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                        ELSE
                           IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                              BC_U_S(BCV,M) = VEL
                              BC_V_S(BCV,M) = ZERO
                              BC_W_S(BCV,M) = ZERO
                           ELSE
                              BC_U_S(BCV,M) = -VEL
                           ENDIF
                        ENDIF
                     CASE ('S')
                        IF (BC_V_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                           IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                               .NOT.COMPARE((-VEL),BC_V_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_V_s', M, BC_V_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                           IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                               .NOT.COMPARE(VEL,BC_V_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_V_s', M, BC_V_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                        ELSE
                           IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                              BC_U_S(BCV,M) = ZERO
                              BC_V_S(BCV,M) = -VEL
                              BC_W_S(BCV,M) = ZERO
                           ELSE
                              BC_V_S(BCV,M) = VEL
                           ENDIF
                        ENDIF
                     CASE ('N')
                        IF (BC_V_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                           IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                               .NOT.COMPARE(VEL,BC_V_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_V_s', M, BC_V_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                           IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                               .NOT.COMPARE((-VEL),BC_V_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_V_s', M, BC_V_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                        ELSE
                           IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                              BC_U_S(BCV,M) = ZERO
                              BC_V_S(BCV,M) = VEL
                              BC_W_S(BCV,M) = ZERO
                           ELSE
                              BC_V_S(BCV,M) = -VEL
                           ENDIF
                        ENDIF
                     CASE ('B')
                        IF (BC_W_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                           IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                               .NOT.COMPARE((-VEL),BC_W_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_W_s', M, BC_W_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                           IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                               .NOT.COMPARE(VEL,BC_W_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_W_s', M, BC_W_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                        ELSE
                           IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                              BC_U_S(BCV,M) = ZERO
                              BC_V_S(BCV,M) = ZERO
                              BC_W_S(BCV,M) = -VEL
                           ELSE
                              BC_W_S(BCV,M) = VEL
                           ENDIF
                        ENDIF
                     CASE ('T')
                        IF (BC_W_S(BCV,M) /= UNDEFINED.AND.DO_VEL_CHECK) THEN
                           IF (BC_TYPE_ENUM(BCV)==MASS_INFLOW .AND. &
                               .NOT.COMPARE(VEL,BC_W_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 VEL, 'BC_W_s', M, BC_W_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                           IF (BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .AND. &
                              .NOT.COMPARE((-VEL),BC_W_S(BCV,M))) THEN
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1300) BCV, &
                                 (-VEL), 'BC_W_s', M, BC_W_S(BCV,M)
                              call mfix_exit(myPE)
                           ENDIF
                        ELSE
                           IF (BC_TYPE_ENUM(BCV) == MASS_INFLOW) THEN
                              BC_U_S(BCV,M) = ZERO
                              BC_V_S(BCV,M) = ZERO
                              BC_W_S(BCV,M) = VEL
                           ELSE
                              BC_W_S(BCV,M) = -VEL
                           ENDIF
                        ENDIF
                     END SELECT    ! end select (trim(bc_plane(bcv))
                  ENDIF   ! end if (bc_volflow_s(bcv,m) /= undefined)
! end solids volumetric flow conversion to velocity
! ----------------------------------------------------------------<<<

               ENDDO   ! end do m = 1,smax
               ENDIF   ! end if (.not.discrete_element)

            ENDIF   ! end if (bc_type_enum(bcv)==mass_inflow or mass_outflow)
         ENDIF   ! end if (bc_defined(bcv))
      ENDDO   ! end do (bcv =1,dimension_bc)


      IF (CONVERTED .AND. (NO_I .OR. NO_J .OR. NO_K) &
          .AND. DMP_LOG)WRITE (UNIT_LOG, 1500)

      RETURN

 1000 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_g) = ',G14.7,/1X,70('*')/)
 1010 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g or',/' a nonzero value for ',A,&
         ' should be specified',/1X,70('*')/)
 1020 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g',/' should be specified',/1X,70('*')/)
 1100 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,') = ',G14.7,/1X,70('*')/)
 1200 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_s',I1,') = ',G14.7,/1X,70('*')/)
 1250 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Non-zero vol. or mass flow specified with BC_ROP_s',&
         I1,' = 0.',/1X,70('*')/)
 1260 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' BC_ROP_s',I1,' not specified',/1X,70('*')/)
 1300 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,I1,') = ',G14.7,/1X,70('*')/)
 1500 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/&
         ' Message: Some volumetric or mass flow rates have been',/&
         '   converted to velocity values.  In 2D simulations ensure',/&
         '   that the third (unused) dimension is correctly specified;',/&
         '   e.g. in axisymmetric cylindrical coordinates ZLENGTH = 2*Pi'/1X,70&
         ('*')/)

 1401 FORMAT(//1X,70('*')/' From: FLOW_TO_VEL',/,' Error 1401:',       &
         ' Solids phase ',I2,' failed sanity check in BC region ',I3,  &
         '. ',/' Please check mfix.dat file.',/1X,70('*')//)


      END SUBROUTINE FLOW_TO_VEL
