!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_OUTFLOW                                         !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message concerning specification  !
! of bc_ep_g + bc_rop_s at a outflow boundary (and pressure inflow)    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_OUTFLOW(M_TOT, BCV)

! Modules
! --------------------------------------------------------------------//
      use bc, only: bc_ep_g, bc_rop_s
      use param1, only: one, undefined, zero
      use physprop, only: ro_s0
      use run, only: solids_model
      use toleranc, only: compare
      use error_manager

      IMPLICIT NONE

! Dummy arguments
! --------------------------------------------------------------------//
! loop/variable indices
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT

! Local variables
! --------------------------------------------------------------------//
      INTEGER :: M
      DOUBLE PRECISION :: SUM_EP
      LOGICAL :: FLAG_WARNING

      FLAG_WARNING = .TRUE.
      CALL INIT_ERR_MSG("CHECK_BC_OUTFLOW")

! if bc_ep_g is defined at the outflow boundary, then the sum of ep_g
! and ep_s at the boundary may not equal one given the code in the
! subroutine set_outflow (see code for details).
! therefore if bc_ep_g and/or bc_rop_s are defined, perform possible
! data consistency checks and, when appropriate, provide the user with
! a warning about their chosen settings.

      IF (BC_EP_G(BCV) /= UNDEFINED) THEN

         SUM_EP = BC_EP_G(BCV)
         DO M = 1, M_TOT

            IF(SOLIDS_MODEL(M) /= 'TFM' .AND. FLAG_WARNING) THEN
               WRITE(ERR_MSG, 1101) trim(iVar('BC_EP_g',BCV))
               CALL FLUSH_ERR_MSG
               FLAG_WARNING = .FALSE.
            ENDIF

            IF(BC_ROP_S(BCV,M) == UNDEFINED) THEN

               IF(BC_EP_G(BCV) == ONE) THEN
! what does it mean to force the bulk density to zero at the
! boundary? (does this value matter anyway?)
                  BC_ROP_S(BCV,M) = ZERO
               ELSEIF(M_TOT == 1 ) THEN
                  BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S0(M)
               ELSE
! bc_ep_g is defined but some bc_rop_s(m) are undefined.
! in this case, ep_p in the outflow boundary will be based on the user
! defined value of bc_ep_g, while rop_s would become based on the
! value in the adjacent fluid cell. consequently, no check ensures
! the result is consistent with a requirement for ep_g+ep_s=1.
                  WRITE(ERR_MSG, 1102) trim(iVar('BC_EP_g',BCV))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1102 FORMAT('Warning 1102: Volume fraction may not sum to one when ',/&
         A,' is defined.')
               ENDIF
            ENDIF  ! end if(bc_rop_s(bcv,m) == undefined)

! by this point bc_rop_s should either be defined or mfix exited
! therefore we can check that sum of void fraction and solids volume
! fractions
            SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M)
        ENDDO

! now verify that the volume fractions sum to one.
        IF(.NOT.COMPARE(SUM_EP,ONE)) THEN
           WRITE(ERR_MSG,1103) BCV, trim(iVal(SUM_EP))
           CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
       ENDIF

! bc_ep_g is not defined but check if any bc_rop_s are defined
      ELSE

         SUM_EP = ZERO
         DO M = 1, M_TOT
            IF(BC_ROP_S(BCV,M) /= UNDEFINED) THEN
               IF(SOLIDS_MODEL(M) /= 'TFM') THEN
                  WRITE(ERR_MSG, 1101) trim(iVar('BC_ROP_s',BCV,M))
                  CALL FLUSH_ERR_MSG
               ENDIF
               SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M)
            ENDIF
         ENDDO

! verify that the sum of any specified volume fractions is not greater
! than one
         IF(SUM_EP > ONE) THEN
            WRITE(ERR_MSG,1103) BCV, trim(iVal(SUM_EP))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1101 FORMAT('Warning 1101: ',A,' should not be specified for ', &
         'outflow BCs',/'with DEM/PIC runs except for a mass outflow ',&
         'boundary with specified ',/ 'flow rate(s). In this case ',&
         'volume fraction data is used for ',/ 'conversion to ',&
         'velocity(s). However, the solids volume fraction data ',/&
         'is effectively disregarded and it is the solids velocity ',&
         'that is ',/'used to direct any solids at the boundary.')

 1103 FORMAT('Error 1103: Illegal boundary condition region: ',I3,'. ',&
         'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
         ')',/'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_P_OUTFLOW                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on bc                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_P_OUTFLOW(M_TOT, BCV)

! Modules
! --------------------------------------------------------------------//
      USE param1, only: UNDEFINED
      USE param1, only: ZERO
      use physprop, only: RO_g0
      use bc, only: BC_P_g
      use error_manager
      IMPLICIT NONE

! Dummy arguments
! --------------------------------------------------------------------//
! loop/variable indices
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
! --------------------------------------------------------------------//

      CALL INIT_ERR_MSG("CHECK_BC_P_OUTFLOW")

      IF (BC_P_G(BCV) == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_P_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN
         WRITE(ERR_MSG, 1100) BCV, trim(iVal(BC_P_G(BCV)))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Pressure must be greater than zero for ',    &
         'compressible flow',/3x,'BC_P_g(',I3,') = ',A,/'Please ',     &
         'correct the mfix.dat file.')

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_P_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_OUTFLOW                                    !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
! Comments:                                                            !
!     The velocities at the outflow face are fixed and the momentum    !
!     equations are not solved in the outflow cells. Since the flow    !
!     is out of the domain none of the other scalars should need to    !
!     be specified (e.g., mass fractions, void fraction, etc.,).       !
!     Such values will become defined according to their adjacent      !
!     fluid cell                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_OUTFLOW(M_TOT, BCV)

! Modules
! --------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_dt_0, bc_massflow_g, bc_volflow_g
      use bc, only: bc_massflow_s, bc_volflow_s
      use bc, only: bc_ep_g, bc_rop_s
      use bc, only: bc_p_g, bc_t_g
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use physprop, only: ro_g0
      use param1, only: undefined, zero
      use error_manager
      IMPLICIT NONE

! Dummy arguments
! --------------------------------------------------------------------//
! loop/variable indices
      INTEGER, intent(in) :: BCV
      INTEGER, intent(in) :: M_TOT
! Local variables
! --------------------------------------------------------------------//
      INTEGER :: M


      CALL INIT_ERR_MSG("CHECK_BC_MASS_OUTFLOW")

      IF(BC_DT_0(BCV) == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('BC_DT_0',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .OR. &
         BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
         IF (BC_EP_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG,1101) trim(iVar('BC_EP_G',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1101 FORMAT('Error 1101: Invalid mass outflow boundary condition: ', /&
         'BC_MASSFLOW_G and/or BC_VOLFLOW_G are DEFINED but ',&
         A,' is not ',/'Please correct the mfix.dat file.')
      ENDIF

      DO M = 1, M_TOT
         IF(BC_MASSFLOW_S(BCV,M) /= UNDEFINED .OR. &
            BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1102) trim(iVar('BC_MASSFLOW_S',BCV,M)), &
               trim(iVar('BC_VOLFLOW_S',BCV,M))
 1102 FORMAT('Warning 1102: ', A,' and/or ', A,' have been defined',/&
         'at a mass outflow boundary. A specified solids flow ',&
         'rate may not be ',/'physically achievable depending on the ',&
         'system and simulation ',/'setup.')

             IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN
                WRITE(ERR_MSG,1103) trim(iVar('BC_ROP_S',BCV,M))
                CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
             ENDIF
 1103 FORMAT('Error 1103: Invalid mass outflow boundary condition: ', /&
         'BC_MASSFLOW_S and/or BC_VOLFLOW_S are DEFINED but ',&
         A,' is not ',/'Please correct the mfix.dat file.')

         ENDIF
      ENDDO

! This check probably needs changed.
      IF(RO_G0 == UNDEFINED .AND. (BC_P_G(BCV) == UNDEFINED .OR.       &
         BC_T_G(BCV) == UNDEFINED) .AND.BC_MASSFLOW_G(BCV) /= ZERO) THEN

         IF(BC_PLANE(BCV)=='W' .OR. BC_PLANE(BCV)=='E') THEN
            IF(BC_U_G(BCV) /= ZERO) THEN
               WRITE(ERR_MSG, 1100) BCV, 'BC_U_g'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSEIF(BC_PLANE(BCV)=='N' .OR. BC_PLANE(BCV)=='S') THEN
            IF(BC_V_G(BCV) /= ZERO) THEN
               WRITE(ERR_MSG, 1100) BCV, 'BC_V_g'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSEIF (BC_PLANE(BCV)=='T' .OR. BC_PLANE(BCV)=='B') THEN
            IF(BC_W_G(BCV) /= ZERO) THEN
               WRITE(ERR_MSG, 1100)  BCV, 'BC_W_g'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDIF   ! end if/else (ro_g0 /=undefined)

 1100 FORMAT('Error 1100: Invalid mass outflow boundary condition: ',  &
         I3,/'RO_g0, BC_P_g, and BC_T_g are UNDEFINED and ',A,' is ',  &
         'non-zero',/'Please correct the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_MASS_OUTFLOW
