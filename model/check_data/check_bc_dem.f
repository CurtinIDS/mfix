!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! minimum amount of geometry data.                                     !
!                                                                      !
! Subroutine: CHECK_BC_DEM                                             !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_DEM(M_TOT)

! Global Variables:
!---------------------------------------------------------------------//
! User specified BC
      use bc
! User specified: BC geometry
      use bc, only: BC_EP_s
! Use specified flag for ignoring PO BC for discrete solids
      USE bc, only: BC_PO_APPLY_TO_DES
! Solids phase identifier
      use run, only: SOLIDS_MODEL
! Number of DEM inlet/outlet BCs detected.
      use des_bc, only: DEM_BCMI, DEM_BCMO
!
      use des_bc, only: DEM_BCMI_MAP
      use des_bc, only: DEM_BCMO_MAP
! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_BC
! Parameter constants
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Passed Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases.
      INTEGER, INTENT(in) :: M_TOT

! Local Variables:
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: BCV, M
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_DEM")

! Initialize
      DEM_BCMI = 0
      DEM_BCMO = 0

! Loop over all BCs looking for DEM solids inlets/outlets
      DO BCV = 1, DIMENSION_BC

         SELECT CASE (BC_TYPE_ENUM(BCV))

! Determine the number of mass inlets that contain DEM solids.
         CASE (MASS_INFLOW)
            M_LP: DO M=1,M_TOT
               IF(SOLIDS_MODEL(M)=='DEM' .AND.                         &
                  BC_EP_s(BCV,M) > ZERO) THEN
                  DEM_BCMI = DEM_BCMI + 1
                  DEM_BCMI_MAP(DEM_BCMI) = BCV
                  EXIT M_LP
               ENDIF
            ENDDO M_LP

! Count the number of pressure outflows.
         CASE (P_OUTFLOW,MASS_OUTFLOW)
            IF(BC_PO_APPLY_TO_DES(BCV)) then
               DEM_BCMO = DEM_BCMO + 1
               DEM_BCMO_MAP(DEM_BCMO) = BCV
            ENDIF

! Flag CG_MI as an error if DEM solids are present.
         CASE (CG_MI)
            DO M=1,M_TOT
               IF(SOLIDS_MODEL(M)=='DEM') THEN
                  IF(BC_EP_s(BCV,M) /= UNDEFINED .AND.                 &
                     BC_EP_s(BCV,M) > ZERO) THEN
                     WRITE(ERR_MSG,1100) trim(iVar('BC_TYPE',BCV)),    &
                        'GC_MI'
                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  ENDIF
               ENDIF
            ENDDO

         CASE (CG_PO)
            WRITE(ERR_MSG,1100) trim(iVar('BC_TYPE',BCV)), 'GC_PO'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         CASE (OUTFLOW, P_INFLOW)
            WRITE(ERR_MSG,1100) trim(iVar('BC_TYPE',BCV)),             &
               BC_TYPE_ENUM(BCV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         END SELECT

      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT('Error 1100: Unsupported boundary condition specified ',  &
         'with',/'DEM simulation: ',A,' = ',A,/'Please correct the ',&
         'mfix.dat file.')

      END SUBROUTINE CHECK_BC_DEM
