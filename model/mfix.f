! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

      PROGRAM MFIX

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE DISCRETELEMENT, ONLY: DES_CONTINUUM_COUPLED, DISCRETE_ELEMENT
      USE ITERATE, ONLY: CONVERGED, DIVERGED, ADJUSTDT
      USE ITERATE, ONLY: ITERATE_INIT, DO_ITERATION, POST_ITERATE
      USE ITERATE, ONLY: LOG_CONVERGED, LOG_DIVERGED, NIT, MAX_NIT
      USE MAIN, ONLY: ADD_COMMAND_LINE_ARGUMENT, INITIALIZE, FINALIZE, EXIT_SIGNAL
      USE RUN, ONLY:  DT, IER, DEM_SOLIDS, PIC_SOLIDS, STEADY_STATE, TIME, TSTOP
      USE STEP, ONLY: TIME_STEP_INIT, TIME_STEP_END

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      INTEGER :: II
      CHARACTER(LEN=80) :: tmp

      DO II=1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(II,tmp)
         CALL ADD_COMMAND_LINE_ARGUMENT(tmp)
      ENDDO

! Initialize the simulation
      CALL INITIALIZE

! Time march loop.

      IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_COUPLED) THEN

! Uncoupled discrete element simulations
         IF (DEM_SOLIDS) CALL DES_TIME_MARCH
         IF (PIC_SOLIDS) CALL PIC_TIME_MARCH

      ELSE

! Transient or steady state simulation
         DO WHILE (TIME + 0.1d0*DT < TSTOP .AND. .NOT. EXIT_SIGNAL)
            CALL TIME_STEP_INIT
            DO
               CALL ITERATE_INIT
               DO WHILE (NIT<MAX_NIT .AND. .NOT.(CONVERGED.OR.DIVERGED))
                  NIT = NIT + 1
                  CALL DO_ITERATION
               ENDDO

               CALL POST_ITERATE

               IF(STEADY_STATE) EXIT
               IF(.NOT.ADJUSTDT()) EXIT
            ENDDO
            CALL TIME_STEP_END
            IF (STEADY_STATE) EXIT
         ENDDO

      ENDIF
      CALL FINALIZE

      STOP
      END PROGRAM MFIX
