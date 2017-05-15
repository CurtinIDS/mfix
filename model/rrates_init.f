!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES_INIT(IER)                                       C
!  Purpose: Initialize reaction rate arrays                            C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RRATES_INIT()

! Global Variables:
!---------------------------------------------------------------------//
      use energy,     only : HOR_g, HOR_s
      use param1,     only : ZERO
      use rxns,       only : R_gp, RoX_gc, SUM_R_g
      use rxns,       only : R_sp, RoX_sc, SUM_R_s
      use rxns,       only : R_PHASE, NO_OF_RXNS, USE_RRATES, RRATE
      use stiff_chem, only : STIFF_CHEMISTRY

      implicit none

! Local Variables:
!----------------------------------------------------------------------!
! NONE

! Flag for invoking reaction rate calculations.
      RRATE = .FALSE.
      IF(NO_OF_RXNS > 0)  RRATE = .TRUE.   ! automated mass balance
      IF(USE_RRATES)      RRATE = .TRUE.   ! legacy hook
      IF(STIFF_CHEMISTRY) RRATE = .FALSE.  ! stiff chemistry solver

! Gas phase source terms:
      R_gp    = ZERO  ! Rate of species formation
      RoX_gc  = ZERO  ! Rate of species consumption (divided X_g)
      SUM_R_G = ZERO  ! Net rate of gas formation/consumption
      HOR_G   = ZERO  ! Heat of reaction

! Solids phase source terms:
      R_sp    = ZERO  ! Rate of species formation
      RoX_sc  = ZERO  ! Rate of species consumption (divided X_s)
      SUM_R_S = ZERO  ! Net rate of solids formation/consumption
      HOR_S   = ZERO  ! Heat of reaction

! Interphase mass transfer.
      R_PHASE = ZERO

      RETURN
      END
