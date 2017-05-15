!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_PIC                                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Calculate the particle stress.                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_PIC


! Global Variables:
!---------------------------------------------------------------------//
! Flag to use Snider's particle stress model
      use mfix_pic, only: MPPIC_SOLID_STRESS_SNIDER
! Particle stress
      use mfix_pic, only: PIC_P_S

! Module procedures:
!---------------------------------------------------------------------//
      use sendrecv, only: SEND_RECV

      IMPLICIT NONE

!......................................................................!


      IF(MPPIC_SOLID_STRESS_SNIDER) THEN
         CALL CALC_PS_PIC_SNIDER
      ELSE
         CALL CALC_PS_PIC_GARG
      ENDIF

      CALL SEND_RECV(PIC_P_S,1)

      RETURN
      END SUBROUTINE CALC_PS_PIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_PIC                                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Evaluate the particle stress model of Snider.              !
!                                                                      !
!  REF: D.M. Snider, "Three-Dimensional Multiphase Particle-in-Cell    !
!     Model for Dense Particle Flows," Journal of Computational        !
!     Physics, Vol. 170, No. 2, pp. 523-549, 2001.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_PIC_SNIDER

! Global Variables:
!---------------------------------------------------------------------//
! Model parameters for Snider particle stress model
      use mfix_pic, only: FRIC_EXP_PIC
      use mfix_pic, only: PSFAC_FRIC_PIC
      use mfix_pic, only: FRIC_NON_SING_FAC
! Calculated particle stress
      use mfix_pic, only: PIC_P_S
! Fluid phase volume fraction
      use fldvar, only: EP_G
! Fluid volume fraction at close-pack
      use constant, only: EP_STAR
! Domain bounds
      use compar, only: IJKSTART3, IJKEND3
! Double precision parameters
      use param1, only: ONE

! Module procedures:
!---------------------------------------------------------------------//
      use functions, only: FLUID_AT
      use functions, only: WEST_OF, EAST_OF
      use functions, only: SOUTH_OF, NORTH_OF
      use functions, only: BOTTOM_OF, TOP_OF

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: IJK
! Volume fraction of cell, modified for wall cells.
      DOUBLE PRECISION :: lEPg
!......................................................................!

      DO IJK = IJKSTART3, IJKEND3

         IF(FLUID_AT(IJK)) THEN
            lEPg = EP_G(IJK)
         ELSE

! Set the volume fraction in the wall to close pack.
            lEPg = EP_STAR
 
! Use the lowest value across all adjacent fluid cells. This is to keep
! cells below close pack from pushing parcels through the walls.
!            lIJK = EAST_OF(IJK)
!            IF(FLUID_AT(lIJK))  lEPg = min(lEPg, EP_G(lIJK))
!!            lIJK = WEST_OF(IJK)
!            IF(FLUID_AT(lIJK))  lEPg = min(lEPg, EP_G(lIJK))
!            lIJK = NORTH_OF(IJK)
!            IF(FLUID_AT(lIJK))  lEPg = min(lEPg, EP_G(lIJK))
!!            lIJK = SOUTH_OF(IJK)
!            IF(FLUID_AT(lIJK))  lEPg = min(lEPg, EP_G(lIJK))
!            IF(DO_K) THEN
!               lIJK = TOP_OF(IJK)
!               IF(FLUID_AT(lIJK))  lEPg = min(lEPg, EP_G(lIJK))
!               lIJK = BOTTOM_OF(IJK)
!               IF(FLUID_AT(lIJK))  lEPg = min(lEPg, EP_G(lIJK))
!            ENDIF
         ENDIF

! Particle stress :: Snider (Eq 33)
         PIC_P_S(IJK,1) = PSFAC_FRIC_PIC *((ONE - lEPg)**FRIC_EXP_PIC)/&
            MAX(lEPg - EP_STAR, FRIC_NON_SING_FAC*lEPg)
      ENDDO



      RETURN
      END SUBROUTINE CALC_PS_PIC_SNIDER


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_PIC                                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Evaluate the particle stress as a coloring function:       !
!     X=0.0 :: cells below packing limit                               !
!     X=EPg :: cells above packing limit                               !
!     X=1.0 :: wall cells (far avobe max packing)                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_PIC_GARG


! Global Variables:
!---------------------------------------------------------------------//
! Calculated particle stress
      use mfix_pic, only: PIC_P_S
! Resulting particle stress force
      use mfix_pic, only: PS_FORCE_PIC
! Fluid phase volume fraction
      use fldvar, only: EP_G
! Fluid volume fraction at close-pack
      use constant, only: EP_STAR
! Domain bounds
      use compar, only: IJKSTART3, IJKEND3
! Double precision parameters
      use param1, only: ZERO, ONE

! Module procedures:
!---------------------------------------------------------------------//
      use functions, only: FLUID_AT


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: IJK
!......................................................................!


! The Garg model uses a coloring function approach. 
      DO IJK = IJKSTART3, IJKEND3
         PS_FORCE_PIC(:,IJK) = ZERO
         IF(FLUID_AT(IJK)) THEN
            IF(EP_G(IJK) < EP_STAR) THEN
               PIC_P_S(IJK,1) = (ONE - EP_G(IJK))
            ELSE
               PIC_P_S(IJK,1) = ZERO
            ENDIF
         ELSE
            PIC_P_S(IJK,1) = ONE
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_PS_PIC_GARG
