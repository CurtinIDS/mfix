!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VAVG_U_g                                               C
!  Purpose: Volume average U_g                                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-APR-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION VAVG_U_G ()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bc
      USE compar
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE utilities, ONLY: mfix_isnan
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          IJK
!
!                      Integral of U_g*EP_g for entire volume
      DOUBLE PRECISION SUM_U_g
!
!                      Total volume of computational cells
      DOUBLE PRECISION SUM_VOL

!  Integrate the velocity values for the whole domain,
!
      SUM_U_G = ZERO
      SUM_VOL = ZERO


!!!!$omp   parallel do private(IJK) reduction(+:SUM_VOL,SUM_U_G)
      DO IJK = IJKSTART3, IJKEND3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK), J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN
            SUM_VOL = SUM_VOL + VOL_U(IJK)
            SUM_U_G = SUM_U_G + U_G(IJK)*EP_G(IJK)*VOL_U(IJK)
         ENDIF
      END DO

      CALL GLOBAL_ALL_SUM(SUM_VOL)
      CALL GLOBAL_ALL_SUM(SUM_U_G)

      VAVG_U_G = SUM_U_G/SUM_VOL
!
! uncomment the following lines to enable trapping NaN's.
!      IF( mfix_isnan(VAVG_U_G) ) THEN
!        write(*,*) VAVG_U_G,  ' NaN being caught in VAVG_U_G.f '
!        AUTOMATIC_RESTART = .TRUE.
!      ENDIF
!
      RETURN
      END FUNCTION VAVG_U_G



      DOUBLE PRECISION FUNCTION VAVG_Flux_U_G ()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE run
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      USE compar
      USE mpi_utility
      USE mflux
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          IJK
!
!                      Integral of U_g*ROP_g*Area
      DOUBLE PRECISION SUM_U_g
!
!                      Total area of computational cells
      DOUBLE PRECISION SUM_AREA

!  Integrate the velocity values for the whole domain,
!
      SUM_U_G = ZERO
      SUM_AREA = ZERO


!!!!$omp   parallel do private(IJK) reduction(+:SUM_AREA,SUM_U_G)
      DO IJK = IJKSTART3, IJKEND3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK), J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN
           IF(.NOT.ADDED_MASS) THEN
              SUM_U_G = SUM_U_G + Flux_gE(IJK)
           ELSE
              SUM_U_G = SUM_U_G + Flux_gSE(IJK)
           ENDIF
           SUM_AREA = SUM_AREA + AYZ(IJK)
         ENDIF
      END DO

      CALL GLOBAL_ALL_SUM(SUM_AREA)
      CALL GLOBAL_ALL_SUM(SUM_U_G)

      VAVG_Flux_U_G = SUM_U_G/SUM_AREA
!
      RETURN
      END FUNCTION VAVG_Flux_U_G

