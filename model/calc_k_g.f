!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_K_g                                                C
!  Purpose: Calculate the effective conductivity of fluid phase        C
!                                                                      C
!                                                                      C
!  Comments:                                                           C
!  This routine will not be called if k_g0 is defined                  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_K_G()

! Modules
!---------------------------------------------------------------------//
      USE param1, only: undefined
      USE physprop, only: k_g0, k_g
      USE sendrecv, only: send_recv
! invoke user defined quantity
      USE usr_prop, only: usr_kg, calc_usr_prop
      USE usr_prop, only: gas_conductivity
      IMPLICIT NONE
!---------------------------------------------------------------------//

      IF (USR_Kg) THEN
         CALL CALC_USR_PROP(Gas_Conductivity,lm=0)
      ELSEIF (K_g0 == UNDEFINED) THEN
! unncessary check but included for clarity
         CALL CALC_DEFAULT_Kg
      ENDIF

      CALL send_recv(K_G, 2)

      RETURN
      END SUBROUTINE CALC_K_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute the default value for gas conductivity where the   C
!  gas phase is assumed to be air                                      C
!  Author:M. Syamlal                                  Date: 24-APR-96  C
!                                                                      C
!  Literature/Document References:                                     C
!  Bird, Stewart, and Lightfoot (1960) --                              C
!    Temperature dependence from formula 8.3-12 on p. 255 and          C
!    conductivity value at 300 K from p. 263                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DEFAULT_Kg

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      USE fldvar, only: T_g
      USE functions, only: fluid_at
      USE param1, only: zero
      USE physprop, only: K_g
      USE run, only: units
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
!---------------------------------------------------------------------//

!!$omp parallel do private(ijk) &
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
! Gas conductivity (air) in cal/(s.cm.K)
            K_G(IJK) = 6.02D-5*SQRT(T_G(IJK)/300.D0)
         ELSE
            K_G(IJK) = ZERO
         ENDIF

! 1 cal = 4.183925D0 J
         IF (UNITS == 'SI') K_G(IJK) = 418.3925D0*K_G(IJK)      !J/s.m.K

      ENDDO

      RETURN
      END SUBROUTINE CALC_DEFAULT_Kg
