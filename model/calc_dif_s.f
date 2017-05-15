!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DIF_S                                              C
!  Purpose: Calculate the effective diffusivity of solids phases       C
!                                                                      C
!                                                                      C
!  Comments:                                                           C
!  This routine will not be called if dif_s0(M) is defined             C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DIF_S(M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: undefined
      USE physprop, only: dif_s0, dif_s
      USE sendrecv, only: send_recv
! invoke user defined quantity
      USE usr_prop, only: usr_difs, calc_usr_prop
      USE usr_prop, only: solids_diffusivity
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER, INTENT(IN) :: M

!---------------------------------------------------------------------//

      IF (USR_Difs(M)) THEN
         CALL CALC_USR_PROP(Solids_Diffusivity,lm=M)
      ELSEIF (Dif_s0(M) == UNDEFINED) THEN
         CALL CALC_DEFAULT_DIF_SOLIDS(M)
      ENDIF

      CALL SEND_RECV(DIF_S, 2)

      RETURN
      END SUBROUTINE CALC_DIF_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DEFAULT_DIF_SOLIDS                                 C
!  Purpose: Compute the default value of each solids phases effective  C
!  diffusivity. Because species are not considered to diffuse in a     C
!  solids the default diffusivsity is zero.                            C
!                                                                      C
!  Author:M. Syamlal                                  Date: 13-EFB-98  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DEFAULT_DIF_SOLIDS(M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      USE fldvar, only: rop_s
      USE param1, only: zero
      USE physprop, only: nmax, dif_s
      USE run, only: units
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: N
! binary diffusivity coefficient
      DOUBLE PRECISION :: Dab
!---------------------------------------------------------------------//

      Dab = ZERO       !cm^2/s
      IF(UNITS == 'SI') Dab = Dab*0.0001D0   !m^2/s

!!$omp  parallel do private(n,ijk) &
!!$omp& schedule(dynamic,chunk_size)
      DO N = 1, NMAX(M)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               DIF_S(IJK,M,N) = ROP_S(IJK,M)*Dab
            ELSE
               DIF_S(IJK,M,N) = ZERO
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_DEFAULT_DIF_SOLIDS



