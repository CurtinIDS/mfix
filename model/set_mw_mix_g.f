!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_MW_MIX_g                                            C
!  Purpose: calculate gas mixture molecular weights                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-OCT-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2, X_g                                  C
!  Variables modified: MW_MIX_g                                        C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_MW_MIX_G

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE geometry
      USE fldvar
      USE constant
      USE indices
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: IJK
!-----------------------------------------------

      IF (MW_AVG /= UNDEFINED) RETURN

!!$omp parallel do private(ijk) &
!!$omp schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
! calculate mw_mix_g in all fluid and flow boundary cells
! set_bc0 will have already defined mw_mix_g in MI and PI boundary cells
! (redundant-remove in set_bc0?)
         IF (.NOT.WALL_AT(IJK)) MW_MIX_G(IJK) = &
            CALC_MW(X_G,DIMENSION_3,IJK,NMAX(0),MW_G)
      ENDDO

      RETURN
      END SUBROUTINE SET_MW_MIX_G


