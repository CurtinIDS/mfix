!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: VF_gs_X                                                 !
!  Author: M. Syamlal                                 Date: 20-MAY-96  !
!                                                                      !
!  Purpose: Calculate the average drag coefficient at i+1/2, j, k and  !
!           multiply with u-momentum cell volume.                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE VF_GS_X(VXF_GS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE discretelement
      USE drag
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Volume x Drag
      DOUBLE PRECISION, INTENT(OUT) :: VxF_gs(DIMENSION_3, DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, IJK, IJKE
! Index of continuum solids phases
      INTEGER :: M
!-----------------------------------------------

      DO M = 1, SMAX
         DO IJK = IJKSTART3, IJKEND3
            IF(IP_AT_E(IJK)) THEN
               VXF_GS(IJK,M) = ZERO
            ELSE
               I = I_OF(IJK)
               IJKE = EAST_OF(IJK)
               VXF_GS(IJK,M) = VOL_U(IJK) * &
                  AVG_X(F_GS(IJK,M),F_GS(IJKE,M),I)
            ENDIF
         ENDDO
      ENDDO

! Calculate the combined effect for all discrete solids.
      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         DO IJK = IJKSTART3, IJKEND3
            IF(IP_AT_E(IJK)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               I = I_OF(IJK)
               IJKE = EAST_OF(IJK)
               VXF_GDS(IJK) = VOL_U(IJK) *                             &
                  AVG_X(F_GDS(IJK),F_GDS(IJKE),I)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE VF_GS_X

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_SS_X                                                 C
!  Purpose: Calculate the average Solid-Solid drag coefficient at      C
!           i+1/2, j, k and multiply with u-momentum cell volume.      C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE VF_SS_X(VXF_SS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE drag
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Volume x Drag
      DOUBLE PRECISION, INTENT(OUT) :: VxF_SS(DIMENSION_3, DIMENSION_LM)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, IJK, IJKE
! Index of continuum solids phases
      INTEGER :: L, M, LM
!-----------------------------------------------

! initialize every call
      VXF_SS(:,:) = ZERO

      DO M = 1, MMAX
         DO L = 1, MMAX
            LM = FUNLM(L,M)
            IF (L .NE. M) THEN
               DO IJK = ijkstart3, ijkend3
                  IF (.NOT.IP_AT_E(IJK)) THEN
                     I = I_OF(IJK)
                     IJKE = EAST_OF(IJK)
                     VXF_SS(IJK,LM) = AVG_X(F_SS(IJK,LM),F_SS(IJKE,LM),I)*VOL_U(IJK)
                  ELSE     !Impermeable wall
                     VXF_SS(IJK,LM) = ZERO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO   ! end do loop (l=1,mmax)
      ENDDO   ! end do loop (m=1,mmax

      IF (DES_CONTINUUM_HYBRID) THEN
         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF (IP_AT_E(IJK)) THEN
                  VXF_SDS(IJK,M) = ZERO
               ELSE
                  I = I_OF(IJK)
                  IJKE = EAST_OF(IJK)
                  VXF_SDS(IJK,M) = AVG_X(F_SDS(IJK,M),F_SDS(IJKE,M),I)*VOL_U(IJK)
               ENDIF
            ENDDO   ! end do loop (dm=1,des_mmax)
         ENDDO   ! end do loop (m=1,mmax)
      ENDIF   ! end if (discrete_element and des_continuum_hybrid)

      RETURN
      END SUBROUTINE VF_SS_X

