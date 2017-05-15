!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_gs_Y                                                 C
!  Purpose: Calculate the average drag coefficient at i, j+1/2, k and  C
!           multiply with v-momentum cell volume.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE VF_GS_Y(VXF_GS)

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
! Volume x Drag
      DOUBLE PRECISION, INTENT(OUT) :: VxF_gs(DIMENSION_3, DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: J, IJK, IJKN
! Index of continuum solids phases
      INTEGER :: M
!-----------------------------------------------


      DO M = 1, SMAX
         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_N(IJK)) THEN
               VXF_GS(IJK,M) = ZERO
            ELSE
               J = J_OF(IJK)
               IJKN = NORTH_OF(IJK)
               VXF_GS(IJK,M) = VOL_V(IJK) *                            &
                  AVG_Y(F_GS(IJK,M),F_GS(IJKN,M),J)
            ENDIF
         ENDDO
      ENDDO


      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         DO IJK = IJKSTART3, IJKEND3
            IF(IP_AT_N(IJK)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               J = J_OF(IJK)
               IJKN = NORTH_OF(IJK)
               VXF_GDS(IJK) = VOL_V(IJK) *                             &
                  AVG_Y(F_GDS(IJK),F_GDS(IJKN),J)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE VF_GS_Y


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_SS_Y                                                 C
!  Purpose: Calculate the average Solid-Solid drag coefficient at      C
!           i, j+1/2, k and multiply with V-momentum cell volume.      C
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

      SUBROUTINE VF_SS_Y(VXF_SS)

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
      INTEGER :: J, IJK, IJKN
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
                  IF (.NOT.IP_AT_N(IJK)) THEN
                     J = J_OF(IJK)
                     IJKN = NORTH_OF(IJK)
                     VXF_SS(IJK,LM) = AVG_Y(F_SS(IJK,LM),F_SS(IJKN,LM),J)*VOL_V(IJK)
                  ELSE     !Impermeable wall
                     VXF_SS(IJK,LM) = ZERO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO   ! end do loop (l=1,mmax)
      ENDDO   ! end do loop (m=1,mmax)

      IF (DES_CONTINUUM_HYBRID) THEN
! initialize every call
         DO M = 1, MMAX
            DO IJK = ijkstart3, ijkend3
               IF (IP_AT_N(IJK)) THEN
                  VXF_SDS(IJK,M) = ZERO
               ELSE
                  J = J_OF(IJK)
                  IJKN = NORTH_OF(IJK)
                  VXF_SDS(IJK,M) = AVG_Y(F_SDS(IJK,M),F_SDS(IJKN,M),J)*VOL_V(IJK)
               ENDIF
            ENDDO      ! end do loop (ijk=ijkstart3,ijkend3)
         ENDDO   ! end do loop (m=1,mmax)
      ENDIF   ! end if (discrete_element and des_continuum_hybrid)

      RETURN
      END SUBROUTINE VF_SS_Y


