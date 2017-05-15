! TODO:
!   p_star calculation should be based on the sum of volume fractions of
!   close-packed solids.

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_P_star                                             C
!  Purpose: Calculate P_star in cells where solids continuity is       C
!     solved                                                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified: P_STAR,                                         C
!     if yu_standish or fedors_landel: ep_star_array,                  C
!                                      ep_g_blend_start,               C
!                                      ep_g_blend_end                  C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_P_STAR(EP_G, P_STAR)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE pgcor
      USE pscor
      USE ur_facs
      USE residual
      USE compar
      USE run
      USE visc_s
      USE solids_pressure
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EP_g(DIMENSION_3)
! Solids pressure
      DOUBLE PRECISION, INTENT(INOUT) :: P_star(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!!   HPF$ align P_star(:) with TT(:)
!!   HPF$ align EP_g(:) with TT(:)

! Indices
      INTEGER :: IJK
! Blend factor
      DOUBLE PRECISION :: blend
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION :: CALC_EP_STAR
!-----------------------------------------------

!!$omp parallel do private(ijk)
!!   HPF$ independent

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN


            IF (YU_STANDISH .OR. FEDORS_LANDEL) THEN
! if Yu_Standish or Fedors_Landel correlations are used, then
! ep_star_array is modified. this is the only time ep_star_array is
! modified (see set_constprop).  (sof Nov-16-2005)
               EP_star_array(ijk) = calc_ep_star(ijk)

! now the values of ep_g_blend_start and ep_g_blend_end need to be
! reassigned based on the new values of ep_star_array
               IF(BLENDING_STRESS.AND.TANH_BLEND) THEN
                  ep_g_blend_start(ijk) = ep_star_array(ijk) * 0.99d0
                  ep_g_blend_end(ijk)   = ep_star_array(ijk) * 1.01d0
               ELSEIF(BLENDING_STRESS.AND.SIGM_BLEND) THEN
                  ep_g_blend_start(ijk) = EP_star_array(ijk) * 0.97d0
                  ep_g_blend_end(ijk) = EP_star_array(ijk) * 1.01d0
               ELSE
                  ep_g_blend_start(ijk) = ep_star_array(ijk)
                  ep_g_blend_end(ijk)   = ep_star_array(ijk)
               ENDIF
            ENDIF

            IF (EP_G(IJK) < EP_g_blend_end(ijk)) THEN
               P_STAR(IJK) = NEG_H(EP_G(IJK),EP_g_blend_end(ijk))
               IF(BLENDING_STRESS) THEN
                  blend =  blend_function(IJK)
                  P_STAR (IJK) = (1.0d0-blend) * P_STAR (IJK)
               ENDIF
            ELSE
               P_STAR(IJK) = ZERO
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_P_STAR


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  FUNCTION: CALC_ep_star                                              C
!  Purpose: calculate the local value of maximum packing               C
!                                                                      C
!  Author: D. Gera and M. Syamlal                     Date: 31-DEC-02  C
!  Reviewer:                                          Date:            C
!  Modified: S. Benyahia                              Date: 02-May-05  C
!                                                                      C
!  Literature/Document References:                                     C
!    A.B. Yu and N. Standish. Powder Tech, 52 (1987) 233-241           C
!    R.F. Fedors and R.F. Landel. Powder Tech, 23 (1979) 225-231       C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      DOUBLE PRECISION FUNCTION CALC_ep_star(IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE toleranc
      USE compar
      USE run
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! IJK index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J

!      DOUBLE PRECISION :: xbar

! start sof modifications (02-May-05)
! maximum packing for the mixture
       DOUBLE PRECISION :: P_IT(MMAX)
! true maximum packing for the mixture
       DOUBLE PRECISION :: EPs_max_local
! maximum packing fraction for a binary mixture
       DOUBLE PRECISION :: P_IJ(MMAX, MMAX)
! particle diameter ratio
       DOUBLE PRECISION :: R_IJ(MMAX, MMAX)
! fractional solids volume corresponding to P_IJ
       DOUBLE PRECISION :: X_IJ(MMAX, MMAX)
! fractional solids volume in a mixture
! this is Xj in eq. 22 of Yu-Standish
       DOUBLE PRECISION :: COMP_X_I(MMAX), SUM_LOCAL
! local aliases for particle diameter, solids volume fraction and the
! maximum solids volume fraction which are used to rearrange solids
! phases from coarsest to finest
       DOUBLE PRECISION :: DP_TMP(MMAX), EPs_TMP(MMAX), &
                           EPs_max_TMP(MMAX), old_value
!-----------------------------------------------


      IF (CALL_DQMOM) THEN
! sort particles to start from coarsest to finest particles
! assigning values to local aliases
         DO I = 1, SMAX
            DP_TMP(I) = D_P(IJK,I)
            EPs_TMP(I) = EP_s(IJK,I)
            EPs_max_TMP(I) = ep_s_max(I)
         ENDDO
! sorting particles from coarse to fine
         DO I = 1, SMAX
            DO J = I , SMAX
! check if phase J is larger than phase I
               IF(DP_TMP(I) < DP_TMP(J)) THEN
! temporarily store phase i diameter
                  old_value = DP_TMP(I)
! overwrite phase i diameter with smaller phase j diameter
                  DP_TMP(I) = DP_TMP(J)
! overwrite phase j diameter with stoired phase i diameter
                  DP_TMP(J) = old_value

                  old_value = EPs_TMP(I)
                  EPs_TMP(I) = EPs_TMP(J)
                  EPs_TMP(J) = old_value

                  old_value = EPs_max_TMP(I)
                  EPs_max_TMP(I) = EPs_max_TMP(J)
                  EPs_max_TMP(J) = old_value
               ENDIF
            ENDDO
         ENDDO

      ELSE  ! not dqmom

! assigning values to local aliases
         DO I = 1, SMAX
            DP_TMP(I) = D_P(IJK,M_MAX(I))
            EPs_TMP(I) = EP_s(IJK,M_MAX(I))
            EPs_max_TMP(I) = ep_s_max(M_MAX(I))
         ENDDO
      ENDIF   ! end if/else (call_dqmom)


! this is the way the algorithm was written by Yu and Standish (sof).
! compute equations 25 in Yu-Standish
! (this is also needed by Fedors_Landel)
      DO I = 1, SMAX
         SUM_LOCAL = ZERO
         DO J = 1, SMAX
            IF(I .GE. J) THEN
               R_IJ(I,J) = DP_TMP(I)/DP_TMP(J)
            ELSE
               R_IJ(I,J) = DP_TMP(J)/DP_TMP(I)
            ENDIF
            SUM_LOCAL = SUM_LOCAL + EPs_TMP(J)
         ENDDO   ! end do (j=1,smax)

         IF(SUM_LOCAL > DIL_EP_s) THEN
! fractional solids volume see eq. 20
            COMP_X_I(I) = EPs_TMP(I)/SUM_LOCAL
         ELSE
! return first phase ep_s_max in case very dilute
            CALC_EP_star = ONE - EPs_max_TMP(1)
            RETURN
         ENDIF
      ENDDO   ! end do (i=1,smax)

! Begin YU_STANDISH section
! ---------------------------------------------------------------->>>
      IF(YU_STANDISH) THEN
! compute equation 23-24 in Yu-Standish
         DO I = 1, SMAX
            DO J = 1, SMAX
               IF(R_IJ(I,J) .LE. 0.741d0) THEN
                  IF(J .LT. I) THEN
                     X_IJ(I,J) = (ONE - R_IJ(I,J)*R_IJ(I,J))/&
                                 (2.0d0 -  EPs_max_TMP(I))
                  ELSE
                     X_IJ(I,J) = ONE - (ONE - R_IJ(I,J)*R_IJ(I,J))/&
                                 (2.0d0 -  EPs_max_TMP(I))
                 ENDIF
                 P_IJ(I, J) = EPs_max_TMP(I) + EPs_max_TMP(I)*&
                    (ONE-EPs_max_TMP(I)) * (ONE - 2.35d0*R_IJ(I,J) + &
                    1.35d0*R_IJ(I,J)*R_IJ(I,J) )
               ELSE
                  P_IJ(I, J) = EPs_max_TMP(I)
               ENDIF
            ENDDO   ! end do (j=1,smax)
         ENDDO   ! end do (i=1,smax)

! Compute equation 22
         EPs_max_local = ONE
         DO I = 1, SMAX
            SUM_LOCAL = ZERO

            IF(I .GE. 2) THEN
               DO J = 1, (I-1)
                  IF(P_IJ(I,J) == EPs_max_TMP(I)) THEN
                     SUM_LOCAL = SUM_LOCAL
                  ELSE
                     SUM_LOCAL = SUM_LOCAL + (ONE - EPs_max_TMP(I)/&
                        P_IJ(I,J))*COMP_X_I(J)/X_IJ(I,J)
                  ENDIF
               ENDDO
            ENDIF

            IF((I+1) .LE. SMAX) THEN
               DO J = (I+1), SMAX
                  IF( P_IJ(I, J) == EPs_max_TMP(I) ) THEN
                     SUM_LOCAL = SUM_LOCAL
                  ELSE
                     SUM_LOCAL = SUM_LOCAL + (ONE - EPs_max_TMP(I)/&
                        P_IJ(I, J))*COMP_X_I(J)/X_IJ(I, J)
                  ENDIF
               ENDDO
            ENDIF

            IF (SUM_LOCAL .NE. ZERO) THEN
               P_IT(I) = EPs_max_TMP(I)/(ONE - SUM_LOCAL)
            ELSE
! do nothing if particles have same diameter
               P_IT(I) = ONE
            ENDIF

            EPs_max_local = MIN(P_IT(I), EPs_max_local)
         ENDDO   ! end do (i=1,smax)

! for the case of all phases having same diameter

         IF (EPs_max_local == ONE) EPs_max_local = EPs_max_TMP(1)
         CALC_EP_star = ONE - EPs_max_local
! end YU_STANDISH section
! ----------------------------------------------------------------<<<

! Part implemented by Dinesh for binary mixture, uncomment to use (Sof)

!       if ((EP_s(IJK,1)+EP_s(IJK,2)) .NE. ZERO) THEN
!          xbar = EP_s(IJK,1)/(EP_s(IJK,1)+EP_s(IJK,2))

!          if (xbar .LE. ep_s_max_ratio(1,2)) THEN
!             CALC_EP_star =MAX(0.36d0, (ONE-(((ep_s_max(1)-ep_s_max(2))+&
!              (ONE-d_p_ratio(1,2))*(ONE-ep_s_max(1))*ep_s_max(2))*(ep_s_max(1)+&
!              (ONE-ep_s_max(1)) *ep_s_max(2))*xbar/ep_s_max(1)+ep_s_max(2))))
!          else
!             CALC_EP_star =MAX(0.36d0, (ONE-((ONE -d_p_ratio(1,2))*(ep_s_max(1)&
!              +(ONE-ep_s_max(1))*ep_s_max(2))*(ONE -xbar) +ep_s_max(1))))
!          end if
!       else
!          CALC_EP_star = ONE - MIN(ep_s_max(1), ep_s_max(2)) !corrected by sof
!       end if

! Use the code (below) instead of the above commented code because the
! phases were not rearranged and I didn't want to modify it (sof)
! If you don't understand what's going on, contact me: sof@fluent.com

! In the case of binary mixture (Fedors-Landel empirical correlation)
! ---------------------------------------------------------------->>>
      ELSEIF(FEDORS_LANDEL) THEN

         IF(COMP_X_I(1) .LE. (EPs_max_TMP(1)/(EPs_max_TMP(1)+ &
                             (ONE - EPs_max_TMP(1))*EPs_max_TMP(2))) ) THEN

            CALC_EP_star = (EPs_max_TMP(1) - EPs_max_TMP(2) + &
               (1 - sqrt(R_IJ(2,1))) * (ONE - EPs_max_TMP(1)) * &
               EPs_max_TMP(2) )*&
               (EPs_max_TMP(1) + (ONE - EPs_max_TMP(1)) * &
               EPs_max_TMP(2)) * COMP_X_I(1)/EPs_max_TMP(1) + &
               EPs_max_TMP(2)
         ELSE
            CALC_EP_star = (ONE-sqrt(R_IJ(2,1))) * (EPs_max_TMP(1)+&
               (ONE-EPs_max_TMP(1)) * EPs_max_TMP(2)) * &
               (ONE - COMP_X_I(1)) + EPs_max_TMP(1)
         ENDIF
! this is gas volume fraction at packing
         CALC_EP_star = ONE - CALC_EP_star
       ENDIF ! for Yu_Standish and Fedors_Landel correlations
! end FEDORS_LANDEL correlation
! ----------------------------------------------------------------<<<

      RETURN
      END FUNCTION CALC_ep_star
