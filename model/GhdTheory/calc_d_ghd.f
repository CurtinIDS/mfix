!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_D_ghd_e(A_m, VxF_gs, d_e, IER)                    C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- East                                C
!                                                                      C
!  Note:  MFIX convention: center coeff is negative, hence:            C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The EAST face area is AYZ                                    C
!         Modifed for GHD theory                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_E calculation in accounting for   C
!           the averaged Solid-Solid drag and multi-solid-particle dragC
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
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
      SUBROUTINE CALC_D_ghd_E(A_M, VXF_GS, D_E)
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA                                        !S. Dartevelle, LANL, Feb.2004
!          Usual Indices
      INTEGER           M, L, I, IJK, IJKE                         !S. Dartevelle, LANL, Feb.2004
!          Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION  EPSA(DIMENSION_M)                          !S. Dartevelle, LANL, Feb.2004
!          ratio of drag and A0 or A_solid and other sum of Solid-Solid drag
      DOUBLE PRECISION  other_ratio_1, FOA1                        !S. Dartevelle, LANL, Feb.2004
!          sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION  SUM_VXF_GS                                 !S. Dartevelle, LANL, Feb.2004
!          What phase momentum equation is activated?
      LOGICAL  Pass1, Pass2                                        !S. Dartevelle, LANL, Feb.2004
!          numerator needed for solid phase M time EPSA
      DOUBLE PRECISION  numeratorxEP(DIMENSION_M)                  !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  denominator(DIMENSION_M)                   !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  other_denominator(DIMENSION_M)             !S. Dartevelle, LANL, Feb.2004
!          total solids volume fraction
      DOUBLE PRECISION  EPStmp
!-----------------------------------------------

   Pass1 = .FALSE.     !initialization
   Pass2 = .FALSE.
   M = MMAX
         if (MOMENTUM_X_EQ(0) .AND. MOMENTUM_X_EQ(M)) then
           Pass1 = .TRUE.    !we have at least one solid phase X-momentum equation
                             !with the gas phase X-momentum equation
         elseif (MOMENTUM_X_EQ(M)) then
           Pass2 = .TRUE.    !we have at least one solid phase X-momentum equation
                             !but the gas phase X-momentum is not solved
         endif


 IF (Pass1) THEN    !Gas and at least one solid phase X-momentum equation
!!!$omp   parallel do private(I,IJK, IJKE, EPGA, EPSA,EPStmp, numeratorxEP, &
!!!$omp&  M, L, Lp, LpL, LM, other_ratio_1, FOA1, &
!!!$omp&  SUM_VXF_GS, other_denominator, denominator ),&
!!!$omp&  schedule(static)
  DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
        DO M= 0, MMAX
         D_E(IJK,M) = ZERO
        END DO
     ELSE
        I = I_OF(IJK)
        IJKE = EAST_OF(IJK)
        EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

        SUM_VXF_GS = ZERO
        EPSA(MMAX) = ZERO
        DO M= 1, SMAX
          EPSA(M) = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
          EPSA(MMAX) = EPSA(MMAX) + EPSA(M)
          SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
        END DO

        M = MMAX
        other_ratio_1 = ( VXF_GS(IJK,M)* (-A_M(IJK,0,M)) /&
                                      ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SMALL_NUMBER )&
                        )

        IF (MODEL_B) THEN   !Model B
          !Linking velocity correction coefficient to pressure - GAS Phase
          if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
             D_E(IJK,0) = P_SCALE*AYZ(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
          else
             D_E(IJK,0) = ZERO
          endif
          !Linking velocity correction coefficient to pressure - SOLID Phase
          M = MMAX
            if ( MOMENTUM_X_EQ(M) .AND. (  ((-A_M(IJK,0,M))>SMALL_NUMBER) .OR. &
                                                     (VXF_GS(IJK,M)>SMALL_NUMBER) ) )   then
               D_E(IJK,M) = D_E(IJK,0)*(&
                                VXF_GS(IJK,M)/((-A_M(IJK,0,M))+VXF_GS(IJK,M))&
                            )
            else
               D_E(IJK,M) = ZERO
            endif
        ELSE                !Model A
          FOA1 = ZERO
          M = MMAX
            FOA1 = FOA1 + (EPSA(M)*VXF_GS(IJK,M)/&
                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SMALL_NUMBER )&
                          )
            other_denominator(M) = VXF_GS(IJK,M)*( (-A_M(IJK,0,0))/&
                                                   ((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)&
                                                 )
            numeratorxEP(M) = ZERO
            denominator(M)  = ZERO
          !Linking velocity correction coefficient to pressure - GAS Phase
          if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
            D_E(IJK,0) = P_SCALE*AYZ(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
          else
            D_E(IJK,0) = ZERO
          endif
          !Linking velocity correction coefficient to pressure - SOLID Phase
          M = MMAX
            if ( MOMENTUM_X_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER)        .OR. &
                                                (other_denominator(M)>SMALL_NUMBER) ) )     then
              D_E(IJK,M) = P_SCALE*AYZ(IJK)*(&
                        ( EPSA(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                        ( (-A_M(IJK,0,M))+other_denominator(M) )&
                                                   )
            else
              D_E(IJK,M) = ZERO
            endif
        ENDIF    !end of Model_B/Model_A if then condition
     ENDIF
  ENDDO

 ELSE IF (MOMENTUM_X_EQ(0)) THEN    !the solid X-momentum equations are not solved
                                    !only gas phase X-momentum
!!!$omp   parallel do private(IJK, I, IJKE, EPGA, M, SUM_VXF_GS), &
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
       D_E(IJK,0) = ZERO
     ELSE
       I = I_OF(IJK)
       IJKE = EAST_OF(IJK)
       EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)


       SUM_VXF_GS = VXF_GS(IJK,MMAX)              !Gas - All Solids VolxDrag summation

       IF ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR. (SUM_VXF_GS>SMALL_NUMBER) ) THEN
         IF (MODEL_B) THEN
           D_E(IJK,0) = P_SCALE*AYZ(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ELSE
           D_E(IJK,0) = P_SCALE*AYZ(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ENDIF
       ELSE
         D_E(IJK,0) = ZERO
       ENDIF
     ENDIF
   END DO

 ELSE IF (Pass2) THEN    !at least one solid phase X-momentum equation is solved
                         !but the gas phase X-momentum is not solved
!!!$omp    parallel do private(IJK, I, IJKE, EPSA,EPStmp, L, Lp, M, &
!!!$omp&   numeratorxEP, denominator), &
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_E(IJK,M) = ZERO
       END DO
     ELSE
       I = I_OF(IJK)
       IJKE = EAST_OF(IJK)

       M = MMAX
       EPStmp = ZERO
       DO L=1,SMAX
          EPStmp = EPStmp+ AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I)
       ENDDO
       EPSA(M) = EPStmp

       !Linking velocity correction coefficient to pressure - SOLID Phase (Model_A only)
       M = MMAX
         if ( MOMENTUM_X_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER) .OR. &
                                        (VXF_GS(IJK,M)>SMALL_NUMBER) ) ) then
           D_E(IJK,M) = P_SCALE*AYZ(IJK)*( EPSA(M) )/&
                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M) )
         else
           D_E(IJK,M) = ZERO
         endif
     ENDIF
   END DO

 ENDIF

 RETURN
 END SUBROUTINE CALC_D_ghd_E


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_D_ghd_n(A_m, VxF_gs, d_n, IER)                C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- North                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Note:  MFIX convention: center coeff is negative, hence             C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The NORTH face area is AXZ                                   C
!         Modifed for GHD theory                                       C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_N calculation in accounting for   C
!                       the averaged Solid-Solid drag                  C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
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
      SUBROUTINE CALC_D_ghd_N(A_M, VXF_GS, D_N)
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA                                        !S. Dartevelle, LANL, Feb.2004
!          Usual Indices
      INTEGER           M, L, J, IJK, IJKN                         !S. Dartevelle, LANL, Feb.2004
!          Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION  EPSA(DIMENSION_M)                          !S. Dartevelle, LANL, Feb.2004
!          ratio of drag and A0 or A_solid and other sum of Solid-Solid drag
      DOUBLE PRECISION  other_ratio_1, FOA1                        !S. Dartevelle, LANL, Feb.2004
!          sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION  SUM_VXF_GS                                 !S. Dartevelle, LANL, Feb.2004
!          What phase momentum equation is activated?
      LOGICAL  Pass1, Pass2                                        !S. Dartevelle, LANL, Feb.2004
!          numerator needed for solid phase M time EPSA
      DOUBLE PRECISION  numeratorxEP(DIMENSION_M)                  !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  denominator(DIMENSION_M)                   !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  other_denominator(DIMENSION_M)             !S. Dartevelle, LANL, Feb.2004
!          total solids volume fraction
      DOUBLE PRECISION  EPStmp
!-----------------------------------------------

   Pass1 = .FALSE.     !initialization
   Pass2 = .FALSE.
   M = MMAX
     if (MOMENTUM_Y_EQ(0) .AND. MOMENTUM_Y_EQ(M)) then
       Pass1 = .TRUE.    !we have at least one solid phase Y-momentum equation
                         !with the gas phase Y-momentum equation
     elseif (MOMENTUM_Y_EQ(M)) then
       Pass2 = .TRUE.    !we have at least one solid phase Y-momentum equation
                         !but the gas phase Y-momentum is not solved
     endif


 IF (Pass1) THEN    !Gas and at least one solid phases Y-momentum equation
!!!$omp   parallel do private(J,IJK, IJKN, EPGA, EPSA, EPStmp, numeratorxEP, &
!!!$omp&  M, L, Lp, LpL, LM, other_ratio_1, FOA1, &
!!!$omp&  SUM_VXF_GS, other_denominator, denominator ),&
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN   !impermeable
       DO M= 0, MMAX
         D_N(IJK,M) = ZERO
       END DO
     ELSE
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)
       EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

       SUM_VXF_GS = ZERO
       EPSA(MMAX) = ZERO
       DO M= 1, SMAX
         EPSA(M) = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
         EPSA(MMAX) = EPSA(MMAX) + EPSA(M)
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
       END DO

       M= MMAX
       other_ratio_1 = ( VXF_GS(IJK,M)* (-A_M(IJK,0,M)) /&
                                      ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SMALL_NUMBER )&
                      )

       IF (MODEL_B) THEN   !Model B
          !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
           D_N(IJK,0) = P_SCALE*AXZ(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
         else
           D_N(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         M = MMAX
           if ( MOMENTUM_Y_EQ(M) .AND. (  ((-A_M(IJK,0,M))>SMALL_NUMBER) .OR. &
                                                     (VXF_GS(IJK,M)>SMALL_NUMBER) ) )   then
             D_N(IJK,M) = D_N(IJK,0)*(&
                                       VXF_GS(IJK,M)/((-A_M(IJK,0,M))+VXF_GS(IJK,M))&
                                               )
           else
             D_N(IJK,M) = ZERO
           endif
       ELSE                !Model A
         FOA1 = ZERO
         M = MMAX
           FOA1 = FOA1 + (EPSA(M)*VXF_GS(IJK,M)/&
                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SMALL_NUMBER )&
                          )
           other_denominator(M) = VXF_GS(IJK,M)*( (-A_M(IJK,0,0))/&
                                                   ((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)&
                                                 )
           numeratorxEP(M) = ZERO
           denominator(M)  = ZERO
         !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
           D_N(IJK,0) = P_SCALE*AXZ(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
         else
           D_N(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         M = MMAX
           if ( MOMENTUM_Y_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER)        .OR. &
                                                (other_denominator(M)>SMALL_NUMBER) ) )     then
             D_N(IJK,M) = P_SCALE*AXZ(IJK)*(&
                        ( EPSA(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                        ( (-A_M(IJK,0,M))+other_denominator(M) )&
                                                   )
           else
             D_N(IJK,M) = ZERO
           endif
       ENDIF    !end of Model_B/Model_A if then condition
     ENDIF
   END DO

 ELSE IF (MOMENTUM_Y_EQ(0)) THEN    !the solid Y-momentum equations are not solved
                                    !only the gas phase
!!!$omp   parallel do private(IJK, J, IJKN, EPGA, M, SUM_VXF_GS), &
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN   !impermeable
       D_N(IJK,0) = ZERO
     ELSE
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)
       EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)


       SUM_VXF_GS = VXF_GS(IJK,MMAX)              !Gas - All Solids VolxDrag summation

       IF ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR. (SUM_VXF_GS>SMALL_NUMBER) ) THEN
         IF (MODEL_B) THEN
           D_N(IJK,0) = P_SCALE*AXZ(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ELSE
           D_N(IJK,0) = P_SCALE*AXZ(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ENDIF
       ELSE
         D_N(IJK,0) = ZERO
       ENDIF
     ENDIF
   END DO

 ELSE IF (Pass2) THEN    !at least one solid phase momentum Y-equation is solved
                         !but the gas phase Y-momentum is not solved
!!!$omp    parallel do private(IJK, J, IJKN, EPSA, EPStmp, L, Lp, M, &
!!!$omp&   numeratorxEP, denominator), &
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_N(IJK,M) = ZERO
       END DO
     ELSE
       J = J_OF(IJK)
       IJKN = NORTH_OF(IJK)

       M = MMAX
       EPStmp = ZERO
       DO L=1,SMAX
         EPStmp = EPStmp+ AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J)
       ENDDO
       EPSA(M) = EPStmp

       !Linking velocity correction coefficient to pressure - SOLID Phase (Model_A only)
       M = MMAX
         if ( MOMENTUM_Y_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER) .OR. &
                                        (VXF_GS(IJK,M)>SMALL_NUMBER) ) ) then
           D_N(IJK,M) = P_SCALE*AXZ(IJK)*( EPSA(M) )/&
                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M) )
         else
           D_N(IJK,M) = ZERO
         endif
     ENDIF
   END DO

 ENDIF

 RETURN
 END SUBROUTINE CALC_D_ghd_N


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_D_ghd_t(A_m, VxF_gs, d_t, IER)                C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- Top                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Note:  MFIX convention: center coeff is negative, hence             C
!                            (-A_M(IJK,0,M)) > or = 0                  C
!         The TOP face area is AXY                                     C
!         Modifed for GHD theory                                       C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_T calculation in accounting for   C
!                       the averaged Solid-Solid drag                  C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: more accurate formulas in the D_s(M) for Model_A           C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 17-MAR-04  C
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
      SUBROUTINE CALC_D_ghd_T(A_M, VXF_GS, D_T)
!
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA                                        !S. Dartevelle, LANL, Feb.2004
!          Usual Indices
      INTEGER           M, L, K, IJK, IJKT                         !S. Dartevelle, LANL, Feb.2004
!          Average solid volume fraction at momentum cell centers
      DOUBLE PRECISION  EPSA(DIMENSION_M)                          !S. Dartevelle, LANL, Feb.2004
!          ratio of drag and A0 or A_solid and other sum of Solid-Solid drag
      DOUBLE PRECISION  other_ratio_1, FOA1                        !S. Dartevelle, LANL, Feb.2004
!          sum of All Solid-Gas VolxDrag
      DOUBLE PRECISION  SUM_VXF_GS                                 !S. Dartevelle, LANL, Feb.2004
!          What phase momentum equation is activated?
      LOGICAL  Pass1, Pass2                                        !S. Dartevelle, LANL, Feb.2004
!          numerator needed for solid phase M time EPSA
      DOUBLE PRECISION  numeratorxEP(DIMENSION_M)                  !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  denominator(DIMENSION_M)                   !S. Dartevelle, LANL, Feb.2004
!          denominator needed for solid phase M
      DOUBLE PRECISION  other_denominator(DIMENSION_M)             !S. Dartevelle, LANL, Feb.2004
!          tmp variable for total solids volume fraction
      DOUBLE PRECISION  EPStmp
!-----------------------------------------------

   Pass1 = .FALSE.     !initialization
   Pass2 = .FALSE.
   M = MMAX
     if (MOMENTUM_Z_EQ(0) .AND. MOMENTUM_Z_EQ(M)) then
       Pass1 = .TRUE.    !we have at least one solid phase Z-momentum equation
                         !with the gas phase Z-momentum equation
     elseif (MOMENTUM_Z_EQ(M)) then
       Pass2 = .TRUE.    !we have at least one solid phase Z-momentum equation
                         !but the gas phase Z-momentum is not solved
     endif


 IF (Pass1) THEN    !Gas and at least one solid phases Z-momentum equation
!!!$omp   parallel do private(K,IJK, IJKT, EPGA, EPSA, EPStmp, numeratorxEP, &
!!!$omp&  M, L, Lp, LpL, LM, other_ratio_1, FOA1,  &
!!!$omp&  SUM_VXF_GS, other_denominator, denominator ),&
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN   !impermeable
       DO M= 0, MMAX
         D_T(IJK,M) = ZERO
       END DO
     ELSE
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)
       EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

       SUM_VXF_GS = ZERO
       EPSA(MMAX) = ZERO
       DO M= 1, SMAX
         EPSA(M) = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
         EPSA(MMAX) = EPSA(MMAX) + EPSA(M)
         SUM_VXF_GS = SUM_VXF_GS + VXF_GS(IJK,M)              !Gas - All Solids VolxDrag summation
       END DO

       M = MMAX
       other_ratio_1 = ( VXF_GS(IJK,M)* (-A_M(IJK,0,M)) /&
                                      ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SMALL_NUMBER )&
                      )

       IF (MODEL_B) THEN   !Model B
         !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
           D_T(IJK,0) = P_SCALE*AXY(IJK)/( (-A_M(IJK,0,0))+other_ratio_1 )
         else
           D_T(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         M = MMAX
           if ( MOMENTUM_Z_EQ(M) .AND. (  ((-A_M(IJK,0,M))>SMALL_NUMBER) .OR. &
                                                     (VXF_GS(IJK,M)>SMALL_NUMBER) ) )   then
             D_T(IJK,M) = D_T(IJK,0)*(&
                                       VXF_GS(IJK,M)/((-A_M(IJK,0,M))+VXF_GS(IJK,M))&
                                               )
           else
             D_T(IJK,M) = ZERO
           endif
       ELSE                !Model A
         FOA1 = ZERO
         M = MMAX
           FOA1 = FOA1 + (EPSA(M)*VXF_GS(IJK,M)/&
                             ( (-A_M(IJK,0,M))+VXF_GS(IJK,M)+SMALL_NUMBER )&
                          )
           other_denominator(M) = VXF_GS(IJK,M)*( (-A_M(IJK,0,0))/&
                                                   ((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)&
                                                 )
           numeratorxEP(M) = ZERO
           denominator(M)  = ZERO
         !Linking velocity correction coefficient to pressure - GAS Phase
         if ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR.  (other_ratio_1>SMALL_NUMBER) ) then
           D_T(IJK,0) = P_SCALE*AXY(IJK)*(EPGA+FOA1)/( (-A_M(IJK,0,0))+other_ratio_1 )
         else
           D_T(IJK,0) = ZERO
         endif
         !Linking velocity correction coefficient to pressure - SOLID Phase
         M = MMAX
           if ( MOMENTUM_Z_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER)        .OR. &
                                                (other_denominator(M)>SMALL_NUMBER) ) )     then
             D_T(IJK,M) = P_SCALE*AXY(IJK)*(&
                        ( EPSA(M) + (VXF_GS(IJK,M)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS+SMALL_NUMBER)) )/&
                        ( (-A_M(IJK,0,M))+other_denominator(M) )&
                                                   )
           else
             D_T(IJK,M) = ZERO
           endif
       ENDIF    !end of Model_B/Model_A if then condition
     ENDIF
   END DO

 ELSE IF (MOMENTUM_Z_EQ(0)) THEN    !the solid Z-momentum equations are not solved
                                    !only the gas phase Z-momentum is solved
!!!$omp   parallel do private(IJK, K, IJKT, EPGA, M, SUM_VXF_GS), &
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN   !impermeable
       D_T(IJK,0) = ZERO
     ELSE
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)
       EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)

       SUM_VXF_GS = VXF_GS(IJK,MMAX)              !Gas - All Solids VolxDrag summation

       IF ( ((-A_M(IJK,0,0))>SMALL_NUMBER) .OR. (SUM_VXF_GS>SMALL_NUMBER) ) THEN
         IF (MODEL_B) THEN
           D_T(IJK,0) = P_SCALE*AXY(IJK)/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ELSE
           D_T(IJK,0) = P_SCALE*AXY(IJK)*EPGA/((-A_M(IJK,0,0))+SUM_VXF_GS)
         ENDIF
       ELSE
         D_T(IJK,0) = ZERO
       ENDIF
     ENDIF
   END DO

 ELSE IF (Pass2) THEN    !at least one solid phase momentum Z-equation is solved
                         !but the gas phase Z-momentum is not solved
!!!$omp    parallel do private(IJK, K, IJKT, EPSA, EPStmp, L, Lp, M, &
!!!$omp&   numeratorxEP, denominator), &
!!!$omp&  schedule(static)
   DO IJK = ijkstart3, ijkend3
     IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK) .OR. MODEL_B) THEN
       DO M= 1, MMAX
         D_T(IJK,M) = ZERO
       END DO
     ELSE
       K = K_OF(IJK)
       IJKT = TOP_OF(IJK)

       M = MMAX
       EPStmp = ZERO
       DO L=1,SMAX
         EPStmp = EPStmp+ AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)
       ENDDO
       EPSA(M) = EPStmp

       !Linking velocity correction coefficient to pressure - SOLID Phase (Model_A only)
       M = MMAX
         if ( MOMENTUM_Z_EQ(M) .AND. ( (-A_M(IJK,0,M)>SMALL_NUMBER) .OR. &
                                        (VXF_GS(IJK,M)>SMALL_NUMBER) ) ) then
            D_T(IJK,M) = P_SCALE*AXY(IJK)*( EPSA(M) )/&
                                          ( (-A_M(IJK,0,M))+VXF_GS(IJK,M) )
         else
           D_T(IJK,M) = ZERO
         endif
     ENDIF
   END DO

 ENDIF

 RETURN
 END SUBROUTINE CALC_D_ghd_T

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
