!``````````````````````````````````````````````````````````````````````!
!  Module: rdf                                                         !
!                                                                      !
!  Calculate radial distribution functions.                            !
!  Note that routines G_0AVG, G_0, and DG_0DNU need to be modified to  !
!  effect a change in the radial distribution function g_0. The old    !
!  routine G_0EP has been replaced with G_0CS, which used only for     !
!  Carnahan-Starling g_0.                                              !
!......................................................................!

MODULE rdf

CONTAINS

!``````````````````````````````````````````````````````````````````````!
!  Function: G_0AVG                                                    !
!  Author: M. Syamlal                                 Date: 05-JAN-05  !
!                                                                      !
!......................................................................!
      DOUBLE PRECISION FUNCTION G_0AVG (IJK1, IJK2, DIR, L, M1, M2)

      USE compar
      USE constant
      use discretelement, only: des_mmax
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE visc_s
      USE run
      USE toleranc

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Cell indices for neighboring cells (passed variable)
      INTEGER, INTENT(IN) :: IJK1, IJK2
! Direction (X, Y, or Z) (passed variable)
      CHARACTER, INTENT(IN) :: DIR
! Direction index (i, j, or k)
      INTEGER, INTENT(IN) :: L
! Solids phase index-1
      INTEGER, INTENT(IN) :: M1
! Solids phase index-2
      INTEGER, INTENT(IN) :: M2

! Local Variables:
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER :: Mx
! Average solids volume fraction
      DOUBLE PRECISION :: EPS
! Average void fraction
      DOUBLE PRECISION :: EPg, EPg_STAR_AVG
! Sum over m (EP_sm/D_pm)
      DOUBLE PRECISION :: EPSoDP
! Sum over EP_s
      DOUBLE PRECISION :: SUM_EPS
! Average number density of solids phase mm
      DOUBLE PRECISION :: NU_MM
! Volume of a single particle of solids phase mm
      DOUBLE PRECISION :: VOLP
! Quantity employed in rdf calculation
      DOUBLE PRECISION :: XI
! Average D_P for phase M1 and M2
      DOUBLE PRECISION :: DP_AVG_M1, DP_AVG_M2, DP_AVG

      IF(IJK1 == IJK2)THEN
         G_0AVG = G_0(IJK1, M1, M2)
      ELSE
         SELECT CASE(RDF_TYPE_ENUM)

! Lebowitz, J.L. (1964) The Physical Review, A133, 895-899
!---------------------------------------------------------------------//
         CASE(LEBOWITZ)
            DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
            DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
            EPSoDP = ZERO

            DO Mx = 1, DES_MMAX+MMAX
               EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
               EPSoDP = EPSoDP + 2d0*EPS / (D_p(IJK1,Mx)+D_p(IJK2,Mx))
            ENDDO

            EPg = AVG_XYZ(EP_g(IJK1), EP_g(IJK2), DIR, L)
            G_0AVG = ONE / EPg + 3.0d0*EPSoDP*DP_AVG_M1*DP_AVG_M2 /    &
               (EPg*EPg *(DP_AVG_M1 + DP_AVG_M2))


! REF: ?
!---------------------------------------------------------------------//
         CASE(MODIFIED_LEBOWITZ)
            DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
            DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
            EPSoDP = ZERO
            SUM_EPS = ZERO

            DO Mx = 1, DES_MMAX+MMAX
               EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
               DP_AVG = HALF*( D_p(IJK1,Mx)+D_p(IJK2,Mx) )
               EPSoDP = EPSoDP + EPS/DP_AVG
               SUM_EPS = SUM_EPS + EPS
            ENDDO

            EPg_STAR_AVG = AVG_XYZ(EP_star_array(IJK1),                &
               EP_star_array(IJK2), DIR, L)


! Prevent G_0 from becoming negative when the solids volume fraction
! approachs a maximum packing.  This may occur during non-converged
! iterations when the local solids volume fraction exceeds the maximum.
            IF(SUM_EPS >= (ONE-EPg_STAR_AVG))                          &
               SUM_EPS = SUM_EPS - DIL_EP_s


            G_0AVG = (ONE/(ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) + 3.0d0 *  &
               ((DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2))*EPSoDP



! Mansoori, GA, Carnahan N.F., Starling, K.E. Leland, T.W. (1971).
! The Journal of Chemical Physics, Vol. 54:1523-1525.
!---------------------------------------------------------------------//
! Extended Carnahan & Starling see Garzo & Dufty (1999) for details.
         CASE(MANSOORI)

            DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
            DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
            SUM_EPS = ZERO
            XI = ZERO

            DO Mx = 1, DES_MMAX+MMAX
               EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
               SUM_EPS = SUM_EPS + EPS

               DP_AVG = HALF*( D_p(IJK1,Mx)+D_p(IJK2,Mx) )
               VOLP = (PI/6.0D0)*DP_AVG**3

               IF(DP_AVG > ZERO) THEN
                  NU_MM = EPS/VOLP
                  XI = XI + NU_MM*DP_AVG*DP_AVG
               ENDIF
            ENDDO

            XI = (PI/6.0D0)*XI

            G_0AVG = (ONE/(ONE-SUM_EPS)) + (3.0D0)*&
               ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) )*        &
               ( XI/((ONE-SUM_EPS)*(ONE-SUM_EPS)) ) + (2.0D0) *        &
               ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) *       &
               ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) *       &
               ( (XI*XI)/((ONE-SUM_EPS)*(ONE-SUM_EPS)*(ONE-SUM_EPS)))


! van Wachem, B.G.M., Schouten, J.C., van den Bleek, C.M., Krishna, R.
! and Sinclair, J. L. (2001). AIChE Journal 47:1035–1051.
!---------------------------------------------------------------------//
         CASE(MODIFIED_MANSOORI)

            DP_AVG_M1 = HALF*(D_p(IJK1,M1)+D_p(IJK2,M1))
            DP_AVG_M2 = HALF*(D_p(IJK1,M2)+D_p(IJK2,M2))
            SUM_EPS = ZERO
            XI = ZERO

            DO Mx = 1, DES_MMAX+MMAX
               EPS = AVG_XYZ(EP_s(IJK1, Mx), EP_s(IJK2, Mx), DIR, L)
               SUM_EPS = SUM_EPS + EPS

               DP_AVG = HALF*( D_p(IJK1,Mx)+D_p(IJK2,Mx) )
               VOLP = (PI/6.0D0)*DP_AVG**3

               IF (DP_AVG > ZERO) THEN
                  NU_MM = EPS/VOLP
                  XI = XI + NU_MM*DP_AVG*DP_AVG
               ENDIF
            ENDDO

            XI = (PI/6.0D0)*XI

            EPg_STAR_AVG = AVG_XYZ(EP_star_array(IJK1),                &
               EP_star_array(IJK2), DIR, L)

            IF(SUM_EPS >= (ONE-EPg_STAR_AVG))                          &
               SUM_EPS = SUM_EPS - DIL_EP_s

            G_0AVG = (ONE/(ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) + (3.0D0)* &
               ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) )*        &
               ( XI/((ONE-SUM_EPS/(ONE-EPg_STAR_AVG))*                 &
               (ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) ) + (2.0D0) *         &
               ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) *       &
               ( (DP_AVG_M1*DP_AVG_M2)/(DP_AVG_M1+DP_AVG_M2) ) *       &
               ( (XI*XI)/((ONE-SUM_EPS/(ONE-EPg_STAR_AVG))*            &
               (ONE-SUM_EPS/(ONE-EPg_STAR_AVG))*                       &
               (ONE-SUM_EPS/(ONE-EPg_STAR_AVG))) )


! Carnahan, N.F. and Starling K.E. (1969).
! The Journal of Chemical Physics, Vol. 51(2):635-636.
         CASE(CARNAHAN_STARLING)
!---------------------------------------------------------------------//

! Do not use when there are more than one granular phase.
            EPS = AVG_XYZ(EP_s(IJK1, M1), EP_s(IJK2, M1), DIR, L)
            G_0AVG = G_0CS(EPS)
         END SELECT

      ENDIF

      RETURN
      END FUNCTION G_0AVG


!``````````````````````````````````````````````````````````````````````!
!  Function: G_0 (IJK, M1, M2)                                         !
!  Author: M. Syamlal                                 Date: 16-MAR-92  !
!                                                                      !
!  Purpose: Calculate radial distribution function at contact for a    !
!           mixture of spheres of different diameters                  !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Lebowitz, J.L., "Exact solution of generalized Percus-Yevick      !
!      equation for a mixture of hard spheres," The Physical Review,   !
!      A133, 895-899 (1964).                                           !
!                                                                      !
!    Iddir, Y.H., "Modeling of the multiphase mixture of particles     !
!         using the kinetic theory approach," PhD Thesis, Illinois     !
!         Institute of Technology, Chicago, Illinois, 2004:            !
!         chapter 2, equations 2-49 through 2-52.                      !
!                                                                      !
!    Mansoori et al. (1971)                                            !
!       This RDF expression is equivalent to that cited by Jenkins and !
!          Mancini (1987) & Garzo and Dufty (1999).                    !
!......................................................................!
      DOUBLE PRECISION FUNCTION G_0 (IJK, M1, M2)

      USE compar
      USE constant
      use discretelement, only: des_mmax
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE toleranc
      USE visc_s

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Solids phase index-1 (passed variable)
      INTEGER, INTENT(IN) :: M1
! Solids phase index-2 (passed variable)
      INTEGER, INTENT(IN) :: M2
! Fluid cell index
      INTEGER, INTENT(IN) :: IJK

! Local Variables:
!---------------------------------------------------------------------//
! Local solids phase index
      INTEGER :: Mx, MM
! Average solids volume fraction
      DOUBLE PRECISION :: EPS
! Average void fraction
      DOUBLE PRECISION :: EPg
! Sum over m (EP_sm/D_pm)
      DOUBLE PRECISION :: EPSoDP
! Sum over EP_s
      DOUBLE PRECISION :: SUM_EPS
! Number density of solids phase mm
      DOUBLE PRECISION :: NU_MM
! Volume of a single particle of solids phase mm
      DOUBLE PRECISION :: VOLP
! Quantity employed in rdf calculation
      DOUBLE PRECISION :: XI
!---------------------------------------------------------------------//

      SUM_EPS = ZERO
      EPg = EP_G(IJK)
      DO MM = 1, DES_MMAX+MMAX
          EPS = EP_s(IJK, MM)
          SUM_EPS = SUM_EPS + EPS
      END DO

      SELECT CASE(RDF_TYPE_ENUM)

! Lebowitz, J.L. (1964) The Physical Review, A133, 895-899
!---------------------------------------------------------------------//
      CASE(LEBOWITZ)

         EPSoDP = ZERO
         DO Mx = 1, DES_MMAX+MMAX
            EPS = EP_s(IJK, Mx)
            EPSoDP = EPSoDP + EPS / D_p(IJK,Mx)
         ENDDO
         EPg = EP_g(IJK)

         G_0 = ONE/EPg + 3.0d0 * EPSoDP * D_p(IJK,M1) * D_p(IJK,M2) /  &
            (EPg*EPg *(D_p(IJK,M1) + D_p(IJK,M2)))

! REF: ?
!---------------------------------------------------------------------//
      CASE(MODIFIED_LEBOWITZ)

         EPSoDP = ZERO
         SUM_EPS = ZERO

         DO MM = 1, DES_MMAX+MMAX
            EPS = EP_s(IJK, MM)
            EPSoDP = EPSoDP + (EPS/D_p(IJK,MM))
            SUM_EPS = SUM_EPS + EPS
         END DO

! Prevent G_0 from becoming negative when the solids volume fraction
! approachs a maximum packing.  This may occur during non-converged
! iterations when the local solids volume fraction exceeds the maximum.
         IF(SUM_EPS >= (ONE-EP_star_array(IJK)))                       &
            SUM_EPS = SUM_EPS - DIL_EP_s

         G_0 = (ONE/(ONE-SUM_EPS/(ONE-EP_star_array(IJK)) )) + 3.0d0 * &
            ((D_p(IJK,M1)*D_p(IJK,M2))/(D_p(IJK,M1)+D_p(IJK,M2)))*EPSoDP

! Mansoori, GA, Carnahan N.F., Starling, K.E. Leland, T.W. (1971).
! The Journal of Chemical Physics, Vol. 54:1523-1525.
!---------------------------------------------------------------------//
! Extended Carnahan & Starling see Garzo & Dufty (1999) for details.
      CASE(MANSOORI)

         SUM_EPS = ZERO
         XI = ZERO

         DO MM = 1, DES_MMAX+MMAX
            EPS = EP_s(IJK, MM)
            SUM_EPS = SUM_EPS + EPS
            VOLP = (PI/6.0D0)*D_P(IJK,MM)**3.0

            IF (D_P(IJK,MM) > ZERO) THEN
               NU_MM = EPS/VOLP
               XI = XI + NU_MM*D_P(IJK,MM)*D_P(IJK,MM)
            ENDIF
         ENDDO
         XI = (PI/6.0D0)*XI

         G_0 = (ONE/(ONE-SUM_EPS)) + (3.0D0)*&
            ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) )*   &
            ( XI/((ONE-SUM_EPS)*(ONE-SUM_EPS)) ) + (2.0D0) *           &
            ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) *  &
            ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) *  &
            ( (XI*XI)/((ONE-SUM_EPS)*(ONE-SUM_EPS)*(ONE-SUM_EPS)) )

! van Wachem, B.G.M., Schouten, J.C., van den Bleek, C.M., Krishna, R.
! and Sinclair, J. L. (2001). AIChE Journal 47:1035–1051.
!---------------------------------------------------------------------//
         CASE(MODIFIED_MANSOORI)

         SUM_EPS = ZERO
         XI = ZERO

         DO MM = 1, DES_MMAX+MMAX
            EPS = EP_s(IJK, MM)
            SUM_EPS = SUM_EPS + EPS
            VOLP = (PI/6.0D0)*D_P(IJK,MM)**3.0

            IF (D_P(IJK,MM) > ZERO) THEN
               NU_MM = EPS/VOLP
               XI = XI + NU_MM*D_P(IJK,MM)*D_P(IJK,MM)
            ENDIF
         ENDDO
         XI = (PI/6.0D0)*XI

         IF(SUM_EPS >= (ONE-EP_star_array(IJK)) )                      &
            SUM_EPS = SUM_EPS - DIL_EP_s

         G_0 = (ONE/(ONE-SUM_EPS/(ONE-EP_star_array(IJK)))) + (3.0D0)* &
            ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) )*   &
            ( XI/((ONE-SUM_EPS/(ONE-EP_star_array(IJK)))*              &
            (ONE-SUM_EPS/(ONE-EP_star_array(IJK)))) ) + (2.0D0) *      &
            ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) *  &
            ( (D_P(IJK,M1)*D_P(IJK,M2))/(D_P(IJK,M1)+D_P(IJK,M2)) ) *  &
            ( (XI*XI)/((ONE-SUM_EPS/(ONE-EP_star_array(IJK)))*         &
            (ONE-SUM_EPS/(ONE-EP_star_array(IJK)))*                    &
            (ONE-SUM_EPS/(ONE-EP_star_array(IJK)))) )

! Carnahan, N.F. and Starling K.E. (1969).
! The Journal of Chemical Physics, Vol. 51(2):635-636.
      CASE(CARNAHAN_STARLING)
!---------------------------------------------------------------------//
         G_0 = G_0CS(EP_S(IJK,M1))

      END SELECT

      RETURN
      END FUNCTION G_0


!``````````````````````````````````````````````````````````````````````!
!  Module name: DG_0DNU (EPs)                                          !
!  Author: K. Agrawal                                 Date: 16-FEB-98  !
!                                                                      !
!  Purpose: Calculate derivative of radial distribution function at    !
!           w.r.t granular volume fraction                             !
!......................................................................!
      DOUBLE PRECISION FUNCTION DG_0DNU (EPS)

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: ZERO, ONE
      USE physprop, only: MMAX

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: EPS

      DG_0DNU = ZERO

! Carnahan-Starling derivative of (G0) wrt EP_s:
! This value is only needed for monodisper simulatoins.
      IF(MMAX == 1) DG_0DNU = (2.5D0-EPS)/(ONE - EPS)**4

      RETURN
      END FUNCTION DG_0DNU

!``````````````````````````````````````````````````````````````````````!
!  Module name: G_0CS(EPs)                                             !
!                                                                      !
!  Purpose: Carnahan-Starling radial distribution function.            !
!......................................................................!
      DOUBLE PRECISION FUNCTION G_0CS (EPS)

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: ONE

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Solids volume fraction
      DOUBLE PRECISION :: EPS

      G_0CS = (ONE-0.5D0*EPS)/(ONE - EPS)**3

      RETURN
      END FUNCTION G_0CS

!``````````````````````````````````````````````````````````````````````!
!  Function: AVG_XYZ                                                   !
!                                                                      !
!  Purpose: Calculate the arithmetic average in X, Y, or Z direction.  !
!......................................................................!
      DOUBLE PRECISION FUNCTION AVG_XYZ (V1, V2, DIR, L)

      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      USE functions
      USE fun_avg

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Direction (X, Y, or Z) (passed variable)
      CHARACTER, INTENT(IN) :: DIR
! Direction index (i, j, or k)
      INTEGER, INTENT(IN) :: L
! The variables to be averaged
      DOUBLE PRECISION, INTENT(IN) :: V1, V2


! Local Variables:
!---------------------------------------------------------------------//
! Fluid cell index

      IF(DIR == 'X')THEN
        AVG_XYZ = AVG_X(V1, V2, L)

      ELSEIF(DIR == 'Y')THEN
        AVG_XYZ = AVG_Y(V1, V2, L)

      ELSEIF(DIR == 'Z')THEN
        AVG_XYZ = AVG_Z(V1, V2, L)

      ELSE
        CALL WRITE_ERROR('AVG_XYZ', 'Unkown direction', 1)
      ENDIF

      RETURN
      END FUNCTION AVG_XYZ

    END MODULE rdf
