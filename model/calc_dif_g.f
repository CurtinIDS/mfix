!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DIF_g                                              C
!  Purpose: Calculate the effective diffusivity of fluid phase         C
!                                                                      C
!                                                                      C
!  Comments:                                                           C
!  This routine will not be called if dif_g is defined                 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DIF_G()

! Modules
!---------------------------------------------------------------------//
      USE param1, only: undefined
      USE physprop, only: dif_g0, dif_g
      USE sendrecv, only: send_recv
! invoke user defined quantity
      USE usr_prop, only: usr_difg, calc_usr_prop
      USE usr_prop, only: gas_diffusivity
      IMPLICIT NONE
!---------------------------------------------------------------------//

      IF (USR_Difg) THEN
         CALL CALC_USR_PROP(Gas_Diffusivity,lm=0)
      ELSEIF (Dif_g0 == UNDEFINED) THEN
! unncessary check but included for clarity
         CALL CALC_DEFAULT_DIF_GAS
      ENDIF

      CALL send_recv(DIF_G, 2)

      RETURN
      END SUBROUTINE CALC_DIF_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DEFAULT_DIF_GAS                                    C
!  Purpose: Compute the default value for diffusivity of each gas      C
!  species; each species is assigned the same value with the base      C
!  value corresponding to CO2 in N2 at T=298K and P~=1atm.             C
!                                                                      C
!  Author:M. Syamlal                                  Date: 13-FEB-98  C
!                                                                      C
!  Revision: Include dilute mixture approximation for calculation of   C
!  multicomponent diffusion coefficients                               C
!  Author:N. Reuge                                    Date: 11-APR-07  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Fuller relation - to account for influence of gas temperature       C
!     and pressure                                                     C
!  Curtiss-Hirschfelder, Wilke & Blanc - dilute mixture approximation  C
!     for multicomponent diffusion. Valid if the mass fraction of the  C
!     carrier species > 0.9                                            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DEFAULT_DIF_GAS

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use fldvar, only: ROP_g, T_g, P_g, X_g
      use functions, only: fluid_at
      use param1, only: zero, one
      use physprop, only: NMAX, Dif_g
      use scales, only: unscale_pressure
      use toleranc, only: zero_x_gs
      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
! Species index
      INTEGER :: N, L
! Binary diffusion coefficient
      DOUBLE PRECISION :: Dab(NMAX(0),NMAX(0))
! Reference temperature and pressure for each species diffusion
! coefficient
      DOUBLE PRECISION :: Tg_ref(NMAX(0))
      DOUBLE PRECISION :: Pg_ref(NMAX(0))
! Intermediate calculation to determine weighted average diffusion
! coefficient
      DOUBLE PRECISION :: lDab, Sum_XgLoDab, Sum_XgL, lXgN
!---------------------------------------------------------------------//

      CALL SET_BINARY_DAB_GAS(Dab, Tg_ref, Pg_ref)

!!$omp  parallel do private(ijk) &
!!$omp& schedule(dynamic,chunk_size)
      DO N = 1, NMAX(0)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN

               IF (NMAX(0) == 1) THEN
! Only 1 species present in gas phase
                  lDab = Dab(N,N)

               ELSE
! Dilute mixture approximation for multi-component diffusion

                  SUM_XgLoDab = ZERO
                  SUM_XgL = ZERO
                  DO L = 1, NMAX(0)
                     IF (L /= N) THEN
                        SUM_XgLoDab = SUM_XgLoDab+&
                                     X_g(IJK,L)/Dab(N,L)
! Sum_XgL = 1-XgN
                        SUM_XgL = SUM_XgL + X_g(IJK,L)
                     ENDIF
                  ENDDO
! it is possible for lXgN to evaluate <0
                  lXgN = ONE-SUM_XgL

                  IF (lXgN > ZERO_X_GS .AND. &
                     SUM_XgL > ZERO_X_GS) THEN
! i.e. when cell is not only species N
! If this criteria is too strict (i.e. when XgN->1), then this section
! may be evaluated when the other XgN do not carry significant value
! (i.e. are effectively zero). As a result, lDab may be evaluated and
! become unrealistically large. This may happen when happen when
! 1-XgN != sum_XgL. Generally, sum_XgL should equal 1-XgN but they may
! differ due to the numerical evaluation of each Xg. If XgN->1, then
! use both sum_XgL and 1-XgN to determine whether calculations should
! proceed. Given the noted limitation of this approximation then an
! additional criteria should probably be added to only evaluate when
! sum_xgL <0.1.
                     IF (SUM_XgLoDab > ZERO) THEN
! for numerical reasons use Sum_XgL rather than ONE-X_g(IJK,N)
                        lDab = SUM_XgL/SUM_XgLoDab
                     ELSE
! this should not occur...
                        lDab = Dab(N,N)
                     ENDIF
                  ELSE
! Address case when the mass fraction of the Nth species is nearly 1.
                     lDab = Dab(N,N)
                  ENDIF
               ENDIF

! Influence of gas temperature and gas pressure from Fuller relation
               DIF_G(IJK,N) = ROP_G(IJK)*lDab* &
                           (T_g(IJK)/Tg_ref(N))**1.75* &
                           Pg_ref(N)/UNSCALE_PRESSURE(P_g(IJK))
            ELSE
               DIF_G(IJK,N) = ZERO
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_DEFAULT_DIF_GAS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BINARY_DAB_GAS                                      C
!  Purpose: Set binary diffusion coefficients for all cross species    C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Bird, Stewart and lightfoot (1960)                                  C
!  Reid, Prausnitz and Poling (1987)                                   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BINARY_DAB_GAS(lDab, lTg_ref, lPg_ref)

! Modules
!---------------------------------------------------------------------//
      use physprop, only: NMAX
      use run, only: units
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Binary diffusion coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDab(NMAX(0),NMAX(0))
      DOUBLE PRECISION, INTENT(OUT) :: lTg_ref(NMAX(0))
      DOUBLE PRECISION, INTENT(OUT) :: lPg_ref(NMAX(0))

! Local Variables
!---------------------------------------------------------------------//
! Species index
      INTEGER :: N, L
!---------------------------------------------------------------------//


! Default gas diffusion coefficient
! Bird, Stewart, and Lightfoot (1960) -- CO2--N2 at 298.2 K
      DO N = 1,NMAX(0)
         DO L = N,NMAX(0)
!         DO L = N+1, NMAX(0)
!            IF (N /= L) THEN
! assign the diaganol case for when nmax = 1 or when only species N present
! in given cell.
               lDab(N,L) = 0.165D0       !cm^2/s
               lDab(L,N) = lDab(N,L)
!            ENDIF
         ENDDO
         lTg_ref(N) = 298.2d0
         lPg_ref(N) = 1.01D6             !dyne
      ENDDO

! Gas diffusion coefficients for 3 species system: SiH4, H2 & N2
! Calculated using relation derived from Chapman and Enskog's theory
! of gases - Reid, Prausnitz and Poling (1987)
! Binary diffusion coefficient SiH4/H2 at 873 K
!      lDab(1,2) = 3.78      ! cm^2/s
!      lDab(2,1) = lDab(1,2)
! Binary diffusion coefficient SiH4/N2 at 873 K
!      lDab(1,3) = 1.02      ! cm^2/s
!      lDab(3,1) = lDab(1,3)
! Binary diffusion coefficient H2/N2 at 873 K
!      lDab(2,3) = 4.52      ! cm^2/s
!      lDab(3,2) = lDab(2,3)
!      DO N = 1,NMAX(0)
!         lTg_ref(N) = 873.0
!         lPg_ref(N) = 1.01e6         ! dyne
!      ENDDO

      IF(UNITS == 'SI') THEN
         DO N = 1,NMAX(0)
            DO L = N,NMAX(0)
!               IF (N /= L) THEN
! assign the diaganol case nmax = 1 or when only species N present
! in given cell.
                  lDab(N,L) = lDab(N,L)*0.0001D0   !m^2/s
                  lDab(L,N) = lDab(N,L)
!               ENDIF
            ENDDO
            lPg_ref(N) = lPg_ref(N)/10.D0          !Pa
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SET_BINARY_DAB_GAS
