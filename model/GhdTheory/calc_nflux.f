!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Module name: CALC_NFLUX(IER)                                        C
!  Purpose: Calculate the convection fluxes based on total species     C
!           number density for use with GHD theory only.               C
!           Adapted from routine calc_mflux.                           C
!                                                                      C
!  Author: S. Benyahia                                Date: 31-DEC-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!      Christine Hrenya handnotes and my own derivation                C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_NFLUX()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE physprop
      USE constant
      USE indices
      USE mflux
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!
!                      Total Number density
      DOUBLE PRECISION Ni_E(DIMENSION_3), Ni_N(DIMENSION_3), Ni_T(DIMENSION_3)
!
!                      Indices
      INTEGER          IJK, IMJK, IJMK, IJKM, L
!
!                      Particle mass
      DOUBLE PRECISION Mp_L
!
!  First compute total number density at faces
!!!$omp  parallel do private( IJK, IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
        Ni_E(IJK) = ZERO
        Ni_N(IJK) = ZERO
        Ni_T(IJK) = ZERO
!
        IF (FLUID_AT(IJK)) THEN
           DO L = 1, SMAX
             MP_L = RO_S(IJK,L) * PI*D_P(IJK,L)**3/6d0
             Ni_E(IJK) = Ni_E(IJK) + ROP_sE(IJK,L) / MP_L
             IMJK = IM_OF(IJK)
             IF (.NOT.FLUID_AT(IMJK)) THEN
                Ni_E(IMJK) = Ni_E(IMJK) + ROP_sE(IMJK,L) / MP_L
             ENDIF
             Ni_N(IJK) = Ni_N(IJK) + ROP_sN(IJK,L) / MP_L
             IJMK = JM_OF(IJK)
             IF (.NOT.FLUID_AT(IJMK)) THEN
               Ni_N(IJMK) = Ni_N(IJMK) + ROP_sN(IJMK,L) / MP_L
             ENDIF
!
             IF (DO_K) THEN
               Ni_T(IJK) = Ni_T(IJK) + ROP_sT(IJK,L) / MP_L
               IJKM = KM_OF(IJK)
                 IF (.NOT.FLUID_AT(IJKM)) THEN
                   Ni_T(IJKM) = Ni_T(IJKM) + ROP_sT(IJKM,L) / MP_L
                 ENDIF
             ENDIF
           ENDDO
        ENDIF
      ENDDO

!
!  Interpolate the face value of density for calculating the convection fluxes
!!!$omp  parallel do private( IJK, IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN
!
!         East face (i+1/2, j, k)
            Flux_nE(IJK) = Ni_E(IJK)*AYZ(IJK)*U_s(IJK,MMAX)
!
!         North face (i, j+1/2, k)
            Flux_nN(IJK) = Ni_N(IJK)*AXZ(IJK)*V_s(IJK,MMAX)
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN
               Flux_nT(IJK) = Ni_T(IJK)*AXY(IJK)*W_s(IJK,MMAX)
            ENDIF
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               Flux_nE(IMJK) = Ni_E(IMJK)*AYZ(IMJK)*U_s(IMJK,MMAX)
            ENDIF
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLUID_AT(IJMK)) THEN
              Flux_nN(IJMK) = Ni_N(IJMK)*AXZ(IJMK)*V_s(IJMK,MMAX)
            ENDIF
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                 Flux_nT(IJKM) = Ni_T(IJKM)*AXY(IJKM)*W_s(IJKM,MMAX)
               ENDIF
            ENDIF
         ENDIF
      END DO

      RETURN
      END SUBROUTINE CALC_NFLUX

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
