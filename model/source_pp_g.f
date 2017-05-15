!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_Pp_g                                             C
!  Purpose: Determine source terms for pressure correction equation.   C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_PP_G(A_M, B_M, B_MMAX)

! Modules
!---------------------------------------------------------------------//
      USE bc, ONLY: IJK_P_G
      USE compar, ONLY: IJKSTART3, IJKEND3
      USE cutcell, ONLY: CARTESIAN_GRID, A_UPG_E, A_VPG_N, A_WPG_T
      USE fldvar, ONLY: U_G, V_G, W_G, U_S, V_S, W_S, ROP_G, ROP_S
      USE fldvar, only: ROP_GO, ROP_SO, RO_G
      USE geometry, ONLY: VOL
      USE param, only: dimension_3, dimension_m
      USE param, only: east, west, south, north, top, bottom
      USE param1, ONLY: SMALL_NUMBER, ONE, ZERO, UNDEFINED
      USE pgcor, ONLY: D_E, D_N, D_T
      USE physprop, ONLY: CLOSE_PACKED, MMAX, RO_G0
      USE run, ONLY: SHEAR, ODT, UNDEFINED_I
      USE rxns, ONLY: SUM_R_G, SUM_R_S
      USE vshear, ONLY: VSH
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)

! Local Variables
!---------------------------------------------------------------------//
! solids phase index
      INTEGER :: M
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! terms of bm expression
      DOUBLE PRECISION bma, bme, bmw, bmn, bms, bmt, bmb, bmr
! error message
      CHARACTER(LEN=80) :: LINE(1)
!---------------------------------------------------------------------//
! loezos
! update to true velocity
      IF (SHEAR) THEN
!!$omp parallel do private(IJK)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               V_G(IJK)=V_G(IJK)+VSH(IJK)
            ENDIF
         ENDDO
      ENDIF

! Calculate convection-diffusion fluxes through each of the faces

!$omp parallel default(none) &
!$omp          private(IJK, IMJK, IJMK, IJKM, M, bma,bme,bmw,bmn,&
!$omp                  bms,bmt,bmb,bmr,line)  &
!$omp          shared(ijkstart3,ijkend3,cartesian_grid,rop_g,rop_go,&
!$omp                 rop_s,rop_so,vol,odt,u_g,v_g,w_g,u_s,v_s,w_s,b_m, &
!$omp                 b_mmax,d_e,d_n,d_t,a_m,a_upg_e,a_vpg_n,a_wpg_t,&
!$omp                 mmax,close_packed,sum_r_s,sum_r_g,ro_g0)
!$omp do
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            bma = (ROP_G(IJK)-ROP_GO(IJK))*VOL(IJK)*ODT
            bme = A_M(IJK,east,0)*U_G(IJK)
            bmw = A_M(IJK,west,0)*U_G(IMJK)
            bmn = A_M(IJK,north,0)*V_G(IJK)
            bms = A_M(IJK,south,0)*V_G(IJMK)
            bmt = A_M(IJK,top,0)*W_G(IJK)
            bmb = A_M(IJK,bottom,0)*W_G(IJKM)
            bmr = SUM_R_G(IJK)*VOL(IJK)
            B_M(IJK,0) = -((-(bma + bme - bmw + &
                              bmn - bms + bmt - bmb ))+ bmr )
            B_MMAX(IJK,0) = max(abs(bma),abs(bme),abs(bmw),&
               abs(bmn),abs(bms),abs(bmt),abs(bmb),abs(bmr) )

            A_M(IJK,east,0) = A_M(IJK,east,0)*D_E(IJK,0)
            A_M(IJK,west,0) = A_M(IJK,west,0)*D_E(IMJK,0)
            A_M(IJK,north,0) = A_M(IJK,north,0)*D_N(IJK,0)
            A_M(IJK,south,0) = A_M(IJK,south,0)*D_N(IJMK,0)
            A_M(IJK,top,0) = A_M(IJK,top,0)*D_T(IJK,0)
            A_M(IJK,bottom,0) = A_M(IJK,bottom,0)*D_T(IJKM,0)

            IF(CARTESIAN_GRID) THEN
               A_M(IJK,east,0) = A_M(IJK,east,0) * A_UPG_E(IJK)
               A_M(IJK,west,0) = A_M(IJK,west,0) * A_UPG_E(IMJK)
               A_M(IJK,north,0) = A_M(IJK,north,0) * A_VPG_N(IJK)
               A_M(IJK,south,0) = A_M(IJK,south,0) * A_VPG_N(IJMK)
               A_M(IJK,top,0) = A_M(IJK,top,0) * A_WPG_T(IJK)
               A_M(IJK,bottom,0) = A_M(IJK,bottom,0) * A_WPG_T(IJKM)
            ENDIF

! solids phase terms are needed for those solids phases that do not
! become close packed. these terms are used to essentially make the
! matrix equation for gas pressure correction a mixture pressure
! correction equation by adding solids phases that do not become
! close packed
            DO M = 1, MMAX
               IF (.NOT.CLOSE_PACKED(M)) THEN
                  B_M(IJK,0) = B_M(IJK,0) - &
                     ((-((ROP_S(IJK,M)-ROP_SO(IJK,M))*VOL(IJK)*ODT+&
                     A_M(IJK,east,M)*U_S(IJK,M)-A_M(IJK,west,M)*U_S(IMJK,M)+&
                     A_M(IJK,north,M)*V_S(IJK,M)-A_M(IJK,south,M)*V_S(IJMK,M)+&
                     A_M(IJK,top,M)*W_S(IJK,M)-A_M(IJK,bottom,M)*W_S(IJKM,M)))+&
                     SUM_R_S(IJK,M)*VOL(IJK))

                  IF(.NOT.CARTESIAN_GRID) THEN
                     A_M(IJK,east,0) = A_M(IJK,east,0) + A_M(IJK,east,M)*D_E(IJK,M)
                     A_M(IJK,west,0) = A_M(IJK,west,0) + A_M(IJK,west,M)*D_E(IMJK,M)
                     A_M(IJK,north,0) = A_M(IJK,north,0) + A_M(IJK,north,M)*D_N(IJK,M)
                     A_M(IJK,south,0) = A_M(IJK,south,0) + A_M(IJK,south,M)*D_N(IJMK,M)
                     A_M(IJK,top,0) = A_M(IJK,top,0) + A_M(IJK,top,M)*D_T(IJK,M)
                     A_M(IJK,bottom,0) = A_M(IJK,bottom,0) + A_M(IJK,bottom,M)*D_T(IJKM,M)
                  ELSE
                     A_M(IJK,east,0) = A_M(IJK,east,0) + A_M(IJK,east,M)*D_E(IJK,M)  * A_UPG_E(IJK)
                     A_M(IJK,west,0) = A_M(IJK,west,0) + A_M(IJK,west,M)*D_E(IMJK,M) * A_UPG_E(IMJK)
                     A_M(IJK,north,0) = A_M(IJK,north,0) + A_M(IJK,north,M)*D_N(IJK,M)  * A_VPG_N(IJK)
                     A_M(IJK,south,0) = A_M(IJK,south,0) + A_M(IJK,south,M)*D_N(IJMK,M) * A_VPG_N(IJMK)
                     A_M(IJK,top,0) = A_M(IJK,top,0) + A_M(IJK,top,M)*D_T(IJK,M)  * A_WPG_T(IJK)
                     A_M(IJK,bottom,0) = A_M(IJK,bottom,0) + A_M(IJK,bottom,M)*D_T(IJKM,M) * A_WPG_T(IJKM)
                  ENDIF
               ENDIF
            ENDDO

            A_M(IJK,0,0) = -(A_M(IJK,east,0)+A_M(IJK,west,0)+&
                             A_M(IJK,north,0)+A_M(IJK,south,0)+&
                             A_M(IJK,top,0)+A_M(IJK,bottom,0))

            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN
                  A_M(IJK,0,0) = -ONE
                  B_M(IJK,0) = ZERO
               ELSEIF (RO_G0 .NE. UNDEFINED) THEN !This is an error only in incompressible flow
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK,0)
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1)
!!$omp             end critical
               ENDIF
            ENDIF


         ELSE   ! if/else branch .not.fluid_at(ijk)
! set the value (correction) in all wall and flow boundary cells to zero
! note the matrix coefficients and source vector should already be zero
! from the initialization of A and B but the following ensures the zero
! value.
            A_M(IJK,east,0) = ZERO
            A_M(IJK,west,0) = ZERO
            A_M(IJK,north,0) = ZERO
            A_M(IJK,south,0) = ZERO
            A_M(IJK,top,0) = ZERO
            A_M(IJK,bottom,0) = ZERO
            A_M(IJK,0,0) = -ONE
            B_M(IJK,0) = ZERO
         ENDIF   ! end if/else branch fluid_at(ijk)
      ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel


! loezos
      IF (SHEAR) THEN
!!$omp parallel do private(IJK)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               V_G(IJK)=V_G(IJK)-VSH(IJK)
            ENDIF
         ENDDO
      ENDIF

! Modification for compressibility 
      IF (RO_G0 == UNDEFINED) CALL COMPRESSIBLE_PP_G(A_M)

! Remove the asymmetry in matrix caused by the pressure outlet or inlet
! boundaries.  Because the P' at such boundaries is zero we may set the
! coefficient in the neighboring fluid cell to zero without affecting
! the linear equation set.
!!$omp    parallel do                                                     &
!!$omp&   private(IJK,IMJK, IPJK, IJMK, IJPK, IJKM, IJKP)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK)
! Cutting the neighbor link between fluid cell and adjacent p_flow_at cell
            if(p_flow_at(imjk)) A_m(IJK, west, 0) = ZERO
            if(p_flow_at(ipjk)) A_m(IJK, east, 0) = ZERO
            if(p_flow_at(ijmk)) A_m(IJK, south, 0) = ZERO
            if(p_flow_at(ijpk)) A_m(IJK, north, 0) = ZERO
            if(p_flow_at(ijkm)) A_m(IJK, bottom, 0) = ZERO
            if(p_flow_at(ijkp)) A_m(IJK, top, 0) = ZERO
         ENDIF
      ENDDO

! Specify P' to zero for incompressible flows. Check set_bc0
! for details on selection of IJK_P_g.
      IF (IJK_P_G /= UNDEFINED_I) THEN
         B_M(IJK_P_G,0) = ZERO
         A_M(IJK_P_G,:,0) = ZERO
         A_M(IJK_P_G,0,0) = -ONE
      ENDIF

      RETURN

CONTAINS

      INCLUDE 'functions.inc'

END SUBROUTINE SOURCE_PP_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: COMPRESSIBLE_PP_G                                       C
!  Purpose: Correction for incompressible flows.                       C
!                                                                      C
!  Notes: A HS_CORRECT option is available as a better approximation   C
!  for high speed flows because it considers density changes in the    C
!  neighboring C cells. However, the code runs faster without it for   C
!  low speed flows. The gas phase mass balance cannot be maintained    C
!  to machine precision with this approximation.                       C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE COMPRESSIBLE_PP_G(A_M)

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use eos, only: droodp_g
      use fldvar, only: ep_g, u_g, v_g, w_g, rop_g, ro_g, p_g
      use functions, only: fluid_at, im_of, jm_of, km_of
      use functions, only: east_of, west_of, north_of, south_of
      use functions, only: top_of, bottom_of
      use geometry, only: vol, ayz, axz, axy, do_k
      use param, only: dimension_3, dimension_m
      use param, only: east, west, south, north, top, bottom
      use param1, only: one
      use run, only: discretize, odt
      use ur_facs, only: ur_fac
      use xsi
      IMPLICIT NONE

! Parameters
!---------------------------------------------------------------------//
      LOGICAL, PARAMETER :: HS_CORRECT = .FALSE.

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! loezos: used for including shearing
      integer :: incr
! temporary use of global arrays:
! xsi_array: convection weighting factors
      DOUBLE PRECISION :: XSI_e(DIMENSION_3)
      DOUBLE PRECISION :: XSI_n(DIMENSION_3)
      DOUBLE PRECISION :: XSI_t(DIMENSION_3)
! under relaxation factor for pressure
      double precision :: fac
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
!---------------------------------------------------------------------//

! since p_g = p_g* + ur_fac * pp_g
      fac = UR_FAC(1)  

     IF (.NOT.HS_CORRECT) THEN
!!$omp    parallel do  &
!!$omp&   private(IJK)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
! account for change in density w.r.t. pressure
            A_M(IJK,0,0) = A_M(IJK,0,0) - ur_fac(1)*&
               DROODP_G(RO_G(IJK),P_G(IJK))*EP_G(IJK)*VOL(IJK)*ODT
         ENDIF
      ENDDO
      ENDIF

      IF (HS_CORRECT) THEN
! since p_g = p_g* + ur_fac * pp_g
! loezos
      incr=0
      CALL CALC_XSI(DISCRETIZE(1),ROP_G,U_G,V_G,W_G,XSI_E,XSI_N,XSI_T,incr)

!!$omp    parallel do  &
!!$omp&   private(IJK,I,J,K, &
!!$omp&            IMJK,IJMK,IJKM,IJKE,IJKW,IJKN,IJKS,IJKT,IJKB)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKW = WEST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKT = TOP_OF(IJK)
            IJKB = BOTTOM_OF(IJK)
            A_M(IJK,0,0) = A_M(IJK,0,0) - fac*DROODP_G(RO_G(IJK),P_G(IJK))*&
               EP_G(IJK)*((ONE - XSI_E(IJK))*U_G(IJK)*AYZ(IJK)-&
                                 XSI_E(IMJK)*U_G(IMJK)*AYZ(IMJK)+&
                             (ONE-XSI_N(IJK))*V_G(IJK)*AXZ(IJK)-&
                                  XSI_N(IJMK)*V_G(IJMK)*AXZ(IJMK))

            A_M(IJK,east,0) = A_M(IJK,east,0) - EP_G(IJKE)*fac*&
               DROODP_G(RO_G(IJKE),P_G(IJKE))*XSI_E(IJK)*U_G(IJK)*AYZ(IJK)
            A_M(IJK,west,0) = A_M(IJK,west,0) + EP_G(IJKW)*fac*&
               DROODP_G(RO_G(IJKW),P_G(IJKW))*(ONE - XSI_E(IMJK))*U_G(IMJK)*AYZ(IMJK)
            A_M(IJK,north,0) = A_M(IJK,north,0) - EP_G(IJKN)*fac*&
               DROODP_G(RO_G(IJKN),P_G(IJKN))*XSI_N(IJK)*V_G(IJK)*AXZ(IJK)
            A_M(IJK,south,0) = A_M(IJK,south,0) + EP_G(IJKS)*fac*&
               DROODP_G(RO_G(IJKS),P_G(IJKS))*(ONE - XSI_N(IJMK))*V_G(IJMK)*AXZ(IJMK)
            IF (DO_K) THEN
               A_M(IJK,0,0) = A_M(IJK,0,0) - fac*DROODP_G(RO_G(IJK),P_G(IJK))*&
                  EP_G(IJK)*((ONE - XSI_T(IJK))*W_G(IJK)*AXY(IJK)-&
                                    XSI_T(IJKM)*W_G(IJKM)*AXY(IJKM))
               A_M(IJK,top,0) = A_M(IJK,top,0) - EP_G(IJKT)*fac*&
                  DROODP_G(RO_G(IJKT),P_G(IJKT))*XSI_T(IJK)*W_G(IJK)*AXY(IJK)
               A_M(IJK,bottom,0) = A_M(IJK,bottom,0) + EP_G(IJKB)*fac*&
                  DROODP_G(RO_G(IJKB),P_G(IJKB))*(ONE - XSI_T(IJKM))*W_G(IJKM)*AXY(IJKM)
            ENDIF

         ENDIF   !end if (fluid_at(ijk))
      ENDDO    ! end do (ijk=ijkstart3,ijkend3)
      ENDIF   ! end if (hs_correct)

      RETURN
      END SUBROUTINE COMPRESSIBLE_PP_G


