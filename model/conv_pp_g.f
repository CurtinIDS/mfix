!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_Pp_g                                               C
!  Purpose: Determine convection terms for Pressure correction         C
!           equation.                                                  C
!                                                                      C
!  Notes: The off-diagonal coefficients calculated here must be        C
!         positive. The center coefficient and the source vector are   C
!         negative.                                                    C
!         Multiplication with factors d_e, d_n, and d_t are carried    C
!         out in source_pp_g.  Constant pressure boundaries are        C
!         handled by holding the Pp_g at the boundaries zero.  For     C
!         specified mass flow boundaries (part of) a's are calculated  C
!         here since b is calculated from a's in source_pp_g.  After   C
!         calculating b, a's are multiplied by d and at the flow       C
!         boundaries are set to zero.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: to use face densities calculated in CONV_ROP               C
!  Author: M. Syamlal                                 Date: 1-JUN-05   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_PP_G(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE run
      USE parallel
      USE physprop
      USE geometry
      USE indices
      USE pgcor
      USE compar
      USE mflux
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
      INTEGER :: M
! local value of A_m
      DOUBLE PRECISION :: am
!-----------------------------------------------


! Calculate convection fluxes through each of the faces

!$omp  parallel do default(none) private( IJK, IPJK, IJPK, IJKP, M, AM, IMJK, IJMK, IJKM) &
!$omp              shared(ijkstart3, ijkend3, rop_ge, rop_gn, rop_gt, axz, axy, ayz, do_k, a_m)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)

! East face (i+1/2, j, k)
            AM = ROP_GE(IJK)*AYZ(IJK)
            A_M(IJK,east,0) = AM
            A_M(IPJK,west,0) = AM

! North face (i, j+1/2, k)
            AM = ROP_GN(IJK)*AXZ(IJK)
            A_M(IJK,north,0) = AM
            A_M(IJPK,south,0) = AM

! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               AM = ROP_GT(IJK)*AXY(IJK)
               A_M(IJK,top,0) = AM
               A_M(IJKP,bottom,0) = AM
            ENDIF

! West face (i-1/2, j, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               AM = ROP_GE(IMJK)*AYZ(IMJK)
               A_M(IJK,west,0) = AM
            ENDIF

! South face (i, j-1/2, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLUID_AT(IJMK)) THEN
               AM = ROP_GN(IJMK)*AXZ(IJMK)
               A_M(IJK,south,0) = AM
            ENDIF

! Bottom face (i, j, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                  AM = ROP_GT(IJKM)*AXY(IJKM)
                  A_M(IJK,bottom,0) = AM
               ENDIF
            ENDIF
         ENDIF   ! end if (fluid_at(ijk))
      ENDDO   ! end do ijk=ijkstart3,ijkend3


! solids phase terms are needed for those solids phases that do not
! become close packed. these terms are used to essentially make the
! matrix equation for gas pressure correction a mixture pressure
! correction equation by adding solids phases that do not become
! close packed
      DO M = 1, MMAX
         IF (.NOT.CLOSE_PACKED(M)) THEN
            DO IJK = ijkstart3, ijkend3
               IF (FLUID_AT(IJK)) THEN
                  IPJK = IP_OF(IJK)
                  IJPK = JP_OF(IJK)
                  IJKP = KP_OF(IJK)

! East face (i+1/2, j, k)
                  AM = ROP_SE(IJK,M)*AYZ(IJK)
                  A_M(IJK,east,M) = AM
                  A_M(IPJK,west,M) = AM

! North face (i, j+1/2, k)
                  AM = ROP_SN(IJK,M)*AXZ(IJK)
                  A_M(IJK,north,M) = AM
                  A_M(IJPK,south,M) = AM

! Top face (i, j, k+1/2)
                  IF (DO_K) THEN
                     AM = ROP_ST(IJK,M)*AXY(IJK)
                     A_M(IJK,top,M) = AM
                     A_M(IJKP,bottom,M) = AM
                  ENDIF

! West face (i-1/2, j, k)
                  IMJK = IM_OF(IJK)
                  IF (.NOT.FLUID_AT(IMJK)) THEN
                     AM = ROP_SE(IMJK,M)*AYZ(IMJK)
                     A_M(IJK,west,M) = AM
                  ENDIF

! South face (i, j-1/2, k)
                  IJMK = JM_OF(IJK)
                  IF (.NOT.FLUID_AT(IJMK)) THEN
                     AM = ROP_SN(IJMK,M)*AXZ(IJMK)
                     A_M(IJK,south,M) = AM
                  ENDIF

! Bottom face (i, j, k-1/2)
                  IF (DO_K) THEN
                     IJKM = KM_OF(IJK)
                     IF (.NOT.FLUID_AT(IJKM)) THEN
                        AM = ROP_ST(IJKM,M)*AXY(IJKM)
                        A_M(IJK,bottom,M) = AM
                     ENDIF
                  ENDIF

               ENDIF   ! end if(fluid_at(ijk))
            ENDDO   ! end do ijk=ijkstart3,ijkend3
         ENDIF   ! end if(.not.close_packed)
      ENDDO   ! end do (m=1,mmax)

      RETURN
      END SUBROUTINE CONV_PP_G
