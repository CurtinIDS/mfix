!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_phi                                              C
!  Purpose: Determine source terms for phi eq. The terms               C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative. The off-diagonal     C
!     coefficients are positive. S_p must be positive.                 C
!                                                                      C
!  Note: This routine is now restricted to Non-Negative scalers when   C
!     using deferred correction.                                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_PHI(S_P, S_C, EP, PHI, M, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE scales
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Source term on LHS.  Must be positive.
      DOUBLE PRECISION, INTENT(IN) :: S_p(DIMENSION_3)
! Source term on RHS
      DOUBLE PRECISION, INTENT(IN) :: S_C(DIMENSION_3)
! Phase volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EP(DIMENSION_3)
! Phi
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3)
! phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM
!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN

! dilute flow
            IF (EP(IJK) <= DIL_EP_S) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

! use the average boundary cell values to compute phi (sof, Aug 23 2005)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IF (EP(WEST_OF(IJK)) > DIL_EP_S .AND. &
                   .NOT.IS_AT_E(IMJK)) A_M(IJK,west,M) = ONE
               IF (EP(EAST_OF(IJK)) > DIL_EP_S .AND. &
                   .NOT.IS_AT_E(IJK)) A_M(IJK,east,M) = ONE
               IF (EP(SOUTH_OF(IJK)) > DIL_EP_S .AND. &
                   .NOT.IS_AT_N(IJMK)) A_M(IJK,south,M) = ONE
               IF (EP(NORTH_OF(IJK)) > DIL_EP_S .AND. &
                   .NOT.IS_AT_N(IJK)) A_M(IJK,north,M) = ONE
               IF(.NOT. NO_K) THEN
                  IF (EP(BOTTOM_OF(IJK)) > DIL_EP_S .AND. &
                     .NOT.IS_AT_T(IJKM)) A_M(IJK,bottom,M) = ONE
                  IF (EP(TOP_OF(IJK)) > DIL_EP_S .AND. &
                     .NOT.IS_AT_T(IJK)) A_M(IJK,top,M) = ONE
               ENDIF

               IF((A_M(IJK,west,M)+A_M(IJK,east,M)+&
                  A_M(IJK,south,M)+A_M(IJK,north,M)+ &
                  A_M(IJK,bottom,M)+A_M(IJK,top,M)) == ZERO) THEN
                  B_M(IJK,M) = -PHI(IJK)
               ELSE
                  A_M(IJK,0,M) = -(A_M(IJK,east,M) + A_M(IJK,west,M) +&
                     A_M(IJK,north,M) + A_M(IJK,south,M) + &
                     A_M(IJK,top,M)+A_M(IJK,bottom,M))
               ENDIF

! Normal case
            ELSE

! Collect the terms
               A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+&
                                A_M(IJK,north,M)+A_M(IJK,south,M)+&
                                A_M(IJK,top,M)+A_M(IJK,bottom,M)+S_P(IJK))

! B_m and A_m are corrected in case deferred corrections computes B_m > S_c
! see CONV_DIF_PHI_DC.
               IF(B_M(IJK,M) < S_C(IJK) .OR. PHI(IJK) == ZERO) THEN
                  B_M(IJK,M) = -S_C(IJK)+B_M(IJK,M)
               ELSE ! disable ELSE statememt if PHI can be negative
                  A_M(IJK,0,M) = A_M(IJK,0,M) - B_M(IJK,M)/PHI(IJK)
                  B_M(IJK,M) = -S_C(IJK)
               ENDIF

            ENDIF   ! end if/else (ep_g(ijk)<=dil_ep_s)
         ENDIF    ! end if (fluid_at(ijk))
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)


      RETURN
      END SUBROUTINE SOURCE_PHI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_PHI                                        C
!                                                                      C
!  Purpose: Adds point sources to the scalar equations.                C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_PHI(PHI, PS_PHI, PS_FLOW,  &
         M, A_M, B_M)

      use compar
      use geometry
      use indices
      use physprop
      use ps
      use run

! To be removed upon complete integration of point source routines.
      use bc
      use usr
      use functions
      use param, only: dimension_3, dimension_m
      use param1, only: one, small_number

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------

! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)


! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(IN) :: PS_PHI(DIMENSION_BC)

! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(IN) :: PS_FLOW(DIMENSION_BC)

      INTEGER, intent(in) :: M

! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Vector b_m
      DOUBLE PRECISION, intent(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------

! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV

! terms of bm expression
      DOUBLE PRECISION pSource

!-----------------------------------------------

      PSV_LP: do PSV = 1, DIMENSION_PS

         if(.NOT.PS_DEFINED(PSV)) cycle PSV_LP
         if(abs(PS_FLOW(PSV)) < small_number) cycle PSV_LP

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(ijk)) cycle

            if(A_M(IJK,0,M) == -ONE .AND. B_M(IJK,M) == -PHI(IJK)) then
               B_M(IJK,M) = -PS_PHI(PSV)
            else

! Calculate the mass flow rate for this cell. (mass/time)
               pSource = PS_FLOW(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

               A_M(IJK,0,M) = A_M(IJK,0,M) - pSource
               B_M(IJK,M) = B_M(IJK,M) - PS_PHI(PSV) * pSource

            endif

         enddo
         enddo
         enddo

      enddo PSV_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_PHI
