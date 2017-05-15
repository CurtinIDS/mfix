!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SET_BC0                                             C
!  Purpose: This module does the initial setting of boundary           C
!           conditions for cut cells only                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, BC_DT_0, TIME, BC_Jet_g0,C
!                        BC_K_b, BC_K_t, BC_J_s, BC_J_n, BC_I_w,       C
!                        BC_I_e, BC_PLANE, BC_EP_g, BC_P_g, BC_T_g,    C
!                        BC_T_s,  BC_U_g, BC_V_g, BC_W_g,              C
!                        MMAX, BC_ROP_s, BC_U_s, BC_V_s, BC_W_s        C
!  Variables modified: BC_TIME, BC_V_g, I, J, K, IJK, EP_g, P_g, T_g,  C
!                      T_s, U_g, V_g, W_g, ROP_s, U_s, V_s, W_s,       C
!                      M                                               C
!                                                                      C
!  Local variables: L, IJK1, IJK2, IJK3                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_SET_BC0
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bc
      USE compar
      USE cutcell
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE quadric
      USE run
      USE scalars
      USE scales
      USE sendrecv
      USE toleranc

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
!                      Local index for boundary condition
      INTEGER          L
!
!                      indices
      INTEGER          IJK, M, NN ,IJKW,IJKS,IJKB
!
!----------------------------------------------

      INTEGER, DIMENSION(8) :: ACCEPTABLE_DEFAULT_WALL=-1
      LOGICAL :: GLOBAL_CORNER

!
!  Define global corners as acceptable default walls
!  These cells should never be used
!

      IF(.NOT.RE_INDEXING.AND.NumPEs==1) THEN

      ACCEPTABLE_DEFAULT_WALL(1) = FUNIJK(IMIN3,JMIN3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(2) = FUNIJK(IMAX3,JMIN3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(3) = FUNIJK(IMIN3,JMAX3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(4) = FUNIJK(IMAX3,JMAX3,KMIN3)
      ACCEPTABLE_DEFAULT_WALL(5) = FUNIJK(IMIN3,JMIN3,KMAX3)
      ACCEPTABLE_DEFAULT_WALL(6) = FUNIJK(IMAX3,JMIN3,KMAX3)
      ACCEPTABLE_DEFAULT_WALL(7) = FUNIJK(IMIN3,JMAX3,KMAX3)
      ACCEPTABLE_DEFAULT_WALL(8) = FUNIJK(IMAX3,JMAX3,KMAX3)

      ENDIF

!      DO N = 1,8
!         print*,'acceptable default=',ACCEPTABLE_DEFAULT_WALL(N)
!      ENDDO


      DO IJK = ijkstart3, ijkend3

         L = BC_ID(IJK)

         IF(L>0) THEN
            IF(BC_TYPE_ENUM(L)==CG_PO) THEN

               P_STAR(IJK) = ZERO
               P_G(IJK) = SCALE_PRESSURE(BC_P_G(L))
   !
               IF (BC_EP_G(L) /= UNDEFINED) EP_G(IJK) = BC_EP_G(L)
               IF (BC_T_G(L) /= UNDEFINED) then
                  T_G(IJK) = BC_T_G(L)
               ELSE
                  T_g(IJK) = TMIN
               ENDIF

               NN = 1
               IF (NMAX(0) > 0) THEN
                  WHERE (BC_X_G(L,:NMAX(0)) /= UNDEFINED) X_G(IJK,:&
                         NMAX(0)) = BC_X_G(L,:NMAX(0))
                  NN = NMAX(0) + 1
               ENDIF

               IF (NScalar > 0) THEN
                  WHERE (BC_Scalar(L,:NScalar) /= UNDEFINED)&
                  Scalar(IJK,:NScalar) = BC_Scalar(L,:NScalar)
               ENDIF

               IF (K_Epsilon) THEN
                  IF (BC_K_Turb_G(L) /= UNDEFINED) K_Turb_G(IJK) = BC_K_Turb_G(L)
                  IF (BC_E_Turb_G(L) /= UNDEFINED) E_Turb_G(IJK) = BC_E_Turb_G(L)
               ENDIF

               DO M = 1, MMAX
                  IF (BC_ROP_S(L,M) /= UNDEFINED) ROP_S(IJK,M) = BC_ROP_S(L,M)
                  IF(BC_T_S(L,M)/=UNDEFINED)T_S(IJK,M)=BC_T_S(L,M)
                  IF (BC_THETA_M(L,M) /= UNDEFINED) THETA_M(IJK,M)= BC_THETA_M(L,M)
                  NN = 1
                  IF (NMAX(M) > 0) THEN
                     WHERE (BC_X_S(L,M,:NMAX(M)) /= UNDEFINED) X_S(&
                            IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M))
                     NN = NMAX(M) + 1
                  ENDIF
               END DO

            ELSEIF(BC_TYPE_ENUM(L)==CG_MI) THEN

               P_STAR(IJK) = ZERO
   !
               EP_G(IJK) = BC_EP_G(L)
               P_G(IJK) = SCALE_PRESSURE(BC_P_G(L))
               T_G(IJK) = BC_T_G(L)

               IF (NMAX(0) > 0) THEN
                  X_G(IJK,:NMAX(0)) = BC_X_G(L,:NMAX(0))
               ENDIF

               IF (NScalar > 0) THEN
                  Scalar(IJK,:NScalar) = BC_Scalar(L,:NScalar)
               ENDIF

               IF (K_Epsilon) THEN
                  K_Turb_G(IJK) = BC_K_Turb_G(L)
                  E_Turb_G(IJK) = BC_E_Turb_G(L)
               ENDIF

               DO M = 1, MMAX
                  ROP_S(IJK,M) = BC_ROP_S(L,M)
                  T_S(IJK,M) = BC_T_S(L,M)
                  THETA_M(IJK,M) = BC_THETA_M(L,M)

                  IF (NMAX(M) > 0) THEN
                     X_S(IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M))
                  ENDIF

               END DO

               IF(BC_U_g(L)/=UNDEFINED) THEN
                  U_G(IJK) =  BC_U_g(L)
               ELSE
                  U_G(IJK) =  BC_VELMAG_g(L)*NORMAL_S(IJK,1)
               ENDIF


               IF(BC_V_g(L)/=UNDEFINED) THEN
                  V_G(IJK) =  BC_V_g(L)
               ELSE
                  V_G(IJK) =  BC_VELMAG_g(L)*NORMAL_S(IJK,2)
               ENDIF

               IF(BC_W_g(L)/=UNDEFINED) THEN
                  W_G(IJK) =  BC_W_g(L)
               ELSE
                  W_G(IJK) =  BC_VELMAG_g(L)*NORMAL_S(IJK,3)
               ENDIF

               IJKW = WEST_OF(IJK)
               IJKS = SOUTH_OF(IJK)
               IJKB = BOTTOM_OF(IJK)

               IF(FLUID_AT(IJKW)) THEN
                  IF(BC_U_g(L)/=UNDEFINED) THEN
                     U_G(IJKW) =  BC_U_g(L)
                  ELSE
                     U_G(IJKW) =  BC_VELMAG_g(L)*NORMAL_S(IJK,1)
                  ENDIF
               ENDIF

               IF(FLUID_AT(IJKS)) THEN
                  IF(BC_V_g(L)/=UNDEFINED) THEN
                     V_G(IJKS) =  BC_V_g(L)
                  ELSE
                     V_G(IJKS) =  BC_VELMAG_g(L)*NORMAL_S(IJK,2)
                  ENDIF
               ENDIF

               IF(FLUID_AT(IJKB)) THEN
                  IF(BC_W_g(L)/=UNDEFINED) THEN
                     W_G(IJKB) =  BC_W_g(L)
                  ELSE
                     W_G(IJKB) =  BC_VELMAG_g(L)*NORMAL_S(IJK,3)
                  ENDIF
               ENDIF

   !
               M = 1

               DO M=1,MMAX

                  IF(BC_U_s(L,M)/=UNDEFINED) THEN
                     U_S(IJK,M) =  BC_U_S(L,M)
                  ELSE
                     U_S(IJK,M) =  BC_VELMAG_S(L,M)*NORMAL_S(IJK,1)
                  ENDIF

                  IF(BC_V_S(L,M)/=UNDEFINED) THEN
                     V_S(IJK,M) =  BC_V_S(L,M)
                  ELSE
                     V_S(IJK,M) =  BC_VELMAG_S(L,M)*NORMAL_S(IJK,2)
                  ENDIF

                  IF(BC_W_S(L,M)/=UNDEFINED) THEN
                     W_S(IJK,M) =  BC_W_S(L,M)
                  ELSE
                     W_S(IJK,M) =  BC_VELMAG_S(L,M)*NORMAL_S(IJK,3)
                  ENDIF

                  IJKW = WEST_OF(IJK)
                  IJKS = SOUTH_OF(IJK)
                  IJKB = BOTTOM_OF(IJK)

                  IF(FLUID_AT(IJKW)) THEN
                     IF(BC_U_S(L,M)/=UNDEFINED) THEN
                        U_S(IJKW,M) =  BC_U_S(L,M)
                     ELSE
                        U_S(IJKW,M) =  BC_VELMAG_S(L,M)*NORMAL_S(IJK,1)
                     ENDIF
                  ENDIF

                  IF(FLUID_AT(IJKS)) THEN
                     IF(BC_V_S(L,M)/=UNDEFINED) THEN
                        V_S(IJKS,M) =  BC_V_S(L,M)
                     ELSE
                        V_S(IJKS,M) =  BC_VELMAG_S(L,M)*NORMAL_S(IJK,2)
                     ENDIF
                  ENDIF

                  IF(FLUID_AT(IJKB)) THEN
                     IF(BC_W_S(L,M)/=UNDEFINED) THEN
                        W_S(IJKB,M) =  BC_W_S(L,M)
                     ELSE
                        W_S(IJKB,M) =  BC_VELMAG_S(L,M)*NORMAL_S(IJK,3)
                     ENDIF
                  ENDIF

               ENDDO

   !            IF (MMAX > 0) THEN
   !               U_S(IJK,:MMAX) = BC_U_S(L,:MMAX)
   !               V_S(IJK,:MMAX) = BC_V_S(L,:MMAX)
   !               W_S(IJK,:MMAX) = BC_W_S(L,:MMAX)
   !
   !               IF(FLUID_AT(IJKW)) THEN
   !                  U_S(IJKW,:MMAX) = BC_U_S(L,:MMAX)
   !               ENDIF
   !
   !               IF(FLUID_AT(IJKS)) THEN
   !                  V_S(IJKS,:MMAX) = BC_V_S(L,:MMAX)
   !               ENDIF
   !
   !               IF(FLUID_AT(IJKB)) THEN
   !                  W_S(IJKB,:MMAX) = BC_W_S(L,:MMAX)
   !               ENDIF
   !               M = MMAX + 1
   !            ENDIF


               IF (MW_AVG == UNDEFINED) THEN
                  MW_MIX_G(IJK) = CALC_MW(X_G,DIMENSION_3,IJK,NMAX(0),MW_G)
               ELSE
                  MW_MIX_G(IJK) = MW_AVG
               ENDIF
               IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G&
                  (IJK),P_G(IJK),T_G(IJK))
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)

            ENDIF
         ENDIF

         IF(DEFAULT_WALL_AT(IJK)) THEN

!            print*,'Default_wall_at IJK=',IJK,I_OF(IJK),J_OF(IJK),K_OF(IJK)

            GLOBAL_CORNER = .FALSE.
            DO NN = 1,8
               IF(IJK==ACCEPTABLE_DEFAULT_WALL(NN)) GLOBAL_CORNER = .TRUE.
            ENDDO

            IF(.NOT.GLOBAL_CORNER.AND..NOT.BLOCKED_CELL_AT(IJK)) THEN

               ICBC_FLAG(IJK)(2:3) = 'CG'

               IF((MyPE == PE_IO).AND.PRINT_WARNINGS) THEN
                  WRITE(*,*) 'WARNING: DEFAULT WALL DETECTED AT I,J,K = ',I_OF(IJK),J_OF(IJK),K_OF(IJK) ,BLOCKED_CELL_AT(IJK)
                  WRITE(*,*) '         WHEN USING CARTESIAN GRID CUT-CELL FEATURE.'
                  WRITE(*,*) '         DEFAULT WALLS ARE NOT ALLOWED WITH CUT-CELLS.'
                  WRITE(*,*) '         THE DEFAULT WALL WAS REMOVED ALONG THIS CELL.'
                  WRITE(*,*) ''
               ENDIF
!               CALL MFIX_EXIT(MYPE)

            ENDIF

         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE CG_SET_BC0

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization: FLAG_G , FLUID_AT_G
!// 360 Check if i,j,k resides on current processor
