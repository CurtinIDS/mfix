!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW                                             C
!  Purpose: Set specified outflow bc for pressure outflow,             C
!  mass outflow, outflow and now also pressure inflow bc               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!                                                                      C
!  Comments:                                                           C
!  If the outflow boundary is on the W, S or B side of the domain and  C
!  the component of velocity through the plane is defined then this    C
!  routine will NOT modify it (i.e., its value from the momentum       C
!  solver is maintained).                                              C
!                                                                      C
!  In general it would seem this routine does little more than what    C
!  is already done within the momentum routines in terms of velocity.  C
!  The normal component is either 1) untouched if the outflow is on    C
!  the W, S, B side of the domain or 2) is set to value of the         C
!  adjacent fluid cell if it is on the E, N, T side (similar to the    C
!  respective momentum bc routines). The primary addition here is      C
!  that the tangential components of a bc cell are set to that of      C
!  the adjacent fluid cell. Note the tangential components are not     C
!  explicitly handled in the momentum _BC_ routines; instead their     C
!  values are based on solution of the momentum equation which is      C
!  replaced here                                                       C
!                                                                      C
!  Several routines are called which perform the following tasks:      C
!  set_outflow_misc - several derived quantities are set in the        C
!      boundary                                                        C
!  set_outflow_ep - the void/volume fraction and bulk densities are    C
!      set in the boundary                                             C
!  set_outflow_fluxes - convective fluxes are set in the boundary      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_OUTFLOW(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e

      use fldvar, only: rop_g, rop_s
      use fldvar, only: u_g, v_g, w_g
      use fldvar, only: u_s, v_s, w_s
      use physprop, only: mmax

      use functions, only: im_of, ip_of, jm_of, jp_of, km_of, kp_of
      use functions, only: fluid_at
      use functions, only: is_on_mype_plus2layers
      use functions, only: funijk
      use compar, only: dead_cell_at

      use param, only: dimension_m
      use param1, only: undefined, zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, M
! index for boundary cell
      INTEGER :: IJK
! index for a fluid cell adjacent to the boundary cell
      INTEGER :: FIJK
! local value for normal component of gas and solids velocity defined
! such that
      DOUBLE PRECISION :: RVEL_G, RVEL_S(DIMENSION_M)
!---------------------------------------------------------------------//

! Loop over the range of boundary cells
      DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)
! Check if current i,j,k resides on this PE
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF(DEAD_CELL_AT(I,J,K)) CYCLE
               IJK = FUNIJK(I,J,K)

! Fluid cell at West
! --------------------------------------------------------------------//
               IF (FLUID_AT(IM_OF(IJK))) THEN
                  FIJK = IM_OF(IJK)
                  RVEL_G = U_G(FIJK)
                  IF (MMAX>0) RVEL_S(:MMAX) = U_S(FIJK,:MMAX)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)
                  IF (BC_TYPE_ENUM(BCV)==P_INFLOW) &
                     CALL SET_PINOUTFLOW(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! Set the boundary cell value of the normal component of velocity
! according to the value in the adjacent fluid cell. Note the value
! of the boundary velocity is a scaled version of the value of the
! adjacent fluid velocity based on the concentration ratio of the fluid
! cell to the boundary cell.
! - For the gas phase, this ratio is most likely 1 except for
!   compressible cases with a PO/PI boundary where P_g of the boundary
!   is set and may differ from the value of the adjacent fluid cell.
! - For the solids phase this seems unnecessary..? Differences may arise
!   if bc_rop_s is set...
                  IF (ROP_G(IJK) > ZERO) THEN
                    U_G(IJK) = ROP_G(FIJK)*U_G(FIJK)/ROP_G(IJK)
                  ELSE
                     U_G(IJK) = ZERO
                  ENDIF

! the tangential components are not explicitly handled in the boundary
! condition routines of the corresponding momentum equation
                  V_G(IJK) = V_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  IF (MMAX > 0) THEN
                     WHERE (ROP_S(IJK,:MMAX) > ZERO)
                       U_S(IJK,:MMAX) = ROP_S(FIJK,:MMAX)*&
                           U_S(FIJK,:MMAX)/ROP_S(IJK,:MMAX)
                     ELSEWHERE
                        U_S(IJK,:MMAX) = ZERO
                     END WHERE
                     V_S(IJK,:MMAX) = V_S(FIJK,:MMAX)
                     W_S(IJK,:MMAX) = W_S(FIJK,:MMAX)
                  ENDIF

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! end if (fluid_at(im_of(ijk)))


! Fluid cell at East
! --------------------------------------------------------------------//
               IF (FLUID_AT(IP_OF(IJK))) THEN
                  FIJK = IP_OF(IJK)
! define normal component such that it is positive when exiting the
! domain
                  RVEL_G = -U_G(IJK)
                  IF (MMAX >0) RVEL_S(:MMAX) = -U_S(IJK,:MMAX)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)
                  IF (BC_TYPE_ENUM(BCV)==P_INFLOW) &
                     CALL SET_PINOUTFLOW(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! provide an initial value for the velocity component through the domain
! otherwise its present value (from solution of the corresponding
! momentum eqn) is kept. values for the velocity components in the off
! directions are modified (needed for PO or O boundaries but not MO or
! PI as velocities should be fully specified by this point)
                  IF (U_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        U_G(IJK) = ROP_G(FIJK)*U_G(FIJK)/ROP_G(IJK)
                     ELSE
                        U_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  V_G(IJK) = V_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  DO M = 1, MMAX
                     IF (U_S(IJK,M) == UNDEFINED) THEN
                        IF (ROP_S(IJK,M) > ZERO) THEN
                           U_S(IJK,M) = ROP_S(FIJK,M)*&
                              U_S(FIJK,M)/ROP_S(IJK,M)
                        ELSE
                           U_S(IJK,M) = ZERO
                        ENDIF
                     ENDIF
                     V_S(IJK,M) = V_S(FIJK,M)
                     W_S(IJK,M) = W_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! end if (fluid_at(ip_of(ijk)))


! Fluid cell at South
! --------------------------------------------------------------------//
               IF (FLUID_AT(JM_OF(IJK))) THEN
                  FIJK = JM_OF(IJK)
                  RVEL_G = V_G(FIJK)
                  IF(MMAX>0) RVEL_S(:MMAX) = V_S(FIJK,:MMAX)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)
                  IF (BC_TYPE_ENUM(BCV)==P_INFLOW) &
                     CALL SET_PINOUTFLOW(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (ROP_G(IJK) > ZERO) THEN
                     V_G(IJK) = ROP_G(FIJK)*V_G(FIJK)/ROP_G(IJK)
                  ELSE
                     V_G(IJK) = ZERO
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  IF (MMAX > 0) THEN
                      WHERE (ROP_S(IJK,:MMAX) > ZERO)
                         V_S(IJK,:MMAX) = ROP_S(FIJK,:MMAX)*&
                            V_S(FIJK,:MMAX)/ROP_S(IJK,:MMAX)
                      ELSEWHERE
                         V_S(IJK,:MMAX) = ZERO
                      END WHERE
                      U_S(IJK,:MMAX) = U_S(FIJK,:MMAX)
                      W_S(IJK,:MMAX) = W_S(FIJK,:MMAX)
                  ENDIF

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! end if (fluid_at(jm_of(ijk)))


! Fluid cell at North
! --------------------------------------------------------------------//
               IF (FLUID_AT(JP_OF(IJK))) THEN
                  FIJK = JP_OF(IJK)
                  RVEL_G = -V_G(IJK)
                  IF (MMAX>0) RVEL_S(:MMAX) = -V_S(IJK,:MMAX)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)
                  IF (BC_TYPE_ENUM(BCV)==P_INFLOW) &
                     CALL SET_PINOUTFLOW(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (V_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        V_G(IJK) = ROP_G(FIJK)*V_G(FIJK)/ROP_G(IJK)
                     ELSE
                        V_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  DO M = 1, MMAX
                     IF (V_S(IJK,M) == UNDEFINED) THEN
                        IF (ROP_S(IJK,M) > ZERO) THEN
                           V_S(IJK,M) = ROP_S(FIJK,M)*&
                              V_S(FIJK,M)/ROP_S(IJK,M)
                        ELSE
                           V_S(IJK,M) = ZERO
                        ENDIF
                     ENDIF
                     U_S(IJK,M) = U_S(FIJK,M)
                     W_S(IJK,M) = W_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! if (fluid_at(jp_of(ijk)))


! Fluid cell at Bottom
! --------------------------------------------------------------------//
               IF (FLUID_AT(KM_OF(IJK))) THEN
                  FIJK = KM_OF(IJK)
                  RVEL_G = W_G(FIJK)
                  IF (MMAX>0) RVEL_S(:MMAX) = W_S(FIJK,:MMAX)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)
                  IF (BC_TYPE_ENUM(BCV)==P_INFLOW) &
                     CALL SET_PINOUTFLOW(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (ROP_G(IJK) > ZERO) THEN
                     W_G(IJK) = ROP_G(FIJK)*W_G(FIJK)/ROP_G(IJK)
                  ELSE
                     W_G(IJK) = ZERO
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  V_G(IJK) = V_G(FIJK)

                  IF (MMAX > 0) THEN
                    WHERE (ROP_S(IJK,:MMAX) > ZERO)
                        W_S(IJK,:MMAX) = ROP_S(FIJK,:MMAX)*&
                           W_S(FIJK,:MMAX)/ROP_S(IJK,:MMAX)
                     ELSEWHERE
                        W_S(IJK,:MMAX) = ZERO
                     END WHERE
                     U_S(IJK,:MMAX) = U_S(FIJK,:MMAX)
                     V_S(IJK,:MMAX) = V_S(FIJK,:MMAX)
                  ENDIF

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! if (fluid_at(km_of(ijk)))


! Fluid cell at Top
! --------------------------------------------------------------------//
               IF (FLUID_AT(KP_OF(IJK))) THEN
                  FIJK = KP_OF(IJK)
                  RVEL_G = -W_G(IJK)
                  IF (MMAX>0) RVEL_S(:MMAX) = -W_S(IJK,:MMAX)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)
                  IF (BC_TYPE_ENUM(BCV)==P_INFLOW) &
                     CALL SET_PINOUTFLOW(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (W_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        W_G(IJK) = ROP_G(FIJK)*W_G(FIJK)/ROP_G(IJK)
                     ELSE
                        W_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  V_G(IJK) = V_G(FIJK)

                  DO M = 1, MMAX
                     IF (W_S(IJK,M) == UNDEFINED) THEN
                        IF (ROP_S(IJK,M) > ZERO) THEN
                           W_S(IJK,M) = ROP_S(FIJK,M)*&
                              W_S(FIJK,M)/ROP_S(IJK,M)
                        ELSE
                           W_S(IJK,M) = ZERO
                        ENDIF
                     ENDIF
                     U_S(IJK,M) = U_S(FIJK,M)
                     V_S(IJK,M) = V_S(FIJK,M)
                  ENDDO

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! if (fluid_at(kp_of(ijk)))

            ENDDO   ! end do (i=i1,i2)
         ENDDO   ! end do (j=j1,j2)
      ENDDO   ! end do (k=k1,k2)

      RETURN
      END SUBROUTINE SET_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW_MISC                                        C
!  Purpose: Set the value of certain variables in the specified        C
!  outflow boundary cell that would not otherwise be set according     C
!  to their value in the adjacent fluid cell.                          C
!                                                                      C
!  READ BEFORE MODIFYING!!                                             C
!  For many of the scalar field variables (e.g., T_g, T_s, X_g, X_s,   C
!  k_turb_g, e_turb_g, theta_m and scalar) their corresponding         C
!  governing equation solver routine sets its own value in the         C
!  boundary cell (see bc_phi). So DO NOT set their value here unless   C
!  a clear explanation to the need is also provided; more likely it    C
!  should simply be dealt with in the governing equation routine.      C
!                                                                      C
!  This routine sets the value of quantites that are derived or with   C
!  unique solver routines (e.g., P_star, P_s, MW_MIX_g, P_g) since no  C
!  other routine will.                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_OUTFLOW_MISC(BCV, IJK, FIJK)

! Global variables
!---------------------------------------------------------------------//
      use bc
      use run, only: kt_type_enum, ghd_2007
      use fldvar, only: p_g, ro_g, T_g
      use fldvar, only: p_s, p_star
      use physprop, only: smax, mmax
      use physprop, only: ro_g0, mw_mix_g
      use eos, only: EOSG

! Global parameters
!---------------------------------------------------------------------//
      use param1, only: undefined
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
!---------------------------------------------------------------------//

      IF (BC_TYPE_ENUM(BCV) /= P_OUTFLOW .AND. &
          BC_TYPE_ENUM(BCV) /= P_INFLOW) P_G(IJK) = P_G(FIJK)

      MW_MIX_G(IJK) = MW_MIX_G(FIJK)

! T_g (bc_phi), P_g (above depending on bc), and MW_MIX_G (above) have
! now been defined at IJK
      IF (RO_G0 == UNDEFINED) RO_G(IJK) = &
         EOSG(MW_MIX_G(IJK),P_G(IJK),T_G(IJK))

      IF (SMAX >0) THEN
         P_STAR(IJK) = P_STAR(FIJK)
         P_S(IJK,:SMAX) = P_S(FIJK,:SMAX)
      ENDIF

      IF (KT_TYPE_ENUM == GHD_2007 .AND. MMAX>0) &
         P_S(IJK,MMAX) = P_S(FIJK,MMAX)

      RETURN
      END SUBROUTINE SET_OUTFLOW_MISC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW_EP                                          C
!  Purpose: Set the volume fraction/bulk density (i.e., ROP_s, EP_g    C
!  ROP_S) in the specified outflow boundary cell that would not        C
!  otherwise be set according to their value in the adjacent fluid     C
!  cell.                                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! Global variables
!---------------------------------------------------------------------//
      use bc
      use run, only: kt_type_enum, ghd_2007
      use physprop, only: smax, mmax
      use fldvar, only: rop_g, ro_g, ep_g
      use fldvar, only: rop_s, ep_s
      use discretelement, only: discrete_element, des_mmax
! Global parameters
!---------------------------------------------------------------------//
      use param, only: dimension_m
      use param1, only: undefined, zero, one
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
! the gas or solids velocity in the fluid cell adjacent to the boundary
! cell dot with the outward normal of that bc plane; defines the gas or
! solids velocity component normal to the bc plane as positive when it
! is flowing into the bc cell from the fluid cell. so, for example, for
! an outflow on the eastern boundary this is the u component of velocity
! while for an outflow on the western boundary this is the -u component,
! etc.
      DOUBLE PRECISION, INTENT(IN) :: RVEL_G
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMENSION_M) :: RVEL_S

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: M
! sum of solids phases volume fractions
      DOUBLE PRECISION :: SUM_EPs
! sum of solids phases bulk densities
      DOUBLE PRECISION :: SUM_ROPS
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
!---------------------------------------------------------------------//

! initializing summation quantities
      SUM_ROPS = ZERO
      SUM_EPS = ZERO

      DO M = 1, SMAX

         IF(BC_TYPE_ENUM(BCV) == P_INFLOW) THEN
            ROP_S(IJK,M) = ROP_S(FIJK,M)
         ELSE
! the outflow type bc do not permit 're-entering' solids, in which
! case the solids are removed
            IF (RVEL_S(M) >= ZERO) THEN
! solids are leaving the domain
               ROP_S(IJK,M) = ROP_S(FIJK,M)
            ELSE
! solids cannot enter the domain at an outflow cell
               ROP_S(IJK,M) = ZERO
            ENDIF
         ENDIF

! if bc_rop_s is defined, set value of rop_s in the ijk boundary cell
! according to user definition
         IF(BC_ROP_S(BCV,M)/=UNDEFINED) ROP_S(IJK,M)=BC_ROP_S(BCV,M)

! add to total solids phase bulk density and solids volume fraction
         SUM_ROPS = SUM_ROPS + ROP_S(IJK,M)
         SUM_EPS = SUM_EPS + EP_S(IJK,M)
      ENDDO   ! end do (m=1,smax)

      IF (KT_TYPE_ENUM == GHD_2007) ROP_S(IJK,MMAX) = SUM_ROPS

! this section must be skipped until after the initial setup of the
! discrete element portion of the simulation (set_bc1 is called once
! before the initial dem setup).
      IF (DISCRETE_ELEMENT .AND. .NOT.FIRST_PASS) THEN
         FIRST_PASS = .FALSE.
         DO M = MMAX+1, DES_MMAX+MMAX
! unlike in the two fluid model, in the discrete element model it is
! possible to actually calculate the bulk density in a flow boundary
! cell. Currently, however, such calculations are not strictly enforced.
! therefore use the bulk density of the adjacent fluid cell
            ROP_S(IJK,M) = ROP_S(FIJK,M)
            SUM_ROPS = SUM_ROPS + ROP_S(IJK,M)
            SUM_EPS = SUM_EPS + EP_S(IJK,M)
         ENDDO
      ENDIF

! if bc_ep_g undefined, set ep_g accordingly (based on flow condition
! or based on bc_rop_s). if bc_ep_g is defined its set value will be
! maintained (from set_bc0).
      IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE - SUM_EPS

! now that ep_g in the boundary cell is known, define the bulk density
! of the gas phase in the boundary cell
      ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)

      RETURN
      END SUBROUTINE SET_OUTFLOW_EP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Update convective fluxes....                               C
!  Set the value of the convective fluxes in the specified boundary    C
!  cell according to their value in the adjacent fluid cell.           C
!                                                                      C
!  Comment/concern:                                                    C
!  Should these be assigned in the same method as the velocity? Note   C
!  if bc_plane is W, S, B then the normal component of velocity may be C
!  assigned a zero value as opposed to value of its neighboring fluid  C
!  cell. This routine would seem to introduce some inconsistency       C
!  between velocity and flux at boundary.                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_OUTFLOW_FLUXES(IJK,FIJK)

! Global variables
!---------------------------------------------------------------------//
      use physprop, only: mmax
      use run, only: kt_type_enum, ghd_2007
      use run, only: added_mass
      use mflux, only: flux_ge, flux_gn, flux_gt
      use mflux, only: flux_se, flux_sn, flux_st
      use mflux, only: flux_ne, flux_nn, flux_nt
      use mflux, only: flux_gse, flux_gsn, flux_gst
      use mflux, only: flux_sse, flux_ssn, flux_sst

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
!---------------------------------------------------------------------//

      Flux_gE(IJK) = Flux_gE(FIJK)
      Flux_gN(IJK) = Flux_gN(FIJK)
      Flux_gT(IJK) = Flux_gT(FIJK)
      IF (MMAX >0) THEN
         Flux_sE(IJK,:MMAX) = Flux_sE(FIJK,:MMAX)
         Flux_sN(IJK,:MMAX) = Flux_sN(FIJK,:MMAX)
         Flux_sT(IJK,:MMAX) = Flux_sT(FIJK,:MMAX)
      ENDIF

      IF(ADDED_MASS) THEN
        Flux_gSE(IJK) = Flux_gSE(FIJK)
        Flux_gSN(IJK) = Flux_gSN(FIJK)
        Flux_gST(IJK) = Flux_gST(FIJK)
        IF (MMAX >0) THEN
           Flux_sSE(IJK) = Flux_sSE(FIJK)
           Flux_sSN(IJK) = Flux_sSN(FIJK)
           Flux_sST(IJK) = Flux_sST(FIJK)
         ENDIF
      ENDIF

      IF (KT_TYPE_ENUM == GHD_2007) THEN
         Flux_nE(IJK) = Flux_nE(FIJK)
         Flux_nN(IJK) = Flux_nN(FIJK)
         Flux_nT(IJK) = Flux_nT(FIJK)
      ENDIF

      RETURN
      END SUBROUTINE SET_OUTFLOW_FLUXES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Purpose: Specify the field variable in the bc cell depending on the  C
! flow direction. If flow is into domain apply user specified BC. If   C
! the flow is out of the domain follow outflow type approach wherein   C
! the value of the adjacent fluid cell is applied to the bc cell.      C
!                                                                      C
! WARNING: this routine only serves to specify the field variables     C
! whose governing solver routine invokes bc_phi! do not insert any     C
! other settings here (such as derived quantities) as it is likely     C
! not appropriate!!                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_PINOUTFLOW(BCV,IJK,FIJK,RVEL_G,RVEL_S)

! Global variables
!---------------------------------------------------------------------//
      use bc, only: bc_t_g, bc_x_g
      use bc, only: bc_scalar, bc_k_turb_g, bc_e_turb_g
      use bc, only: bc_t_s, bc_x_s, bc_theta_m
      use fldvar, only: t_g, t_s, theta_m
      use fldvar, only: x_g, x_s, scalar
      use fldvar, only: k_turb_g, e_turb_g
! needed because of ghd theory
      use fldvar, only: rop_s, ro_s, d_p
      use run, only: kt_type_enum, ghd_2007
      use run, only: k_epsilon
      use physprop, only: nmax, smax, mmax
      use scalars, only: nscalar, phase4scalar

! Global parameters
!---------------------------------------------------------------------//
      use param, only: dimension_m
      use param1, only: zero
      use constant, only: pi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! index for a fluid cell adjacent to the boundary cell
      INTEGER, INTENT(IN) :: FIJK
! gas and solids velocity defined to be positive when flow is leaving
! the domain (outflow!)
      DOUBLE PRECISION, INTENT(IN) :: RVEL_G, RVEL_S(DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! number density
      DOUBLE PRECISION :: Nm, NTot
! indices
      INTEGER :: M, N, Ms
!---------------------------------------------------------------------//

! address gas phase
      IF (RVEL_G < 0) THEN   ! inflow (apply BC)
         T_G(IJK) = BC_T_G(BCV)
         IF (NMAX(0) > 0) &
           X_G(IJK,:NMAX(0)) = BC_X_G(BCV,:NMAX(0))
         IF (K_Epsilon) THEN
           K_Turb_G(IJK) = BC_K_Turb_G(BCV)
           E_Turb_G(IJK) = BC_E_Turb_G(BCV)
         ENDIF
         IF (NScalar > 0) THEN
            DO N = 1,NScalar
               Ms = Phase4Scalar(N)
               IF (Ms == 0) &
                  Scalar(IJK, N) = BC_Scalar(BCV,N)
            ENDDO
         ENDIF
      ELSE   ! outflow (apply neighbor fluid)
         T_G(IJK) = T_G(FIJK)
         IF (NMAX(0) > 0) &
            X_G(IJK,:NMAX(0)) = X_G(FIJK,:NMAX(0))
         IF(K_Epsilon) THEN
            K_Turb_G(IJK) = K_Turb_G(FIJK)
            E_Turb_G(IJK) = E_Turb_G(FIJK)
         ENDIF
         IF (NScalar >0) THEN
            DO N = 1, NScalar
               Ms = Phase4Scalar(N)
               IF (Ms == 0) &
                  Scalar(IJK, N) = Scalar(FIJK, N)
            ENDDO
         ENDIF
      ENDIF   ! end gas phase

! address solids phases
      DO M = 1, SMAX
         IF (RVEL_S(M) < 0) THEN   ! inflow (apply BC)
            T_S(IJK,M) = BC_T_S(BCV,M)
            THETA_M(IJK,M) = BC_THETA_M(BCV,M)
            IF (NMAX(M) > 0) &
               X_S(IJK,M,:NMAX(M)) = BC_X_S(BCV,M,:NMAX(M))
            IF (NScalar > 0) THEN
               DO N = 1,NScalar
                  Ms = Phase4Scalar(N)
                  IF (Ms == M) &
                     Scalar(IJK, N) = BC_Scalar(BCV,N)
               ENDDO
            ENDIF
         ELSE   ! outflow (apply neighbor fluid)
            T_S(IJK,M) = T_S(FIJK,M)
            THETA_M(IJK,M) =  THETA_M(FIJK,M)
            IF (NMAX(M) > 0) &
               X_S(IJK,M,:NMAX(M)) = X_S(FIJK,M,:NMAX(M))
            IF (NScalar > 0) THEN
               DO N = 1,NScalar
                  Ms = Phase4Scalar(N)
                  IF (Ms == M) &
                     Scalar(IJK, N) = Scalar(FIJK,N)
               ENDDO
            ENDIF
         ENDIF
      ENDDO


! derived BC field quantity (since it has no explicit user defined BC)
      IF(KT_TYPE_ENUM == GHD_2007) THEN
         nTOT = zero
         DO M = 1, SMAX
             nM = ROP_S(IJK,M)*6d0/ &
                 (PI*D_p(IJK,M)**3*RO_S(IJK,M))
             nTOT = nTOT + nM
             THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) + &
                nM*THETA_M(IJK,M)
         ENDDO
         IF (NTOT > zero) &
            THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) / nTOT
      ENDIF

      RETURN
      END SUBROUTINE SET_PINOUTFLOW
