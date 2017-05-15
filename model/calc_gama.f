!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_GAMA                                               C
!                                                                      C
!  Purpose: Calculate gas-solids heat transfer coefficients            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_GAMA(M)

! Modules
!---------------------------------------------------------------------//
      USE sendrecv, only: send_recv
      USE energy, only: gama_gs
      USE usr_prop, only: usr_gama, calc_usr_prop
      USE usr_prop, only: gassolids_heattransfer
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M
!---------------------------------------------------------------------//

      IF (USR_GAMA(M)) THEN
         CALL CALC_USR_PROP(GasSolids_HeatTransfer,lm=M)
      ELSE
         CALL CALC_DEFAULT_GAMA_GS(M)
      ENDIF
      call send_recv(GAMA_GS,2)

      RETURN
      END SUBROUTINE CALC_GAMA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DEFAULT_GAMA_GS                                    C
!  Purpose: Compute the default value for the gas-solids heat transfer C
!  coefficients.                                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-JUL-92  C
!                                                                      C
!  Literature/Document References:                                     C
!  Gunn, D. J., (1978), International Journal of Heat and Mass         C
!    Transfer, Vol. 21, p 467-476.                                     C
!  Bird, Stewart, and Lightfoot (1960, p.663)                          C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DEFAULT_GAMA_GS(M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE energy, only: gama_gs
      USE fldvar, only: u_g, v_g, w_g, u_s, v_s, w_s
      USE fldvar, only: ep_g, ro_g, d_p, ep_s
      Use fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      USE functions, only: fluidorp_flow_at
      use functions, only: im_of, jm_of, km_of
      use indices, only: i_of
      USE param1, only: zero, half, one, small_number, large_number
      USE physprop, only: mu_g, k_g, c_pg
      USE rxns, only: r_phase
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, IJK, IMJK, IJMK, IJKM
! Cell center value of U_g, Vg, and Wg 
      DOUBLE PRECISION :: UGC, VGC, WGC
! Cell center value of U_sm, V_sm, W_sm
      DOUBLE PRECISION :: USCM, VSCM, WSCM
! Gas-solids relative velocity
      DOUBLE PRECISION :: VREL
! function of Prandtl number, Pr^(1/3)
      DOUBLE PRECISION :: Pr1o3
! Reynolds number, Re
      DOUBLE PRECISION :: Re
! EP_g^2
      DOUBLE PRECISION :: EP_g2
! a factor
      DOUBLE PRECISION :: FAC
! index for storing interphase mass transfer coefficients in the upper
! triangle of the matrix. 
      INTEGER :: LM
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF (FLUIDorP_FLOW_AT(IJK)) THEN
            I = I_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            EP_G2 = EP_G(IJK)*EP_G(IJK)

! Calculate Prandtl number to the 1/3 power
            if(K_G(IJK) > ZERO) then
              PR1O3 = (C_PG(IJK)*MU_G(IJK)/K_G(IJK))**(1.D0/3.D0)
            else
              PR1O3 = LARGE_NUMBER
            endif

! Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! Calculate the magnitude of gas-solids relative velocity
            VREL=SQRT((UGC-USCM)**2+(VGC-VSCM)**2+(WGC-WSCM)**2)

            if(MU_G(IJK) > ZERO)then
               RE = EP_G(IJK)*D_P(IJK,M)*VREL*RO_G(IJK)/MU_G(IJK)
            else
               RE = LARGE_NUMBER
            endif

! Calculate gas-solids heat transfer coefficient (Gunn 1978)
            GAMA_GS(IJK,M) = ((7.D0 - 10.D0*EP_G(IJK)+5.D0*EP_G2)*&
                              (ONE+0.7D0*RE**0.2D0*PR1O3)+&
                              (1.33D0 - 2.4D0*EP_G(IJK)+1.2D0*EP_G2)*&
                              RE**0.7D0*PR1O3)*(K_G(IJK)/D_P(IJK,M))*&
                              (6.D0*EP_S(IJK,M)/D_P(IJK,M))

! Correct the heat transfer coefficient for transpiration
! Bird, Stewart, and Lightfoot (1960, p.663)
            IF (GAMA_GS(IJK,M) > SMALL_NUMBER) THEN
! Only the effect of gas-solids heat transfer is accounted for so find
! corresponding lm index for gas-solids mass transfer
               LM = 1+ (M-1)*M/2
               FAC = R_PHASE(IJK,LM)*C_PG(IJK)/GAMA_GS(IJK,M)
               IF (ABS(FAC) < 0.1D0) THEN
                  GAMA_GS(IJK,M) = GAMA_GS(IJK,M)/&
                                   (ONE+FAC/2.D0+FAC*FAC/6.D0)
               ELSE
                  IF (R_PHASE(IJK,LM) > ZERO) THEN
                     GAMA_GS(IJK,M) = R_PHASE(IJK,LM)*C_PG(IJK)*&
                                      EXP((-FAC))/(ONE - EXP((-FAC)))
                  ELSE
                     GAMA_GS(IJK,M) = R_PHASE(IJK,LM)*C_PG(IJK)/&
                                      (EXP(FAC) - ONE)
                  ENDIF
               ENDIF
            ENDIF   ! end if gama_gs(ijk,M) > 0

         ENDIF   ! end if (fluidorp_flow_at(ijk)
      ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)


      RETURN
      END SUBROUTINE CALC_DEFAULT_GAMA_GS
