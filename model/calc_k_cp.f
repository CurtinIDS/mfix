! TO DO:
! 1. Kcp needs to be defined for each solids phase (?).
! 2. The part of Kcp from P_star should be based on the
!    sum of EP_s of close-packed solids.
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine:  CALC_K_cp                                              C
!  Purpose: Calculate and store dPodEp_s                               C
!                                                                      C
!  Notes: MCP must be defined to call this routine.                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-FEB-97   C
!  Reviewer:                                          Date:            C
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

      SUBROUTINE CALC_K_cp(Kcp)

!-----------------------------------------------
! Modules
!----------------------------------------------
      USE compar
      USE constant
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE pscor
      USE rdf
      USE run
      USE sendrecv
      USE solids_pressure
      USE trace
      USE visc_s
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! dPodEP_s
      DOUBLE PRECISION :: Kcp(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: IJK, M
! Other variables
      DOUBLE PRECISION :: Pc, DPcoDEPs, Mu, Mu_b, Mu_zeta, ZETA
      DOUBLE PRECISION :: F2, DF2oDEPs, Pf, Pfmax, N_Pff
! Blend Factor
      Double Precision :: blend
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION :: DZETAoDEPs
!-----------------------------------------------

! initializing
      KCP(:) = ZERO

      IF (MCP == UNDEFINED_I) THEN
! this error should be caught earlier in the routines so that this
! branch should never be entered
         RETURN
      ELSE
! the lowest solids phase index of those solids phases that can close
! pack (i.e. close_packed=T) and the index of the solids phase that is
! used to form the solids correction equation.
         M = MCP
      ENDIF


! by definition M must be close_packed (this is a redundant check)
      IF(CLOSE_PACKED(M)) THEN

!!$omp    parallel do default(shared)                  &
!!$omp    private( IJK, DPcoDEPS, Pc,                  &
!!$omp             Mu, Mu_b, Mu_zeta, ZETA, N_Pff,     &
!!$omp             F2, DF2oDEPs, Pf, Pfmax, blend )

         DO IJK = ijkstart3, ijkend3
            IF(.NOT.WALL_AT(IJK))THEN

               IF (FRICTION) THEN

                  IF ((ONE-EP_G(IJK)).GT.EPS_f_min) THEN
! if friction and sufficiently packed to invoke friction
! ---------------------------------------------------------------->>>

                     IF ((ONE-EP_G(IJK)).GT.(ONE-ep_star_array(ijk))) THEN

! Linearized form of Pc; this is more stable and provides continuous
! function.
                       DPcoDEPS = (to_SI*Fr)*((delta**5)*&
                          (2d0*(ONE-ep_star_array(IJK)-delta) - &
                          2d0*eps_f_min)+&
                          ((ONE-ep_star_array(ijk)-delta)-eps_f_min)*&
                          (5*delta**4))/(delta**10)

                       Pc = (to_SI*Fr)*&
                          ( ((ONE-ep_star_array(IJK)-delta)-&
                          EPS_f_min)**N_Pc )/(delta**D_Pc)
                       Pc = Pc + DPcoDEPS*( (ONE-EP_G(IJK))+delta-&
                          (ONE-ep_star_array(IJK)))

!                        Pc = 1d25*( ((ONE-EP_G(IJK)) - &
!                           (ONE-ep_star_array(ijk)))**10d0)
!                        DPcoDEPS = 1d26*(((ONE-EP_G(IJK)) - &
!                           (ONE-ep_star_array(ijk)))**9d0)
                     ELSE
                        Pc = Fr*(((ONE-EP_G(IJK)) - EPS_f_min)**N_Pc)/&
                           (((ONE-ep_star_array(ijk)) - (ONE-EP_G(IJK))&
                           + SMALL_NUMBER)**D_Pc)
                        DPcoDEPs =&
                           Fr*(((ONE-EP_G(IJK)) - EPS_f_min)**(N_Pc - ONE))&
                           *(N_Pc*((ONE-ep_star_array(ijk)) - (ONE-EP_G(IJK)))&
                           +D_Pc*((ONE-EP_G(IJK)) - EPS_f_min))&
                           / (((ONE-ep_star_array(ijk)) - (ONE-EP_G(IJK)) + &
                           SMALL_NUMBER)**(D_Pc + ONE))
                     ENDIF

                     Mu = (5d0*DSQRT(Pi*Theta_m(IJK,M))&
!QX
                        *D_p(IJK,M)*RO_S(IJK,M))/96d0

                     Mu_b = (256d0*Mu*EP_s(IJK,M)*EP_s(IJK,M)&
                        *G_0(IJK,M,M))/(5d0*Pi)

                     IF (SAVAGE.EQ.1) THEN
                        Mu_zeta =&
                           ((2d0+ALPHA)/3d0)*((Mu/(Eta*(2d0-Eta)*&
                           G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                           G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                           EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b*Eta))
!QX
                        ZETA = ((48d0*Eta*(1d0-Eta)*RO_S(IJK,M)*EP_s(IJK,M)*&
                            EP_s(IJK,M)*G_0(IJK,M,M)*&
                            (Theta_m(IJK,M)**1.5d0))/&
                            (SQRT_Pi*D_p(IJK,M)*2d0*Mu_zeta))**0.5d0

                     ELSEIF (SAVAGE.EQ.0) THEN
                        ZETA = (SMALL_NUMBER +&
                             trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                             trD_s_C(IJK,M))/3.d0))**0.5d0

                     ELSE
                        ZETA = ((Theta_m(IJK,M)/(D_p(IJK,M)*D_p(IJK,M))) +&
                             (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                             trD_s_C(IJK,M))/3.d0)))**0.5d0

                     ENDIF

                     IF (trD_s_C(IJK,M) .GE. ZERO) THEN
                        N_Pff = DSQRT(3d0)/(2d0*Sin_Phi) !dilatation
                     ELSE
                        N_Pff = N_Pf !compaction
                     ENDIF

                     IF ((trD_s_C(IJK,M)/(ZETA*N_Pff*DSQRT(2d0)*&
                          Sin_Phi)) .GT. 1d0) THEN
                        F2 = 0d0
                        DF2oDEPs = ZERO
                     ELSEIF(trD_s_C(IJK,M) == ZERO) THEN
                        F2 = ONE
                        DF2oDEPs = ZERO
                     ELSE
                        F2 = (1d0 - (trD_s_C(IJK,M)/(ZETA*N_Pff*&
                           DSQRT(2d0)*Sin_Phi)))**(N_Pff-1d0)
                        IF (SAVAGE.EQ.1) THEN
                           DF2oDEPs = (N_Pff-1d0)*(F2**(N_Pff-2d0))*&
                              trD_s_C(IJK,M)&
                              *DZETAoDEPs(EP_s(IJK,M), IJK, M)&
                              / (ZETA*ZETA*&
                              N_Pff*DSQRT(2d0)*Sin_Phi)
                        ELSE
                           DF2oDEPs=ZERO
                        ENDIF

                        Pf = Pc*F2
                        Pfmax = Pc*((N_Pf/(N_Pf-1d0))**(N_Pf-1d0))

                        IF (Pf> Pfmax) THEN
                           F2 = (N_Pf/(N_Pf-1d0))**(N_Pf-1d0)
                           DF2oDEPS = ZERO
                        ENDIF
                     ENDIF

! Contributions to Kcp(IJK) from kinetic theory have been left out in
! the expressions below as they cause convergence problems at low solids
! volume fraction
                     Kcp(IJK) = F2*DPcoDEPS + Pc*DF2oDEPS

                  ELSE
! the solids are not sufficiently packed to invoke friction model
                     Kcp(IJK) = ZERO

                  ENDIF   ! end if/else branch (one-ep_g(ijk)>eps_f_min)
! end if friction and sufficiently packed to invoke friction
! ----------------------------------------------------------------<<<

               ELSE ! FRICTION = .FALSE.


                  IF(EP_g(IJK) .LT. ep_g_blend_end(ijk)) THEN
! not friction but the solids are packed so that the plastic pressure
! model is invoked
! ---------------------------------------------------------------->>>

                     Kcp(IJK) = dPodEP_s(EP_s(IJK, M),ep_g_blend_end(ijk))

                     IF(BLENDING_STRESS) THEN
                        blend =  blend_function(IJK)
                        Kcp(IJK) = (1.0d0-blend) * Kcp(IJK)
                     ENDIF
                  ELSE
! the solids are not sufficiently packed to invoke the plastic stress
! model
                     Kcp(IJK) = ZERO
                  ENDIF
! end if not friction but sufficiently packed to invoke plastic pressure
! ----------------------------------------------------------------<<<

               ENDIF   ! end if/else branch if(friction)

            ELSE    ! else branch of if (.not.wall_at(ijk))
               Kcp(IJK) = ZERO
            ENDIF   ! end if/else (.not.wall_at(ijk))

         ENDDO  ! do ijk=ijkstart3,ijkend3
      ENDIF   ! end if (close_packed(m))

      CALL send_recv(Kcp, 2)

      RETURN
      END


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: DZETAoDEPs (EPs,IJK,M)                                    C
!  Purpose: Calculate derivative of zeta                               C
!           w.r.t granular volume fraction                             C
!                                                                      C
!  Author: A. Srivastava                              Date: 8-JUNE-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      DOUBLE PRECISION FUNCTION DZETAoDEPs(EPs, IJK, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE pscor
      USE rdf
      USE run
      USE trace
      USE visc_s
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! solids volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPs
! indices
      INTEGER, INTENT(IN) :: IJK,M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local radial distribution function
      DOUBLE PRECISION :: g0
! Other variables
      DOUBLE PRECISION :: Mu, Mu_b,  DEPs2G_0oDEPs, F1, DF1oDEPs
!-----------------------------------------------

      g0 = G_0(IJK, M, M)
!QX
      Mu = (5d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(IJK,M)*RO_S(IJK,M))/96d0

      Mu_b = (256d0*Mu*EPs*EPs*g0/(5d0*Pi))

      DEPs2G_0oDEPs = EPs*EPs*DG_0DNU(EPs) + 2d0*EPs*g0

      F1 = ((2d0+ALPHA)/3d0)*((2*Mu/(Eta*(2d0-Eta)*&
           g0))*(1d0+1.6d0*Eta*EPs*g0)*(1d0+1.6d0*Eta*(3d0*Eta-2d0)&
           *EPs*g0)+(1.2d0*Mu_b*Eta))

      DF1oDEPs = ((2d0+ALPHA)/3d0)*((2*Mu/(Eta*(2d0-Eta))*&
         ((-DG_0DNU(EPs)/(g0*g0)) + (1.6d0*Eta*(3d0*Eta-1d0))&
        + (64d0*Eta*Eta*(3d0*Eta-2d0)*DEPs2G_0oDEPs/25d0))) +&
        3.2d0*Eta*RO_S(IJK,M)*D_p(IJK,M)*((Theta_m(IJK,M)/Pi)**0.5d0)&
          *DEPs2G_0oDEPs)

      DZETAoDEPs = 0.5d0*((48d0*Eta*(1d0-Eta)*RO_S(IJK,M)*F1*&
              (Theta_m(IJK,M)**1.5d0)/&
           (SQRT_Pi*D_p(IJK,M)*EPs*EPs*g0))**0.5d0)*&
           (F1*DEPs2G_0oDEPs - EPs*Eps*g0*DF1oDEPs)/(F1*F1)


      RETURN
      END

