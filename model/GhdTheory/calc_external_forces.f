!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: calc_external_forces(IER)                              C
!  Purpose: Calculate the 3-components of external forces at cell      C
!           faces                                                      C
!                                                                      C
!  Author: S. Benyahia                              Date: 13-MAY-09    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:   JoiX  JoiY   JoiZ                             C
!                                                                      C
!     Local variables: all terms in mass flux                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_EXTERNAL_FORCES()
!
!-----------------------------------------------
!     Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE ghdtheory
      USE physprop
      USE run
      USE constant
      USE drag
      USE bc
      use scales
      use bodyforce
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
!                      Index
      INTEGER          IJK, I, J, K
!
!                      Index
      INTEGER          IJKE, IJKN, IJKT
!
!                      Solids phase
      INTEGER          M ,IM
!
!     number densities to compute del(Nj)
      DOUBLE PRECISION NjC, NjE, NjN, NjT
!
!     mass, volume of species
      DOUBLE PRECISION Mj, Vj
!
!     drag force on a particle
      DOUBLE PRECISION dragFc, dragFe, dragFn, dragFt
!
!     pressure terms in mass mobility
      DOUBLE PRECISION PGE, PGN, PGT, SDPx, SDPy, SDPz

!     off-diagonal terms for HYS drag relation
      DOUBLE PRECISION  avgDragx, avgDragy, avgDragz
!
!-----------------------------------------------
!     Function subroutines
!-----------------------------------------------

      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)


          IF ( FLUID_AT(IJK) ) THEN

               IJKE = EAST_OF(IJK)
               IJKN = NORTH_OF(IJK)
               IJKT = TOP_OF(IJK)

! pressure term (no need to have it inside smax loop)
               PGE = P_G(IJKE)
               PGN = P_G(IJKN)
               PGT = P_G(IJKT)
               IF (CYCLIC_X_PD) THEN
                 IF (CYCLIC_AT_E(IJK)) PGE = P_G(IJKE) - DELP_X
               ENDIF
               IF (CYCLIC_Y_PD) THEN
                 IF (CYCLIC_AT_N(IJK)) PGN = P_G(IJKN) - DELP_Y
               ENDIF
               IF (CYCLIC_Z_PD) THEN
                 IF (CYCLIC_AT_T(IJK)) PGT = P_G(IJKT) - DELP_Z
               ENDIF
               SDPx = -P_SCALE*(PGE - P_G(IJK))  * oDX_E(I)
               SDPy = -P_SCALE*(PGN - P_G(IJK))  * oDY_N(J)
               SDPz = -P_SCALE*(PGT - P_G(IJK))  * (oX_E(I)*oDZ_T(K))

               DO M = 1, SMAX
                 Vj = (PI/6.d0)*D_P(IJK,M)**3 ! particle volume
                 Mj = Vj * RO_S(IJK,M)            ! particle mass

                 NjC = ROP_s(IJK,M) / Mj
                 NjE = ROP_S(IJKE,M) / Mj
                 NjN = ROP_S(IJKN,M) / Mj
                 NjT = ROP_S(IJKT,M) / Mj

! drag force on a particle in -x -y -z directions
                 dragFc = zero
                 dragFe = zero
                 dragFn = zero
                 dragFt = zero
                 if(NjC > zero) dragFc = F_GS(IJK ,M)/NjC
                 if(NjE > zero) dragFe = F_GS(IJKE,M)/NjE
                 if(NjN > zero) dragFn = F_GS(IJKN,M)/NjN
                 if(NjT > zero) dragFt = F_GS(IJKT,M)/NjT

                 dragFxflux(IJK,M) = AVG_X(dragFc,dragFe,I) * (U_g(IJK) - U_s(IJK,M))
                 dragFyflux(IJK,M) = AVG_Y(dragFc,dragFn,J) * (V_g(IJK) - V_s(IJK,M))
                 dragFzflux(IJK,M) = AVG_Z(dragFc,dragFt,K) * (W_g(IJK) - W_s(IJK,M))

                 dragFx(IJK,M) = AVG_X(dragFc,dragFe,I) * (U_g(IJK))
                 dragFy(IJK,M) = AVG_Y(dragFc,dragFn,J) * (V_g(IJK))
                 dragFz(IJK,M) = AVG_Z(dragFc,dragFt,K) * (W_g(IJK))

                 beta_cell_X(IJK,M)=AVG_X(dragFc,dragFe,I)
                 beta_cell_Y(IJK,M)=AVG_Y(dragFc,dragFn,J)
                 beta_cell_Z(IJK,M)=AVG_Z(dragFc,dragFt,K)

                 IF(DRAG_TYPE_ENUM == HYS)THEN
                    DO IM=1,SMAX
                       IF(IM /= M)THEN

! HYS additional drag force on a particle in -x -y -z directions
                          dragFc = zero
                          dragFe = zero
                          dragFn = zero
                          dragFt = zero
                          if(NjC > zero) dragFc = beta_ij(IJK,M,IM) /NjC
                          if(NjE > zero) dragFe = beta_ij(IJKE,M,IM)/NjE
                          if(NjN > zero) dragFn = beta_ij(IJKN,M,IM)/NjN
                          if(NjT > zero) dragFt = beta_ij(IJKT,M,IM)/NjT

                          avgDragx = AVG_X(dragFc,dragFe,I)
                          dragFx(IJK,M) = dragFx(IJK,M) - avgDragx*(U_g(IJK))
                          dragFxflux(IJK,M) = dragFxflux(IJK,M) - avgDragx*(U_g(IJK) - U_s(IJK,IM))
                          beta_ij_cell_X(IJK,M,IM)= AVG_X(dragFc,dragFe,I)

                          avgDragy = AVG_Y(dragFc,dragFn,J)
                          dragFy(IJK,M) = dragFy(IJK,M) - avgDragy*(V_g(IJK))
                          dragFyflux(IJK,M) = dragFyflux(IJK,M) - avgDragy*(V_g(IJK) - V_s(IJK,IM))
                          beta_ij_cell_Y(IJK,M,IM)=AVG_Y(dragFc,dragFn,J)

                          avgDragz = AVG_Z(dragFc,dragFt,K)
                          dragFz(IJK,M) = dragFz(IJK,M) - avgDragz*(W_g(IJK))
                          dragFzflux(IJK,M) = dragFzflux(IJK,M) - avgDragz*(W_g(IJK) - W_s(IJK,IM))
                          beta_ij_cell_Z(IJK,M,IM)=AVG_Z(dragFc,dragFt,K)

                       ENDIF
                    ENDDO
                 ENDIF


                 FiXvel(IJK,M) =  (Mj * BFX_S(IJK,M)+dragFx(IJK,M) +Vj*SDPx)
                 FiYvel(IJK,M) =  (Mj * BFY_S(IJK,M)+dragFy(IJK,M) +Vj*SDPy)
                 FiZvel(IJK,M) =  (Mj * BFZ_S(IJK,M)+dragFz(IJK,M) +Vj*SDPz)

                 FiX(IJK,M) =  (Mj * BFX_S(IJK,M)+dragFxflux(IJK,M) +Vj*SDPx)
                 FiY(IJK,M) =  (Mj * BFY_S(IJK,M)+dragFyflux(IJK,M) +Vj*SDPy)
                 FiZ(IJK,M) =  (Mj * BFZ_S(IJK,M)+dragFzflux(IJK,M) +Vj*SDPz)

                 FiMinusDragX(IJK,M) =  (Mj * BFX_S(IJK,M) + Vj*SDPx)
                 FiMinusDragY(IJK,M) =  (Mj * BFY_S(IJK,M) + Vj*SDPy)
                 FiMinusDragZ(IJK,M) =  (Mj * BFZ_S(IJK,M) + Vj*SDPz)
               ENDDO
          ENDIF     ! Fluid_at
 200  CONTINUE     ! outer IJK loop

      RETURN
      END SUBROUTINE CALC_EXTERNAL_FORCES
