!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: B_m_P_star_e(IJK, IER)                                 C
!  Purpose: Determine p_star when there is an interface at east        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-AUG-96   C
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
!
      SUBROUTINE B_M_P_STAR_E(B_M, IJK)
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
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE bodyforce
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKE
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average density
      DOUBLE PRECISION ROPGA
!
!                      LHS (similar to A_m)
      DOUBLE PRECISION A
!
!                      RHS (similar to b_m)
      DOUBLE PRECISION BB
!
!                      b_m
      DOUBLE PRECISION b_m
!
!                      sum of ep_s
      DOUBLE PRECISION Eps
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sdp, Sdps
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION V0, Vmt, Vbf
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IJKE = EAST_OF(IJK)
!
      A = ZERO
      BB = ZERO
      EPS = ZERO
!
      DO M = 1, MMAX
         IF (CLOSE_PACKED(M)) THEN
            EPGA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
!
!         Surface forces
!
!         Pressure term
            SDP = -P_SCALE*EPGA*(P_G(IJKE)-P_G(IJK))*AYZ(IJK)
            SDPS = -EPGA*(P_S(IJKE,M)-P_S(IJK,M))*AYZ(IJK)
!
!           Shear stress terms
!
!         Volumetric forces
            ROPGA = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)
!
!         Previous time step
            V0 = AVG_X(ROP_SO(IJK,M),ROP_SO(IJKE,M),I)*ODT
!
!         Interphase mass transfer
            VMT = AVG_X(SUM_R_S(IJK,M),SUM_R_S(IJKE,M),I)
!
!         Body force
            VBF = ROPGA*BFX_S(IJK,M)
!
!         Collect the terms
            A = A - ((V0 + ZMAX(VMT))*VOL_U(IJK))
            BB=BB-(SDP+SDPS+((V0+ZMAX((-VMT)))*U_SO(IJK,M)+VBF)*VOL_U(IJK))
            EPS = EPS + EPGA
         ENDIF
      END DO
      B_M = -((BB - A)/(EPS*AYZ(IJK))-P_STAR(IJK))
!
      RETURN
      END SUBROUTINE B_M_P_STAR_E
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: B_m_P_star_n(IJK, IER)                                 C
!  Purpose: Determine p_star when there is an interface at north       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-AUG-96   C
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
!
      SUBROUTINE B_M_P_STAR_N(B_M, IJK)
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
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE bodyforce
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKN
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average density
      DOUBLE PRECISION ROPGA
!
!                      LHS (similar to A_m)
      DOUBLE PRECISION A
!
!                      RHS (similar to b_m)
      DOUBLE PRECISION BB
!
!                      b_m
      DOUBLE PRECISION b_m
!
!                      sum of ep_s
      DOUBLE PRECISION Eps
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sdp, Sdps
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION V0, Vmt, Vbf
!
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IJKN = NORTH_OF(IJK)
!
      A = ZERO
      BB = ZERO
      EPS = ZERO
!
      DO M = 1, MMAX
         IF (CLOSE_PACKED(M)) THEN
            EPGA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
!
!         Surface forces
!
!         Pressure term
            SDP = -P_SCALE*EPGA*(P_G(IJKN)-P_G(IJK))*AXZ(IJK)
            SDPS = -EPGA*(P_S(IJKN,M)-P_S(IJK,M))*AXZ(IJK)
!
!           Shear stress terms
!
!         Volumetric forces
            ROPGA = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)
!
!         Previous time step
            V0 = AVG_Y(ROP_SO(IJK,M),ROP_SO(IJKN,M),J)*ODT
!
!         Interphase mass transfer
            VMT = AVG_Y(SUM_R_S(IJK,M),SUM_R_S(IJKN,M),J)
!
!         Body force
            VBF = ROPGA*BFY_S(IJK,M)
!
!         Collect the terms
            A = A - ((V0 + ZMAX(VMT))*VOL_V(IJK))
            BB=BB-(SDP+SDPS+((V0+ZMAX((-VMT)))*V_SO(IJK,M)+VBF)*VOL_V(IJK))
            EPS = EPS + EPGA
         ENDIF
      END DO
      B_M = -((BB - A)/(EPS*AXZ(IJK))-P_STAR(IJK))
!
      RETURN
      END SUBROUTINE B_M_P_STAR_N
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: B_m_P_star_t(IJK, IER)                                 C
!  Purpose: Determine p_star when there is an interface at top         C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-AUG-96   C
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
!
      SUBROUTINE B_M_P_STAR_T(B_M, IJK)
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
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE bodyforce
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKT
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average density
      DOUBLE PRECISION ROPGA
!
!                      LHS (similar to A_m)
      DOUBLE PRECISION A
!
!                      RHS (similar to b_m)
      DOUBLE PRECISION BB
!
!                      b_m
      DOUBLE PRECISION b_m
!
!                      sum of ep_s
      DOUBLE PRECISION Eps
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sdp, Sdps
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION V0, Vmt, Vbf
!
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IJKT = TOP_OF(IJK)
!
      A = ZERO
      BB = ZERO
      EPS = ZERO
!
      DO M = 1, MMAX
         IF (CLOSE_PACKED(M)) THEN
            EPGA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
!
!         Surface forces
!
!         Pressure term
            SDP = -P_SCALE*EPGA*(P_G(IJKT)-P_G(IJK))*AXY(IJK)
            SDPS = -EPGA*(P_S(IJKT,M)-P_S(IJK,M))*AXY(IJK)
!
!           Shear stress terms
!
!         Volumetric forces
            ROPGA = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)
!
!         Previous time step
            V0 = AVG_Z(ROP_SO(IJK,M),ROP_SO(IJKT,M),K)*ODT
!
!         Interphase mass transfer
            VMT = AVG_Z(SUM_R_S(IJK,M),SUM_R_S(IJKT,M),K)
!
!         Body force
            VBF = ROPGA*BFZ_S(IJK,M)
!
!         Collect the terms
            A = A - ((V0 + ZMAX(VMT))*VOL_W(IJK))
            BB=BB-(SDP+SDPS+((V0+ZMAX((-VMT)))*W_SO(IJK,M)+VBF)*VOL_W(IJK))
            EPS = EPS + EPGA
         ENDIF
      END DO
      B_M = -((BB - A)/(EPS*AXY(IJK))-P_STAR(IJK))
!
      RETURN
      END SUBROUTINE B_M_P_STAR_T
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: B_m_P_star_w(b_m, IJK, IER)                            C
!  Purpose: Determine p_star when there is an interface at west        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-AUG-96   C
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
!
      SUBROUTINE B_M_P_STAR_W(B_M, IJK)
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
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE bodyforce
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          IM, J, K, IJK, IJKW, IMJK
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average density
      DOUBLE PRECISION ROPGA
!
!                      LHS (similar to A_m)
      DOUBLE PRECISION A
!
!                      RHS (similar to b_m)
      DOUBLE PRECISION BB
!
!                      b_m
      DOUBLE PRECISION b_m
!
!                      sum of ep_s
      DOUBLE PRECISION Eps
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sdp, Sdps
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION V0, Vmt, Vbf
!
!-----------------------------------------------

      IM = IM1(I_OF(IJK))
      J = J_OF(IJK)
      K = K_OF(IJK)
      IJKW = WEST_OF(IJK)
      IMJK = IM_OF(IJK)
!
      A = ZERO
      BB = ZERO
      EPS = ZERO
!
      DO M = 1, MMAX
         IF (CLOSE_PACKED(M)) THEN
            EPGA = AVG_X(EP_S(IJKW,M),EP_S(IJK,M),IM)
!
!         Surface forces
!
!         Pressure term
            SDP = -P_SCALE*EPGA*(P_G(IJK)-P_G(IJKW))*AYZ(IMJK)
            SDPS = -EPGA*(P_S(IJK,M)-P_S(IJKW,M))*AYZ(IMJK)
!
!           Shear stress terms
!
!         Volumetric forces
            ROPGA = AVG_X(ROP_S(IJKW,M),ROP_S(IJK,M),IM)
!
!         Previous time step
            V0 = AVG_X(ROP_SO(IJKW,M),ROP_SO(IJK,M),IM)*ODT
!
!         Interphase mass transfer
            VMT = AVG_X(SUM_R_S(IJKW,M),SUM_R_S(IJK,M),IM)
!
!         Body force
            VBF = ROPGA*BFX_S(IJKW,M)
!
!         Collect the terms
            A = A - ((V0 + ZMAX(VMT))*VOL_U(IJKW))
            BB=BB-(SDP+SDPS+((V0+ZMAX((-VMT)))*U_SO(IMJK,M)+VBF)*VOL_U(IJKW))
            EPS = EPS + EPGA
         ENDIF
      END DO
      B_M = -(((-(BB - A)/(EPS*AYZ(IMJK))))-P_STAR(IJK))
!
      RETURN
      END SUBROUTINE B_M_P_STAR_W
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: B_m_P_star_s(b_m, IJK, IER)                            C
!  Purpose: Determine p_star when there is an interface at south       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-AUG-96   C
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
!
      SUBROUTINE B_M_P_STAR_S(B_M, IJK)
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
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE bodyforce
      USE fun_avg
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
!                      Indices
      INTEGER          I, JM, K, IJK, IJKS, IJMK
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average density
      DOUBLE PRECISION ROPGA
!
!                      LHS (similar to A_m)
      DOUBLE PRECISION A
!
!                      RHS (similar to b_m)
      DOUBLE PRECISION BB
!
!                      b_m
      DOUBLE PRECISION b_m
!
!                      sum of ep_s
      DOUBLE PRECISION Eps
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sdp, Sdps
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION V0, Vmt, Vbf
!-----------------------------------------------

      I = I_OF(IJK)
      JM = JM1(J_OF(IJK))
      K = K_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJMK = JM_OF(IJK)
!
      A = ZERO
      BB = ZERO
      EPS = ZERO
!
      DO M = 1, MMAX
         IF (CLOSE_PACKED(M)) THEN
            EPGA = AVG_Y(EP_S(IJKS,M),EP_S(IJK,M),JM)
!
!         Surface forces
!
!         Pressure term
            SDP = -P_SCALE*EPGA*(P_G(IJK)-P_G(IJKS))*AXZ(IJK)
            SDPS = -EPGA*(P_S(IJK,M)-P_S(IJKS,M))*AXZ(IJK)
!
!           Shear stress terms
!
!         Volumetric forces
            ROPGA = AVG_Y(ROP_S(IJKS,M),ROP_S(IJK,M),JM)
!
!         Previous time step
            V0 = AVG_Y(ROP_SO(IJKS,M),ROP_SO(IJK,M),JM)*ODT
!
!         Interphase mass transfer
            VMT = AVG_Y(SUM_R_S(IJKS,M),SUM_R_S(IJK,M),JM)
!
!         Body force
            VBF = ROPGA*BFY_S(IJKS,M)
!
!         Collect the terms
            A = A - ((V0 + ZMAX(VMT))*VOL_V(IJKS))
            BB=BB-(SDP+SDPS+((V0+ZMAX((-VMT)))*V_SO(IJMK,M)+VBF)*VOL_V(IJKS))
            EPS = EPS + EPGA
         ENDIF
      END DO
      B_M = -(((-(BB - A)/(EPS*AXZ(IJK))))-P_STAR(IJK))
!
      RETURN
      END SUBROUTINE B_M_P_STAR_S
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: B_m_P_star_b(b_m, IJK, IER)                            C
!  Purpose: Determine p_star when there is an interface at bottom      C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-AUG-96   C
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
!
      SUBROUTINE B_M_P_STAR_B(B_M, IJK)
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
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE rxns
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE bodyforce
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, KM, IJK, IJKB, IJKM
!
!                      Phase index
      INTEGER          M
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average density
      DOUBLE PRECISION ROPGA
!
!                      LHS (similar to A_m)
      DOUBLE PRECISION A
!
!                      RHS (similar to b_m)
      DOUBLE PRECISION BB
!
!                      b_m
      DOUBLE PRECISION b_m
!
!                      sum of ep_s
      DOUBLE PRECISION Eps
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sdp, Sdps
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION V0, Vmt, Vbf
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      KM = KM1(K_OF(IJK))
      IJKB = BOTTOM_OF(IJK)
      IJKM = KM_OF(IJK)
!
      A = ZERO
      BB = ZERO
      EPS = ZERO
!
      DO M = 1, MMAX
         IF (CLOSE_PACKED(M)) THEN
            EPGA = AVG_Z(EP_S(IJKB,M),EP_S(IJK,M),KM)
!
!         Surface forces
!
!         Pressure term
            SDP = -P_SCALE*EPGA*(P_G(IJK)-P_G(IJKB))*AXY(IJK)
            SDPS = -EPGA*(P_S(IJK,M)-P_S(IJKB,M))*AXY(IJK)
!
!           Shear stress terms
!
!         Volumetric forces
            ROPGA = AVG_Z(ROP_S(IJKB,M),ROP_S(IJK,M),KM)
!
!         Previous time step
            V0 = AVG_Z(ROP_SO(IJKB,M),ROP_SO(IJK,M),KM)*ODT
!
!         Interphase mass transfer
            VMT = AVG_Z(SUM_R_S(IJKB,M),SUM_R_S(IJK,M),KM)
!
!         Body force
            VBF = ROPGA*BFZ_S(IJKB,M)
!
!         Collect the terms
            A = A - ((V0 + ZMAX(VMT))*VOL_W(IJKB))
            BB=BB-(SDP+SDPS+((V0+ZMAX((-VMT)))*W_SO(IJKM,M)+VBF)*VOL_W(IJKB))
            EPS = EPS + EPGA
         ENDIF
      END DO
      B_M = -(((-(BB - A)/(EPS*AXY(IJK))))-P_STAR(IJK))
!
      RETURN
      END SUBROUTINE B_M_P_STAR_B

