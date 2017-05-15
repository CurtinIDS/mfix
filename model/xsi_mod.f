!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XSI                                                    C
!  Purpose: Determine convection weighting factors for higher order    C
!  discretization.                                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-MAR-97   C
!                                                                      C
!  Comments: contains following functions subroutines & functions:     C
!    calc_xsi, cxs, dw,                                                C
!    xsi_func                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE XSI
      IMPLICIT NONE

      CONTAINS


      SUBROUTINE CALC_XSI(DISCR, PHI, U, V, W, XSI_E, XSI_N, XSI_T, &
                          incr)

! Modules
!---------------------------------------------------------------------//
      USE chischeme, only: chi_flag
      USE chischeme, only: chi_e, chi_n, chi_t

      USE compar, only: ijkstart3, ijkend3

      USE discretization, only: phi_c_of
      USE discretization, only: superbee
      USE discretization, only: chi_smart, smart
      USE discretization, only: ultra_quick
      USE discretization, only: quickest
      USE discretization, only: chi_muscl, muscl
      USE discretization, only: vanleer
      USE discretization, only: minmod
      USE discretization, only: central_scheme

      USE functions, only: east_of, west_of, north_of, south_of
      USE functions, only: top_of, bottom_of

      USE geometry, only: do_k
      USE geometry, only: odx, odx_e, ody, ody_n, odz, odz_t
      USE geometry, only: ox

      USE indices, only: i_of, j_of, k_of
      USE indices, only: im1, ip1, jm1, jp1, km1, kp1

      USE param, only: dimension_3

      USE param1, only: zero

      USE run, only: shear, dt

      USE sendrecv, only: send_recv

      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: ival, flush_err_msg
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! discretization method
      INTEGER, INTENT(IN) :: DISCR
! convected quantity
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Convection weighting factors
      DOUBLE PRECISION, INTENT(OUT) :: XSI_e(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_n(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_t(DIMENSION_3)
! shear indicator
      INTEGER, INTENT(IN) :: incr

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IJKC, IJKD, IJKU, I, J, K
!
      DOUBLE PRECISION :: PHI_C
! down wind factor
      DOUBLE PRECISION :: dwf
! Courant number
      DOUBLE PRECISION :: cf
! cell widths for QUICKEST
      DOUBLE PRECISION :: oDXc, oDXuc, oDYc, oDYuc, oDZc, oDZuc
!---------------------------------------------------------------------//

      IF (SHEAR) THEN
! calculate XSI_E, XSI_N, XSI_T when periodic shear BCs are used
         call CXS(incr, DISCR, U, V, W, PHI, XSI_E, XSI_N, XSI_T)
      ELSE

       SELECT CASE (DISCR)                    !first order upwinding
       CASE (:1)

!$omp    parallel do default(none) private(IJK) shared(ijkstart3, ijkend3, u, v, w, xsi_e, xsi_n, xsi_t, do_k)
          DO IJK = ijkstart3, ijkend3
             XSI_E(IJK) = XSI_func(U(IJK),ZERO)
             XSI_N(IJK) = XSI_func(V(IJK),ZERO)
             IF (DO_K) XSI_T(IJK) = XSI_func(W(IJK),ZERO)
          ENDDO


       CASE (2)                               !Superbee

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF)
          DO IJK = ijkstart3, ijkend3

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = SUPERBEE(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = SUPERBEE(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = SUPERBEE(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (3)                               !SMART

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF)
          DO IJK = ijkstart3, ijkend3

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             if(Chi_flag) then
                DWF = Chi_SMART(PHI_C, CHI_e(IJK))
             else
                DWF = SMART(PHI_C)
             endif
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             if(Chi_flag) then
                DWF = Chi_SMART(PHI_C, CHI_n(IJK))
             else
                DWF = SMART(PHI_C)
             endif
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                if(Chi_flag) then
                   DWF = Chi_SMART(PHI_C, CHI_t(IJK))
                else
                   DWF = SMART(PHI_C)
                endif
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (4)                               !ULTRA-QUICK

!!!$omp    parallel do private(IJK, I,J,K, IJKC,IJKD,IJKU, PHI_C,DWF,CF)
          DO IJK = ijkstart3, ijkend3

             I = I_OF(IJK)
             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(U(IJK))*DT*ODX_E(I)
             DWF = ULTRA_QUICK(PHI_C,CF)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             J = J_OF(IJK)
             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(V(IJK))*DT*ODY_N(J)
             DWF = ULTRA_QUICK(PHI_C,CF)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                K = K_OF(IJK)
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                CF = ABS(W(IJK))*DT*OX(I)*ODZ_T(K)
                DWF = ULTRA_QUICK(PHI_C,CF)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (5)                               !QUICKEST

!!!$omp    parallel do &
!!!$omp&   private(IJK,I,J,K, IJKC,IJKD,IJKU, &
!!!$omp&           ODXC,ODXUC, PHI_C,CF,DWF, &
!!!$omp&           ODYC,ODYUC,  ODZC,ODZUC )
          DO IJK = ijkstart3, ijkend3

             I = I_OF(IJK)
             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
                ODXC = ODX(I)
                ODXUC = ODX_E(IM1(I))
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
                ODXC = ODX(IP1(I))
                ODXUC = ODX_E(IP1(I))
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(U(IJK))*DT*ODX_E(I)
             DWF = QUICKEST(PHI_C,CF,ODXC,ODXUC,ODX_E(I))
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             J = J_OF(IJK)
             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
                ODYC = ODY(J)
                ODYUC = ODY_N(JM1(J))
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
                ODYC = ODY(JP1(J))
                ODYUC = ODY_N(JP1(J))
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             CF = ABS(V(IJK))*DT*ODY_N(J)
             DWF = QUICKEST(PHI_C,CF,ODYC,ODYUC,ODY_N(J))
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                K = K_OF(IJK)
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                   ODZC = ODZ(K)
                   ODZUC = ODZ_T(KM1(K))
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                   ODZC = ODZ(KP1(K))
                   ODZUC = ODZ_T(KP1(K))
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                CF = ABS(W(IJK))*DT*OX(I)*ODZ_T(K)
                DWF = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ_T(K))
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (6)                               !MUSCL

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF )
          DO IJK = ijkstart3, ijkend3

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             if(Chi_flag) then
                DWF = Chi_MUSCL(PHI_C, CHI_e(IJK))
             else
                DWF = MUSCL(PHI_C)
             endif
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             if(Chi_flag) then
                DWF = Chi_MUSCL(PHI_C, CHI_n(IJK))
             else
                DWF = MUSCL(PHI_C)
             endif
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                if(Chi_flag) then
                   DWF = Chi_MUSCL(PHI_C, CHI_t(IJK))
                else
                   DWF = MUSCL(PHI_C)
                endif
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (7)                               !Van Leer

!!!$omp    parallel do private( IJK, IJKC,IJKD,IJKU,  PHI_C,DWF )
          DO IJK = ijkstart3, ijkend3

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = VANLEER(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = VANLEER(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = VANLEER(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (8)                               !Minmod

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF )
          DO IJK = ijkstart3, ijkend3

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = MINMOD(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = MINMOD(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF

                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = MINMOD(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE (9)                               ! Central

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF)
          DO IJK = ijkstart3, ijkend3

             IF (U(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = EAST_OF(IJK)
                IJKU = WEST_OF(IJKC)
             ELSE
                IJKC = EAST_OF(IJK)
                IJKD = IJK
                IJKU = EAST_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_E(IJK) = XSI_func(U(IJK),DWF)

             IF (V(IJK) >= ZERO) THEN
                IJKC = IJK
                IJKD = NORTH_OF(IJK)
                IJKU = SOUTH_OF(IJKC)
             ELSE
                IJKC = NORTH_OF(IJK)
                IJKD = IJK
                IJKU = NORTH_OF(IJKC)
             ENDIF
             PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
             DWF = CENTRAL_SCHEME(PHI_C)
             XSI_N(IJK) = XSI_func(V(IJK),DWF)

             IF (DO_K) THEN
                IF (W(IJK) >= ZERO) THEN
                   IJKC = IJK
                   IJKD = TOP_OF(IJK)
                   IJKU = BOTTOM_OF(IJKC)
                ELSE
                   IJKC = TOP_OF(IJK)
                   IJKD = IJK
                   IJKU = TOP_OF(IJKC)
                ENDIF
                PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD))
                DWF = CENTRAL_SCHEME(PHI_C)
                XSI_T(IJK) = XSI_func(W(IJK),DWF)
             ENDIF
          ENDDO


       CASE DEFAULT                           !Error
! should never hit this
          CALL INIT_ERR_MSG("CALC_XSI")
          WRITE(ERR_MSG, 1100) IVAL(DISCR)
 1100 FORMAT('ERROR 1100: Invalid DISCRETIZE= ', A,' The check_data ',&
         'routines should',/, 'have already caught this error and ',&
         'prevented the simulation from ',/,'running. Please notify ',&
         'the MFIX developers.')
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
          CALL FINL_ERR_MSG

       END SELECT

      ENDIF   ! end if/else shear

      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)

      RETURN
      END SUBROUTINE CALC_XSI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CXS                                                     C
!  Purpose: Calculates XSI_E, XSI_N,XSI_T when using periodic shear    C
!  BC's.  Uses true velocities.                                        C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CXS(INCR, DISCR, U, V, WW, PHI, XSI_E, XSI_N, XSI_T)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE vshear
      USE compar
      USE sendrecv
      USE functions
      IMPLICIT NONE

! Dummy argments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: incr
      INTEGER, INTENT(IN) :: DISCR
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: WW(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: XSI_T(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: VV(DIMENSION_3)
      DOUBLE PRECISION :: SRT
      DOUBLE PRECISION :: DWFE,DWFN,DWFT
      DOUBLE PRECISION :: PHICU, PHIDU, PHIUU
      DOUBLE PRECISION :: PHICV, PHIDV, PHIUV
      DOUBLE PRECISION :: PHICW, PHIDW, PHIUW
      INTEGER :: IJK, IJKC, IJKD, IJKU, I
!---------------------------------------------------------------------//

      IF (INCR .eq. 2) THEN                   !V momentum balance
         SRT=(2d0*V_sh/XLENGTH)

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU,DWFE,DWFN,DWFT,&
!!!$omp&        PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV,PHICW,PHIDW,PHIUW,I)&
!!!$omp&  shared(DISCR)
         DO IJK = ijkstart3, ijkend3

            I=I_OF(IJK)
            VV(IJK)=V(IJK)+VSH(IJK)

            IF (U(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = EAST_OF(IJK)
               IJKU = WEST_OF(IJKC)
               PHICU=PHI(IJKC)
               PHIDU=PHI(IJKD)+SRT*1d0/oDX_E(I)
               PHIUU=PHI(IJKU)-SRT*1d0/oDX_E(IM1(I))
            ELSE
               IJKC = EAST_OF(IJK)
               IJKD = IJK
               IJKU = EAST_OF(IJKC)
               PHICU=PHI(IJKC)+SRT*1d0/oDX_E(I)
               PHIDU=PHI(IJKD)
               PHIUU=PHI(IJKU)+SRT*1d0/oDX_E(I)&
                    +SRT*1d0/oDX_E(IP1(I))
            ENDIF

            IF (VV(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = NORTH_OF(IJK)
               IJKU= SOUTH_OF(IJKC)
               PHICV=PHI(IJKC)
               PHIDV=PHI(IJKD)
               PHIUV=PHI(IJKU)
            ELSE
               IJKC= NORTH_OF(IJK)
               IJKD= IJK
               IJKU= NORTH_OF(IJKC)
               PHICV=PHI(IJKC)
               PHIDV=PHI(IJKD)
               PHIUV=PHI(IJKU)
            ENDIF

            IF (DO_K) THEN
               IF (WW(IJK) >= ZERO) THEN
                  IJKC = IJK
                  IJKD = TOP_OF(IJK)
                  IJKU = BOTTOM_OF(IJKC)
               ELSE
                  IJKC = TOP_OF(IJK)
                  IJKD = IJK
                  IJKU = TOP_OF(IJKC)
               ENDIF
               PHICW=PHI(IJKC)
               PHIDW=PHI(IJKD)
               PHIUW=PHI(IJKU)
            ELSE
               PHICW=0d0
               PHIDW=0d0
               PHIUW=0d0
               DWFT=0d0
            ENDIF

            CALL DW(U(IJK), VV(IJK), WW(IJK), IJK, PHICU, PHIDU, PHIUU, &
                    PHICV, PHIDV, PHIUV, PHICW, PHIDW, PHIUW, &
                    DWFE, DWFN, DWFT, DISCR)

            XSI_E(IJK) = XSI_func(U(IJK),DWFE)
            XSI_N(IJK) = XSI_func(VV(IJK),DWFN)
            XSI_T(IJK) = XSI_func(WW(IJK),DWFT)

            VV(IJK)=V(IJK)-VSH(IJK)
         ENDDO


      ELSEIF (INCR .eq. 1) THEN                  !u momentum balance

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU,DWFE,DWFN,DWFT,&
!!!$omp&        PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV,PHICW,PHIDW,PHIUW,I)&
!!!$omp&  shared(DISCR)
         DO IJK = ijkstart3, ijkend3

            VV(IJK)=V(IJK)+VSHE(IJK)
            IF (U(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = EAST_OF(IJK)
               IJKU = WEST_OF(IJKC)
               PHICU=PHI(IJKC)
               PHIDU=PHI(IJKD)
               PHIUU=PHI(IJKU)
            ELSE
               IJKC = EAST_OF(IJK)
               IJKD = IJK
               IJKU = EAST_OF(IJKC)
               PHICU=PHI(IJKC)
               PHIDU=PHI(IJKD)
               PHIUU=PHI(IJKU)
            ENDIF

            IF (VV(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = NORTH_OF(IJK)
               IJKU= SOUTH_OF(IJKC)
               PHICV=PHI(IJKC)
               PHIDV=PHI(IJKD)
               PHIUV=PHI(IJKU)
            ELSE
               IJKC= NORTH_OF(IJK)
               IJKD= IJK
               IJKU= NORTH_OF(IJKC)
               PHICV=PHI(IJKC)
               PHIDV=PHI(IJKD)
               PHIUV=PHI(IJKU)
            ENDIF

            IF (DO_K) THEN
               IF (WW(IJK) >= ZERO) THEN
                  IJKC = IJK
                  IJKD = TOP_OF(IJK)
                  IJKU = BOTTOM_OF(IJKC)
               ELSE
                  IJKC = TOP_OF(IJK)
                  IJKD = IJK
                  IJKU = TOP_OF(IJKC)
               ENDIF
               PHICW=PHI(IJKC)
               PHIDW=PHI(IJKD)
               PHIUW=PHI(IJKU)
            ELSE
               PHICW=0d0
               PHIDW=0d0
               PHIUW=0d0
               DWFT=0d0
            ENDIF

            CALL DW(U(IJK), VV(IJK), WW(IJK), IJK, PHICU, PHIDU, PHIUU, &
                    PHICV, PHIDV, PHIUV, PHICW, PHIDW, PHIUW, &
                    DWFE, DWFN, DWFT, DISCR)

            XSI_E(IJK) = XSI_func(U(IJK),DWFE)
            XSI_N(IJK) = XSI_func(VV(IJK),DWFN)
            XSI_T(IJK) = XSI_func(WW(IJK),DWFT)

            VV(IJK)=V(IJK)-VSHE(IJK)
         ENDDO


      ELSEIF (INCR .eq. 0) THEN                  !scalars and w momentum

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU,DWFE,DWFN,DWFT,&
!!!$omp&        PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV,PHICW,PHIDW,PHIUW,I)&
!!!$omp&  shared(DISCR)
         DO IJK = ijkstart3, ijkend3

            VV(IJK)=V(IJK)+VSH(IJK)
            IF (U(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = EAST_OF(IJK)
               IJKU = WEST_OF(IJKC)
               PHICU=PHI(IJKC)
               PHIDU=PHI(IJKD)
               PHIUU=PHI(IJKU)
            ELSE
               IJKC = EAST_OF(IJK)
               IJKD = IJK
               IJKU = EAST_OF(IJKC)
               PHICU=PHI(IJKC)
               PHIDU=PHI(IJKD)
               PHIUU=PHI(IJKU)
            ENDIF

            IF (VV(IJK) >= ZERO) THEN
               IJKC = IJK
               IJKD = NORTH_OF(IJK)
               IJKU= SOUTH_OF(IJKC)
               PHICV=PHI(IJKC)
               PHIDV=PHI(IJKD)
               PHIUV=PHI(IJKU)
            ELSE
               IJKC= NORTH_OF(IJK)
               IJKD= IJK
               IJKU= NORTH_OF(IJKC)
               PHICV=PHI(IJKC)
               PHIDV=PHI(IJKD)
               PHIUV=PHI(IJKU)
            ENDIF

            IF (DO_K) THEN
               IF (WW(IJK) >= ZERO) THEN
                  IJKC = IJK
                  IJKD = TOP_OF(IJK)
                  IJKU = BOTTOM_OF(IJKC)
               ELSE
                  IJKC = TOP_OF(IJK)
                  IJKD = IJK
                  IJKU = TOP_OF(IJKC)
               ENDIF
               PHICW=PHI(IJKC)
               PHIDW=PHI(IJKD)
               PHIUW=PHI(IJKU)
            ELSE
               PHICW=0d0
               PHIDW=0d0
               PHIUW=0d0
               DWFT=0d0
            ENDIF

            CALL DW(U(IJK), VV(IJK), WW(IJK), IJK, PHICU, PHIDU, PHIUU, &
                    PHICV, PHIDV, PHIUV, PHICW, PHIDW, PHIUW, &
                    DWFE, DWFN, DWFT, DISCR)

            XSI_E(IJK) = XSI_func(U(IJK),DWFE)
            XSI_N(IJK) = XSI_func(VV(IJK),DWFN)
            XSI_T(IJK) = XSI_func(WW(IJK),DWFT)

            VV(IJK)=V(IJK)-VSH(IJK)
         ENDDO

      ELSE
         write(*,*) 'INCR ERROR'
      ENDIF

      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)
      RETURN
      END SUBROUTINE CXS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DW                                                      C
!  Purpose: Calculates DWFs for various discretization methods when    C
!  period shear BCs are used                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DW(UU, VV, WW, IJK, PHICU, PHIDU, PHIUU, &
                    PHICV, PHIDV, PHIUV, &
                    PHICW, PHIDW, PHIUW, &
                    DWFE, DWFN, DWFT, DISCR)

! Modules
!---------------------------------------------------------------------//
      USE compar
      USE discretization
      USE exit, only: mfix_exit
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE run
      USE sendrecv
      USE vshear
      IMPLICIT NONE

! Dummy argments
!---------------------------------------------------------------------//
!
      DOUBLE PRECISION, INTENT(IN) :: UU, VV, WW
!
      INTEGER, INTENT(IN) :: IJK
!
      DOUBLE PRECISION, INTENT(IN) :: PHICU, PHIDU, PHIUU
      DOUBLE PRECISION, INTENT(IN) :: PHICV, PHIDV, PHIUV
      DOUBLE PRECISION, INTENT(IN) :: PHICW, PHIDW, PHIUW
!
      DOUBLE PRECISION, INTENT(OUT) :: DWFE, DWFN, DWFT
! discretization method
      INTEGER, INTENT(IN) :: DISCR

! Local variables
!---------------------------------------------------------------------//
!
      DOUBLE PRECISION :: PHI_C
! cell widths for QUICKEST
      DOUBLE PRECISION oDXc, oDXuc, oDYc, oDYuc, oDZc, oDZuc
      DOUBLE PRECISION CF, DZUC
! Indices
      INTEGER :: I, J, K
!---------------------------------------------------------------------//

      SELECT CASE (DISCR)                        !first order upwinding
      CASE (:1)
         DWFE=0d0
         DWFN=0d0
         DWFT=0d0


      CASE (2)                                   !Superbee
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         DWFE = SUPERBEE(PHI_C)
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         DWFN = SUPERBEE(PHI_C)
         IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            DWFT = SUPERBEE(PHI_C)
         ENDIF


      CASE (3)                               !SMART
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         DWFE = SMART(PHI_C)
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         DWFN = SMART(PHI_C)
         IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            DWFT = SMART(PHI_C)
         ENDIF


      CASE (4)                               !ULTRA-QUICK
         I=I_OF(IJK)
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         CF = ABS(UU)*DT*ODX_E(I)
         DWFE = ULTRA_QUICK(PHI_C,CF)

         J=J_OF(IJK)
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         CF = ABS(VV)*DT*ODY_N(J)
         DWFN = ULTRA_QUICK(PHI_C,CF)

         IF (DO_K) THEN
            K=K_OF(IJK)
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            CF = ABS(WW)*DT*OX(I)*ODZ_T(K)
            DWFT = ULTRA_QUICK(PHI_C,CF)
         ENDIF


      CASE(5)                               !QUICKEST
         I=I_OF(IJK)
         IF (UU >= ZERO) THEN
            ODXC = ODX(I)
            ODXUC = ODX_E(IM1(I))
         ELSE
            ODXC = ODX(IP1(I))
            ODXUC = ODX_E(IP1(I))
         ENDIF
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         CF = ABS(UU)*DT*ODX_E(I)
         DWFE = QUICKEST(PHI_C,CF,ODXC,ODXUC,ODX_E(I))

         J=J_OF(IJK)
         IF (VV >= ZERO) THEN
            ODYC = ODY(J)
            ODYUC = ODY_N(JM1(J))
         ELSE
            ODYC = ODY(JP1(J))
            ODYUC = ODY_N(JP1(J))
         ENDIF
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         CF = ABS(VV)*DT*ODY_N(J)
         DWFN = QUICKEST(PHI_C,CF,ODYC,ODYUC,ODY_N(J))

         IF (DO_K) THEN
            K=K_OF(IJK)
            IF (WW >= ZERO) THEN
               ODZC = ODZ(K)
               ODZUC = ODZ_T(KM1(K))
            ELSE
               ODZC = ODZ(KP1(K))
               DZUC = ODZ_T(KP1(K))
            ENDIF
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            CF = ABS(WW)*DT*OX(I)*ODZ_T(K)
            DWFT = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ_T(K))
         ENDIF


      CASE (6)                               !MUSCL
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         DWFE = MUSCL(PHI_C)
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         DWFN = MUSCL(PHI_C)
         IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            DWFT = MUSCL(PHI_C)
         ENDIF

      CASE (7)                               !Van Leer
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         DWFE = VANLEER(PHI_C)
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         DWFN = VANLEER(PHI_C)
         IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            DWFT = VANLEER(PHI_C)
         ENDIF


      CASE (8)                               !Minmod
         PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU)
         DWFE = MINMOD(PHI_C)
         PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)
         DWFN = MINMOD(PHI_C)
         IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            DWFT = MINMOD(PHI_C)
         ENDIF


      CASE DEFAULT                           !Error
         WRITE (*,*) 'DISCRETIZE = ', DISCR, ' not supported.'
         CALL MFIX_EXIT(myPE)

      END SELECT
      RETURN
      END SUBROUTINE DW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: Xsi_func                                                  C
!  Purpose: Special function for xsi that should be similar to:        C
!      xsi(v,dw) = merge( v >= 0, dwf, one-dwf)                        C
!                                                                      C
!  Slight difference when v is exactly zero but may not be             C
!  significant.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION XSI_func(XXXv,XXXdwf)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXv, XXXdwf
      XSI_func = (sign(1d0,(-XXXv))+1d0)/(2d0) + &
         sign(1d0,XXXv)*XXXdwf
      END FUNCTION XSI_func

      END MODULE XSI
