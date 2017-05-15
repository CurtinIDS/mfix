!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_OUTFLOW(BCV, I1, I2, J1, J2, K1, K2)               C
!  Purpose: Set specified pressure outflow bc for a specified range of C
!           cells                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX                                          C
!                                                                      C
!  Variables modified: I, J, K, RO_g, ROP_g,                           C
!                      EP_g, C
!                      T_g, T_s,  M, ROP_s, U_g, U_s, V_g, V_s,  C
!                      W_g, W_s,
!                                                                      C
!  Local variables: IJK, LFLUID                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_SET_OUTFLOW
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bc
      USE compar        !//d
      USE cutcell
      USE eos, ONLY: EOSG
      USE fldvar
      USE geometry
      USE indices
      USE mflux
      USE param
      USE param1
      USE physprop
      USE quadric
      USE run
      USE scalars
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      indices
      INTEGER          I, J, K, M, NN
!
!                      Local index for boundary cell
      INTEGER          IJK
!
!                      Boundary condition number
      INTEGER          BCV
!
!                      Locall index for a fluid cell near the boundary cell
      INTEGER          LFLUID

      INTEGER          IJKW,IJKWW,IJKS,IJKSS,IJKB

      INTEGER :: BCT1,BCT2,BCT3,BCT4

      LOGICAL :: TEST1,TEST2
!-----------------------------------------------
!
!      print*,'top of cg_set_outflow'
      DO IJK = IJKSTART3, IJKEND3
      IF(INTERIOR_CELL_AT(IJK)) THEN


         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

!//SP Check if current i,j,k resides on this PE
         IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
!        IJK = FUNIJK(I,J,K)
!
! Fluid cell at West
!
!      print*,'west'
          BCV = BC_U_ID(IJK)


          IF(BCV>0)THEN

            IF (BC_TYPE_ENUM(BCV) == CG_PO) THEN

               IF (FLUID_AT(IM_OF(IJK))) THEN
                  LFLUID = IM_OF(IJK)
!            print*,'west treatment:IJK,LFLUID=',IJK,LFLUID
!            read(*,*)
!                  IF (U_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
!                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID)
                     T_G(IJK) = T_G(LFLUID)
                     NN = 1
                     IF (NMAX(0) > 0) THEN
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0))
                        NN = NMAX(0) + 1
                     ENDIF
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK))
!                  ENDIF
                  P_STAR(IJK) = P_STAR(LFLUID)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  DO NN = 1, NScalar
                     M = Phase4Scalar(NN)
                     IF(M == 0)Then
!                      IF (U_G(LFLUID)>=ZERO) THEN
                          Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                      ENDIF
                     Else
!                      IF (U_s(LFLUID, M)>=ZERO) THEN
                          Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                     ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
!                    IF (U_G(LFLUID) >= ZERO) THEN
                     K_Turb_G(IJK) = K_Turb_G(LFLUID)
                     E_Turb_G(IJK) = E_Turb_G(LFLUID)
!                   ENDIF
                  ENDIF


                  DO M = 1, MMAX
                     P_S(IJK,M) = P_S(LFLUID,M)
!                     IF (U_S(LFLUID,M) >= ZERO) THEN
                        ROP_S(IJK,M) = ROP_S(LFLUID,M)
                        T_S(IJK,M) = T_S(LFLUID,M)
!                     ELSE
!                        ROP_S(IJK,M) = ZERO
!                     ENDIF
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
!
                     NN = 1
                     IF (NMAX(M) > 0) THEN
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M))
                        NN = NMAX(M) + 1
                     ENDIF
                  END DO
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!                  IF (ROP_G(IJK) > ZERO) THEN
!                     U_G(IJK) = ROP_G(LFLUID)*U_G(LFLUID)/ROP_G(IJK)
!                  ELSE
!                     U_G(IJK) = ZERO
!                  ENDIF
!                  V_G(IJK) = V_G(LFLUID)
!                  W_G(IJK) = W_G(LFLUID)
!                  Flux_gE(IJK) = Flux_gE(LFLUID)
!                 Flux_gN(IJK) = Flux_gN(LFLUID)
!                 Flux_gT(IJK) = Flux_gT(LFLUID)

!                  IF (MMAX > 0) THEN
!                     WHERE (ROP_S(IJK,:MMAX) > ZERO)
!                        U_S(IJK,:MMAX) = ROP_S(LFLUID,:MMAX)*U_S(LFLUID,:MMAX)/&
!                           ROP_S(IJK,:MMAX)
!                     ELSEWHERE
!                        U_S(IJK,:MMAX) = ZERO
!                    END WHERE
!                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX)
!                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX)
!                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
!                    Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
!                    Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX)
!                  ENDIF
               ENDIF


            ENDIF
         ENDIF

         IJKW = WEST_OF(IJK)
         IJKWW = WEST_OF(IJKW)

         BCT1=BLANK
         BCT2=BLANK
         BCT3=BLANK
         BCT4=BLANK
         BCV = BC_ID(IJK)
         IF(BCV>0) BCT1 = BC_TYPE_ENUM(BCV)
         BCV = BC_U_ID(IJK)
         IF(BCV>0) BCT2 = BC_TYPE_ENUM(BCV)
         BCV = BC_U_ID(IJKW)
         IF(BCV>0) BCT3 = BC_TYPE_ENUM(BCV)
         BCV = BC_ID(IJKW)
         IF(BCV>0) BCT4 = BC_TYPE_ENUM(BCV)

         TEST1= ((BCT1 == CG_PO).AND.(BCT2 /= CG_PO).AND.(BCT3 == CG_PO))

         TEST2 = (BCT2 == CG_PO).AND.(BCT4 /= CG_PO)

!         TEST2 = (FLAG(IJK) == 11)


         IF(TEST1) THEN
            BCV = BC_ID(IJK)

!         IF(TEST1.OR.TEST2) THEN

!            IF(FLUID_AT(IJKW)) THEN

!            ELSEIF(FLUID_AT(IJKWW)) THEN
!               IJKW = IJKWW
!            ENDIF


!         IF((BCT1 == 'CG_PO').AND.(BCT2 /= 'CG_PO').AND.(BCT3 == 'CG_PO')) THEN



!         IF(BCT2 == 'CG_PO')


!            print*,'IJK,IJKW=',IJK,IJKW
!            read(*,*)

            T_G(IJK) = T_G(IJKW)
            NN = 1
            IF (NMAX(0) > 0) THEN
               X_G(IJK,:NMAX(0)) = X_G(IJKW,:NMAX(0))
               NN = NMAX(0) + 1
            ENDIF
            MW_MIX_G(IJK) = MW_MIX_G(IJKW)
            RO_G(IJK) = RO_G(IJKW)
            P_STAR(IJK) = P_STAR(IJKW)
            EP_G(IJK) = EP_G(IJKW)
            ROP_G(IJK) = ROP_G(IJKW)
            DO NN = 1, NScalar
               M = Phase4Scalar(NN)
               IF(M == 0)Then
                  Scalar(IJK, NN) = Scalar(IJKW, NN)
               ELSE
                  Scalar(IJK, NN) = Scalar(IJKW, NN)
               ENDIF
            ENDDO
            IF(K_Epsilon) THEN
               K_Turb_G(IJK) = K_Turb_G(IJKW)
               E_Turb_G(IJK) = E_Turb_G(IJKW)
            ENDIF
!            DO M = 1, MMAX
!               P_S(IJK,M) = P_S(IJKW,M)
!               ROP_S(IJK,M) = ROP_S(IJKW,M)
!               T_S(IJK,M) = T_S(IJKW,M)
!               N = 1
!               IF (NMAX(M) > 0) THEN
!                 X_S(IJK,M,:NMAX(M)) = X_S(IJKW,M,:NMAX(M))
!                 N = NMAX(M) + 1
!               ENDIF
!            ENDDO

            DO M = 1, MMAX
               P_S(IJK,M) = P_S(IJKW,M)
!               IF (U_S(IJKW,M) >= ZERO) THEN
                  ROP_S(IJK,M) = ROP_S(IJKW,M)
                  T_S(IJK,M) = T_S(IJKW,M)
!               ELSE
!                  ROP_S(IJK,M) = ZERO
!               ENDIF
!
               IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
               IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
!
               NN = 1
               IF (NMAX(M) > 0) THEN
                  X_S(IJK,M,:NMAX(M)) = X_S(IJKW,M,:NMAX(M))
                  NN = NMAX(M) + 1
               ENDIF
            END DO
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!            IF (ROP_G(IJK) > ZERO) THEN
!               U_G(IJK) = ROP_G(IJKW)*U_G(IJKW)/ROP_G(IJK)
!            ELSE
!               U_G(IJK) = ZERO
!            ENDIF
!            V_G(IJK) = V_G(IJKW)
!            W_G(IJK) = W_G(IJKW)
!            Flux_gE(IJK) = Flux_gE(IJKW)
!            Flux_gN(IJK) = Flux_gN(IJKW)
!            Flux_gT(IJK) = Flux_gT(IJKW)

!            IF (MMAX > 0) THEN
!               WHERE (ROP_S(IJK,:MMAX) > ZERO)
!                  U_S(IJK,:MMAX) = ROP_S(IJKW,:MMAX)*U_S(IJKW,:MMAX)/&
!                     ROP_S(IJK,:MMAX)
!               ELSEWHERE
!                  U_S(IJK,:MMAX) = ZERO
!               END WHERE
!               V_S(IJK,:MMAX) = V_S(IJKW,:MMAX)
!               W_S(IJK,:MMAX) = W_S(IJKW,:MMAX)
!               Flux_sE(IJK,:MMAX) = Flux_sE(IJKW,:MMAX)
!              Flux_sN(IJK,:MMAX) = Flux_sN(IJKW,:MMAX)
!              Flux_sT(IJK,:MMAX) = Flux_sT(IJKW,:MMAX)
!            ENDIF



!            CYCLE

         ENDIF




!
! Fluid cell at East
!      print*,'east'

          BCV = BC_U_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE_ENUM(BCV) == CG_PO) THEN


               IF (FLUID_AT(IP_OF(IJK))) THEN
                  LFLUID = IP_OF(IJK)
!            print*,'east treatment:IJK,LFLUID=',IJK,LFLUID
!            read(*,*)
!                  IF (U_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
!                     IF (BC_TYPE_ENUM(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID)
                     T_G(IJK) = T_G(LFLUID)
                     NN = 1
                     IF (NMAX(0) > 0) THEN
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0))
                        NN = NMAX(0) + 1
                     ENDIF
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK))
!                  ENDIF
                  P_STAR(IJK) = P_STAR(LFLUID)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  DO NN = 1, NScalar
                     M = Phase4Scalar(NN)
                     IF(M == 0)Then
!                      IF (U_G(LFLUID) <= ZERO) THEN
                          Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                      ENDIF
                     Else
!                      IF (U_s(LFLUID, M) <= ZERO) THEN
                          Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                      ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
!                    IF (U_G(LFLUID) <= ZERO) THEN
                        K_Turb_G(IJK) = K_Turb_G(LFLUID)
                        E_Turb_G(IJK) = E_Turb_G(LFLUID)
!                    ENDIF
                  ENDIF

                  DO M = 1, MMAX
                     P_S(IJK,M) = P_S(LFLUID,M)
!                     IF (U_S(IJK,M) <= ZERO) THEN
                        ROP_S(IJK,M) = ROP_S(LFLUID,M)
                        T_S(IJK,M) = T_S(LFLUID,M)
!                     ELSE
!                        ROP_S(IJK,M) = ZERO
!                     ENDIF
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
                     NN = 1
                     IF (NMAX(M) > 0) THEN
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M))
                        NN = NMAX(M) + 1
                     ENDIF
                  END DO
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!                  IF (U_G(IJK) == UNDEFINED) THEN
!                     IF (ROP_G(IJK) > ZERO) THEN
!                        U_G(IJK) = ROP_G(LFLUID)*U_G(LFLUID)/ROP_G(IJK)
!                     ELSE
!                        U_G(IJK) = ZERO
!                     ENDIF
!                  ENDIF
!                  V_G(IJK) = V_G(LFLUID)
!                  W_G(IJK) = W_G(LFLUID)
!                  Flux_gE(IJK) = Flux_gE(LFLUID)
!                  Flux_gN(IJK) = Flux_gN(LFLUID)
!                  Flux_gT(IJK) = Flux_gT(LFLUID)

!                  DO M = 1, MMAX
!                     IF (U_S(IJK,M) == UNDEFINED) THEN
!                        IF (ROP_S(IJK,M) > ZERO) THEN
!                           U_S(IJK,M) = ROP_S(LFLUID,M)*U_S(LFLUID,M)/ROP_S(IJK&
!                              ,M)
!                        ELSE
!                           U_S(IJK,M) = ZERO
!                        ENDIF
!                     ENDIF
!                     V_S(IJK,M) = V_S(LFLUID,M)
!                     W_S(IJK,M) = W_S(LFLUID,M)
!                     Flux_sE(IJK,M) = Flux_sE(LFLUID,M)
!                    Flux_sN(IJK,M) = Flux_sN(LFLUID,M)
!                    Flux_sT(IJK,M) = Flux_sT(LFLUID,M)
!                  END DO
               ENDIF

            ENDIF
         ENDIF

!
! Fluid cell at South
!
!      print*,'south'
          BCV = BC_V_ID(IJK)
!      print*,ijk,I,J,K,bcv,BC_TYPE_ENUM(BCV)
          IF(BCV>0)THEN

            IF (BC_TYPE_ENUM(BCV) == CG_PO) THEN


               IF (FLUID_AT(JM_OF(IJK))) THEN
                  LFLUID = JM_OF(IJK)
!            print*,'south treatment:IJK,LFLUID=',IJK,LFLUID
!            read(*,*)


!                  IF (V_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
!                     IF (BC_TYPE_ENUM(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID)
                     T_G(IJK) = T_G(LFLUID)
                     NN = 1
                     IF (NMAX(0) > 0) THEN
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0))
                        NN = NMAX(0) + 1
                     ENDIF
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK))
!                  ENDIF
                  P_STAR(IJK) = P_STAR(LFLUID)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  DO NN = 1, NScalar
                     M = Phase4Scalar(NN)
                     IF(M == 0)Then
!                      IF (V_G(LFLUID) >= ZERO) THEN
                          Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                      ENDIF
                     Else
!                      IF (V_s(LFLUID, M) >= ZERO) THEN
                       Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                     ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
!                    IF (V_G(LFLUID) >= ZERO) THEN
                        K_Turb_G(IJK) = K_Turb_G(LFLUID)
                        E_Turb_G(IJK) = E_Turb_G(LFLUID)
!                   ENDIF
                  ENDIF

                  DO M = 1, MMAX
                     P_S(IJK,M) = P_S(LFLUID,M)
!                     IF (V_S(LFLUID,M) >= 0.) THEN
                        ROP_S(IJK,M) = ROP_S(LFLUID,M)
                        T_S(IJK,M) = T_S(LFLUID,M)
!                     ELSE
!                        ROP_S(IJK,M) = ZERO
!                     ENDIF
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
                     NN = 1
                     IF (NMAX(M) > 0) THEN
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M))
                        NN = NMAX(M) + 1
                     ENDIF
                  END DO
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!                  U_G(IJK) = U_G(LFLUID)
!                  IF (ROP_G(IJK) > ZERO) THEN
!                     V_G(IJK) = ROP_G(LFLUID)*V_G(LFLUID)/ROP_G(IJK)
!                  ELSE
!                     V_G(IJK) = ZERO
!                  ENDIF
!                  W_G(IJK) = W_G(LFLUID)
!                  Flux_gE(IJK) = Flux_gE(LFLUID)
!                  Flux_gN(IJK) = Flux_gN(LFLUID)
!                  Flux_gT(IJK) = Flux_gT(LFLUID)

!                  IF (MMAX > 0) THEN
!                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX)
!                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX)
!                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX)
!                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
!                    Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
!                    Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX)
!                  ENDIF
               ENDIF

            ENDIF
         ENDIF

         IJKS = SOUTH_OF(IJK)
         IJKSS = SOUTH_OF(IJKS)
!         print*,'ijk,ijks=',ijk,ijks,JM_OF(IJK)
         BCT1=BLANK
         BCT2=BLANK
         BCT3=BLANK
         BCV = BC_ID(IJK)
!         print*,'bcv=',bcv
         IF(BCV>0) BCT1 = BC_TYPE_ENUM(BCV)
         BCV = BC_V_ID(IJK)
!         print*,'bcv=',bcv
         IF(BCV>0) BCT2 = BC_TYPE_ENUM(BCV)
         BCV = BC_V_ID(IJKS)
!         print*,'bcv=',bcv
         IF(BCV>0) BCT3 = BC_TYPE_ENUM(BCV)

!         IF((BCT1 == 'CG_PO').AND.(BCT2 /= 'CG_PO').AND.(BCT3 == 'CG_PO')) THEN



         TEST1 = (BCT1 == CG_PO).AND.(BCT2 /= CG_PO).AND.(BCT3 == CG_PO)
         TEST2 = (FLAG(IJK) == 11).AND.(BCT2 == CG_PO)

         TEST2 = (FLAG(IJK) == 11)

         IF(TEST1) THEN

            BCV = BC_ID(IJK)

!         IF(TEST1.OR.TEST2) THEN

!            IF(FLUID_AT(IJKS)) THEN

!            ELSEIF(FLUID_AT(IJKSS)) THEN
!               IJKS = IJKSS
!            ENDIF


!            print*,'IJK,IJKS=',IJK,IJKS
!            read(*,*)

            T_G(IJK) = T_G(IJKS)
            NN = 1
            IF (NMAX(0) > 0) THEN
               X_G(IJK,:NMAX(0)) = X_G(IJKS,:NMAX(0))
               NN = NMAX(0) + 1
            ENDIF
            MW_MIX_G(IJK) = MW_MIX_G(IJKS)
            RO_G(IJK) = RO_G(IJKS)
            P_STAR(IJK) = P_STAR(IJKS)
            EP_G(IJK) = EP_G(IJKS)
            ROP_G(IJK) = ROP_G(IJKS)
            DO NN = 1, NScalar
               M = Phase4Scalar(NN)
               IF(M == 0)Then
                  Scalar(IJK, NN) = Scalar(IJKS, NN)
               ELSE
                  Scalar(IJK, NN) = Scalar(IJKS, NN)
               ENDIF
            ENDDO
            IF(K_Epsilon) THEN
               K_Turb_G(IJK) = K_Turb_G(IJKS)
               E_Turb_G(IJK) = E_Turb_G(IJKS)
            ENDIF
!            DO M = 1, MMAX
!               P_S(IJK,M) = P_S(IJKS,M)
!               ROP_S(IJK,M) = ROP_S(IJKS,M)
!               T_S(IJK,M) = T_S(IJKS,M)
!               N = 1
!               IF (NMAX(M) > 0) THEN
!                 X_S(IJK,M,:NMAX(M)) = X_S(IJKS,M,:NMAX(M))
!                 N = NMAX(M) + 1
!               ENDIF
!            ENDDO

            DO M = 1, MMAX
               P_S(IJK,M) = P_S(IJKS,M)
!               IF (V_S(IJKS,M) >= 0.) THEN
                  ROP_S(IJK,M) = ROP_S(IJKS,M)
                  T_S(IJK,M) = T_S(IJKS,M)
!               ELSE
!                  ROP_S(IJK,M) = ZERO
!               ENDIF
!
               IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
               IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
               NN = 1
               IF (NMAX(M) > 0) THEN
                  X_S(IJK,M,:NMAX(M)) = X_S(IJKS,M,:NMAX(M))
                  NN = NMAX(M) + 1
               ENDIF
            END DO
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!            U_G(IJK) = U_G(IJKS)
!            IF (ROP_G(IJK) > ZERO) THEN
!               V_G(IJK) = ROP_G(IJKS)*V_G(IJKS)/ROP_G(IJK)
!            ELSE
!               V_G(IJK) = ZERO
!            ENDIF
!            W_G(IJK) = W_G(IJKS)
!            Flux_gE(IJK) = Flux_gE(IJKS)
!           Flux_gN(IJK) = Flux_gN(IJKS)
!           Flux_gT(IJK) = Flux_gT(IJKS)

!            IF (MMAX > 0) THEN
!               U_S(IJK,:MMAX) = U_S(IJKS,:MMAX)
!               V_S(IJK,:MMAX) = V_S(IJKS,:MMAX)
!               W_S(IJK,:MMAX) = W_S(IJKS,:MMAX)
!               Flux_sE(IJK,:MMAX) = Flux_sE(IJKS,:MMAX)
!              Flux_sN(IJK,:MMAX) = Flux_sN(IJKS,:MMAX)
!              Flux_sT(IJK,:MMAX) = Flux_sT(IJKS,:MMAX)
!            ENDIF




!            CYCLE

         ENDIF



!
! Fluid cell at North
!

!      print*,'north'
          BCV = BC_V_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE_ENUM(BCV) == CG_PO) THEN


               IF (FLUID_AT(JP_OF(IJK))) THEN
                  LFLUID = JP_OF(IJK)

!            print*,'north treatment:IJK,LFLUID=',IJK,LFLUID
!            read(*,*)

!                  IF (V_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
!                     IF (BC_TYPE_ENUM(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID)
                     T_G(IJK) = T_G(LFLUID)
                     NN = 1
                     IF (NMAX(0) > 0) THEN
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0))
                        NN = NMAX(0) + 1
                     ENDIF
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK))
!                  ENDIF
                  P_STAR(IJK) = P_STAR(LFLUID)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  DO NN = 1, NScalar
                     M = Phase4Scalar(NN)
                     IF(M == 0)Then
!                      IF (V_G(LFLUID) <= ZERO) THEN
                         Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                     ENDIF
                     Else
!                      IF (V_s(LFLUID, M) <= ZERO) THEN
                         Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                     ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
!                    IF (V_G(LFLUID) <= ZERO) THEN
                       K_Turb_G(IJK) = K_Turb_G(LFLUID)
                       E_Turb_G(IJK) = E_Turb_G(LFLUID)
!                   ENDIF
                  ENDIF

                  DO M = 1, MMAX
                     P_S(IJK,M) = P_S(LFLUID,M)
!                     IF (V_S(IJK,M) <= ZERO) THEN
                        ROP_S(IJK,M) = ROP_S(LFLUID,M)
                        T_S(IJK,M) = T_S(LFLUID,M)
!                     ELSE
!                        ROP_S(IJK,M) = ZERO
!                     ENDIF
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
                     NN = 1
                     IF (NMAX(M) > 0) THEN
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M))
                        NN = NMAX(M) + 1
                     ENDIF
                  END DO
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!                  U_G(IJK) = U_G(LFLUID)
!                  IF (V_G(IJK) == UNDEFINED) THEN
!                     IF (ROP_G(IJK) > ZERO) THEN
!                        V_G(IJK) = ROP_G(LFLUID)*V_G(LFLUID)/ROP_G(IJK)
!                     ELSE
!                        V_G(IJK) = ZERO
!                     ENDIF
!                  ENDIF
!                  W_G(IJK) = W_G(LFLUID)
!                  Flux_gE(IJK) = Flux_gE(LFLUID)
!                 Flux_gN(IJK) = Flux_gN(LFLUID)
!                 Flux_gT(IJK) = Flux_gT(LFLUID)

!                  IF (MMAX > 0) THEN
!                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX)
!                     WHERE (V_S(IJK,:MMAX) == UNDEFINED) V_S(IJK,:MMAX) = V_S(&
!                        LFLUID,:MMAX)
!                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX)
!                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
!                    Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
!                    Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX)
!                  ENDIF
               ENDIF

            ENDIF
         ENDIF

!
! Fluid cell at Bottom
!
!      print*,'bottom'

          BCV = BC_W_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE_ENUM(BCV) == CG_PO) THEN

               IF (FLUID_AT(KM_OF(IJK))) THEN
                  LFLUID = KM_OF(IJK)
!                  IF (W_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
!                     IF (BC_TYPE_ENUM(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID)
                     T_G(IJK) = T_G(LFLUID)
                     NN = 1
                     IF (NMAX(0) > 0) THEN
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0))
                        NN = NMAX(0) + 1
                     ENDIF
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK))
!                  ENDIF
                  P_STAR(IJK) = P_STAR(LFLUID)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  DO NN = 1, NScalar
                     M = Phase4Scalar(NN)
                     IF(M == 0)Then
!                      IF (W_G(LFLUID) >= ZERO) THEN
                         Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                      ENDIF
                     Else
!                      IF (W_s(LFLUID, M) >= ZERO) THEN
                         Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                     ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
!                    IF (W_G(LFLUID) >= ZERO) THEN
                       K_Turb_G(IJK) = K_Turb_G(LFLUID)
                       E_Turb_G(IJK) = E_Turb_G(LFLUID)
!                   ENDIF
                  ENDIF

                  DO M = 1, MMAX
                     P_S(IJK,M) = P_S(LFLUID,M)
!                     IF (W_S(LFLUID,M) >= 0.) THEN
                        ROP_S(IJK,M) = ROP_S(LFLUID,M)
                        T_S(IJK,M) = T_S(LFLUID,M)
!                     ELSE
                        ROP_S(IJK,M) = ZERO
!                     ENDIF
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
                     NN = 1
                     IF (NMAX(M) > 0) THEN
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M))
                        NN = NMAX(M) + 1
                     ENDIF
                  END DO
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!                  U_G(IJK) = U_G(LFLUID)
!                  V_G(IJK) = V_G(LFLUID)
!                  IF (ROP_G(IJK) > ZERO) THEN
!                     W_G(IJK) = ROP_G(LFLUID)*W_G(LFLUID)/ROP_G(IJK)
!                  ELSE
!                     W_G(IJK) = ZERO
!                  ENDIF
!                  Flux_gE(IJK) = Flux_gE(LFLUID)
!                 Flux_gN(IJK) = Flux_gN(LFLUID)
!                 Flux_gT(IJK) = Flux_gT(LFLUID)

!                  IF (MMAX > 0) THEN
!                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX)
!                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX)
!                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX)
!                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
!                    Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
!                    Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX)
!                  ENDIF
               ENDIF


            ENDIF
         ENDIF

         IJKB = BOTTOM_OF(IJK)
         BCT1=BLANK
         BCT2=BLANK
         BCT3=BLANK
         BCV = BC_ID(IJK)
         IF(BCV>0) BCT1 = BC_TYPE_ENUM(BCV)
         BCV = BC_W_ID(IJK)
         IF(BCV>0) BCT2 = BC_TYPE_ENUM(BCV)
         BCV = BC_W_ID(IJKB)
         IF(BCV>0) BCT3 = BC_TYPE_ENUM(BCV)

         IF((BCT1 == CG_PO).AND.(BCT2 /= CG_PO).AND.(BCT3 == CG_PO)) THEN
!            print*,'IJK,IJKB=',IJK,IJKB
!            read(*,*)

            BCV = BC_ID(IJK)

            T_G(IJK) = T_G(IJKB)
            NN = 1
            IF (NMAX(0) > 0) THEN
               X_G(IJK,:NMAX(0)) = X_G(IJKB,:NMAX(0))
               NN = NMAX(0) + 1
            ENDIF
            MW_MIX_G(IJK) = MW_MIX_G(IJKB)
            RO_G(IJK) = RO_G(IJKB)
            P_STAR(IJK) = P_STAR(IJKB)
            EP_G(IJK) = EP_G(IJKB)
            ROP_G(IJK) = ROP_G(IJKB)
            DO NN = 1, NScalar
               M = Phase4Scalar(NN)
               IF(M == 0)Then
                  Scalar(IJK, NN) = Scalar(IJKB, NN)
               ELSE
                  Scalar(IJK, NN) = Scalar(IJKB, NN)
               ENDIF
            ENDDO
            IF(K_Epsilon) THEN
               K_Turb_G(IJK) = K_Turb_G(IJKB)
               E_Turb_G(IJK) = E_Turb_G(IJKB)
            ENDIF
!            DO M = 1, MMAX
!               P_S(IJK,M) = P_S(IJKB,M)
!               ROP_S(IJK,M) = ROP_S(IJKB,M)
!               T_S(IJK,M) = T_S(IJKB,M)
!               N = 1
!               IF (NMAX(M) > 0) THEN
!                 X_S(IJK,M,:NMAX(M)) = X_S(IJKB,M,:NMAX(M))
!                 N = NMAX(M) + 1
!               ENDIF
!            ENDDO

            DO M = 1, MMAX
               P_S(IJK,M) = P_S(IJKB,M)
!               IF (W_S(IJKB,M) >= 0.) THEN
                  ROP_S(IJK,M) = ROP_S(IJKB,M)
                  T_S(IJK,M) = T_S(IJKB,M)
!               ELSE
!                  ROP_S(IJK,M) = ZERO
!               ENDIF
!
               IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
               IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
               NN = 1
               IF (NMAX(M) > 0) THEN
                  X_S(IJK,M,:NMAX(M)) = X_S(IJKB,M,:NMAX(M))
                  NN = NMAX(M) + 1
               ENDIF
            END DO
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!            U_G(IJK) = U_G(IJKB)
!            V_G(IJK) = V_G(IJKB)
!            IF (ROP_G(IJK) > ZERO) THEN
!               W_G(IJK) = ROP_G(IJKB)*W_G(IJKB)/ROP_G(IJK)
!            ELSE
!               W_G(IJK) = ZERO
!            ENDIF
!            Flux_gE(IJK) = Flux_gE(IJKB)
!           Flux_gN(IJK) = Flux_gN(IJKB)
!           Flux_gT(IJK) = Flux_gT(IJKB)

!            IF (MMAX > 0) THEN
!               U_S(IJK,:MMAX) = U_S(IJKB,:MMAX)
!               V_S(IJK,:MMAX) = V_S(IJKB,:MMAX)
!               W_S(IJK,:MMAX) = W_S(IJKB,:MMAX)
!               Flux_sE(IJK,:MMAX) = Flux_sE(IJKB,:MMAX)
!              Flux_sN(IJK,:MMAX) = Flux_sN(IJKB,:MMAX)
!              Flux_sT(IJK,:MMAX) = Flux_sT(IJKB,:MMAX)
!            ENDIF

!            CYCLE

         ENDIF


!
! Fluid cell at Top
!
!      print*,'top'
          BCV = BC_W_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE_ENUM(BCV) == CG_PO) THEN


               IF (FLUID_AT(KP_OF(IJK))) THEN
                  LFLUID = KP_OF(IJK)
!                  IF (W_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN
!                     IF (BC_TYPE_ENUM(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID)
                     T_G(IJK) = T_G(LFLUID)
                     NN = 1
                     IF (NMAX(0) > 0) THEN
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0))
                        NN = NMAX(0) + 1
                     ENDIF
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK))
!                  ENDIF
                  P_STAR(IJK) = P_STAR(LFLUID)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  DO NN = 1, NScalar
                     M = Phase4Scalar(NN)
                     IF(M == 0)Then
!                      IF (W_G(LFLUID) <= ZERO) THEN
                         Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                     ENDIF
                     Else
!                      IF (W_s(LFLUID, M) <= ZERO) THEN
                         Scalar(IJK, NN) = Scalar(LFLUID, NN)
!                      ENDIF
                     Endif
                  END DO

                  IF(K_Epsilon) THEN
!                    IF (W_G(LFLUID) <= ZERO) THEN
                       K_Turb_G(IJK) = K_Turb_G(LFLUID)
                       E_Turb_G(IJK) = E_Turb_G(LFLUID)
!                   ENDIF
                  ENDIF

                  DO M = 1, MMAX
                     P_S(IJK,M) = P_S(LFLUID,M)
!                     IF (W_S(IJK,M) <= ZERO) THEN
                        ROP_S(IJK,M) = ROP_S(LFLUID,M)
                        T_S(IJK,M) = T_S(LFLUID,M)
!                     ELSE
!                        ROP_S(IJK,M) = ZERO
!                     ENDIF
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M)
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M)
                     NN = 1
                     IF (NMAX(M) > 0) THEN
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M))
                        NN = NMAX(M) + 1
                     ENDIF
                  END DO
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
!                  U_G(IJK) = U_G(LFLUID)
!                  V_G(IJK) = V_G(LFLUID)
!                  IF (W_G(IJK) == UNDEFINED) THEN
!                     IF (ROP_G(IJK) > ZERO) THEN
!                        W_G(IJK) = ROP_G(LFLUID)*W_G(LFLUID)/ROP_G(IJK)
!                     ELSE
!                        W_G(IJK) = ZERO
!                     ENDIF
!                  ENDIF
!                  Flux_gE(IJK) = Flux_gE(LFLUID)
!                 Flux_gN(IJK) = Flux_gN(LFLUID)
!                 Flux_gT(IJK) = Flux_gT(LFLUID)

!                  IF (MMAX > 0) THEN
!                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX)
!                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX)
!                     WHERE (W_S(IJK,:MMAX) == UNDEFINED) W_S(IJK,:MMAX) = W_S(&
!                        LFLUID,:MMAX)
!                     Flux_sE(IJK,:MMAX) = Flux_sE(LFLUID,:MMAX)
!                    Flux_sN(IJK,:MMAX) = Flux_sN(LFLUID,:MMAX)
!                    Flux_sT(IJK,:MMAX) = Flux_sT(LFLUID,:MMAX)
!                  ENDIF
               ENDIF

            ENDIF
         ENDIF

      ENDIF
      END DO

!      print*,'bottom of cg_set_outflow'
      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CG_SET_OUTFLOW
