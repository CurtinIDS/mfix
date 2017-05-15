! TODO:
! 1. Account for solids-solids transfer terms for scalars.  Will need
!    to pass the exchange coefficients to this routine (unlike momentum
!    eqs)
! 2. repetitive code may be eliminated by defining and using functions


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PARTIAL_ELIM_S                                          C
!  Purpose: Do partial elimination for scalar quantities               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments:                                                           C
!     Currently does not evaluate partial elimination terms between    C
!     unlike solids phases.  The exchange coefficient must be hard     C
!     coded explicitly for this case - it is not a passed dummy        C
!     argument                                                         C
!                                                                      C
!  Literature/Document References:                                     C
!     Spalding, D.B., 1980, Numerical Computation of Multi-phase       C
!        fluid flow and heat transfer, Recent Advances in Numerical    C
!        Methods in Fluids, C. Taylor et al., eds, Pineridge Press     C
!     Syamlal, M., 1998, MFIX Documentation: Numerical Technique,      C
!        Technical Note, DOE/MC-31346-5824, NTIS/DE98002029, EG&G      C
!        Technical Services of West Virginia, Inc., Morgantown, WV.    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTIAL_ELIM_S(VAR_G, VAR_S, VXF, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE physprop
      USE indices
      USE compar
      USE drag
      USE fldvar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! gas phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_g(DIMENSION_3)
! solids phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_s(DIMENSION_3, DIMENSION_M)
! Volume x gas-solids transfer coefficient
      DOUBLE PRECISION, INTENT(IN) :: VxF(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT
      INTEGER :: L, M, LP
!
      DOUBLE PRECISION :: SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
! a0, b0 etc.
      DOUBLE PRECISION :: a(0:DIMENSION_M), BB(0:DIMENSION_M),&
                          F(0:DIMENSION_M,0:DIMENSION_M),&
                          Saxf(0:DIMENSION_M)
!-----------------------------------------------

!!$omp  parallel do private( IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, den) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IJKW = WEST_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            DO M=0, MMAX
               A(M)=A_M(IJK,0,M)
               BB(M)=B_M(IJK,M)

               if (m .ne. 0) then
                  F(M,0)=-VXF(IJK,M)
                  F(0,M)=F(M,0)
               else
                  F(0,0) = ZERO
               endif

               DO L =1, MMAX
                  IF ((L .NE. M) .AND. (M .NE. 0)) THEN
                     F(M,L) = ZERO  !insert solids-solids exchange coefficients here
                  ELSE
                     F(M,L) = ZERO
                  ENDIF
                  F(L,M) = ZERO
               ENDDO

               IF (M == 0 ) THEN
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_G(IJKE)+A_M(IJK,west,M)*VAR_G(IJKW&
                       )+A_M(IJK,north,M)*VAR_G(IJKN)+A_M(IJK,south,M)*VAR_G(IJKS))
               ELSE
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_S(IJKE,M)+A_M(IJK,west,M)*VAR_S(&
                       IJKW,M)+A_M(IJK,north,M)*VAR_S(IJKN,M)+A_M(IJK,south,M)*VAR_S(&
                       IJKS,M))
               ENDIF

               IF (DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  IJKT = TOP_OF(IJK)
                  IF ( M ==0) THEN
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_G(IJKT)+A_M(IJK,bottom,M)*&
                           VAR_G(IJKB))
                  ELSE
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_S(IJKT,M)+A_M(IJK,bottom,M)*&
                           VAR_S(IJKB,M))
                  ENDIF
               ENDIF
            ENDDO

            DO M=0,MMAX
               SUM_A = ZERO
               SUM_B = ZERO

               DO L =0,MMAX
                  SUM_A_LPRIME = ZERO
                  SUM_B_LPRIME = ZERO

                  DO LP=0,MMAX
                     IF ( LP .NE. M) THEN
                        SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
                        IF (LP == 0) THEN
                           SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
                        ELSE
                           SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
                        ENDIF
                     ENDIF
                  ENDDO

                  DEN = A(L) + SUM_A_LPRIME + F(L,M)
                  IF ( DEN .NE. ZERO) THEN
                     SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
                     SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
                             )/DEN
                  ENDIF
               ENDDO

               A_M(IJK,0,M)= SUM_A+A(M)
               B_M(IJK,M)  = SUM_B+BB(M)
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE PARTIAL_ELIM_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PARTIAL_ELIM_IA                                         C
!  Purpose: Do partial elimination for granular temperature coupling   C
!                                                                      C
!                                                                      C
!  Author: J. Galvin                                  Date: 21-MAY-05  C
!  Reviewer: S. Benyahia                              Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!    Iddir, Y.H., Modeling of the multiphase mixture of particles      C
!       using the kinetic theory approach, PhD Thesis, Illinois        C
!       Institute of Technology, Chicago, Illinois, 2004               C
!    Iddir, Y.H., & H. Arastoopour, Modeling of multitype particle     C
!       flow using the kinetic theory approach, AIChE J., Vol 51,      C
!       No 6, June 2005                                                C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTIAL_ELIM_IA(VAR_S, VXTCSS, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE physprop
      USE indices
      USE compar
      USE drag
      USE fldvar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_s(DIMENSION_3, DIMENSION_M)
! Volume x solids-solids transfer coefficient
      DOUBLE PRECISION, INTENT(IN) :: VxTCss (DIMENSION_3, DIMENSION_LM)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT
      INTEGER :: L, M, LP, LM
!
      DOUBLE PRECISION :: SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
! a0, b0 etc.
      DOUBLE PRECISION :: a(DIMENSION_M), BB(DIMENSION_M),&
                          F(DIMENSION_M,DIMENSION_M),&
                          Saxf(DIMENSION_M)
!-----------------------------------------------

!!$omp  parallel do private( IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, den) &
!!$omp&  schedule(static)

      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN
            IJKW = WEST_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)

            F(:,:) = ZERO

            DO M = 1, MMAX
               A(M)=A_M(IJK,0,M)
               BB(M)=B_M(IJK,M)

               DO L = 1, MMAX
                  IF ( L .NE. M ) THEN
                     LM = FUNLM(L,M)
                     F(M,L)=-VxTCSS(IJK,LM)
                  ENDIF
                  F(L,M) = F(M,L)
               ENDDO

               SAXF(M) = -(A_M(IJK,east,M)*VAR_S(IJKE,M)+&
                  A_M(IJK,west,M)*VAR_S(IJKW,M)+&
                  A_M(IJK,north,M)*VAR_S(IJKN,M)+&
                  A_M(IJK,south,M)*VAR_S(IJKS,M))

               IF (DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  IJKT = TOP_OF(IJK)
                  SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_S(IJKT,M)+&
                     A_M(IJK,bottom,M)*VAR_S(IJKB,M))
               ENDIF
            ENDDO

            DO M = 1, MMAX
               SUM_A = ZERO
               SUM_B = ZERO

               DO L = 1, MMAX
                  SUM_A_LPRIME = ZERO
                  SUM_B_LPRIME = ZERO

                  DO LP = 1, MMAX
                     IF ( LP .NE. M) THEN
                        SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
                        SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
                     ENDIF
                  ENDDO

                  DEN = A(L) + SUM_A_LPRIME + F(L,M)

                  IF ( DEN .NE. ZERO) THEN
                     SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
                     SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+&
                        SUM_B_LPRIME)/DEN
                  ENDIF
               ENDDO

               A_M(IJK,0,M)= SUM_A+A(M)
               B_M(IJK,M)  = SUM_B+BB(M)
            ENDDO

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE PARTIAL_ELIM_IA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PARTIAL_ELIM_U                                          C
!  Purpose: Do partial elimination for X vector quantities             C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments:                                                           C
!     Unlike gas-solids exchange coefficient (multiplied by volume),   C
!     the solids-solids exchange coefficient is not as a passed dummy  C
!     argument to this routine but simply re-evaluated here.           C
!                                                                      C
!  Literature/Document References:                                     C
!     Spalding, D.B., 1980, Numerical Computation of Multi-phase       C
!        fluid flow and heat transfer, Recent Advances in Numerical    C
!        Methods in Fluids, C. Taylor et al., eds, Pineridge Press     C
!     Syamlal, M., 1998, MFIX Documentation: Numerical Technique,      C
!        Technical Note, DOE/MC-31346-5824, NTIS/DE98002029, EG&G      C
!        Technical Services of West Virginia, Inc., Morgantown, WV.    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTIAL_ELIM_U(VAR_G, VAR_S, VXF, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE physprop
      USE indices
      USE run
      USE compar
      USE drag
      USE fldvar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! gas phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_g(DIMENSION_3)
! solids phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_s(DIMENSION_3, DIMENSION_M)
! Volume x gas-solids transfer coefficient
      DOUBLE PRECISION, INTENT(IN) :: VxF(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, I, IJKE
      INTEGER :: L, M, LP, LM
!
      DOUBLE PRECISION :: SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
! a0, b0 etc.
      DOUBLE PRECISION :: a(0:DIMENSION_M), BB(0:DIMENSION_M), &
                          F(0:DIMENSION_M,0:DIMENSION_M),&
                          Saxf(0:DIMENSION_M)
!-----------------------------------------------

!!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,  &
!!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, den) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_E(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

            F = ZERO
            DO M=0, MMAX
               A(M) = A_M(IJK,0,M)
               BB(M)= B_M(IJK,M)

! locally storing vol x avg gas-solids drag coefficient at momentum cell
! center
               IF (M .NE. 0) THEN
                  F(M,0)=-VXF(IJK,M)
                  F(0,M) = F(M,0)
               ENDIF

! locally storing vol x avg solids-solids drag coefficient at momentum cell
! center
               DO L =1, MMAX
                  IF ((L .NE. M) .AND. (M .NE. 0)) THEN
                     LM = FUNLM(L,M)
                     IF (.NOT.IP_AT_E(IJK)) THEN
                        I = I_OF(IJK)
                        IJKE = EAST_OF(IJK)
                        F(M,L) = -AVG_X(F_SS(IJK,LM),F_SS(IJKE,LM),I)*VOL_U(IJK)
                     ENDIF
                  ENDIF
                  F(L,M)=F(M,L)
               ENDDO

               IF (M == 0) THEN
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_G(IPJK)+A_M(IJK,west,M)*VAR_G(IMJK&
                            )+A_M(IJK,north,M)*VAR_G(IJPK)+A_M(IJK,south,M)*VAR_G(IJMK))
               ELSE
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_S(IPJK,M)+A_M(IJK,west,M)*VAR_S(&
                       IMJK,M)+A_M(IJK,north,M)*VAR_S(IJPK,M)+A_M(IJK,south,M)*VAR_S(&
                       IJMK,M))
               ENDIF

               IF (DO_K) THEN
                  IJKM = KM_OF(IJK)
                  IJKP = KP_OF(IJK)
                  IF (M ==0) THEN
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_G(IJKP)+A_M(IJK,bottom,M)*&
                        VAR_G(IJKM))
                  ELSE
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_S(IJKP,M)+A_M(IJK,bottom,M)*&
                        VAR_S(IJKM,M))
                  ENDIF
               ENDIF
           ENDDO   ! end do loop for M over 0 to MMAX

           DO M=0,MMAX
               IF (MOMENTUM_X_EQ(M)) THEN
                  SUM_A = ZERO
                  SUM_B = ZERO

                  DO L =0,MMAX
! - should only need to loop for L/=M. since the resulting term should
!   evaluate to zero anyway, don't bother to add an explicit conditional.
! - if phase L's momentum equation is not solved then do not evaluate the
!   outer sum for phase L (as phase L has no momentum equation). however,
!   the inner sum may include terms from phase L as all phases still see
!   phase L and therefore experience drag with phase L
                     IF (MOMENTUM_X_EQ(L)) THEN
                        SUM_A_LPRIME = ZERO
                        SUM_B_LPRIME = ZERO

                        DO LP=0,MMAX    ! only need to loop for LP/=L; but
                                        ! evaluates to zero so don't bother
                                        ! with an explicit conditional
                           IF ( LP .NE. M) THEN
                              SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
                              IF (LP == 0) THEN
                                 SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
                              ELSE
                                  SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
                              ENDIF
                           ENDIF
                        ENDDO

                        DEN = A(L) + SUM_A_LPRIME + F(L,M)
                        IF ( DEN .NE. ZERO) THEN
                           SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
                           SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
                                 )/DEN
                        ENDIF
                     ELSE
                        SUM_A = SUM_A + F(L,M)
                        IF (L == 0) THEN
                           SUM_B = SUM_B + F(L,M)*VAR_G(IJK)
                        ELSE
                           SUM_B = SUM_B + F(L,M)*VAR_S(IJK,L)
                        ENDIF
                     ENDIF   ! end if momentum_x_eq(L)
                  ENDDO   !end do loop for L over 0 to mmax

                  A_M(IJK,0,M) = SUM_A+A(M)
                  B_M(IJK,M)  =  SUM_B+BB(M)
               ENDIF    ! end if momentum_x_eq(M)
            ENDDO       ! end do loop for M over 0 to MMAX

         ENDIF      ! end if flow_at_e(ijk)
      ENDDO         ! end do loop for ijk over ijkstart3,ijkend3

      RETURN
      END SUBROUTINE PARTIAL_ELIM_U


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_V                                         C
!  Purpose: Do partial elimination for Y vector quantities             C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments:                                                           C
!     Unlike gas-solids exchange coefficient (multiplied by volume),   C
!     the solids-solids exchange coefficient is not as a passed dummy  C
!     argument to this routine but simply re-evaluated here.           C
!                                                                      C
!  Literature/Document References:                                     C
!     Spalding, D.B., 1980, Numerical Computation of Multi-phase       C
!        fluid flow and heat transfer, Recent Advances in Numerical    C
!        Methods in Fluids, C. Taylor et al., eds, Pineridge Press     C
!     Syamlal, M., 1998, MFIX Documentation: Numerical Technique,      C
!        Technical Note, DOE/MC-31346-5824, NTIS/DE98002029, EG&G      C
!        Technical Services of West Virginia, Inc., Morgantown, WV.    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTIAL_ELIM_V(VAR_G, VAR_S, VXF, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE physprop
      USE indices
      USE run
      USE compar
      USE drag
      USE fldvar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! gas phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_g(DIMENSION_3)
! solids phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_s(DIMENSION_3, DIMENSION_M)
! Volume x gas-solids transfer coefficient
      DOUBLE PRECISION, INTENT(IN) :: VxF(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, IJKN, J
      INTEGER ::  L, M, LP, LM
!
      DOUBLE PRECISION :: SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
! a0, b0 etc.
      DOUBLE PRECISION :: a(0:DIMENSION_M), BB(0:DIMENSION_M),&
                          F(0:DIMENSION_M,0:DIMENSION_M),&
                          Saxf(0:DIMENSION_M)
!-----------------------------------------------

!!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,  &
!!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, DEN) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_N(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

            F = ZERO
            DO M = 0, MMAX
               A(M)=A_M(IJK,0,M)
               BB(M)=B_M(IJK,M)

! locally storing vol x avg gas-solids drag coefficient at momentum cell
! center
               IF (M .NE. 0) THEN
                  F(M,0)=-VXF(IJK,M)
                  F(0,M)=F(M,0)
               ENDIF

! locally storing vol x avg solids-solids drag coefficient at momentum cell
! center
               DO L =1, MMAX
                  IF ((L .NE. M) .AND. (M .NE. 0)) THEN
                     LM = FUNLM(L,M)
                     IF (.NOT.IP_AT_N(IJK)) THEN
                        J = J_OF(IJK)
                        IJKN = NORTH_OF(IJK)
                        F(M,L)=-AVG_Y(F_SS(IJK,LM),F_SS(IJKN,LM),J)*VOL_V(IJK)
                     ENDIF
                  ENDIF
                  F(L,M)=F(M,L)
               ENDDO

               IF (M == 0) THEN
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_G(IPJK)+A_M(IJK,west,M)*VAR_G(IMJK&
                        )+A_M(IJK,north,M)*VAR_G(IJPK)+A_M(IJK,south,M)*VAR_G(IJMK))
               ELSE
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_S(IPJK,M)+A_M(IJK,west,M)*VAR_S(&
                       IMJK,M)+A_M(IJK,north,M)*VAR_S(IJPK,M)+A_M(IJK,south,M)*VAR_S(&
                       IJMK,M))
               ENDIF

               IF (DO_K) THEN
                  IJKM = KM_OF(IJK)
                  IJKP = KP_OF(IJK)
                  IF (M == 0) THEN
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_G(IJKP)+A_M(IJK,bottom,M)*&
                        VAR_G(IJKM))
                  ELSE
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_S(IJKP,M)+A_M(IJK,bottom,M)*&
                        VAR_S(IJKM,M))
                  ENDIF
               ENDIF
            ENDDO     ! end do loop for M over 0 to MMAX


            DO M=0,MMAX
               IF (MOMENTUM_Y_EQ(M)) THEN
                  SUM_A = ZERO
                  SUM_B = ZERO

                  DO L =0,MMAX
                     IF(MOMENTUM_Y_EQ(L)) THEN
                        SUM_A_LPRIME = ZERO
                        SUM_B_LPRIME = ZERO

                        DO LP=0,MMAX
                           IF ( LP .NE. M) THEN
                              SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
                              IF (LP == 0) THEN
                                 SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
                              ELSE
                                 SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
                              ENDIF
                           ENDIF
                        ENDDO

                        DEN = A(L) + SUM_A_LPRIME + F(L,M)
                        IF ( DEN .NE. ZERO) THEN
                           SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
                           SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
                                   )/DEN
                        ENDIF
                     ELSE
                        SUM_A = SUM_A + F(L,M)
                        IF (L == 0) THEN
                           SUM_B = SUM_B + F(L,M)*VAR_G(IJK)
                        ELSE
                           SUM_B = SUM_B + F(L,M)*VAR_S(IJK,L)
                        ENDIF
                     ENDIF   ! end if momentum_y_eq(l)
                  ENDDO   !end do loop for L over 0 to mmax

                  A_M(IJK,0,M)=SUM_A+A(M)
                  B_M(IJK,M) = SUM_B+BB(M)
               ENDIF    ! end if momentum_y_eq(M)
            ENDDO       ! end do loop for M over 0 to MMAX

         ENDIF      ! end if flow_at_n(ijk)
      ENDDO         ! end do loop for ijk over ijkstart3,ijkend3

      RETURN
      END SUBROUTINE PARTIAL_ELIM_V


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PARTIAL_ELIM_W                                          C
!  Purpose: Do partial elimination for Z vector quantities
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments:                                                           C
!     Unlike gas-solids exchange coefficient (multiplied by volume),   C
!     the solids-solids exchange coefficient is not as a passed dummy  C
!     argument to this routine but simply re-evaluated here.           C
!                                                                      C
!  Literature/Document References:                                     C
!     Spalding, D.B., 1980, Numerical Computation of Multi-phase       C
!        fluid flow and heat transfer, Recent Advances in Numerical    C
!        Methods in Fluids, C. Taylor et al., eds, Pineridge Press     C
!     Syamlal, M., 1998, MFIX Documentation: Numerical Technique,      C
!        Technical Note, DOE/MC-31346-5824, NTIS/DE98002029, EG&G      C
!        Technical Services of West Virginia, Inc., Morgantown, WV.    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTIAL_ELIM_W(VAR_G, VAR_S, VXF, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE physprop
      USE indices
      USE run
      USE compar
      USE drag
      USE fldvar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! gas phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_g(DIMENSION_3)
! solids phase variable
      DOUBLE PRECISION, INTENT(IN) :: Var_s(DIMENSION_3, DIMENSION_M)
! Volume x gas-solids transfer coefficient
      DOUBLE PRECISION, INTENT(IN) :: VxF(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, K, IJKT
      INTEGER :: L, M, LP, LM
!
      DOUBLE PRECISION :: SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
!                      a0, b0 etc.
      DOUBLE PRECISION :: a(0:DIMENSION_M), BB(0:DIMENSION_M), &
                          F(0:DIMENSION_M,0:DIMENSION_M),&
                          Saxf(0:DIMENSION_M)
!-----------------------------------------------

!!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, &
!!$omp&  a, bb, F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, DEN) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_T(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

            F = ZERO
            DO M=0, MMAX
               A(M)=A_M(IJK,0,M)
               BB(M)=B_M(IJK,M)

! locally storing vol x avg gas-solids drag coefficient at momentum cell
! center
               IF (M .NE. 0) THEN
                  F(M,0)=-VXF(IJK,M)
                  F(0,M)=F(M,0)
               ENDIF

! locally storing vol x avg solids-solids drag coefficient at momentum cell
! center
               DO L =1, MMAX
                  IF ((L .NE. M) .AND. (M .NE. 0)) THEN
                     LM = FUNLM(L,M)
                     IF (.NOT.IP_AT_T(IJK)) THEN
                        K = K_OF(IJK)
                        IJKT = TOP_OF(IJK)
                        F(M,L) = -AVG_Z(F_SS(IJK,LM),F_SS(IJKT,LM),K)*VOL_W(IJK)
                     ENDIF
                  ENDIF
                  F(L,M)=F(M,L)
               ENDDO

               IF (M == 0) THEN
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_G(IPJK)+A_M(IJK,west,M)*VAR_G(IMJK&
                       )+A_M(IJK,north,M)*VAR_G(IJPK)+A_M(IJK,south,M)*VAR_G(IJMK))
               ELSE
                  SAXF(M) = -(A_M(IJK,east,M)*VAR_S(IPJK,M)+A_M(IJK,west,M)*VAR_S(&
                       IMJK,M)+A_M(IJK,north,M)*VAR_S(IJPK,M)+A_M(IJK,south,M)*VAR_S(&
                       IJMK,M))
               ENDIF

               IF (DO_K) THEN
                  IJKM = KM_OF(IJK)
                  IJKP = KP_OF(IJK)
                  IF (M == 0) THEN
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_G(IJKP)+A_M(IJK,bottom,M)*&
                               VAR_G(IJKM))
                  ELSE
                     SAXF(M) = SAXF(M) - (A_M(IJK,top,M)*VAR_S(IJKP,M)+A_M(IJK,bottom,M)*&
                               VAR_S(IJKM,M))
                  ENDIF
               ENDIF
            ENDDO     ! end do loop for M over 0 to MMAX


            DO M=0,MMAX
               IF (MOMENTUM_Z_EQ(M)) THEN
                  SUM_A = ZERO
                  SUM_B = ZERO

                  DO L =0,MMAX
                     IF (MOMENTUM_Z_EQ(L)) THEN
                        SUM_A_LPRIME = ZERO
                        SUM_B_LPRIME = ZERO

                        DO LP=0,MMAX
                           IF ( LP .NE. M) THEN
                              SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
                              IF (LP == 0) THEN
                                 SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
                              ELSE
                                 SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
                              ENDIF
                           ENDIF
                        ENDDO

                        DEN = A(L) + SUM_A_LPRIME + F(L,M)
                        IF ( DEN .NE. ZERO) THEN
                           SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
                           SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
                           )/DEN
                        ENDIF
                     ELSE
                        SUM_A = SUM_A + F(L,M)
                        IF (L == 0) THEN
                           SUM_B = SUM_B + F(L,M)*VAR_G(IJK)
                        ELSE
                           SUM_B = SUM_B + F(L,M)*VAR_S(IJK,L)
                        ENDIF
                     ENDIF   ! end if momentum_z_eq(L)
                  ENDDO   ! end do loop for L over 0 to mmax

                  A_M(IJK,0,M)=SUM_A+A(M)
                  B_M(IJK,M) = SUM_B+BB(M)
               ENDIF    ! end if momentum_Z_eq(M)
            ENDDO       ! end do loop for M over 0 to MMAX

         ENDIF      ! end if flow_at_t(ijk)
      ENDDO         ! end do loop for ijk over ijkstart3,ijkend3

      RETURN
      END SUBROUTINE PARTIAL_ELIM_W


