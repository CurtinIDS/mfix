!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_VORTICITY                                         C
!  Purpose: Computes the vorticity                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CALC_VORTICITY

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE fldvar
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      INTEGER :: I,J,K,IJK,IP,IM,JP,JM,KP,KM,IJKE,IJKN,IJKT,IJKW,IJKS,IJKB
      DOUBLE PRECISION :: DU_DX,DU_DY,DU_DZ,DV_DX,DV_DY,DV_DZ,DW_DX,DW_DY,DW_DZ
      DOUBLE PRECISION :: OMEGA_X,OMEGA_Y,OMEGA_Z
      DOUBLE PRECISION :: LAMBDA_1,LAMBDA_2,LAMBDA_3
      DOUBLE PRECISION,DIMENSION(3,3) :: OMEGA,SS,AA
      DOUBLE PRECISION,DIMENSION(4) :: POLY

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IM = I - 1
         IP = I + 1

         JM = J - 1
         JP = J + 1

         KM = K - 1
         KP = K + 1


         IJKE = FUNIJK(IP,J,K)
         IJKN = FUNIJK(I,JP,K)
         IJKT = FUNIJK(I,J,KP)

         IJKW = FUNIJK(IM,J,K)
         IJKS = FUNIJK(I,JM,K)
         IJKB = FUNIJK(I,J,KM)



         IF (FLUID_AT(IJK).AND.INTERIOR_CELL_AT(IJK)) THEN

            DU_DX = ZERO
            DU_DY = ZERO
            DU_DZ = ZERO

            DV_DX = ZERO
            DV_DY = ZERO
            DV_DZ = ZERO

            DW_DX = ZERO
            DW_DY = ZERO
            DW_DZ = ZERO

!======================================================================
!  Get Velocity derivatives
!======================================================================

!======================================================================
!  du/dx
!======================================================================

            IF ((FLUID_AT(IJKE)).AND.(FLUID_AT(IJKW)))  THEN

               DU_DX = (U_g(IJKE) - U_g(IJKW)) / (X_U(IJKE) - X_U(IJKW))

            ELSE IF ((FLUID_AT(IJKE)).AND.(.NOT.FLUID_AT(IJKW)))  THEN

               DU_DX = (U_g(IJKE) - U_g(IJK)) / (X_U(IJKE) - X_U(IJK))

            ELSE IF ((FLUID_AT(IJKW)).AND.(.NOT.FLUID_AT(IJKE)))  THEN

               DU_DX = (U_g(IJK) - U_g(IJKW)) / (X_U(IJK) - X_U(IJKW))

            END IF

!======================================================================
!  du/dy
!======================================================================

            IF ((FLUID_AT(IJKN)).AND.(FLUID_AT(IJKS)))  THEN

               DU_DY = (U_g(IJKN) - U_g(IJKS)) / (Y_U(IJKN) - Y_U(IJKS))

            ELSE IF ((FLUID_AT(IJKN)).AND.(.NOT.FLUID_AT(IJKS)))  THEN

               DU_DY = (U_g(IJKN) - U_g(IJK)) / (Y_U(IJKN) - Y_U(IJK))

            ELSE IF ((FLUID_AT(IJKS)).AND.(.NOT.FLUID_AT(IJKN)))  THEN

               DU_DY = (U_g(IJK) - U_g(IJKS)) / (Y_U(IJK) - Y_U(IJKS))

            END IF

!======================================================================
!  du/dz
!======================================================================

            IF(DO_K) THEN
               IF ((FLUID_AT(IJKT)).AND.(FLUID_AT(IJKB)))  THEN

                  DU_DZ = (U_g(IJKT) - U_g(IJKB)) / (Z_U(IJKT) - Z_U(IJKB))

               ELSE IF ((FLUID_AT(IJKT)).AND.(.NOT.FLUID_AT(IJKB)))  THEN

                  DU_DZ = (U_g(IJKT) - U_g(IJK)) / (Z_U(IJKT) - Z_U(IJK))

               ELSE IF ((FLUID_AT(IJKB)).AND.(.NOT.FLUID_AT(IJKT)))  THEN

                  DU_DZ = (U_g(IJK) - U_g(IJKB)) / (Z_U(IJK) - Z_U(IJKB))

               END IF
            ENDIF

!======================================================================
!  dv/dx
!======================================================================

            IF ((FLUID_AT(IJKE)).AND.(FLUID_AT(IJKW)))  THEN

               DV_DX = (V_g(IJKE) - V_g(IJKW)) / (X_V(IJKE) - X_V(IJKW))

            ELSE IF ((FLUID_AT(IJKE)).AND.(.NOT.FLUID_AT(IJKW)))  THEN

               DV_DX = (V_g(IJKE) - V_g(IJK)) / (X_V(IJKE) - X_V(IJK))

            ELSE IF ((FLUID_AT(IJKW)).AND.(.NOT.FLUID_AT(IJKE)))  THEN

               DV_DX = (V_g(IJK) - V_g(IJKW)) / (X_V(IJK) - X_V(IJKW))

            END IF

!======================================================================
!  dv/dy
!======================================================================

            IF ((FLUID_AT(IJKN)).AND.(FLUID_AT(IJKS)))  THEN

               DV_DY = (V_g(IJKN) - V_g(IJKS)) / (Y_V(IJKN) - Y_V(IJKS))

            ELSE IF ((FLUID_AT(IJKN)).AND.(.NOT.FLUID_AT(IJKS)))  THEN

               DV_DY = (V_g(IJKN) - V_g(IJK)) / (Y_V(IJKN) - Y_V(IJK))

            ELSE IF ((FLUID_AT(IJKS)).AND.(.NOT.FLUID_AT(IJKN)))  THEN

               DV_DY = (V_g(IJK) - V_g(IJKS)) / (Y_V(IJK) - Y_V(IJKS))

            END IF

!======================================================================
!  dv/dz
!======================================================================

            IF(DO_K) THEN
               IF ((FLUID_AT(IJKT)).AND.(FLUID_AT(IJKB)))  THEN

                  DV_DZ = (V_g(IJKT) - V_g(IJKB)) / (Z_V(IJKT) - Z_V(IJKB))

               ELSE IF ((FLUID_AT(IJKT)).AND.(.NOT.FLUID_AT(IJKB)))  THEN

                  DV_DZ = (V_g(IJKT) - V_g(IJK)) / (Z_V(IJKT) - Z_V(IJK))

               ELSE IF ((FLUID_AT(IJKB)).AND.(.NOT.FLUID_AT(IJKT)))  THEN

                  DV_DZ = (V_g(IJK) - V_g(IJKB)) / (Z_V(IJK) - Z_V(IJKB))

               END IF
            ENDIF

!======================================================================
!  dw/dx
!======================================================================

            IF(DO_K) THEN
               IF ((FLUID_AT(IJKE)).AND.(FLUID_AT(IJKW)))  THEN

                  DW_DX = (W_g(IJKE) - W_g(IJKW)) / (X_W(IJKE) - X_W(IJKW))

               ELSE IF ((FLUID_AT(IJKE)).AND.(.NOT.FLUID_AT(IJKW)))  THEN

                  DW_DX = (W_g(IJKE) - W_g(IJK)) / (X_W(IJKE) - X_W(IJK))

               ELSE IF ((FLUID_AT(IJKW)).AND.(.NOT.FLUID_AT(IJKE)))  THEN

                  DW_DX = (W_g(IJK) - W_g(IJKW)) / (X_W(IJK) - X_W(IJKW))

               END IF

!======================================================================
!  dw/dy
!======================================================================

               IF ((FLUID_AT(IJKN)).AND.(FLUID_AT(IJKS)))  THEN

                  DW_DY = (W_g(IJKN) - W_g(IJKS)) / (Y_W(IJKN) - Y_W(IJKS))

               ELSE IF ((FLUID_AT(IJKN)).AND.(.NOT.FLUID_AT(IJKS)))  THEN

                  DW_DY = (W_g(IJKN) - W_g(IJK)) / (Y_W(IJKN) - Y_W(IJK))

               ELSE IF ((FLUID_AT(IJKS)).AND.(.NOT.FLUID_AT(IJKN)))  THEN

                  DW_DY = (W_g(IJK) - W_g(IJKS)) / (Y_W(IJK) - Y_W(IJKS))

               END IF

!======================================================================
!  dw/dz
!======================================================================

               IF ((FLUID_AT(IJKT)).AND.(FLUID_AT(IJKB)))  THEN

                  DW_DZ = (W_g(IJKT) - W_g(IJKB)) / (Z_W(IJKT) - Z_W(IJKB))

               ELSE IF ((FLUID_AT(IJKT)).AND.(.NOT.FLUID_AT(IJKB)))  THEN

                  DW_DZ = (W_g(IJKT) - W_g(IJK)) / (Z_W(IJKT) - Z_W(IJK))

               ELSE IF ((FLUID_AT(IJKB)).AND.(.NOT.FLUID_AT(IJKT)))  THEN

                  DW_DZ = (W_g(IJK) - W_g(IJKB)) / (Z_W(IJK) - Z_W(IJKB))

               END IF
            ENDIF

!======================================================================
!  Build Vorticity components and magnitude
!======================================================================

            OMEGA_X = DW_DY - DV_DZ
            OMEGA_Y = DU_DZ - DW_DX
            OMEGA_Z = DV_DX - DU_DY


            VORTICITY(IJK) = DSQRT(OMEGA_X**2 + OMEGA_Y**2 + OMEGA_Z**2)

!======================================================================
!  Build Matrices OMEGA , SS, and AA
!
!  OMEGA_ij = 1/2 *(dui/dxj + duj/dxi)
!  SS_ij    = 1/2 *(dui/dxj - duj/dxi)
!  AA       = OMEGA^2 + SS^2
!======================================================================

            OMEGA(1,1) = DU_DX
            OMEGA(1,2) = HALF * (DU_DY + DV_DX)
            OMEGA(1,3) = HALF * (DU_DZ + DW_DX)

            OMEGA(2,1) = OMEGA(1,2)
            OMEGA(2,2) = DV_DY
            OMEGA(2,3) = HALF * (DV_DZ + DW_DY)

            OMEGA(3,1) = OMEGA(1,3)
            OMEGA(3,2) = OMEGA(2,3)
            OMEGA(3,3) = DW_DZ


            SS(1,1) = ZERO
            SS(1,2) = HALF * (DU_DY - DV_DX)
            SS(1,3) = HALF * (DU_DZ - DW_DX)

            SS(2,1) = - SS(1,2)
            SS(2,2) = ZERO
            SS(2,3) = HALF * (DV_DZ - DW_DY)

            SS(3,1) = - SS(1,3)
            SS(3,2) = - SS(2,3)
            SS(3,3) = ZERO


            AA = MATMUL(OMEGA,OMEGA) + MATMUL(SS,SS)

!======================================================================
!  Build Characteristic polynomial of AA
!======================================================================

            POLY(1) = -  ONE
            POLY(2) =    AA(1,1) + AA(2,2) + AA(3,3)
            POLY(3) =    AA(2,1)*AA(1,2) + AA(3,1)*AA(1,3) + AA(3,2)*AA(2,3) &
                       -( AA(1,1)*AA(2,2) + AA(1,1)*AA(3,3) + AA(2,2)*AA(3,3) )
            POLY(4) =    AA(1,1)*AA(2,2)*AA(3,3) + AA(1,2)*AA(2,3)*AA(3,1) + AA(2,1)*AA(3,2)*AA(1,3) &
                       -( AA(3,1)*AA(2,2)*AA(1,3) + AA(1,2)*AA(2,1)*AA(3,3) + AA(2,3)*AA(3,2)*AA(1,1) )



!======================================================================
!  Find roots of characteristic polynomial = eigenvalues
!  using Bairstow method
!======================================================================

            CALL BAIRSTOW(POLY,LAMBDA_1,LAMBDA_2,LAMBDA_3)

            LAMBDA2(IJK) = LAMBDA_2

         ENDIF

      END DO

      RETURN


      END SUBROUTINE CALC_VORTICITY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BAIRSTOW                                               C
!  Purpose: FIND THE ROOTS OF A POLYNOMIAL USING BAIRSTOW METHOD       C
!           LIMITED TO POLYNOMIALS OF DEGREE 3                         C
!           (USED TO FIND THE EIGENVALUES OF A 3x3 MATRIX)             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 17-Jul-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE BAIRSTOW(A,X1,X2,X3)
      USE param1

      IMPLICIT NONE
      INTEGER :: N,NMAX,ITERMAX,J,ITER

      PARAMETER(NMAX = 4)  ! MAXIMUM NUMBER OF POLYNOMIAL COEFFICIENTS ALLOWED

      DOUBLE PRECISION, DIMENSION(NMAX)::  A,B,C
      DOUBLE PRECISION::  R,S,DELTA,DELTAR,DELTAS,DENOM
      DOUBLE PRECISION:: EPS
      DOUBLE PRECISION:: X1,X2,X3
      DOUBLE PRECISION:: BUFFER1

      LOGICAL :: TEST1, TEST2, TEST3

      N = 3            ! This subroutine is designed for a ploynomial of degree three only
      ITERMAX = 200    ! It is not anticipated that ITERMAX nor EPS would need to
      EPS = 1.0D-6     ! be modified

!=======================================================================
!     POLYNOMIAL REDUCTION
!     ALWAYS START WITH INITIAL GUESS R = - 0.9 AND S = -1.0
!     AND ITERATE
!=======================================================================

      R = -0.9D0
      S = -1.0D0

      ITER = 0

!=======================================================================
!     COMPUTING B'S ANS C'S
!=======================================================================

20       CONTINUE

      B(1) = A(1)
      B(2) = A(2) + R * B(1)
      C(1) = B(1)
      C(2) = B(2) + R * C(1)
      B(3) = A(3) + R * B(2) + S * B(1)
      C(3) = B(3) + R * C(2) + S * C(1)
      B(4) = A(4) + R * B(3) + S * B(2)

!=======================================================================
!     UPDATING VALUES OF R AND S
!=======================================================================

      DENOM = C(N-1) * C(N-1) - C(N)*C(N-2)

      DELTAR = (-B(N)*C(N-1) + B(N+1)*C(N-2) ) / DENOM
      DELTAS = (-C(N-1)*B(N+1) + B(N)*C(N) ) / DENOM

      R = R + DELTAR
      S = S + DELTAS

      ITER = ITER +1

!=======================================================================
!     CHECKING CONVERGENCE
!=======================================================================

      TEST1 = (ABS( B(N) ) > EPS)
      TEST2 = (ABS( B(N+1)) > EPS)
      TEST3 = (ITER <= ITERMAX)

      IF(TEST1.AND.TEST2.AND.TEST3) GOTO 20

      IF (.NOT.TEST3) THEN
!         WRITE(*,*)'ERROR IN BAIRSTOW SUBROUTINE:'
!         WRITE(*,*)'DIVERGENCE... THE PROGRAM WILL BE TERMINATED NOW.'
!         CALL MFIX_EXIT(myPE)
         X1 = UNDEFINED
         X2 = UNDEFINED
         X3 = UNDEFINED
      ENDIF

!=======================================================================
!     NOW, THE ITERATIVE SCHEME HAS CONVERGED
!     THE QUADRATIC FACTOR IS X*X-R*X-S
!     SOLVING AND DISPLAYING THE TWO ROOTS (LIMITED TO REAL ROOTS)
!=======================================================================

      DELTA = R*R + 4.0D0 * S

      IF (DELTA.LT.0.0D0) THEN

!         WRITE(*,*) ' ERROR IN BAIRSTOW SUBROUTINE:TWO COMPLEX ROOTS:'

         X1 = UNDEFINED
         X2 = UNDEFINED

      ELSE

         X1 = ( R + DSQRT(DELTA)) / ( 2.0D0 )
         X2 = ( R - DSQRT(DELTA)) / ( 2.0D0 )

      ENDIF

!=======================================================================
!     IN PREPARATION FOR THE LAST LINEAR FACTOR,
!     SET N = N-2
!     RENAME THE B'S AS A'S
!=======================================================================

      N = N - 2

      DO J=1,N+1
         A(J) = B(J)
      END DO

      IF(DABS(A(1))<1.0D-9) THEN
         X3 = UNDEFINED
      ELSE
         X3 = -A(2)/A(1)
      ENDIF

!=======================================================================
!     SORING OUT ROOTS IN ASCENDING ORDER
!=======================================================================

      IF(X1>X2) THEN
         BUFFER1 = X1
         X1 = X2
         X2 = BUFFER1
      ENDIF

      IF(X2>X3) THEN
         BUFFER1 = X2
         X2 = X3
         X3 = BUFFER1
      ENDIF

      IF(X1>X2) THEN
         BUFFER1 = X1
         X1 = X2
         X2 = BUFFER1
      ENDIF

      RETURN

      END SUBROUTINE BAIRSTOW


