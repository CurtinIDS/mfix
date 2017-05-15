!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_GMRES                                               C
!  Purpose: Solve system of linear system using GMRES method           C
!           generalized minimal residual                               C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_GMRES(VNAME, VNO, VAR, A_M, B_M, &
                    cmethod, TOL, ITMAX, MAX_IT, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE leqsol, only: leq_msolve, leq_matvec
      USE mpi_utility, only: global_all_and
      USE param, only: dimension_3
      USE param1, only: zero

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m
! Vector b_m
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! maximum number of outer iterations
! (currently set to 1 by calling subroutine-solve_lin_eq)
      INTEGER, INTENT(IN) :: MAX_IT
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-------------------------------------------------
! Local Variables
!-------------------------------------------------
      LOGICAL :: IS_BM_ZERO, ALL_IS_BM_ZERO
!-------------------------------------------------

      IS_BM_ZERO = (MAXVAL( ABS(B_M(:)) ) .EQ. ZERO)
      CALL GLOBAL_ALL_AND( IS_BM_ZERO, ALL_IS_BM_ZERO )

      IF (ALL_IS_BM_ZERO) THEN
          VAR(:) = ZERO
          RETURN
      ENDIF

      CALL LEQ_GMRES0( VNAME, VNO, VAR, A_M, B_M,  &
           CMETHOD, TOL, ITMAX, MAX_IT, LEQ_MATVEC, LEQ_MSOLVE, IER )

      RETURN
      END SUBROUTINE LEQ_GMRES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_GMRES0                                              C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_GMRES0(VNAME, VNO, VAR, A_M, B_M,  &
                            CMETHOD, TOL, ITMAX, MAX_IT, &
                            MATVEC, MSOLVE, IER )

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE debug, only: idebug, write_debug
      USE functions, only: funijk
      USE gridmap, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE leqsol, only: dot_product_par
      USE mpi_utility, only: global_all_and, global_all_or, global_all_max, global_all_min
      USE param, only: dimension_3
      USE param1, only: one, zero

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments/procedure
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here-see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3,-3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! maximum number of outer iterations
      INTEGER, INTENT(IN) :: MAX_IT
! error indicator
      INTEGER, INTENT(INOUT) :: IER
! dummy arguments/procedures set as indicated
!     matvec->leq_matvec
! for preconditioner (leq_pc)
!     msolve->leq_msolve
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER JDEBUG
      PARAMETER(JDEBUG=0)
!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION, dimension(:,:), allocatable :: V
      DOUBLE PRECISION H(ITMAX+1,ITMAX)
      DOUBLE PRECISION CS(ITMAX)
      DOUBLE PRECISION SN(ITMAX)
      DOUBLE PRECISION Y(ITMAX)

      DOUBLE PRECISION, dimension(:), allocatable :: R,TEMP,WW

      DOUBLE PRECISION E1(ITMAX+2)
      DOUBLE PRECISION SS(ITMAX+2)

      DOUBLE PRECISION BNRM2, ERROR, ERROR0
      DOUBLE PRECISION NORM_R0, NORM_R, NORM_W
      DOUBLE PRECISION INV_NORM_R, NORM_S, NORM_Y
      DOUBLE PRECISION YII

      INTEGER IJK, II, JJ, KK, I, K
      INTEGER RESTRT, M, ITER, MDIM
      DOUBLE PRECISION DTEMP
      DOUBLE PRECISION MINH, MAXH, CONDH, ALL_MINH, ALL_MAXH
      DOUBLE PRECISION INV_H_IP1_I

      LOGICAL IS_BM_ZERO, IS_ERROR, IS_CONVERGED
      LOGICAL ALL_IS_BM_ZERO, ALL_IS_ERROR, ALL_IS_CONVERGED

      CHARACTER(LEN=40) :: NAME
!-----------------------------------------------


      ! allocating
      allocate(V(DIMENSION_3,ITMAX+1))
      allocate(R(DIMENSION_3))
      allocate(TEMP(DIMENSION_3))
      allocate(WW(DIMENSION_3))

! initializing
      NAME = 'LEQ_GMRES0 ' // TRIM(VNAME)
      RESTRT = ITMAX
      M = RESTRT
      ITER = 0
      IER = 0


! clear arrays
! --------------------------------
      R(:) = ZERO
      TEMP(:) = ZERO

      H(:,:) = ZERO
      CS(:) = ZERO
      E1(:) = ZERO
      E1(1) = ONE
      SS(:) = ZERO
      SN(:) = ZERO
      WW(:) = ZERO
      V(:,:) = ZERO


      BNRM2 = DOT_PRODUCT_PAR( B_M, B_M )
      BNRM2 = sqrt( BNRM2 )

      if (idebug.ge.1) then
         call write_debug(name, 'bnrm2 = ', bnrm2 )
      endif

      IS_BM_ZERO = BNRM2 .EQ. ZERO
      CALL GLOBAL_ALL_AND( IS_BM_ZERO, ALL_IS_BM_ZERO )
      IF (ALL_IS_BM_ZERO) THEN
        BNRM2 = ONE
      ENDIF

! r = M \ (b - A*x)
! error = norm(r) / bnrm2
! --------------------------------

      CALL MATVEC(VNAME, VAR, A_M, R)   ! returns R=A*VAR

!!$omp   parallel do private(ii,jj,kk,ijk)
      DO KK=KSTART3,KEND3
        DO JJ=JSTART3,JEND3
          DO II=ISTART3,IEND3
             IJK = FUNIJK( II,JJ,KK )
             TEMP(IJK) = B_M(IJK) - R(IJK)
          ENDDO
        ENDDO
      ENDDO

! Solve A*R(:) = TEMP(:)
      CALL MSOLVE(VNAME, TEMP, A_M,  R, CMETHOD)   ! returns R

      NORM_R = DOT_PRODUCT_PAR( R, R )
      NORM_R = SQRT( NORM_R )

      ERROR = NORM_R/BNRM2
      ERROR0 = ERROR
      NORM_R0 = NORM_R

      IF (JDEBUG.GE.1) THEN
         CALL WRITE_DEBUG( NAME,  ' INITIAL ERROR ', ERROR0 )
         CALL WRITE_DEBUG( NAME,  ' INITIAL RESIDUAL ', NORM_R0)
      ENDIF


! begin iteration
      DO ITER=1,MAX_IT
! ---------------------------------------------------------------->>>

! r = M \ (b-A*x)
! --------------------------------
         CALL MATVEC( VNAME, VAR, A_M, R )   ! Returns R=A*VAR
!        TEMP(:) = B_M(:) - R(:)

!!$omp    parallel do private(ii,jj,kk,ijk)
         DO KK=KSTART3,KEND3
           DO JJ=JSTART3,JEND3
             DO II=ISTART3,IEND3
                IJK = FUNIJK( II,JJ,KK )
                TEMP(IJK) = B_M(IJK) - R(IJK)
             ENDDO
           ENDDO
         ENDDO

! Solve A*R(:) = TEMP(:)
         CALL MSOLVE(VNAME, TEMP, A_M, R, CMETHOD)   ! returns R

         NORM_R =  DOT_PRODUCT_PAR( R, R )
         NORM_R = SQRT( NORM_R )
         INV_NORM_R = ONE / NORM_R

!         V(:,1) = R(:) * INV_NORM_R
!!$omp    parallel do private(ii,jj,kk,ijk)
         DO KK=KSTART3,KEND3
           DO JJ=JSTART3,JEND3
             DO II=ISTART3,IEND3
                IJK = FUNIJK( II,JJ,KK )
                V(IJK,1) = R(IJK) * INV_NORM_R
             ENDDO
           ENDDO
         ENDDO

         SS(:) = NORM_R * E1(:)


! construct orthonormal basis using Gram-Schmidt
! ---------------------------------------------------------------->>>
         DO I=1,M    ! M->restrt->itmax
! w = M \ (A*V(:,i))
! --------------------------------
            CALL MATVEC(VNAME, V(:,I), A_M, TEMP)   ! returns TEMP=A*V
! Solve A*WW(:) = TEMP(:)
            CALL MSOLVE(VNAME, TEMP, A_M, WW, CMETHOD)   ! returns WW

            DO K=1,I
               DTEMP = DOT_PRODUCT_PAR( WW(:), V(:,K) )
               H(K,I) = DTEMP
!               WW(:) = WW(:) - H(K,I)*V(:,K)

!!$omp          parallel do private(ii,jj,kk,ijk)
               DO KK=KSTART3,KEND3
                 DO JJ=JSTART3,JEND3
                   DO II=ISTART3,IEND3
                      IJK = FUNIJK( II,JJ,KK )
                      WW(IJK) = WW(IJK) - H(K,I)*V(IJK,K)
                   ENDDO
                 ENDDO
               ENDDO
            ENDDO

            NORM_W =  DOT_PRODUCT_PAR( WW(:), WW(:) )
            NORM_W = SQRT( NORM_W )
            H(I+1,I) = NORM_W
!            V(:,I+1) = WW(:) / H(I+1,I)
            INV_H_IP1_I = ONE / H(I+1,I)

!!$omp       parallel do private(ii,jj,kk,ijk)
            DO KK=KSTART3,KEND3
              DO JJ=JSTART3,JEND3
                DO II=ISTART3,IEND3
                   IJK = FUNIJK( II,JJ,KK )
                   V(IJK, I+1) = WW(IJK) * INV_H_IP1_I
                ENDDO
              ENDDO
            ENDDO

! apply Givens rotation
! --------------------------------
            DO K=1,I-1
               DTEMP    =  CS(K)*H(K,I) + SN(K)*H(K+1,I)
               H(K+1,I) = -SN(K)*H(K,I) + CS(K)*H(K+1,I)
               H(K,I)   = DTEMP
            ENDDO

! form i-th rotation matrix approximate residual norm
! --------------------------------
            CALL ROTMAT( H(I,I), H(I+1,I), CS(I), SN(I) )

            DTEMP = CS(I)*SS(I)
            SS(I+1) = -SN(I)*SS(I)
            SS(I) = DTEMP
            H(I,I) = CS(I)*H(I,I) + SN(I)*H(I+1,I)
            H(I+1,I) = ZERO
            ERROR = ABS( SS(I+1) ) / BNRM2

            IS_CONVERGED = (ERROR .LE. TOL*ERROR0)
            CALL GLOBAL_ALL_AND( IS_CONVERGED, ALL_IS_CONVERGED )
            IF (ALL_IS_CONVERGED) THEN
! update approximation and exit
! --------------------------------

! triangular solve with y = H(1:i,1:i) \ s(1:i)
! --------------------------------
               MDIM = I
               DO II=1,MDIM
                 Y(II) = SS(II)
               ENDDO

               DO II=MDIM,1,-1
                  YII = Y(II)/H(II,II)
                  Y(II) = YII
                  DO JJ=1,II-1
                     Y(JJ) = Y(JJ) - H(JJ,II)*YII
                  ENDDO
               ENDDO

! double check
! --------------------------------
               IF (JDEBUG.GE.1) THEN
                  MAXH = ABS(H(1,1))
                  MINH = ABS(H(1,1))
                  DO II=1,MDIM
                     MAXH = MAX( MAXH, ABS(H(II,II)) )
                     MINH = MIN( MINH, ABS(H(II,II)) )
                  ENDDO

                  CALL GLOBAL_ALL_MAX( MAXH, ALL_MAXH )
                  CALL GLOBAL_ALL_MIN( MINH, ALL_MINH )
                  MAXH = ALL_MAXH
                  MINH = ALL_MINH
                  CONDH = MAXH/MINH

                  TEMP(1:MDIM) = SS(1:MDIM) - &
                  MATMUL( H(1:MDIM,1:MDIM ), Y(1:MDIM) )
                  DTEMP = DOT_PRODUCT( TEMP(1:MDIM), TEMP(1:MDIM) )
                  DTEMP = SQRT( DTEMP )

                  NORM_S = DOT_PRODUCT( SS(1:MDIM),SS(1:MDIM) )
                  NORM_S = SQRT( NORM_S )

                  NORM_Y = DOT_PRODUCT( Y(1:MDIM),Y(1:MDIM) )
                  NORM_Y = SQRT( NORM_Y )

                  IS_ERROR = (DTEMP .GT. CONDH*NORM_S)
                  CALL GLOBAL_ALL_OR( IS_ERROR, ALL_IS_ERROR )
                  IF (ALL_IS_ERROR) THEN
                     CALL WRITE_DEBUG(NAME, &
                        'DTEMP, NORM_S ', DTEMP, NORM_S )
                     CALL WRITE_DEBUG(NAME, &
                        'CONDH, NORM_Y ', CONDH, NORM_Y )
                     CALL WRITE_DEBUG(NAME, &
                        '** STOP IN LEQ_GMRES ** ')
                     IER = 999
                     RETURN
                  ENDIF
               ENDIF   ! end if(jdebug>=1)


!               VAR(:)=VAR(:)+MATMUL(V(:,1:I),Y(1:I))
!!$omp          parallel do private(ii,jj,kk,ijk)
               DO KK=KSTART3,KEND3
                 DO JJ=JSTART3,JEND3
                   DO II=ISTART3,IEND3
                      IJK = FUNIJK( II,JJ,KK )
                      VAR(IJK) = VAR(IJK) + &
                         DOT_PRODUCT( V(IJK,1:I), Y(1:I) )
                   ENDDO
                 ENDDO
               ENDDO

               EXIT
            ENDIF   !end if(all_is_converged)
         ENDDO   !end do I=1,m (construct orthonormal basis using Gram-Schmidt)
! ----------------------------------------------------------------<<<

         IS_CONVERGED = ( ERROR .LE. TOL*ERROR0)
         CALL GLOBAL_ALL_AND( IS_CONVERGED, ALL_IS_CONVERGED )
         IF ( ALL_IS_CONVERGED ) THEN
            EXIT
         ENDIF

! update approximations
! --------------------------------

! y = H(1:m,1:m) \ s(1:m)
! x = x + V(:,1:m)*y
! r = M \ (b-A*x)
! --------------------------------
         MDIM = M
         DO II=1,MDIM
            Y(II) = SS(II)
         ENDDO

         DO II=MDIM,1,-1
            YII = Y(II)/H(II,II)
            Y(II) = YII
            DO JJ=1,II-1
               Y(JJ) = Y(JJ) - H(JJ,II)*YII
            ENDDO
         ENDDO

! double check
! --------------------------------
         IF (JDEBUG.GE.1) THEN
            MAXH = ABS(H(1,1))
            MINH = ABS(H(1,1))
            DO II=1,MDIM
              MAXH = MAX( MAXH, ABS(H(II,II)) )
              MINH = MIN( MINH, ABS(H(II,II)) )
            ENDDO
            CONDH = MAXH/MINH

            TEMP(1:MDIM) = SS(1:MDIM) - &
               MATMUL( H(1:MDIM,1:MDIM ), Y(1:MDIM) )

            DTEMP = DOT_PRODUCT( TEMP(1:MDIM), TEMP(1:MDIM) )
            DTEMP = SQRT( DTEMP )

            NORM_S = DOT_PRODUCT( SS(1:MDIM),SS(1:MDIM) )
            NORM_S = SQRT( NORM_S )

            NORM_Y = DOT_PRODUCT( Y(1:MDIM),Y(1:MDIM) )
            NORM_Y = SQRT( NORM_Y )

            IS_ERROR = (DTEMP .GT. CONDH*NORM_S)
            CALL GLOBAL_ALL_OR( IS_ERROR, ALL_IS_ERROR )
            IF (ALL_IS_ERROR) THEN
               CALL WRITE_DEBUG(NAME, &
                  'DTEMP, NORM_S ', DTEMP, NORM_S )
               CALL WRITE_DEBUG(NAME, &
                  'CONDH, NORM_Y ', CONDH, NORM_Y )
               CALL WRITE_DEBUG(NAME, &
                  '** STOP IN LEQ_GMRES ** ')
               IER = 999
               RETURN
            ENDIF
         ENDIF   ! end if(jdebug>=1)

!         VAR(:) = VAR(:) + MATMUL( V(:,1:M), Y(1:M) )
!!$omp    parallel do private(ii,jj,kk,ijk)
         DO KK=KSTART3,KEND3
           DO JJ=JSTART3,JEND3
             DO II=ISTART3,IEND3
                IJK = FUNIJK( II,JJ,KK )
                VAR(IJK) = VAR(IJK) + &
                   DOT_PRODUCT( V(IJK,1:M), Y(1:M) )
             ENDDO
           ENDDO
         ENDDO

         CALL MATVEC(VNAME, VAR, A_M, R)   ! returns R=A*VAR

!         TEMP(:) = B_M(:) - R(:)
!!$omp    parallel do private(ii,jj,kk,ijk)
         DO KK=KSTART3,KEND3
           DO JJ=JSTART3,JEND3
             DO II=ISTART3,IEND3
                IJK = FUNIJK( II,JJ,KK )
                TEMP(IJK) = B_M(IJK) - R(IJK)
             ENDDO
           ENDDO
         ENDDO

! Solve A*R(:)=TEMP(:)
         CALL MSOLVE( VNAME, TEMP, A_M, R, CMETHOD)   ! Returns R
         NORM_R =  DOT_PRODUCT_PAR( R, R )
         NORM_R = SQRT( NORM_R )

         SS(I+1) = NORM_R
         ERROR = SS(I+1) / BNRM2

         IF (JDEBUG.GE.1) THEN
            CALL WRITE_DEBUG(NAME, 'LEQ_GMRES: I, ITER ', I, ITER)
            CALL WRITE_DEBUG(NAME, 'LEQ_GMRES:  ERROR, NORM_R ',  &
               ERROR, NORM_R )
         ENDIF

         IS_CONVERGED = (ERROR .LE. TOL*ERROR0)
         CALL GLOBAL_ALL_AND( IS_CONVERGED, ALL_IS_CONVERGED )
         IF (ALL_IS_CONVERGED)  THEN
            EXIT
         ENDIF

      ENDDO   ! end do iter=1,max_it
! ----------------------------------------------------------------<<<

      IS_ERROR = (ERROR .GT. TOL*ERROR0)
      CALL GLOBAL_ALL_OR( IS_ERROR, ALL_IS_ERROR )
      IF (ALL_IS_ERROR) THEN
         IER = 1
      ENDIF

      IF (JDEBUG.GE.1) THEN
         CALL MATVEC(VNAME, VAR, A_M, R)   ! Returns R=A*VAR

!         R(:) = R(:) - B_M(:)
!!$omp    parallel do private(ii,jj,kk,ijk)
         DO KK=KSTART3,KEND3
           DO JJ=JSTART3,JEND3
             DO II=ISTART3,IEND3
                IJK = FUNIJK( II,JJ,KK )
                R(IJK) = R(IJK) - B_M(IJK)
             ENDDO
           ENDDO
         ENDDO

         NORM_R=DOT_PRODUCT_PAR(R,R)
         NORM_R = SQRT( NORM_R )
         CALL WRITE_DEBUG(NAME, 'ITER, I ', ITER, I )
         CALL WRITE_DEBUG(NAME, 'VNAME = ' // VNAME )
         CALL WRITE_DEBUG(NAME, 'IER  = ', IER )
         CALL WRITE_DEBUG(NAME, 'ERROR0  = ', ERROR0 )
         CALL WRITE_DEBUG(NAME, 'NORM_R0  = ', NORM_R0 )
         CALL WRITE_DEBUG(NAME, 'FINAL ERROR  = ', ERROR )
         CALL WRITE_DEBUG(NAME, 'FINAL RESIDUAL ', NORM_R )
         CALL WRITE_DEBUG(NAME, 'ERR RATIO ', ERROR/ERROR0 )
         CALL WRITE_DEBUG(NAME, 'RESID RATIO ', NORM_R/NORM_R0 )
      ENDIF   ! end if (jdebug>=1)

      deallocate(V)
      deallocate(R)
      deallocate(TEMP)
      deallocate(WW)

      RETURN
      END SUBROUTINE LEQ_GMRES0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ROTMAT(  A, B, C, S )

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      DOUBLE PRECISION A,B,C,S
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION ONE,ZERO
      PARAMETER(ONE=1.0D0,ZERO=0.0D0)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION TEMP
!-----------------------------------------------

      IF (B.EQ.ZERO) THEN
            C = ONE
            S = ZERO
      ELSEIF (ABS(B) .GT. ABS(A)) THEN
           TEMP = A / B
           S = ONE / SQRT( ONE + TEMP*TEMP )
           C = TEMP * S
      ELSE
           TEMP = B / A
           C = ONE / SQRT( ONE + TEMP*TEMP )
           S = TEMP * C
      ENDIF

      RETURN
      END

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: Leq_check                                               C
!  Purpose: Verify boundary nodes in A_m are set correctly             C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      subroutine leq_check( vname, A_m )

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE PARAM
      USE PARAM1
      USE GEOMETRY
      USE INDICES
      USE debug
      USE compar
      USE mpi_utility
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3)
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: VNAME
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      logical, parameter :: do_reset = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer, dimension(-3:3) :: ijktable
      integer :: istartl, iendl, jstartl, jendl, kstartl, kendl, elstart, elend
      integer :: i,j,k,el,   ijk,ijk2,   ii,jj,kk
      integer :: nerror, all_nerror
      logical :: is_in_k, is_in_j, is_in_i, is_in
      logical :: is_bc_k, is_bc_j, is_bc_i, is_bc, is_ok
!-----------------------------------------------

      kstartl = kstart2
      kendl = kend2
      jstartl = jstart2
      jendl = jend2
      istartl = istart2
      iendl = iend2

      if (no_k) then
        kstartl  = 1
        kendl = 1
      endif

      nerror = 0

      do k=kstartl,kendl
      do j=jstartl,jendl
      do i=istartl,iendl
         is_in_k = (kstart1 <= k) .and. (k <= kend1)
         is_in_j = (jstart1 <= j) .and. (j <= jend1)
         is_in_i = (istart1 <= i) .and. (i <= iend1)

         is_in = is_in_k .and. is_in_j .and. is_in_i
         if (is_in) cycle

         ijk = funijk(i,j,k)
         ijktable( -2 ) = jm_of(ijk)
         ijktable( -1 ) = im_of(ijk)
         ijktable(  0 ) = ijk
         ijktable(  1 ) = ip_of(ijk)
         ijktable(  2 ) = jp_of(ijk)

         if (.not. no_k) then
            ijktable( -3 ) = km_of(ijk)
            ijktable(  3 ) = kp_of(ijk)
         endif

         elstart = -3
         elend = 3
         if (no_k) then
            elstart = -2
            elend =  2
         endif

         do el=elstart,elend
            ijk2 = ijktable( el )
            kk = k_of(ijk2)
            jj = j_of(ijk2)
            ii = i_of(ijk2)
            is_bc_k = (kk < kstart2) .or. (kk > kend2)
            is_bc_j = (jj < jstart2) .or. (jj > jend2)
            is_bc_i = (ii < istart2) .or. (ii > iend2)
            is_bc = is_bc_k .or. is_bc_j .or. is_bc_i
            if (is_bc) then
               is_ok = (A_m(ijk,el).eq.zero)
               if (.not.is_ok) then
                  nerror = nerror + 1
                  if (do_reset) A_m(ijk,el) = zero
               endif
            endif
         enddo ! do el
      enddo
      enddo
      enddo

      call global_sum( nerror, all_nerror )
      nerror = all_nerror

      if ((nerror >= 1) .and. (myPE.eq.PE_IO)) then
         if (do_reset) then
              call write_debug( 'leq_check: ' // trim( vname ), &
                        'A_m is reset. nerror = ', nerror )
         else
              call write_debug( 'leq_check: ' // trim( vname ), &
                        'A_m is not reset. nerror = ', nerror )
         endif
      endif

      return
      end subroutine leq_check
