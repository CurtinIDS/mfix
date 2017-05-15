!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_LIN_EQ                                            C
!  Purpose: Interface for linear equation solver                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_LIN_EQ(VNAME, Vno, VAR, A_M, B_M, M, ITMAX,&
                              METHOD, SWEEP, TOL1, PC, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE residual
      USE toleranc
      USE leqsol
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number
!     Note: not really used beyond this subroutine. here it is
!     used for potentially adjusting the tolerances but it is
!     currently disabled code.
!     1 = pressure correction equation
!     2 = solids correction equation or gas/solids continuity
!     3 = gas/solids u-momentum
!     4 = gas/solids v-momentum
!     5 = gas/solids w-momentum
!     6 = temperature
!     7 = species
!     8 = granular temperature
!     9 = scalar, E_Turb_G, k_Turb_G
      INTEGER, INTENT(IN) :: Vno
! variable
!     e.g., pp_g, epp, rop_g, rop_s, u_g, u_s, v_g, v_s, w_g,
!     w_s, T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! phase index
      INTEGER, INTENT(IN) :: M
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! linear equation solver method (generally leq_method)
!     1 = sor
!     2 = bicgstab (default)
!     3 = gmres
!     5 = cg
      INTEGER, INTENT(IN) :: METHOD
! sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=4), INTENT(IN) :: SWEEP
! convergence tolerance for leq solver (leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL1
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) :: PC
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Adjust LEQ tolerance flag
      LOGICAL, PARAMETER :: adjust_leq_tol = .FALSE.
      LOGICAL, PARAMETER :: leq_tol_scheme1 = .FALSE.
! currently only used for gmres routine
      INTEGER, PARAMETER :: MAX_IT = 1
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
      DOUBLE PRECISION :: max_resid_local, tol_resid_max
! convergence tolerance for leq solver
      DOUBLE PRECISION :: TOL
! transpose of septadiaganol matrix A_M
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A_MT
! indices
      INTEGER :: II, IJK
! for constructing local character strings
      CHARACTER(LEN=80) :: LINE0, LINE1
!-----------------------------------------------


! Adjusting the tolerances
! ---------------------------------------------------------------->>>
      IF(adjust_leq_tol) THEN
         max_resid_local = maxval(resid(:,M),1)
         tol_resid_max   = max(TOL_RESID, TOL_RESID_T, TOL_RESID_TH, TOL_RESID_X)
         IF(leq_tol_scheme1.AND.resid(Vno,M).LT.1.0D-1) THEN
            if(Vno.le.5) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID)
            elseif (Vno.eq.6) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_T)
            elseif (Vno.eq.7) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_X)
            elseif (Vno.eq.8) then
               TOL = MAX(TOL1,TOL1*RESID(Vno,M)/TOL_RESID_Th)
            endif
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, resid(Vno,M)
         ELSEIF(max_resid_local.LT.1.0D-1) THEN
            TOL = MAX(TOL1,TOL1*max_resid_local/TOL_RESID_max)
            Write(*,*) 'Adjusting LEQ_Tolerance', Vname, tol, max_resid_local
         ENDIF
      ELSE
        TOL = TOL1
      ENDIF
! ----------------------------------------------------------------<<<


! Solve the linear system of equations
! ---------------------------------------------------------------->>>
      SELECT CASE (METHOD)
      CASE (1)
! SOR: Successive Sver Relaxation method from Templates
        CALL LEQ_SOR (VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), &
                      ITMAX, IER)

      CASE (2)
! BICGSTAB: BIConjugate Gradients STabilized method
         IF(do_transpose) THEN  ! mfix.dat keyword default=false
            allocate( A_mt(-3:3, ijkstart3:ijkend3 ))
!!$omp parallel do private(ijk,ii)
            DO ijk=ijkstart3,ijkend3
               do ii=-3,3
                  A_mt(ii,ijk) = A_m(ijk,ii,M)
               enddo
            ENDDO
            call leq_bicgst(VNAME, VNO, VAR, A_Mt(:,:), B_M(:,M), &
                            SWEEP, TOL, PC, ITMAX, IER)
            deallocate( A_mt )
         ELSE
            call leq_bicgs(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, PC, ITMAX, IER)
         ENDIF


      CASE (3)
! GMRES: A Generalized Minimal RESidual Algorithm
         call leq_gmres(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                        SWEEP, TOL, ITMAX, MAX_IT, IER)

      CASE (4)
! Mix:
         IER = 0
         call leq_bicgs(VNAME,VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP,&
                       TOL, PC, ITMAX, IER)
         IF (IER .eq. -2) THEN
            IER = 0
            print*,'calling leq_gmres', Vname
            call leq_gmres(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M),&
                           SWEEP, TOL, ITMAX, MAX_IT, IER)
         ENDIF


      CASE (5)
! CG: Conjugate Gradients
         call leq_cg(VNAME, VNO, VAR, A_M(:,:,M), B_M(:,M), SWEEP,&
                     TOL, PC, ITMAX, IER)

!     CASE (6) - Disabled
! LSOR: Line Successive Over Relaxation method
!       CALL LEQ_LSOR(VNAME, VAR, A_M(:,:,M), B_M(:,M), ITMAX, IER)


      CASE DEFAULT
         LINE0(1:14) = 'SOLVE_LIN_EQ: '
         LINE0(15:80)= VName
         WRITE(LINE1,'(A, I2, A)') &
             'Error: LEQ_METHOD = ', METHOD, ' is invalid'
         CALL WRITE_ERROR(LINE0, LINE1, 1)
         CALL mfix_exit(myPE)
      END SELECT
! ----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE SOLVE_LIN_EQ
