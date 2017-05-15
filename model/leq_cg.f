!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_CG                                                  C
!  Purpose: Solve system of linear system using CG method              C
!           conjugate gradients                                        C
!                                                                      C
!  Author: S. Pannala                                 Date: 18-JUL-07  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE LEQ_CG(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                        TOL, PC, ITMAX, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE indices
      USE leqsol
      USE funits
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
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3), INTENT(INOUT) :: A_m
! Vector b_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
! Note: this setting only seems to matter when leq_pc='line'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) ::  PC
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-------------------------------------------------

      if(PC.eq.'LINE') then   ! default
         call LEQ_CG0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_CG0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, IER )
      elseif(PC.eq.'NONE') then
         call LEQ_CG0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) &
           'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_CG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_CG0                                                 C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: S. Pannala                                 Date: 18-JUL-07  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_CG0(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                         TOL, ITMAX, MATVEC, MSOLVE, IER )

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE indices
      USE leqsol
      USE functions
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
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(ijkstart3:ijkend3,-3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(ijkstart3:ijkend3)
! Sweep direction of leq solver (leq_sweep)
!     options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
! dummy arguments/procedures set as indicated
!     matvec->leq_matvec
! for preconditioner (leq_pc)
!    'line' msolve->leq_msolve  (default)
!    'diag' msolve->leq_msolve1
!    'none' msolve->leq_msolve0
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 0
      DOUBLE PRECISION, PARAMETER :: ratiotol = 0.2
      logical, parameter :: do_unit_scaling = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: &
                                  R, P, Zvec, Q
      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: &
                                  alpha, beta, rho
      DOUBLE PRECISION :: RxZ, PxQ, oam
      DOUBLE PRECISION :: Rnorm, Rnorm0, TOLMIN
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk
      INTEGER :: iter
!-----------------------------------------------

      is_serial = numPEs.eq.1.and.is_serial

! these scalars should not be necessary to initialize but done as failsafe
      rnorm = ZERO
      rnorm0 = ZERO

! initializing
      alpha(:)  = zero
      beta(:)   = zero
      rho(:)    = zero


! zero out R,Zvec, P and Q
! --------------------------------
      if (use_doloop) then
!!$omp  parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = zero
            Zvec(ijk) = zero
            P(ijk) = zero
            Q(ijk) = zero
         enddo
      else
         R(:) = zero
         Zvec(:) = zero
         P(:) = zero
         Q(:) = zero
      endif

      TOLMIN = EPSILON( one )

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) then
!!$omp parallel do private(ijk,i,j,k,oam,aijmax)
         do k = kstart2,kend2
            do i = istart2,iend2
               do j = jstart2,jend2
                  IJK = funijk(i,j,k)
!                  aijmax = maxval(abs(A_M(ijk,:)) )
!                  if(aijmax.ne.abs(A_M(ijk,0))) &
!                  write(*,*) 'Not positive definite', k,i,j,(A_M(ijk,:))
                  OAM = one/A_M(ijk,0)
                  A_M(IJK,:) = A_M(IJK,:)*OAM
                  B_M(IJK) = B_M(IJK)*OAM
               enddo
            enddo
         enddo
      endif
! ----------------------------------------------------------------<<<


! assume initial guess in Var + some small random number
! r = b - A*x : Line 1
!     call random_number(Xinit(:))
!     if (use_doloop) then
!!$omp   parallel do private(ijk)
!        do ijk=ijkstart3,ijkend3
!           Xinit(ijk) = Var(ijk)*(ONE + (2.0d0*Xinit(ijk)-1.0d0)*1.0d-6)
!        enddo
!     else
!        Xinit(:) = Var(:)* (ONE + (2.0d0*Xinit(:)-1.0d0)*1.0d-6)
!     endif
!     Xinit(:) = Zero

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_cg, initial: ', Vname,' resid ', Rnorm0
      endif


! Compute initial residual, R = b - A*x
! ---------------------------------------------------------------->>>
      call MATVEC(Vname, Var, A_M, R)   ! returns R=A_M*Var

      if (use_doloop) then
!!$omp   parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            R(ijk) = B_m(ijk) - R(ijk)
         enddo
      else
         R(:) = B_m(:) - R(:)
      endif

      if(is_serial) then
         Rnorm0 = zero
         if (use_doloop) then
!!$omp          parallel do private(ijk) reduction(+:Rnorm0)
            do ijk=ijkstart3,ijkend3
               Rnorm0 = Rnorm0 + R(ijk)*R(ijk)
            enddo
         else
            Rnorm0 = dot_product(R,R)
         endif
         Rnorm0 = sqrt( Rnorm0 )
      else
         Rnorm0 = sqrt( dot_product_par( R, R ) )
      endif

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_cg, initial: ', Vname,' resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<


! Main loop : Line 2
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax
! Solve M Zvec(:) = R(:) : Line 3
! --------------------------------
         call MSOLVE( Vname, R, A_m, Zvec, CMETHOD)   ! returns Zvec

! Solve Rho = RxZ : Line 4
! --------------------------------
         if(is_serial) then
            if (use_doloop) then
               RxZ = zero
!!$omp        parallel do private(ijk) reduction(+:RxZ)
               do ijk=ijkstart3,ijkend3
                  RxZ = RxZ + R(ijk) * Zvec(ijk)
               enddo
               rho(i-1) = RxZ
            else
               rho(i-1) = dot_product( R, Zvec )
            endif
         else
            rho(i-1) = dot_product_par( R, Zvec )
         endif                  ! is_serial


         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
               ier = -2
            else
! converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif                  ! rho(i-1).eq.0

         if (i .eq. 1) then
! P_1 = Z_0 : Line 6
! --------------------------------
            if (use_doloop) then
!!$omp        parallel do private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = Zvec(ijk)
               enddo
            else
               P(:) = Zvec(:)
            endif
         else
! beta = rho(i-1)/rho(i-2) : Line 8
! P = Z + beta*P : Line 9
! --------------------------------
            beta(i-1) = ( rho(i-1)/rho(i-2) )
            if (use_doloop) then
!!!$omp        parallel do private(ijk)
               do ijk=ijkstart3,ijkend3
                  P(ijk) = Zvec(ijk) + beta(i-1)* P(ijk)
               enddo
            else
               P(:) = Zvec(:) + beta(i-1)*P(:)
            endif
         endif                  ! i.eq.1

! Q(:) = A*P(:) : Line 10
! --------------------------------
         call MATVEC(Vname, P, A_m, Q)   ! Returns Q = A_m*P

         if(is_serial) then
            if (use_doloop) then
               PxQ = zero
!!$omp         parallel do private(ijk) reduction(+:PxQ)
               do ijk=ijkstart3,ijkend3
                  PxQ = PxQ + P(ijk) * Q(ijk)
               enddo
            else
               PxQ = dot_product( P, Q )
            endif
         else
            PxQ = dot_product_par( P, Q )
         endif                  ! is_serial

!  alpha = rho/PxQ : Line 11
! --------------------------------
         alpha(i) = rho(i-1)/PxQ

! x = x + alpha*p : Line 12
! r = r - alpha*q : Line 13
! --------------------------------
         if (use_doloop) then
!!$omp     parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - alpha(i) * Q(ijk)
               Var(ijk) = Var(ijk) + alpha(i) * P(ijk)
            enddo
         else
            R(:) = R(:) - alpha(i) * Q(:)
            Var(:) = Var(:) + alpha(i) * P(:)
         endif                  ! use_doloop

! Check norm of R(:); if small enough, Exit
! --------------------------------
         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!!$omp       parallel do private(ijk) reduction(+:Rnorm)
               do ijk=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk) * R(ijk)
               enddo
            else
               Rnorm = dot_product( R, R )
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R, R ) )
         endif                  ! is_serial

         if (idebugl.ge.1) then
            print*,'leq_cs, initial: ', Vname,' Vnorm ', Rnorm
            if (myPE.eq.PE_IO) then
               print*,'iter, Rnorm ', iter, Rnorm
               print*,'PxQ, rho(i-1) ', PxQ, rho(i-1)
               print*,'alpha(i), beta(i-1) ', alpha(i), beta(i-1)
            endif
         endif

         isconverged = (Rnorm <= TOL*Rnorm0)

         if (isconverged) then
            iter_tot(vno) = iter_tot(vno) + iter + 1
            EXIT
         endif

! Advance the iteration count
         iter = iter + 1

      enddo
! end of linear solver loop
! ----------------------------------------------------------------<<<


      if (idebugl >= 1) then
         call MATVEC( Vname, Var, A_m, R )
         if (use_doloop) then
!!$omp  parallel do private(ijk)
            do ijk=ijkstart3,ijkend3
               R(ijk) = R(ijk) - B_m(ijk)
            enddo
         else
            R(:) = R(:) - B_m(:)
         endif

         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!!$omp         parallel do private(ijk) reduction(+:Rnorm)
               do ijk=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk) * R(ijk)
               enddo
            else
               Rnorm = dot_product( R,R)
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R,R) )
         endif

         if(myPE.eq.0) print*,'leq_cg: final Rnorm ', Rnorm

         if(myPE.eq.0)  print*,'leq_cg ratio : ', Vname,' ',iter,     &
         ' L-2', Rnorm/Rnorm0
      endif

      IER = 0
      if (.not.isconverged) then
         IER = -1
         iter_tot(vno) = iter_tot(vno) + iter
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      call send_recv(var,2)

      return
      end subroutine LEQ_CG0
