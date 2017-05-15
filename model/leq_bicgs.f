!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine LEQ_BICGS                                                C
!  Purpose: Solve system of linear system using BICGS method           C
!           Biconjugate gradients stabilized                           C
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
! Handan Liu wrote below:               !Jan 22 2013
! The modification is as below:
!       Adding a loop of 2D RSRS sweep and parallelizing for OpenMP.
!       Splitting the existing 3D RSRS loop into two loops for OpenMP
!               due to data dependency
!       Adding openmp directives in all loops in leq_bicgs0.
!
      SUBROUTINE LEQ_BICGS(VNAME, VNO, VAR, A_M, B_m, cmethod, &
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
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: Var
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3,-3:3), INTENT(INOUT) :: A_m
     DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m

! Vector b_m
!      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3), INTENT(INOUT) :: B_m
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
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
! Local Variables
!-------------------------------------------------

      if(PC.eq.'LINE') then   ! default
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, .true., IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, .true., IER )
      elseif(PC.eq.'NONE') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, .false., IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) &
           'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_BICGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_BICGS0                                              C
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

      SUBROUTINE LEQ_BICGS0(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                            TOL, ITMAX, MATVEC, MSOLVE, USE_PC, IER )

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE mpi_utility
      USE sendrecv
      USE indices
      USE leqsol
      USE cutcell
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
!      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(INOUT) :: A_m(ijkstart3:ijkend3,-3:3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3,-3:3), INTENT(INOUT) :: A_m
! Vector b_m
!      DOUBLE PRECISION, INTENT(INOUT) :: B_m(ijkstart3:ijkend3)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: B_m
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! indicate whether to use preconditioner
      LOGICAL, INTENT(IN) :: USE_PC
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
      LOGICAL, PARAMETER :: do_unit_scaling = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION(:), allocatable :: R,Rtilde, Tvec,V
      DOUBLE PRECISION, DIMENSION(:), allocatable, target :: P, P_preconditioned
      DOUBLE PRECISION, DIMENSION(:), allocatable, target :: Svec, Svec_preconditioned

      ! Phat points to either preconditioned value of P, or P itself (to avoid copying for efficiency)
      DOUBLE PRECISION, POINTER :: Phat(:), Shat(:)

      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: &
                        alpha, beta, omega, rho
      DOUBLE PRECISION :: TxS, TxT, RtildexV, &
                          aijmax, oam
      DOUBLE PRECISION :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk
      INTEGER :: iter
      DOUBLE PRECISION, DIMENSION(2) :: TxS_TxT
!-----------------------------------------------

! Initialize the error flag.
      IER=0

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) then

         IF(RE_INDEXING) THEN  ! Loop only over active cells
!$omp parallel do default(shared) private(ijk,oam,aijmax)
            DO IJK = IJKSTART3,IJKEND3
               aijmax = maxval(abs(A_M(ijk,:)) )
               if (aijmax > 0.0)then
                  OAM = one/aijmax
                  A_M(IJK,:) = A_M(IJK,:)*OAM
                  B_M(IJK) = B_M(IJK)*OAM
               else
                  ier = -2
               endif
            ENDDO

         ELSE

!$omp parallel do default(shared) private(ijk,i,j,k,oam,aijmax)
            do k = kstart2,kend2
               do i = istart2,iend2
                  do j = jstart2,jend2
                     IJK = funijk(i,j,k)
                     aijmax = maxval(abs(A_M(ijk,:)) )
                     if (aijmax > 0.0) then
                        OAM = one/aijmax
                        A_M(IJK,:) = A_M(IJK,:)*OAM
                        B_M(IJK) = B_M(IJK)*OAM
                     else
                        ier = -2
                     endif
                  enddo
               enddo
            enddo

         ENDIF

      endif

! A singlular matrix was detected.
      if(IER /= 0) RETURN
! ----------------------------------------------------------------<<<



      allocate(R(DIMENSION_3))
      allocate(Rtilde(DIMENSION_3))
      allocate(P(DIMENSION_3))
      allocate(P_preconditioned(DIMENSION_3))
      allocate(Svec(DIMENSION_3))
      allocate(Svec_preconditioned(DIMENSION_3))
      allocate(Tvec(DIMENSION_3))
      allocate(V(DIMENSION_3))

! these scalars should not be necessary to initialize but done as failsafe
      rnorm = ZERO
      rnorm0 = ZERO
      snorm = ZERO
      pnorm = ZERO

! initializing
      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero

!$omp parallel sections
      R(:) = zero
!$omp section
      Rtilde(:) = zero
!$omp section
      P(:) = zero
!$omp section
      P_preconditioned(:) = zero
!$omp section
      Svec(:) = zero
!$omp section
      Svec_preconditioned(:) = zero
!$omp section
      Tvec(:) = zero
!$omp section
      V(:) = zero
!$omp end parallel sections

      TOLMIN = EPSILON( one )

! Compute initial residual (R = b-A*x) for Ax=b
!    assume initial guess in Var
!    rtilde = r
! ---------------------------------------------------------------->>>
      call MATVEC(Vname, Var, A_M, R)   ! returns R=A*Var

!$omp parallel workshare
      R(:) = B_m(:) - R(:)
!$omp end parallel workshare

      call send_recv(R,nlayers_bicgs)

      Rnorm0 = sqrt( dot_product_par( R, R ) )

! determine an initial guess for the residual = residual + small random
! number (so it could be set to anything). note that since random_number
! is used to supply the guess, this line could potentially be the source
! of small differences between runs.  the random number is shifted below
! between -1 and 1 and then scaled by factor 1.0D-6*Rnorm0
      call random_number(Rtilde(:))

! Shift random number array to be consistent with case when RE_INDEXING is .FALSE.
       IF(RE_INDEXING) CALL SHIFT_DP_ARRAY(Rtilde)

!$omp parallel workshare
       Rtilde(:) = R(:) + (2.0d0*Rtilde(:)-1.0d0)*1.0d-6*Rnorm0
!$omp end parallel workshare

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<


! Main loop
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax

         rho(i-1) = dot_product_par( Rtilde, R )

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
               ier = -2
            else
! Method converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif ! rho(i-1).eq.0

         if (i .eq. 1) then
!$omp parallel workshare
            P(:) = R(:)
!$omp end parallel workshare
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
!$omp parallel workshare
            P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
!$omp end parallel workshare
         endif ! i.eq.1

! Solve A*Phat(:) = P(:)
! V(:) = A*Phat(:)
! --------------------------------
         if (USE_PC) then
            call MSOLVE(Vname, P, A_m, P_preconditioned, CMETHOD) ! returns P_preconditioned
            Phat => P_preconditioned
         else
            Phat => P
         endif

         call MATVEC(Vname, Phat, A_m, V)   ! returns V=A*Phat

         RtildexV = dot_product_par( Rtilde, V )

! compute alpha
! --------------------------------
         alpha(i) = rho(i-1) / RtildexV

! compute Svec
! --------------------------------
!$omp parallel workshare
         Svec(:) = R(:) - alpha(i) * V(:)
!$omp end parallel workshare

! Check norm of Svec(:); if small enough:
! set X(:) = X(:) + alpha(i)*Phat(:) and stop
! --------------------------------
         if(.not.minimize_dotproducts) then
            Snorm = sqrt( dot_product_par( Svec, Svec ) )

            if (Snorm <= TOLMIN) then
!$omp parallel workshare
               Var(:) = Var(:) + alpha(i)*Phat(:)
!$omp end parallel workshare

! Recompute residual norm
! --------------------------------
               if (idebugl >= 1) then
                  call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
!                  Rnorm = sqrt( dot_product_par( Var, Var ) )
!                  print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm

!$omp parallel workshare
                  R(:) = B_m(:) - R(:)
!$omp end parallel workshare

                  Rnorm = sqrt( dot_product_par( R, R ) )
               endif            ! idebugl >= 1
               isConverged = .TRUE.
               EXIT
            endif               ! end if (Snorm <= TOLMIN)
         endif                  ! end if (.not.minimize_dotproducts)

! Solve A*Shat(:) = Svec(:)
! Tvec(:) = A*Shat(:)
! --------------------------------

         if (USE_PC) then
            call MSOLVE(Vname, Svec, A_m, Svec_preconditioned, CMETHOD) ! returns S_preconditioned
            Shat => Svec_preconditioned
         else
            Shat => Svec
         endif

         call MATVEC( Vname, Shat, A_m, Tvec )   ! returns Tvec=A*Shat

         if(.not.minimize_dotproducts) then
!!     $omp parallel sections
            TxS = dot_product_par( Tvec, Svec )
!!     $omp section
            TxT = dot_product_par( Tvec, Tvec )
!!     $omp end parallel sections
         else
            TxS_TxT = dot_product_par2(Tvec, Svec, Tvec, Tvec )
            TxS = TxS_TxT(1)
            TxT = TxS_TxT(2)
         endif

         IF(TxT.eq.Zero) TxT = SMALL_NUMBER

! compute omega
! --------------------------------
         omega(i) = TxS / TxT

! compute new guess for Var
! --------------------------------

!$omp parallel sections
            Var(:) = Var(:) + alpha(i)*Phat(:) + omega(i)*Shat(:)
!$omp section
            R(:) = Svec(:) - omega(i)*Tvec(:)
!$omp end parallel sections

! --------------------------------
         if(.not.minimize_dotproducts.or.(mod(iter,icheck_bicgs).eq.0)) then
            Rnorm = sqrt( dot_product_par(R, R) )

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'iter, Rnorm ', iter, Rnorm, Snorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

!           call mfix_exit(myPE)

! Check convergence; continue if necessary
! for continuation, it is necessary that omega(i) .ne. 0
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) then
               iter_tot(vno) = iter_tot(vno) + iter + 1
               EXIT
            endif
         endif                  ! end if(.not.minimize_dotproducts)

! Advance the iteration count
         iter = iter + 1

      enddo   ! end do i=1,itmax
! end of linear solver loop
! ----------------------------------------------------------------<<<

      if (idebugl >= 1) then
         call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var

!$omp parallel workshare
         R(:) = R(:) - B_m(:)
!$omp end parallel workshare

         Rnorm = sqrt( dot_product_par( R,R) )

         if(myPE.eq.0) print*,'leq_bicgs: final Rnorm ', Rnorm

         if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,  &
         ' L-2', Rnorm/Rnorm0
      endif   ! end if(idebugl >=1)

!      isconverged = (real(Rnorm) <= TOL*Rnorm0);
      if(.NOT.isConverged) isconverged = (real(Rnorm) <= TOL*Rnorm0);
!     write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
      IER = 0
      if (.not.isconverged) then
         IER = -1
         iter_tot(vno) = iter_tot(vno) + iter
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      call send_recv(var,2)

      deallocate(R)
      deallocate(Rtilde)
      deallocate(P)
      deallocate(P_preconditioned)
      deallocate(Svec)
      deallocate(Svec_preconditioned)
      deallocate(Tvec)
      deallocate(V)

      return
      end subroutine LEQ_BICGS0
