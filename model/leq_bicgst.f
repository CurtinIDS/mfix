!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_BICGSt(Vname, Var, A_m, B_m,                       C
!                         cmethod, TOL, ITMAX, IER )
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_BICGSt(VNAME, VNO, VAR, A_M, B_m,  cmethod, TOL, PC, ITMAX,IER)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE exit, only: mfix_exit
      USE funits
      USE geometry
      USE indices
      USE leqsol
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error indicator
      INTEGER ::          IER
!                      maximum number of iterations
      INTEGER ::          ITMAX
!                      variable number
      INTEGER ::          VNO
!                      convergence tolerance
      DOUBLE PRECISION ::  TOL
!                      Preconditioner
      CHARACTER(LEN=4)   ::  PC
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(-3:3,ijkstart3:ijkend3) :: A_m
!                      Vector b_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: B_m
!                      Variable name
      CHARACTER(LEN=*) ::    Vname
!                      Variable
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: Var
!                    sweep direction
      CHARACTER(LEN=*) :: CMETHOD
!
!-------------------------------------------------
      EXTERNAL LEQ_MATVECt, LEQ_MSOLVEt, LEQ_MSOLVE0t, LEQ_MSOLVE1t
!--------------------------------------------------

      if(PC.eq.'LINE') then
         call LEQ_BICGS0t( Vname, Vno, Var, A_m, B_m,                        &
         cmethod, TOL, ITMAX, LEQ_MATVECt, LEQ_MSOLVEt, IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_BICGS0t( Vname, Vno, Var, A_m, B_m,                        &
         cmethod, TOL, ITMAX, LEQ_MATVECt, LEQ_MSOLVE1t, IER )
      elseif(PC.eq.'NONE') then
         call LEQ_BICGS0t( Vname, Vno, Var, A_m, B_m,                        &
         cmethod, TOL, ITMAX, LEQ_MATVECt, LEQ_MSOLVE0t, IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_BICGSt




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_BICGS0(Vname, Var, A_m, B_m,                       C
!                         cmethod, TOL, ITMAX, MATVEC, MSOLVE, IER )   C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_BICGS0t(VNAME, VNO, VAR, A_M, B_m,  cmethod, TOL, ITMAX,  &
                            MATVECt, MSOLVEt, IER )
!-----------------------------------------------
!   M o d u l e s
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
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error indicator
      INTEGER ::          IER
!                      maximum number of iterations
      INTEGER ::          ITMAX
!                      variable number
      INTEGER ::          VNO
!                      convergence tolerance
      DOUBLE PRECISION ::  TOL
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION, DIMENSION(-3:3,ijkstart3:ijkend3) :: A_m
!                      Vector b_m
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: B_m
!                      Variable name
      CHARACTER(LEN=*) ::    Vname
!                      Variable
      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) :: Var
!                    sweep direction
      CHARACTER(LEN=*) :: CMETHOD
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 1
      DOUBLE PRECISION :: ratiotol = 0.2

      DOUBLE PRECISION, DIMENSION(ijkstart3:ijkend3) ::                       &
                                R,Rtilde, P,Phat, Svec, Shat, Tvec,V
      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: alpha,beta,omega,rho
      DOUBLE PRECISION :: TxS, TxT, oam,RtildexV,                   &
                      RtildexR, aijmax, Rnorm=0, Rnorm0, Snorm, TOLMIN
      LOGICAL :: isconverged
      INTEGER :: i, j, k, iter
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
      INTEGER :: ijk2
      DOUBLE PRECISION, DIMENSION(2) :: TxS_TxT

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL  MATVECt, MSOLVEt

      logical, parameter :: do_unit_scaling = .true.

      is_serial = numPEs.eq.1.and.is_serial

      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero

!
!     ---------------------------------------------
!     zero out R,Rtilde, P,Phat, Svec, Shat, Tvec,V
!     ---------------------------------------------
      if (use_doloop) then

!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 201)
!!$omp  parallel do private(ijk2)
         do ijk2=ijkstart3,ijkend3
            R(ijk2) = zero
            Rtilde(ijk2) = zero
            P(ijk2) = zero
            Phat(ijk2) = zero
            Svec(ijk2) = zero
            Shat(ijk2) = zero
            Tvec(ijk2) = zero
            V(ijk2) = zero
         enddo

      else

         R(:) = zero
         Rtilde(:) = zero
         P(:) = zero
         Phat(:) = zero
         Svec(:) = zero
         Shat(:) = zero
         Tvec(:) = zero
         V(:) = zero

      endif


      TOLMIN = EPSILON( one )

      if (do_unit_scaling) then
!
!     Scale matrix to have unit diagonal
!
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 233)
!!$omp parallel do private(ijk2,i,j,k,oam,aijmax) collapse (3)
         do k = kstart2,kend2
            do i = istart2,iend2
               do j = jstart2,jend2
                  IJK2 = funijk(i,j,k)

                  aijmax = maxval(abs(A_M(:,ijk2)) )

                  OAM = one/aijmax

                  A_M(:,IJK2) = A_M(:,IJK2)*OAM

                  B_M(IJK2) = B_M(IJK2)*OAM

               enddo
            enddo
         enddo
      endif

!
!    Compute initial residual, assume initial guess in Var
!    r = b - A*x
!    rtilde = r


      call MATVECt( Vname, Var, A_M, R )


      if (use_doloop) then
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 264)
!!!$omp   parallel do private(ijk2)
         do ijk2=ijkstart3,ijkend3
            R(ijk2) = B_m(ijk2) - R(ijk2)
         enddo
      else
         R(:) = B_m(:) - R(:)
      endif

      if(is_serial) then
         Rnorm0 = zero
         if (use_doloop) then
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 276)
!!!$omp          parallel do private(ijk2) reduction(+:Rnorm0)
            do ijk2=ijkstart3,ijkend3
               Rnorm0 = Rnorm0 + R(ijk2)*R(ijk2)
            enddo
         else
            Rnorm0 = dot_product(R,R)
         endif
         Rnorm0 = sqrt( Rnorm0 )
      else
         Rnorm0 = sqrt( dot_product_par( R, R ) )
      endif

      call random_number(Rtilde(:))

      if (use_doloop) then
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
!!!$omp   parallel do private(ijk2)
         do ijk2=ijkstart3,ijkend3
            Rtilde(ijk2) = R(ijk2) + (2.0d0*Rtilde(ijk2)-1.0d0)*1.0d-6*Rnorm0
         enddo
      else
         Rtilde(:) = R(:) + (2.0d0*Rtilde(:)-1.0d0)*1.0d-6*Rnorm0
      endif

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
      endif
!
!     Main loop
!
      iter = 1
      do i=1,itmax

         if(is_serial) then
            if (use_doloop) then
               RtildexR = zero
!AIKE PFUPGRADE 091409
!!!$omp        parallel do private(ijk2) reduction(+:RtildexR)
               do ijk2=ijkstart3,ijkend3
                  RtildexR = RtildexR + Rtilde(ijk2) * R(ijk2)
               enddo
               rho(i-1) = RtildexR
            else
               rho(i-1) = dot_product( Rtilde, R )
            endif
         else
            rho(i-1) = dot_product_par( Rtilde, R )
         endif ! is_serial

!     print*,'leq_bicgs, initial: ', Vname,' rho(i-1) ', rho(i-1)

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
!     ------------
!     Method fails
!     ------------
!     print*, 'leq_bicgs,',Vname,': rho(i-1) == 0 '
               ier = -2
            else
!     ------------
!     converged.  residual is already zero
!     ------------
               ier = 0
            endif
            call send_recv(var,2)
            return
         endif ! rho(i-1).eq.0

         if (i .eq. 1) then
            if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp        parallel do private(ijk2)
               do ijk2=ijkstart3,ijkend3
                  P(ijk2) = R(ijk2)
               enddo
            else
               P(:) = R(:)
            endif
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
            if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp        parallel do private(ijk2)
               do ijk2=ijkstart3,ijkend3
                  P(ijk2) = R(ijk2) + beta(i-1)*( P(ijk2) - omega(i-1)*V(ijk2) )
               enddo
            else
               P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
            endif
         endif ! i.eq.1

!
!     Solve M Phat(:) = P(:)
!     V(:) = A*Phat(:)
!

         call MSOLVEt( Vname, P, A_m, Phat, CMETHOD)

         call MATVECt( Vname, Phat, A_m, V )

         if(is_serial) then
            if (use_doloop) then
               RtildexV = zero
!AIKE PFUPGRADE 091409
!!!$omp         parallel do private(ijk2) reduction(+:RtildexV)
               do ijk2=ijkstart3,ijkend3
                  RtildexV = RtildexV + Rtilde(ijk2) * V(ijk2)
               enddo
            else
               RtildexV = dot_product( Rtilde, V )
            endif
         else
            RtildexV = dot_product_par( Rtilde, V )
         endif ! is_serial

!     print*,'leq_bicgs, initial: ', Vname,' RtildexV ', RtildexV

         alpha(i) = rho(i-1) / RtildexV

         if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp     parallel do private(ijk2)
            do ijk2=ijkstart3,ijkend3
               Svec(ijk2) = R(ijk2) - alpha(i) * V(ijk2)
            enddo
         else
            Svec(:) = R(:) - alpha(i) * V(:)
         endif ! use_doloop

         if(.not.minimize_dotproducts) then
!
!     Check norm of Svec(:); if small enough:
!     set X(:) = X(:) + alpha(i)*Phat(:) and stop
!
            if(is_serial) then
               if (use_doloop) then
                  Snorm = zero
!AIKE PFUPGRADE 091409
!!!$omp       parallel do private(ijk2) reduction(+:Snorm)
                  do ijk2=ijkstart3,ijkend3
                     Snorm = Snorm + Svec(ijk2) * Svec(ijk2)
                  enddo
               else
                  Snorm = dot_product( Svec, Svec )
               endif
               Snorm = sqrt( Snorm )
            else
               Snorm = sqrt( dot_product_par( Svec, Svec ) )
            endif               ! is_serial
!     print*,'leq_bicgs, initial: ', Vname,' Snorm ', real(Snorm)


            if (Snorm <= TOLMIN) then
               if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp          parallel do private(ijk2)
                  do ijk2=ijkstart3,ijkend3
                     Var(ijk2) = Var(ijk2) + alpha(i)*Phat(ijk2)
                  enddo
               else
                  Var(:) = Var(:) + alpha(i)*Phat(:)
               endif            ! use_doloop

               if (idebugl >= 1) then
!
!     Recompute residual norm
!
                  call MATVECt( Vname, Var, A_m, R )

!     Rnorm = sqrt( dot_product_par( Var, Var ) )
!     print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm

                  if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp          parallel do private(ijk2)
                     do ijk2=ijkstart3,ijkend3
                        R(ijk2) = B_m(ijk2) - R(ijk2)
                     enddo
                  else
                     R(:) = B_m(:) - R(:)
                  endif

                  if(is_serial) then
                     if (use_doloop) then
                        Rnorm = zero
!AIKE PFUPGRADE 091409
!!!$omp            parallel do private(ijk2) reduction(+:Rnorm)
                        do ijk2=ijkstart3,ijkend3
                           Rnorm = Rnorm + R(ijk2)*R(ijk2)
                        enddo
                     else
                        Rnorm =  dot_product( R, R )
                     endif
                     Rnorm = sqrt( Rnorm )
                  else
                     Rnorm = sqrt( dot_product_par( R, R ) )
                  endif
!     print*,'leq_bicgs, initial: ', Vname,' Rnorm ', Rnorm
               endif            ! idebugl >= 1

               EXIT
            endif               ! Snorm <= TOLMIN

         endif                  ! .not.minimize_dotproducts

!
!     Solve M Shat(:) = Svec(:)
!     Tvec(:) = A * Shat(:)
!
         call MSOLVEt( Vname, Svec, A_m, Shat, CMETHOD)

         call MATVECt( Vname, Shat, A_m, Tvec )

         if(is_serial) then
            if (use_doloop) then
               TxS = zero
               TxT = zero
!AIKE PFUPGRADE 091409
!!!$omp  parallel do private(ijk2) reduction(+:TxS,TxT)
               do ijk2=ijkstart3,ijkend3
                  TxS = TxS + Tvec(ijk2)  * Svec(ijk2)
                  TxT = TxT + Tvec(ijk2)  * Tvec(ijk2)
               enddo
            else
               TxS = dot_product( Tvec, Svec )
               TxT = dot_product( Tvec, Tvec )
            endif
         else
            if(.not.minimize_dotproducts) then
               TxS = dot_product_par( Tvec, Svec )
               TxT = dot_product_par( Tvec, Tvec )
            else
               TxS_TxT = dot_product_par2(Tvec, Svec, Tvec, Tvec )
               TxS = TxS_TxT(1)
               TxT = TxS_TxT(2)
            endif
         endif
         IF(TxT.eq.Zero) TxT = SMALL_NUMBER
         omega(i) = TxS / TxT


         if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp    parallel do private(ijk2)
            do ijk2=ijkstart3,ijkend3
               Var(ijk2) = Var(ijk2) +                           &
               alpha(i)*Phat(ijk2) + omega(i)*Shat(ijk2)
               R(ijk2) = Svec(ijk2) - omega(i)*Tvec(ijk2)
            enddo
         else
            Var(:) = Var(:) +                           &
            alpha(i)*Phat(:) + omega(i)*Shat(:)
            R(:) = Svec(:) - omega(i)*Tvec(:)
         endif

         if(.not.minimize_dotproducts.or.(mod(iter,5).eq.0)) then
            if(is_serial) then
               if (use_doloop) then
                  Rnorm = zero
!AIKE PFUPGRADE 091409
!!!$omp       parallel do private(ijk2) reduction(+:Rnorm)
                  do ijk2=ijkstart3,ijkend3
                     Rnorm = Rnorm + R(ijk2) * R(ijk2)
                  enddo
               else
                  Rnorm =  dot_product(R, R )
               endif
               Rnorm = sqrt( Rnorm )
            else
               Rnorm = sqrt( dot_product_par(R, R) )
            endif               ! is_serial

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'iter, Rnorm ', iter, Rnorm, Snorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

!     call mfix_exit(myPE)

!     Check convergence; continue if necessary
!     for continuation, it is necessary that omega(i) .ne. 0
!
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) then
               iter_tot(vno) = iter_tot(vno) + iter + 1
               EXIT
            endif

         endif                  ! .not.minimize_dotproducts

!     Advance the iteration count
         iter = iter + 1

      enddo

      if (idebugl >= 1) then
         call MATVECt( Vname, Var, A_m, R )
         if (use_doloop) then
!AIKE PFUPGRADE 091409
!!!$omp  parallel do private(ijk2)
            do ijk2=ijkstart3,ijkend3
               R(ijk2) = R(ijk2) - B_m(ijk2)
            enddo
         else
            R(:) = R(:) - B_m(:)
         endif

         if(is_serial) then
            if (use_doloop) then
               Rnorm = zero
!AIKE PFUPGRADE 091409
!!!$omp         parallel do private(ijk2) reduction(+:Rnorm)
               do ijk2=ijkstart3,ijkend3
                  Rnorm = Rnorm + R(ijk2) * R(ijk2)
               enddo
            else
               Rnorm = dot_product( R,R)
            endif
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R,R) )
         endif

         if(myPE.eq.0) print*,'leq_bicgs: final Rnorm ', Rnorm

         if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,     &
         ' L-2', Rnorm/Rnorm0
      endif

      isconverged = (real(Rnorm) <= TOL*Rnorm0);
!     write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
      IER = 0
      if (.not.isconverged) then
         iter_tot(vno) = iter_tot(vno) + iter
         IER = -1
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      call send_recv(var,2)

      return
      end subroutine LEQ_BICGS0t


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_ISWEEPt(I, Vname, Var, A_m, B_m )                  C
!  Purpose: Perform line sweep at coordiante I                         C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_ISWEEPt(I,Vname, VAR, A_M, B_M)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE indices
      USE funits
      USE sendrecv
      USE mpi_utility
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          I
!
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3,ijkstart3:ijkend3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER(LEN=*)    Vname

!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NSTART, NEND
      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC,DD,EE,BB
      INTEGER :: INFO, IJK, J, K, IM1JK, IP1JK

      NEND = JEND
      NSTART = JSTART
      K = 1

      DO J=NSTART, NEND
!     IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IJK = FUNIJK(I,J,K)
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)

         DD(J) = A_M(0, IJK)
         CC(J) = A_M(-2, IJK)
         EE(J) = A_M(2, IJK)
         BB(J) = B_M(IJK) -  A_M(-1, IJK) * Var( IM1JK )         &
         -  A_M(1, IJK) * Var( IP1JK )

      END DO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0

!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB,  JEND-JSTART+1, INFO )

      IF (INFO.NE.0) THEN
         RETURN
      ENDIF

      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) =  BB(J)
      ENDDO

      RETURN
      END SUBROUTINE  LEQ_ISWEEPt
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_IKSWEEPt(I, K, Vname, Var, A_m, B_m )              C
!  Purpose: Perform line sweep at coordiante I,K                       C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_IKSWEEPt(I,K,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE funits
      USE indices
      USE sendrecv
      USE mpi_utility
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          I,K
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3,ijkstart3:ijkend3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER(LEN=*)    Vname

!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)
!

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!


!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC,DD,EE, BB
      INTEGER :: NSTART, NEND, INFO, IJK, J,  IM1JK, IP1JK, IJKM1, IJKP1

      NEND = JEND
      NSTART = JSTART

!!!$omp parallel do private(j,ijk,im1jk,ip1jk,ijkm1,ijkp1)
      DO J=NSTART, NEND

!     IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IJK = FUNIJK(I,J,K)
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)
         IJKM1 = KM_OF(IJK)
         IJKP1 = KP_OF(IJK)

         DD(J) = A_M(0, IJK)
         CC(J) = A_M(-2, IJK)
         EE(J) = A_M(2, IJK)
         BB(J) = B_M(IJK) -  A_M(-1, IJK) * Var( IM1JK )         &
         -  A_M(1, IJK) * Var( IP1JK )         &
         -  A_M(-3, IJK) * Var( IJKM1 )         &
         -  A_M(3, IJK) * Var( IJKP1 )

      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB,  JEND-JSTART+1, INFO )

      IF (INFO.NE.0) THEN
         write(*,*) 'leq_iksweep',INFO, myPE
         RETURN
      ENDIF

      DO J=NSTART, NEND

         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(J)

      ENDDO

      RETURN
      END SUBROUTINE  LEQ_IKSWEEPt


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_MATVECt(Vname, Var, A_m, B_m )                     C
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_MATVECt(VNAME, VAR, A_M, Avar )
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      USE mpi_utility
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
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3, ijkstart3:ijkend3)
!
!                      Vector AVar
      DOUBLE PRECISION AVar(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER(LEN=*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)
!                      Variable
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
      INTEGER          I,  J, K, IJK2

      integer :: im1jk,ip1jk, ijm1k,ijp1k, ijkm1,ijkp1

      if (do_k) then

!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 938)
!!$omp    parallel  do &
!!$omp&   private(     &
!!$omp&           ijk2,i,j,k, &
!!$omp&           im1jk,ip1jk,ijm1k,ijp1k,ijkm1,ijkp1) collapse (3)
         do k = kstart,kend
            do i = istart,iend
               do j = jstart,jend

                  IJK2 = funijk(i,j,k)

                  im1jk = im_of(ijk2)
                  ip1jk = ip_of(ijk2)
                  ijm1k = jm_of(ijk2)
                  ijp1k = jp_of(ijk2)
!
                  ijkm1 = km_of(ijk2)
                  ijkp1 = kp_of(ijk2)


                  AVar(ijk2) =      A_m(-3, ijk2) * Var(ijkm1)   &
                  + A_m(-2, ijk2) * Var(ijm1k)   &
                  + A_m(-1, ijk2) * Var(im1jk)   &
                  + A_m(0, ijk2) * Var(ijk2)     &
                  + A_m(1, ijk2) * Var(ip1jk)   &
                  + A_m(2, ijk2) * Var(ijp1k)   &
                  + A_m(3, ijk2) * Var(ijkp1)

               enddo
            enddo
         enddo

      else
         k = 1
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 971)
!!$omp parallel do private(i,j,ijk2,im1jk,ip1jk,ijm1k,ijp1k) collapse (2)
         do i = istart,iend
            do j = jstart,jend


               IJK2 = funijk(i,j,k)


               im1jk = im_of(ijk2)
               ip1jk = ip_of(ijk2)
               ijm1k = jm_of(ijk2)
               ijp1k = jp_of(ijk2)
               AVar(ijk2) =      A_m(-2, ijk2) * Var(ijm1k)   &
               + A_m(-1, ijk2) * Var(im1jk)   &
               + A_m(0, ijk2) * Var(ijk2)     &
               + A_m(1, ijk2) * Var(ip1jk)   &
               + A_m(2, ijk2) * Var(ijp1k)

            enddo
         enddo

      endif

      call send_recv(Avar,nlayers_bicgs)
      return
      END SUBROUTINE LEQ_MATVECt
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_MSOLVE(Vname, B_m, A_m, Var, CMETHOD)
!  Purpose: Successive line over-relaxation method -- Cyclic bc        C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_MSOLVEt(VNAME, B_m, A_M, Var, CMETHOD)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3, ijkstart3:ijkend3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER(LEN=*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)


!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
      INTEGER ::   I, J, K, ITER, NITER, IJK2
      INTEGER ::   I1 , K1 , I2, K2, IK, ISIZE, KSIZE
      INTEGER ::   ICASE

!     CHARACTER(LEN=4), PARAMETER :: CMETHOD = 'II'
!     sweep direction
      CHARACTER(LEN=4) :: CMETHOD
      CHARACTER :: CH
      LOGICAL :: DO_ISWEEP, DO_JSWEEP, DO_KSWEEP
      LOGICAL :: DO_SENDRECV, DO_REDBLACK
      LOGICAL, PARAMETER :: USE_IKLOOP = .FALSE.

      LOGICAL, PARAMETER :: SETGUESS = .TRUE.

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime
!       by Tingwen
!!$      omp_start=omp_get_wtime()

      IF (SETGUESS) THEN

!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 1077)
!!$omp   parallel do private(i,j,k,ijk2) collapse (3)
         do k = kstart3,kend3
            do i = istart3,iend3
               do j = jstart3,jend3

                  IJK2 = funijk(i,j,k)

                  VAR(IJK2) = B_M(IJK2)

               enddo
            enddo
         enddo

      call send_recv(var,nlayers_bicgs)

      ENDIF

      NITER = LEN( CMETHOD )

      DO ITER=1,NITER
!
!     Perform sweeps
!
         CH = CMETHOD( ITER:ITER )
         DO_ISWEEP = (CH .EQ. 'I') .OR. (CH .EQ. 'i')
         DO_JSWEEP = (CH .EQ. 'J') .OR. (CH .EQ. 'j')
         DO_KSWEEP = (CH .EQ. 'K') .OR. (CH .EQ. 'k')
         DO_SENDRECV = (CH .EQ. 'S') .OR. (CH .EQ. 's')
         DO_REDBLACK = (CH .EQ. 'R') .OR. (CH .EQ. 'r')

         IF (NO_K) THEN
! 2D run no need to enable openmp parallel
            IF ( DO_ISWEEP ) THEN
!!!$omp   parallel do private(I)
               DO I=istart,iend,1
                  CALL LEQ_ISWEEPt( I, Vname, Var, A_m, B_m )
               ENDDO
            ENDIF

         ELSE

            IF(DO_REDBLACK) THEN

               i1 = istart
               k1 = kstart
               i2 = iend
               k2 = kend
               isize = i2-i1+1
               ksize = k2-k1+1

               DO icase = 1, 2
!!$omp   parallel do private(K,I,IK)
                  DO IK=icase, ksize*isize, 2
                     if (mod(ik,isize).ne.0) then
                        k = int( ik/isize ) + k1
                     else
                        k = int( ik/isize ) + k1 -1
                     endif
                     i = (ik-1-(k-k1)*isize) + i1 + mod(k,2)
                     if(i.gt.i2) i=i-i2 + i1 -1
                     CALL LEQ_IKSWEEPt( I,K, Vname, Var, A_m, B_m )
                  ENDDO
               ENDDO

            ENDIF
!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
            IF(USE_IKLOOP) THEN

               i1 = istart
               k1 = kstart
               i2 = iend
               k2 = kend
               isize = i2-i1+1
               ksize = k2-k1+1

               IF (DO_ISWEEP) THEN
!!!$omp   parallel do private(K,I,IK)
                  DO IK=1, ksize*isize
                     if (mod(ik,isize).ne.0) then
                        k = int( ik/isize ) + k1
                     else
                        k = int( ik/isize ) + k1 -1
                     endif
                     i = (ik-1-(k-k1)*isize) + i1
                     CALL LEQ_IKSWEEPt( I,K, Vname, Var, A_m, B_m )
                  ENDDO
               ENDIF

               IF (DO_KSWEEP) THEN
!!!$omp   parallel do private(K,I,IK)
                  DO IK=1, ksize*isize
                     if (mod(ik,ksize).ne.0) then
                        i = int( ik/ksize ) + i1
                     else
                        i = int( ik/ksize ) + i1 -1
                     endif
                     k = (ik-1-(i-i1)*ksize) + k1

                     CALL LEQ_IKSWEEPt( I,K, Vname, Var, A_m, B_m )
                  ENDDO
               ENDIF

            ELSE
!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
               IF (DO_ISWEEP) THEN
!!!$omp   parallel do private(K,I)
                  DO K=kstart,kend
                     DO I=istart,iend
                        CALL LEQ_IKSWEEPt( I,K, Vname, Var, A_m, B_m )
                     ENDDO
                  ENDDO
               ENDIF

               IF (DO_KSWEEP) THEN
!!!$omp   parallel do private(K,I)
                  DO I=istart,iend
                     DO K=kstart,kend
                        CALL LEQ_IKSWEEPt( I,K, Vname, Var, A_m, B_m )
                     ENDDO
                  ENDDO
               ENDIF

            ENDIF

            IF (DO_SENDRECV) call send_recv(var,nlayers_bicgs)

         ENDIF

      ENDDO
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'leq_msolvet:',omp_end - omp_start
      RETURN
      END SUBROUTINE LEQ_MSOLVEt


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_JKSWEEPt(J, K, Vname, Var, A_m, B_m )              C
!  Purpose: Perform line sweep at coordiante I,K                       C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_JKSWEEPt(J,K,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE exit, only: mfix_exit
      USE functions
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          J,K
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3, ijkstart3:ijkend3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER(LEN=*)    Vname
!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION (IMAX2) :: CC,DD,EE,BB
      INTEGER :: NN, INFO, IJK, I

      NN = IMAX2

      DO I=1,NN
         IJK = FUNIJK(I,J,K)

         DD(I) = A_M(0, IJK)
         CC(I) = A_M(-1, IJK)
         EE(I) = A_M(1, IJK)
         BB(I) = B_M(IJK)    -  A_M(-2,IJK) * Var( JM_OF(IJK) )         &
         -  A_M(2, IJK) * Var( JP_OF(IJK) )         &
         -  A_M(-3, IJK) * Var( KM_OF(IJK) )         &
         -  A_M(3, IJK) * Var( KP_OF(IJK) )

      ENDDO

      CC(1) = ZERO
      EE(NN) = ZERO
!     DL(1:NEND-1) = CC(2:NEND)
      INFO = 0
      CALL DGTSL( NN, CC, DD, EE, BB, INFO )
!     CALL DGTSV( JEND-JSTART+1, 1, DL, DD, EE, BB,  JEND-JSTART+1, INFO )

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO I=1,NN
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(I)
      ENDDO

      RETURN
      END SUBROUTINE  LEQ_JKSWEEPt

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_IJSWEEPt(I,J, Vname, Var, A_m, B_m )               C
!  Purpose: Perform line sweep at coordiante I,K                       C
!                                                                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_IJSWEEPt(I,J,Vname, VAR, A_M, B_M )

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE exit, only: mfix_exit
      USE functions
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Line position
      INTEGER          I,J
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3, ijkstart3:ijkend3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!                      Variable name
      CHARACTER(LEN=*)    Vname

!
!                      Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION (KMAX2) :: CC,DD,EE,BB
      INTEGER :: NN, INFO, IJK, K

      NN = KMAX2

      DO K=1,NN
         IJK = FUNIJK(I,J,K)

         DD(K) = A_M(0,IJK)
         CC(K) = A_M(-3,IJK)
         EE(K) = A_M(3,IJK)
         BB(K) = B_M(IJK)    -  A_M(-2,IJK) * Var( JM_OF(IJK) )         &
         -  A_M(2,IJK) * Var( JP_OF(IJK) )         &
         -  A_M(-1,IJK) * Var( IM_OF(IJK) )         &
         -  A_M(1,IJK) * Var( IP_OF(IJK) )

      ENDDO

      CC(1) = ZERO
      EE(NN) = ZERO
!     DL(1:NEND-1) = CC(2:NEND)
      INFO = 0
      CALL DGTSL( NN, CC, DD, EE, BB, INFO )
!     CALL DGTSV( JEND-JSTART+1, 1, DL, DD, EE, BB,  JEND-JSTART+1, INFO )

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO K=1,NN
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(K)
      ENDDO

      RETURN
      END SUBROUTINE  LEQ_IJSWEEPt

!-----------------------------------------------
      SUBROUTINE LEQ_MSOLVE0t(VNAME, B_m, A_M, Var, CMETHOD )
!-----------------------------------------------
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE indices
      USE sendrecv

      use parallel
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!     Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3, ijkstart3:ijkend3)
!
!     Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!     Variable name
      CHARACTER(LEN=*)    Vname
!
!     Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)
      integer :: ijk

      CHARACTER(LEN=4) :: CMETHOD

!     do nothing or no preconditioning

      if (use_doloop) then

!!!$omp  parallel do private(ijk)
         do ijk=ijkstart3,ijkend3
            var(ijk) = b_m(ijk)
         enddo
      else
         var(:) = b_m(:)
      endif
      call send_recv(var,nlayers_bicgs)

      return
      end subroutine leq_msolve0t

!-----------------------------------------------
      SUBROUTINE LEQ_msolve1t(VNAME, B_m, A_M, Var, CMETHOD )
!-----------------------------------------------
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE indices
      USE sendrecv
      USE parallel
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!     Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(-3:3, ijkstart3:ijkend3)
!
!     Vector b_m
      DOUBLE PRECISION B_m(ijkstart3:ijkend3)
!
!     Variable name
      CHARACTER(LEN=*)    Vname
!
!     Variable
      DOUBLE PRECISION Var(ijkstart3:ijkend3)

      CHARACTER(LEN=4) :: CMETHOD

      integer :: i,j,k, ijk2

      if (use_doloop) then
!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 1516)
!!!$omp    parallel do private(ijk2)
         do ijk2=ijkstart3,ijkend3
            var(ijk2) = zero
         enddo
      else
         var(:) = ZERO
      endif

!     diagonal scaling

!AIKE PFUPGRADE 091409 Modified ijk to ijk2 to avoid compilation error since PF upgrade
! PGF90-S-0155-ijk may not appear in a PRIVATE clause (leq_bicgst.f: 1526)
!!$omp   parallel do private(i,j,k,ijk2) collapse (3)
      do k=kstart2,kend2
         do i=istart2,iend2
            do j=jstart2,jend2

               ijk2 = funijk( i,j,k )
               var(ijk2) = b_m(ijk2)/A_m(0,ijk2)

            enddo
         enddo
      enddo

      call send_recv(var,nlayers_bicgs)

      return
      end subroutine leq_msolve1t

