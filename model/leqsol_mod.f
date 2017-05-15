MODULE leqsol

  use param, only: DIM_EQS
  use exit, only: mfix_exit

! Automatic adjustment of leq parameters possible (set in iterate after
! the completion of first iteration).
  LOGICAL :: LEQ_ADJUST

! Maximum number of linear equation solver iterations
  INTEGER :: LEQ_IT(DIM_EQS)

! Linear equation solver method
  INTEGER :: LEQ_METHOD(DIM_EQS)

! Total Iterations
  INTEGER :: ITER_TOT(DIM_EQS) = 0

! Linear equation solver sweep direction
  CHARACTER(LEN=4) :: LEQ_SWEEP(DIM_EQS)

! Linear equation solver tolerance
  DOUBLE PRECISION :: LEQ_TOL(DIM_EQS)

! Preconditioner option
  CHARACTER(LEN=4) :: LEQ_PC(DIM_EQS)

! Option to minimize dot products
  LOGICAL :: MINIMIZE_DOTPRODUCTS

! Option to transpose A_m
  LOGICAL :: DO_TRANSPOSE

! Frequency of convergence check in BiCGStab
  INTEGER :: ICHECK_BICGS

! Optimize for massively parallel machine
  LOGICAL :: OPT_PARALLEL

! Linear and non-linear solver statistics
  LOGICAL :: SOLVER_STATISTICS

CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
  SUBROUTINE REPORT_SOLVER_STATS(TNIT, STEPS)

    use error_manager

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: TNIT, STEPS

    INTEGER :: LC

    WRITE(ERR_MSG,1100) iVal(TNIT), iVal(TNIT/STEPS)
    CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

1100 FORMAT(/2x,'Total number of non-linear iterations: ', A,/2x,&
         'Average number per time-step: ',A)

    WRITE(ERR_MSG,1200)
    CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

1200 FORMAT(2x,'|',10('-'),'|',13('-'),'|',14('-'),'|',/&
         2x,'| Equation |  Number of  |  Avg Solves  |',/&
         2x,'|  Number  |   Solves    |   for NIT    |',/&
         2x,'|',10('-'),'|',13('-'),'|',14('-'),'|')

    DO LC = 1, DIM_EQS
       WRITE(ERR_MSG,1201) LC, ITER_TOT(LC), ITER_TOT(LC)/TNIT
       CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
    ENDDO

1201 FORMAT(2x,'|',3x,I3,4x,'|',2x,I9,2x,'|',2x,I10,2x,'|',/ &
         2x,'|',10('-'),'|',13('-'),'|',14('-'),'|')


    RETURN
  END SUBROUTINE REPORT_SOLVER_STATS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MATVEC                                              C
!  Purpose: Compute matrix vector multiplication                       C
!           (for linear equation Ax=b compute Ax                       C
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

  SUBROUTINE LEQ_MATVEC(VNAME, VAR, A_M, Avar)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    USE compar, ONLY: istart, iend, jstart, jend, kstart, kend, IJKSTART3, IJKEND3, nlayers_bicgs, c0, c1, c2, mype
    USE cutcell, ONLY: re_indexing, CARTESIAN_GRID
    USE geometry, ONLY: do_k, use_corecell_loop, CORE_ISTART, CORE_IEND, CORE_JSTART, CORE_JEND, CORE_KSTART, CORE_KEND
    USE indices
    USE param, ONLY: DIMENSION_3
    USE sendrecv, ONLY: send_recv
    IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
!      DOUBLE PRECISION, INTENT(IN) :: Var(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
    DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3)
! Vector AVar
!      DOUBLE PRECISION, INTENT(OUT) :: AVar(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(OUT) :: AVar(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
    INTEGER :: I, J, K, IJK
    integer :: im1jk, ip1jk, ijm1k, ijp1k, ijkm1, ijkp1
    integer :: class, interval
    integer :: j_start(2), j_end(2)
!-----------------------------------------------

    IF(RE_INDEXING) THEN

       DO IJK = IJKSTART3, IJKEND3  ! Loop only over active cells

          im1jk = im_of(ijk)
          ip1jk = ip_of(ijk)
          ijm1k = jm_of(ijk)
          ijp1k = jp_of(ijk)

          AVar(ijk) =      A_m(ijk,-2) * Var(ijm1k)   &
               + A_m(ijk,-1) * Var(im1jk)   &
               + A_m(ijk, 0) * Var(ijk)     &
               + A_m(ijk, 1) * Var(ip1jk)   &
               + A_m(ijk, 2) * Var(ijp1k)

          if (do_k) then
             ijkm1 = km_of(ijk)
             ijkp1 = kp_of(ijk)


             AVar(ijk) =   AVar(ijk) + A_m(ijk,-3) * Var(ijkm1)   &
                  + A_m(ijk, 3) * Var(ijkp1)

          endif

       enddo

    ELSE

          core_istart = istart+2
          core_iend = iend-2

          core_jstart = jstart+2
          core_jend = jend-2

          if (do_k) then
             core_kstart = kstart+2
             core_kend = kend-2
          else
             core_kstart = 1
             core_kend = 1
             kstart = 1
             kend = 1
          endif

          if (USE_CORECELL_LOOP) then

          class = cell_class(funijk(core_istart,core_jstart,core_kstart))

!$omp    parallel do default(none) shared(c0,c1,c2,avar,a_m,var,do_k,increment_for_mp,istart,jstart,kstart,iend,jend,kend,cell_class,core_istart,core_jstart,core_kstart,core_iend,core_jend,core_kend,use_corecell_loop,class) &
!$omp&   private(ijk,i,j,k) collapse (3)
             do k = core_kstart,core_kend
                do i = core_istart,core_iend
                   do j = core_jstart,core_jend
                      ijk = (j + c0 + i*c1 + k*c2)

                      AVar(ijk) = &
                           + A_m(ijk,-2) * Var(ijk+INCREMENT_FOR_MP(3,class))   &
                           + A_m(ijk,-1) * Var(ijk+INCREMENT_FOR_MP(1,class))   &
                           + A_m(ijk, 0) * Var(ijk)     &
                           + A_m(ijk, 1) * Var(ijk+INCREMENT_FOR_MP(2,class))   &
                           + A_m(ijk, 2) * Var(ijk+INCREMENT_FOR_MP(4,class))

                      if (do_k) then
                         AVar(ijk) =  AVar(ijk) + A_m(ijk,-3) * Var(ijk+INCREMENT_FOR_MP(5,class))
                         AVar(ijk) =  AVar(ijk) + A_m(ijk, 3) * Var(ijk+INCREMENT_FOR_MP(6,class))
                      endif
                   enddo
                enddo
             enddo
          endif

          j_start(1) = jstart
          j_end(1) = jend
          j_start(2) = 0 ! no iterations
          j_end(2) = -1  ! no iterations

!$omp    parallel do default(none) shared(c0,c1,c2,avar,a_m,var,do_k,increment_for_mp,istart,jstart,kstart,iend,jend,kend,cell_class,core_istart,core_jstart,core_kstart,core_iend,core_jend,core_kend,use_corecell_loop) &
!$omp&   private(ijk,i,j,k,class,interval) firstprivate(j_start,j_end) collapse (2)
          do k = kstart,kend
             do i = istart,iend

                if  (USE_CORECELL_LOOP) then
                   if (core_istart<= i .and. i <= core_iend .and. core_kstart <= k .and. k<=core_kend) then
                      j_start(1) = jstart
                      j_end(1) = core_jstart-1
                      j_start(2) = core_jend+1
                      j_end(2) = jend
                   else
                      j_start(1) = jstart
                      j_end(1) = jend
                      j_start(2) = 0 ! no iterations
                      j_end(2) = -1  ! no iterations
                   endif
                endif

                do interval=1,2
                   do j = j_start(interval),j_end(interval)
                      ijk = (j + c0 + i*c1 + k*c2)
                      class = cell_class(ijk)

                      AVar(ijk) = &
                           + A_m(ijk,-2) * Var(ijk+INCREMENT_FOR_MP(3,class))   &
                           + A_m(ijk,-1) * Var(ijk+INCREMENT_FOR_MP(1,class))   &
                           + A_m(ijk, 0) * Var(ijk)     &
                           + A_m(ijk, 1) * Var(ijk+INCREMENT_FOR_MP(2,class))   &
                           + A_m(ijk, 2) * Var(ijk+INCREMENT_FOR_MP(4,class))

                      if (do_k) then
                         AVar(ijk) =  AVar(ijk) + A_m(ijk,-3) * Var(ijk+INCREMENT_FOR_MP(5,class))
                         AVar(ijk) =  AVar(ijk) + A_m(ijk, 3) * Var(ijk+INCREMENT_FOR_MP(6,class))
                      endif
                   enddo
                enddo
             enddo
          enddo

       ENDIF ! RE_INDEXING

       call send_recv(Avar,nlayers_bicgs)
    RETURN

  CONTAINS

    INCLUDE 'functions.inc'

  END SUBROUTINE LEQ_MATVEC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE                                              C
!  Purpose:                                                            C
!  Notes: if leq_method is biggs or cg then this subroutine is         C
!         invoked when leq_pc='line'. if leq_method is gmres then      C
!         this subroutine is invoked (leq_pc setting does not matter)  C
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

  SUBROUTINE LEQ_MSOLVE(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
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
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Vector b_m
!      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3)
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
    DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3)
! Variable
!      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
    CHARACTER(LEN=4), INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
    LOGICAL, PARAMETER :: USE_IKLOOP = .FALSE.
    LOGICAL, PARAMETER :: SETGUESS = .TRUE.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
    INTEGER :: ITER, NITER
    INTEGER :: IJK, I , J, K
    INTEGER :: I1, J1, K1, I2, J2, K2, IK, JK, IJ
    INTEGER :: ISIZE, JSIZE, KSIZE
    INTEGER :: ICASE

!     CHARACTER(LEN=4), PARAMETER :: CMETHOD = 'II'
    CHARACTER :: CH
    LOGICAL :: DO_ISWEEP, DO_JSWEEP, DO_KSWEEP
    LOGICAL :: DO_SENDRECV, DO_REDBLACK, DO_ALL

!-----------------------------------------------
!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime
!       by Tingwen
!!$      omp_start=omp_get_wtime()

    IF (SETGUESS) THEN

! !!omp   parallel do private(i,j,k,ijk) collapse(3)
!        do k = kstart3,kend3
!           do i = istart3,iend3
!              do j = jstart3,jend3
!                 IJK = (J + C0 + I*C1 + K*C2)
!                 VAR(IJK) = B_M(IJK)
!              enddo
!           enddo
!        enddo

       VAR(:) = B_M(:)

       call send_recv(var,nlayers_bicgs)
    ENDIF

    NITER = LEN( CMETHOD )

    DO ITER=1,NITER

! Perform sweeps
       CH = CMETHOD( ITER:ITER )
       DO_ISWEEP = (CH .EQ. 'I') .OR. (CH .EQ. 'i')
       DO_JSWEEP = (CH .EQ. 'J') .OR. (CH .EQ. 'j')
       DO_KSWEEP = (CH .EQ. 'K') .OR. (CH .EQ. 'k')
       DO_ALL = (CH .EQ. 'A') .OR. (CH .EQ. 'a')
       DO_REDBLACK = (CH .EQ. 'R') .OR. (CH .EQ. 'r')
       DO_SENDRECV = (CH .EQ. 'S') .OR. (CH .EQ. 's')

       IF (NO_K) THEN   ! two dimensional
! 2D run no need to enable openmp parallel
          IF ( DO_ISWEEP ) THEN
!!$omp   parallel do private(I)
             DO I=istart,iend,1
                CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
             ENDDO
          ENDIF
! ----------------------------------------------------------------<<<
! Handan Liu added 2D RSRS sweep and parallelized this loop on Jan 22 2013:
          IF (DO_REDBLACK) THEN
!$omp parallel do private(I)
             DO I=istart,iend,2
                CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
             ENDDO
!$omp parallel do private(I)
             DO I=istart+1,iend,2
                CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
             ENDDO
          ENDIF
! ---------------------------------------------------------------->>>
       ELSE   ! three dimensional


! do_all true only for leq_pc='asas'
! ---------------------------------------------------------------->>>
          IF(DO_ALL) THEN        ! redblack for all sweeps, not used by default
! JK Loop
! --------------------------------
             j1 = jstart
             k1 = kstart
             j2 = jend
             k2 = kend
             jsize = j2-j1+1
             ksize = k2-k1+1
             DO icase = 1, 2
!!$omp   parallel do private(K,J,JK)
                DO JK=icase, ksize*jsize, 2
                   if (mod(jk,jsize).ne.0) then
                      k = int( jk/jsize ) + k1
                   else
                      k = int( jk/jsize ) + k1 -1
                   endif
                   j = (jk-1-(k-k1)*jsize) + j1 + mod(k,2)
                   if(j.gt.j2) j=j-j2 + j1 -1
                   CALL LEQ_JKSWEEP(J, K, Vname, Var, A_m, B_m)
                ENDDO

             ENDDO
             call send_recv(var,nlayers_bicgs)

! IJ Loop
! --------------------------------
             i1 = istart
             j1 = jstart
             i2 = iend
             j2 = jend
             isize = i2-i1+1
             jsize = j2-j1+1
             DO icase = 1, 2
!!$omp   parallel do private(J,I,IJ)
                DO IJ=icase, jsize*isize, 2
                   if (mod(ij,isize).ne.0) then
                      j = int( ij/isize ) + j1
                   else
                      j = int( ij/isize ) + j1 -1
                   endif
                   i = (ij-1-(j-j1)*isize) + i1 + mod(j,2)
                   if(i.gt.i2) i=i-i2 + i1 -1
                   CALL LEQ_IJSWEEP(I, J, Vname, Var, A_m, B_m)
                ENDDO

             ENDDO
             call send_recv(var,nlayers_bicgs)

! IK Loop
! --------------------------------
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
                   CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                ENDDO

             ENDDO
          ENDIF ! end DO_ALL
! ----------------------------------------------------------------<<<

! do_redblack only true leq_pc='rsrs'
! ---------------------------------------------------------------->>>
          IF(DO_REDBLACK) THEN
!i1 = istart
!k1 = kstart
!i2 = iend
!k2 = kend
!isize = i2-i1+1
!ksize = k2-k1+1
!               DO icase = 1, 2
!!$omp   parallel do private(K,I,IK)
!                  DO IK=icase, ksize*isize, 2
!                     if (mod(ik,isize).ne.0) then
!                        k = int( ik/isize ) + k1
!                     else
!                        k = int( ik/isize ) + k1 -1
!                     endif
!                     i = (ik-1-(k-k1)*isize) + i1 + mod(k,2)
!                     if(i.gt.i2) i=i-i2 + i1 -1
!                     CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
!                  ENDDO
!               ENDDO
!             ELSE
! Handan Liu split above loop for OpenMP at May 22 2013, modified at July 17
!$omp parallel do default(shared) private(I,K) schedule(auto)
             DO k=kstart,kend
                IF(mod(k,2).ne.0)THEN
                   DO I=istart+1,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ELSE
                   DO I=istart,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDIF
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared) private(I,K) schedule(auto)
             DO k=kstart,kend
                IF(mod(k,2).ne.0)THEN
                   DO I=istart,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ELSE
                   DO I=istart+1,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDIF
             ENDDO
!$omp end parallel do

          ENDIF       ! end if(do_redblack)
! ----------------------------------------------------------------<<<

!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
          IF(USE_IKLOOP) THEN
! use_ikloop is currently hard-wired to false (so goto else branch)
! ---------------------------------------------------------------->>>
             i1 = istart
             k1 = kstart
             i2 = iend
             k2 = kend
             isize = i2-i1+1
             ksize = k2-k1+1
             IF (DO_ISWEEP) THEN
!!$omp   parallel do private(K,I,IK)
                DO IK=1, ksize*isize
                   if (mod(ik,isize).ne.0) then
                      k = int( ik/isize ) + k1
                   else
                      k = int( ik/isize ) + k1 -1
                   endif
                   i = (ik-1-(k-k1)*isize) + i1
                   CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                ENDDO
             ENDIF
             IF (DO_KSWEEP) THEN
!!$omp   parallel do private(K,I,IK)
                DO IK=1, ksize*isize
                   if (mod(ik,ksize).ne.0) then
                      i = int( ik/ksize ) + i1
                   else
                      i = int( ik/ksize ) + i1 -1
                   endif
                   k = (ik-1-(i-i1)*ksize) + k1
                   CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                ENDDO
             ENDIF
! ----------------------------------------------------------------<<<
          ELSE   ! else branch of if(use_ikloop)
!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
! ---------------------------------------------------------------->>>
             IF (DO_ISWEEP) THEN
!!$omp   parallel do private(K,I)
                DO K=kstart,kend
                   DO I=istart,iend
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDDO
             ENDIF
             IF (DO_KSWEEP) THEN
!!$omp   parallel do private(K,I)
                DO I=istart,iend
                   DO K=kstart,kend
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDDO
             ENDIF
          ENDIF   ! end if/else (use(ikloop)
! ----------------------------------------------------------------<<<

       ENDIF   ! end if/else (do_k)


! this is called for all settings of leq_pc
       IF (DO_SENDRECV) call send_recv(var,nlayers_bicgs)


    ENDDO   ! end do iter=1,niter
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'leq_msolve:',omp_end - omp_start

    RETURN
  END SUBROUTINE LEQ_MSOLVE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE0                                             C
!  Notes: do nothing or no preconditioning (leq_pc='none')             C
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

  SUBROUTINE LEQ_MSOLVE0(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
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
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Vector b_m
!      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3)
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
    DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3)
! Variable
!      DOUBLE PRECISION, INTENT(OUT) :: Var(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(OUT) :: Var(DIMENSION_3)
! sweep direction
    CHARACTER(LEN=4), INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    integer :: ijk
!-----------------------------------------------

! do nothing or no preconditioning
    if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp  parallel do private(ijk)
       do ijk=ijkstart3,ijkend3
          var(ijk) = b_m(ijk)
       enddo
    else
       var(:) = b_m(:)
    endif
    call send_recv(var,nlayers_bicgs)

    return
  end subroutine leq_msolve0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE1                                             C
!  Notes: diagonal scaling (leq_pc='diag')                             C
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

  SUBROUTINE LEQ_MSOLVE1(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    USE param
    USE param1
    USE geometry
    USE compar
    USE indices
    USE sendrecv
    use parallel
!      USE cutcell, only: RE_INDEXING,INTERIOR_CELL_AT
    USE cutcell
    USE functions
    IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Vector b_m
!      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3)
! Septadiagonal matrix A_m
!      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
    DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3)
! Variable
!      DOUBLE PRECISION, INTENT(OUT) :: Var(ijkstart3:ijkend3)
    DOUBLE PRECISION, INTENT(OUT) :: Var(DIMENSION_3)
! sweep direction
    CHARACTER(LEN=4), INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    integer :: i,j,k, ijk
!-----------------------------------------------

    if (use_doloop) then   ! mfix.dat keyword default=false
!!$omp    parallel do private(ijk)
       do ijk=ijkstart3,ijkend3
          var(ijk) = zero
       enddo
    else
       var(:) = ZERO
    endif

! diagonal scaling
    IF(.NOT.RE_INDEXING) THEN
!$omp   parallel do private(i,j,k,ijk)  collapse (3)
       do k=kstart2,kend2
          do i=istart2,iend2
             do j=jstart2,jend2
                ijk = funijk( i,j,k )
                var(ijk) = b_m(ijk)/A_m(ijk,0)
             enddo
          enddo
       enddo
    ELSE
!$omp   parallel do private(ijk)  collapse (1)
       DO IJK=IJKSTART3,IJKEND3
          var(ijk) = b_m(ijk)/A_m(ijk,0)
       ENDDO
    ENDIF

    call send_recv(var,nlayers_bicgs)

    return
  end subroutine leq_msolve1


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_ISWEEP                                              C
!  Purpose: Perform line sweep at coordinate I                         C
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

      SUBROUTINE LEQ_ISWEEP(I, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
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
! Dummy arguments
!-----------------------------------------------
!  Line position
      INTEGER, INTENT(IN) :: I
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO
      INTEGER :: IJK, J, K, IM1JK, IP1JK
!-----------------------------------------------

      NEND = JEND
      NSTART = JSTART
      K = 1

      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         IM1JK = IM_OF(IJK)
         IP1JK = IP_OF(IJK)
         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IM1JK )  &
                          -  A_M(IJK, 1) * Var( IP1JK )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB, JEND-JSTART+1, INFO )

      IF (INFO.NE.0) THEN
         RETURN
      ENDIF

      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) =  BB(J)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_ISWEEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_IKSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
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

      SUBROUTINE LEQ_IKSWEEP(I, K, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
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
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, K
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(JSTART:JEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO
      INTEGER :: IJK, J, CLASS
!-----------------------------------------------

      NEND = JEND
      NSTART = JSTART

!!$omp parallel do private(j,ijk,im1jk,ip1jk,ijkm1,ijkp1)
      DO J=NSTART, NEND
!         IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))
         IJK = (J + C0 + I*C1 + K*C2)
         CLASS = CELL_CLASS(IJK)
         DD(J) = A_M(IJK,  0)
         CC(J) = A_M(IJK, -2)
         EE(J) = A_M(IJK,  2)
         BB(J) = B_M(IJK) -  A_M(IJK,-1) * Var( IJK+INCREMENT_FOR_MP(1,class) ) &
                          -  A_M(IJK, 1) * Var( IJK+INCREMENT_FOR_MP(2,class) ) &
                          -  A_M(IJK,-3) * Var( IJK+INCREMENT_FOR_MP(5,class) ) &
                          -  A_M(IJK, 3) * Var( IJK+INCREMENT_FOR_MP(6,class) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         write(*,*) 'leq_iksweep',INFO, myPE
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IKSWEEP'
         RETURN
      ENDIF

      DO J=NSTART, NEND
         Var(J + C0 + I*C1 + K*C2) = BB(J)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_IKSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_JKSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
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

      SUBROUTINE LEQ_JKSWEEP(J, K, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE funits
      USE compar
      USE indices
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: J, K
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (ISTART:IEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO, IJK, I
!-----------------------------------------------

      NEND = IEND
      NSTART = ISTART

      DO I=NSTART,NEND
         IJK = FUNIJK(I,J,K)
         DD(I) = A_M(IJK,  0)
         CC(I) = A_M(IJK, -1)
         EE(I) = A_M(IJK,  1)
         BB(I) = B_M(IJK)    -  A_M(IJK,-2) * Var( JM_OF(IJK) ) &
                             -  A_M(IJK, 2) * Var( JP_OF(IJK) ) &
                             -  A_M(IJK,-3) * Var( KM_OF(IJK) ) &
                             -  A_M(IJK, 3) * Var( KP_OF(IJK) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' JKSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO I=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(I)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_JKSWEEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_IJSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
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

      SUBROUTINE LEQ_IJSWEEP(I, J, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE funits
      USE compar
      USE indices
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, J
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(ijkstart3:ijkend3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(ijkstart3:ijkend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(ijkstart3:ijkend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (KSTART:KEND) :: CC, DD, EE, BB
      INTEGER :: NEND, NSTART, INFO, IJK, K
!-----------------------------------------------

      NEND = KEND
      NSTART = KSTART

      DO K=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         DD(K) = A_M(IJK,  0)
         CC(K) = A_M(IJK, -3)
         EE(K) = A_M(IJK,  3)
         BB(K) = B_M(IJK)    -  A_M(IJK,-2) * Var( JM_OF(IJK) ) &
                             -  A_M(IJK, 2) * Var( JP_OF(IJK) ) &
                             -  A_M(IJK,-1) * Var( IM_OF(IJK) ) &
                             -  A_M(IJK, 1) * Var( IP_OF(IJK) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IJSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO K=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(K)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_IJSWEEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  double precision function dot_product_par(r1,r2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    use compar
    use cutcell
    use functions
    use geometry
    use indices
    use parallel, only: is_serial
    use mpi_utility
    use param, only: dimension_3
    implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!      double precision, intent(in), dimension(ijkstart3:ijkend3) :: r1,r2
    double precision, intent(in), dimension(DIMENSION_3) :: r1,r2
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
    logical, parameter :: do_global_sum = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    DOUBLE PRECISION, allocatable, Dimension(:) :: r1_g, r2_g
    double precision :: prod
    integer :: i, j, k, ijk
!-----------------------------------------------

    if (numPEs.eq.1.and.is_serial) then
!!$  if (.false.) then
       dot_product_par = dot_product(r1,r2)
!!$  endif

!!$       dot_product_par = 0
!!$omp parallel do private(ijk) reduction(+:dot_product_par)
!!$       do ijk = 1, DIMENSION_3
!!$          dot_product_par = dot_product_par + r1(ijk)*r2(ijk)
!!$       enddo
       return
    endif

    if(do_global_sum) then
       prod = 0.0d0

       IF(RE_INDEXING) THEN
!         IF(.FALSE.) THEN
! Somehow, looping in this order leads to smaller time step than k,i,j nested loop below ....
          DO IJK = IJKSTART3,IJKEND3
             IF(INTERIOR_CELL_AT(IJK)) prod = prod + r1(ijk)*r2(ijk)
          ENDDO

          call global_all_sum(prod, dot_product_par)

       ELSE

!$omp parallel do private(i,j,k,ijk) reduction(+:prod)  collapse (3)
          do k = kstart1, kend1
             do i = istart1, iend1
                do j = jstart1, jend1
                   ijk = funijk_map_c (i,j,k)
                   prod = prod + r1(ijk)*r2(ijk)
                enddo
             enddo
          enddo

          call global_all_sum(prod, dot_product_par)

       ENDIF

    else
       if(myPE.eq.root) then
          allocate (r1_g(1:ijkmax3))
          allocate (r2_g(1:ijkmax3))
       else
          allocate (r1_g(10))
          allocate (r2_g(10))
       endif
       call gather(r1,r1_g)
       call gather(r2,r2_g)

       if(myPE.eq.root) then
          prod = 0.0d0

!$omp parallel do private(i,j,k,ijk) reduction(+:prod)  collapse (3)
          do k = kmin1, kmax1
             do i = imin1, imax1
                do j = jmin1, jmax1
                   ijk = funijk_gl (imap_c(i),jmap_c(j),kmap_c(k))
!                     ijk = funijk_gl (i,j,k)
                   prod = prod + r1_g(ijk)*r2_g(ijk)
                enddo
             enddo
          enddo

       endif
       call bcast( prod)

       dot_product_par = prod

       deallocate (r1_g)
       deallocate (r2_g)

    endif

  end function dot_product_par


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  function dot_product_par2(r1,r2,r3,r4)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    use mpi_utility
    use geometry
    use compar
    use indices
    use functions
    implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
    double precision, intent(in), dimension(ijkstart3:ijkend3) :: r1,r2,r3,r4
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
    logical, parameter :: do_global_sum = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    DOUBLE PRECISION, allocatable, Dimension(:,:) :: r_temp, rg_temp
    double precision, Dimension(2) :: prod, dot_product_par2
    integer :: i, j, k, ijk
!-----------------------------------------------

    if(do_global_sum) then

       prod(:) = 0.0d0

!$omp parallel do private(i,j,k,ijk) reduction(+:prod)  collapse (3)
       do k = kstart1, kend1
          do i = istart1, iend1
             do j = jstart1, jend1

                IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

                ijk = funijk_map_c (i,j,k)
                prod(1) = prod(1) + r1(ijk)*r2(ijk)
                prod(2) = prod(2) + r3(ijk)*r4(ijk)
             enddo
          enddo
       enddo

       call global_all_sum(prod, dot_product_par2)

    else
       allocate (r_temp(ijkstart3:ijkend3,4))
       r_temp(:,1) = r1
       r_temp(:,2) = r2
       r_temp(:,3) = r3
       r_temp(:,4) = r4

       if(myPE.eq.root) then
          allocate (rg_temp(1:ijkmax3,4))
       else
          allocate (rg_temp(10,4))
       endif
       call gather(r_temp,rg_temp)

       if(myPE.eq.root) then
          prod = 0.0d0
!$omp parallel do private(i,j,k,ijk) reduction(+:prod)  collapse (3)
          do k = kmin1, kmax1
             do i = imin1, imax1
                do j = jmin1, jmax1
                   IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                   ijk = funijk_gl (imap_c(i),jmap_c(j),kmap_c(k))
!                     ijk = funijk_gl (i,j,k)
                   prod(1) = prod(1) + rg_temp(ijk,1)*rg_temp(ijk,2)
                   prod(2) = prod(2) + rg_temp(ijk,3)*rg_temp(ijk,4)
                enddo
             enddo
          enddo
       endif
       call bcast( prod)

       dot_product_par2 = prod

       deallocate (r_temp)
       deallocate (rg_temp)

    endif

  end function dot_product_par2

END MODULE leqsol
