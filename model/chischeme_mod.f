!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Determine factors for Chi-scheme of Darwish and Moukalled  C
!  to ensure consistency of equation sets (e.g., species mass          C
!  fractions add up to 1)                                              C
!                                                                      C
!  To initiate Chi-Scheme                                              C
!     call set_chi( ...)                                               C
!  and to terminate Chi-Scheme                                         C
!     call unset_chi(ier)                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-AUG-03   C
!                                                                      C
!                                                                      C
!  References:                                                         C
!  -Darwish, M. and Moukalled, F., "The Chi-shemes: a new consistent   C
!   high-resolution formulation based on the normalized variable       C
!   methodology," Comput. Methods Appl. Mech. Engrg., 192, 1711-1730   C
!   (2003)                                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE ChiScheme

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_e
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_n
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_t


      LOGICAL :: ChiScheme_allocated = .false.
      LOGICAL :: Chi_flag = .false.


      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE Set_Chi(DISCR, PHI, Nmax, U, V, W)

      USE compar, only: ijkstart3, ijkend3
      USE param, only: dimension_3
      USE param1, only: large_number
      USE sendrecv, only: send_recv

      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: flush_err_msg

! Dummy arguments
!---------------------------------------------------------------------//
! discretization method
      INTEGER, INTENT(IN) :: DISCR
! Second dimension of Phi array
      INTEGER, INTENT(IN) :: NMax
! convected quantity
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3, Nmax)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_e_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_n_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_t_temp
! index
      INTEGER :: IJK, N
!---------------------------------------------------------------------//

      if(.not.ChiScheme_allocated)then
         Allocate( Chi_e(DIMENSION_3) , &
                  Chi_n(DIMENSION_3) , &
                  Chi_t(DIMENSION_3)     )
         ChiScheme_allocated = .true.
      endif

      if(Chi_flag)then
! Error: Chi-Scheme is already active.  This routine cannot be called
! again before unsetting the flag
         CALL INIT_ERR_MSG("SET_CHI")
         WRITE(ERR_MSG, 1102)
 1102 FORMAT('ERROR 1102: Cannot call Set_Chi again, before Unset_chi')
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         CALL FINL_ERR_MSG
      else
! Set Chi_flag to indicate that future calls to calc_Xsi will use
! the Chi-Scheme for discretization
         Chi_e = large_number
         Chi_n = large_number
         Chi_t = large_number
         Chi_flag = .true.
      Endif

      Allocate( Chi_e_temp(DIMENSION_3) , &
                Chi_n_temp(DIMENSION_3) , &
                Chi_t_temp(DIMENSION_3)  )

! Start Chi calculations
      DO N = 1, Nmax
         CALL CALC_CHI(DISCR, PHI(1,N), U, V, W, Chi_e_temp, &
                       Chi_n_temp, Chi_t_temp)

!!!$omp    parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            Chi_e(IJK) = MIN(Chi_e(IJK), Chi_e_temp(IJK))
            Chi_n(IJK) = MIN(Chi_n(IJK), Chi_n_temp(IJK))
            Chi_t(IJK) = MIN(Chi_t(IJK), Chi_t_temp(IJK))
         ENDDO
      ENDDO

      call send_recv(CHI_E,2)
      call send_recv(CHI_N,2)
      call send_recv(CHI_T,2)

      Deallocate( Chi_e_temp , &
                  Chi_n_temp , &
                  Chi_t_temp  )


      RETURN
      END SUBROUTINE Set_Chi


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE Unset_Chi()
      IMPLICIT NONE
      Chi_flag = .false.
      RETURN
      END SUBROUTINE Unset_Chi


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CHI                                                C
!  Purpose: Determine CHI factors for higher order discretization.     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CHI(DISCR, PHI, U, V, W, CHI_E, CHI_N, CHI_T)


! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE discretization, only: chi4smart, chi4muscl, phi_c_of
      USE functions, only: wall_at
      USE functions, only: east_of, west_of, north_of, south_of
      USE functions, only: top_of, bottom_of
      USE geometry, only: do_k
      USE param, only: dimension_3
      USE param1, only: zero
      USE run, only: shear
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
      DOUBLE PRECISION, INTENT(INOUT) :: CHI_e(DIMENSION_3)
      DOUBLE PRECISION, INTENT(INOUT) :: CHI_n(DIMENSION_3)
      DOUBLE PRECISION, INTENT(INOUT) :: CHI_t(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IJKC, IJKD, IJKU
      DOUBLE PRECISION :: PHI_C
!---------------------------------------------------------------------//

! calculate CHI_E,CHI_N,CHI_T when periodic shear BCs are used
      IF (SHEAR) THEN
! this needs implementation...
! note mfix will error before this is hit
!         call CXS(incr, DISCR, U, V, W, PHI, CHI_E, CHI_N, CHI_T)
      ELSE


      SELECT CASE (DISCR)


      CASE (:1)                                  !first order upwinding
!!!$omp    parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            CHI_E(IJK) = ZERO
            CHI_N(IJK) = ZERO
            IF (DO_K) CHI_T(IJK) = ZERO
         ENDDO


!      CASE (2)                                   !Superbee


      CASE (3)                                   !SMART
!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C)
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.WALL_AT(IJK)) THEN ! no need to do these calculations for walls
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
               CHI_E(IJK) = CHI4SMART(PHI_C, PHI(IJKU), PHI(IJKC), PHI(IJKD))

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
               CHI_N(IJK) = CHI4SMART(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))

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
                  CHI_T(IJK) = CHI4SMART(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
               ENDIF
            ELSE
               CHI_E(IJK) = ZERO
               CHI_N(IJK) = ZERO
               CHI_T(IJK) = ZERO
            ENDIF  ! endif (.not. wall_at)
         ENDDO   ! end do ijk


!      CASE (4)                                   !ULTRA-QUICK


!      CASE (5)                                   !QUICKEST


      CASE (6)                                   !MUSCL

!!!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C )
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.WALL_AT(IJK)) THEN ! no need to do these calculations for walls
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
               CHI_E(IJK) = CHI4MUSCL(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))

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
               CHI_N(IJK) = CHI4MUSCL(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))

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
                  CHI_T(IJK) = CHI4MUSCL(PHI_C,PHI(IJKU),PHI(IJKC),PHI(IJKD))
               ENDIF
            ELSE
               CHI_E(IJK) = ZERO
               CHI_N(IJK) = ZERO
               CHI_T(IJK) = ZERO
            ENDIF   ! endif (.not.wall_at)
         ENDDO   ! end do ijk


!      CASE (7)                                   !Van Leer


!      CASE (8)                                   !Minmod


      CASE DEFAULT                               !Error
! should never hit this
          CALL INIT_ERR_MSG("CALC_CHI")
          WRITE(ERR_MSG, 1103) IVAL(DISCR)
 1103 FORMAT('ERROR 1103: Invalid DISCRETIZE= ', A,' when using '&
         'chi_scheme',/'The check_data routines should have already ',&
         'caught this error and ',/,'pevented the simulation from ',&
         'running. Please notify the MFIX ',/,'MFIX developers.')
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
          CALL FINL_ERR_MSG
      END SELECT

      ENDIF

      call send_recv(CHI_E,2)
      call send_recv(CHI_N,2)
      call send_recv(CHI_T,2)

      RETURN
      END SUBROUTINE CALC_CHI

      END MODULE CHIScheme
