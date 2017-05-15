!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLAGS                                               C
!  Purpose: This module assigns a flag to a cell to identify its type. C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: In cylindrical geometry, set the b.c at X=0 to free-slip   C
!  Author: M. Syamlal                                 Date: 09-APR-92  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Initialize FLAG_E, FLAG_N, FLAG_T for IS specifications.   C
!           Change definition of Flags.                                C
!  Author: M. Syamlal                                 Date: 21-OCT-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: Define FLAG using info from ICBC_FLAG                      C
!  Author: M. Syamlal                                 Date: 21-APR-93  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, BC_DEFINED, BC_TYPE,     C
!                        BC_K_b, BC_K_t, BC_J_s, BC_J_n, BC_I_w,       C
!                        IS_K_b, IS_K_t, IS_J_s, IS_J_n, IS_I_w,       C
!                        BC_I_e, NO_I, NO_J, NO_K                      C
!  Variables modified: FLAG, FLAG_E, FLAG_N, FLAG_T                    C
!                                                                      C
!  Local variables: I, J, K, IJK, L, FLAGX                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLAGS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE exit, only: mfix_exit
      USE fldvar
      USE function3
      USE functions
      USE funits
      USE geometry
      USE indices
      USE is
      USE parallel
      USE param
      USE param1
      USE physprop
      USE sendrecv
      USE sendrecv3
      use mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IJK1
! Local DO loop index for b.c. specification
      INTEGER :: L
! Temporary storage for FLAG value
      INTEGER :: FLAGX
      integer, allocatable :: arr1(:)
!-----------------------------------------------

!  Cell flag definitions
!  FLAG  ICBC_FLAG BC_TYPE        Cell type
!  ----- --------- -------        ---------
!   1       .        -            Cell containing gas or solids or both
!  10       p      P_INFLOW       Specified pressure inflow cell
!  11       P      P_OUTFLOW      Specified pressure outflow cell
!  20       I      MASS_INFLOW    Specified mass flux inflow cell
!  21       O      MASS_OUTFLOW   Specified mass flux outflow cell
!  31       o      OUTFLOW        outflow cell
! 100       W      NO_SLIP_WALL   Internal/external wall with no-slip b.c.
! 101       S      FREE_SLIP_WALL Internal/external wall with free-slip
! 102       s      PAR_SLIP_WALL  Internal/external wall with partial-slip b.c.
! 106       c      CYCLIC         Cyclic b.c.
! 107       C      CYCLIC_PD      Cyclic b.c. with pressure drop
! Flag values greater than 100 are considered to be wall cells
! (see function.inc).

! make the wall cells adjacent to flow boundaries free-slip wall to
! avoid unphysical strain rates in fluid cells adjacent to the flow
! boundary
! ---------------------------------------------------------------->>>
!!$omp  parallel do private( IJK) &
!!$omp  schedule(static)
      DO i = istart4, iend4
         DO j = jstart4, jend4
            DO k = kstart4, kend4

              IJK = funijk(i, j, k)
              SELECT CASE (TRIM(ICBC_FLAG(IJK)(1:1)))
                CASE ('p', 'P', 'I', 'O', 'o')

                ijk1 = bound_funijk(i+1, j, k)
                IF(TRIM(ICBC_FLAG(IJK1)(1:1)) == 'W')ICBC_FLAG(IJK1)(1:1)='S'

                ijk1 = bound_funijk(i-1, j, k)
                IF(TRIM(ICBC_FLAG(IJK1)(1:1)) == 'W')ICBC_FLAG(IJK1)(1:1)='S'

                ijk1 = bound_funijk(i, j+1, k)
                IF(TRIM(ICBC_FLAG(IJK1)(1:1)) == 'W')ICBC_FLAG(IJK1)(1:1)='S'

                ijk1 = bound_funijk(i, j-1, k)
                IF(TRIM(ICBC_FLAG(IJK1)(1:1)) == 'W')ICBC_FLAG(IJK1)(1:1)='S'

                ijk1 = bound_funijk(i, j, k+1)
                IF(TRIM(ICBC_FLAG(IJK1)(1:1)) == 'W')ICBC_FLAG(IJK1)(1:1)='S'

                ijk1 = bound_funijk(i, j, k-1)
                IF(TRIM(ICBC_FLAG(IJK1)(1:1)) == 'W')ICBC_FLAG(IJK1)(1:1)='S'
              END SELECT
            ENDDO
          ENDDO
      ENDDO
! ----------------------------------------------------------------<<<

! Define the numerical value of the variable flag for all cells based
! on the corresponding character value of icbc_flag.  By this point the
! icbc_flag has been defined in all cells (see check_data_06 and
! check_data07 -> get_flow_bc and get_wall_bc)
! ---------------------------------------------------------------->>>
!!$omp  parallel do private( IJK) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         SELECT CASE (TRIM(ICBC_FLAG(IJK)(1:1)))
         CASE ('.')
            FLAG(IJK) = 1
         CASE ('p')
            FLAG(IJK) = 10
         CASE ('P')
            FLAG(IJK) = 11
         CASE ('I')
            FLAG(IJK) = 20
         CASE ('O')
            FLAG(IJK) = 21
         CASE ('o')
            FLAG(IJK) = 31
         CASE ('W')
            FLAG(IJK) = 100
         CASE ('S')
            FLAG(IJK) = 101
         CASE ('s')
            FLAG(IJK) = 102
         CASE ('c')
            FLAG(IJK) = 106
         CASE ('C')
            FLAG(IJK) = 107
         CASE DEFAULT

! Access to only one thread at a time
!!$omp       critical
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000) IJK, ICBC_FLAG(IJK)
            call mfix_exit(myPE)
!!$omp       end critical
         END SELECT
! ----------------------------------------------------------------<<<

! Initialize cell face flags.  UNDEFINED_I should be a large +ve value.
         FLAG_E(IJK) = UNDEFINED_I
         FLAG_N(IJK) = UNDEFINED_I
         FLAG_T(IJK) = UNDEFINED_I
      ENDDO



! Setting up flags for the higher-order implementation.
! ---------------------------------------------------------------->>>
      call send_recv(flag)

      DO i = istart3, iend3
         DO j = jstart3, jend3
            DO k = kstart3, kend3
               Flag3(funijk3(i,j,k)) = Flag(funijk(i,j,k))
            ENDDO
         ENDDO
      ENDDO

      DO i = istart4, iend4
         DO j = jstart4, jend4
            DO k = kstart4, kend4
               If(i.eq.istart4.and.istart4.ne.istart3) then
                  Flag3(funijk3(i,j,k)) = Flag3(funijk3(i+1,j,k))
               endif

               If(j.eq.jstart4.and.kstart4.ne.kstart3) then
                  Flag3(funijk3(i,j,k)) = Flag3(funijk3(i,j+1,k))
               endif

               If(k.eq.kstart4.and.kstart4.ne.kstart3) then
                  Flag3(funijk3(i,j,k)) = Flag3(funijk3(i,j,k+1))
               endif

               If(i.eq.iend4.and.iend4.ne.iend3) then
                  Flag3(funijk3(i,j,k)) = Flag3(funijk3(i-1,j,k))
               endif

               If(j.eq.jend4.and.jend4.ne.jend3) then
                  Flag3(funijk3(i,j,k)) = Flag3(funijk3(i,j-1,k))
               endif

               If(k.eq.kend4.and.kend4.ne.kend3) then
                  Flag3(funijk3(i,j,k)) = Flag3(funijk3(i,j,k-1))
               endif

            ENDDO
         ENDDO
      ENDDO

      call send_recv3(flag3)
! ----------------------------------------------------------------<<<


! Set flag_e, flag_n and flag_b to indicate any internal surfaces. If
! the flag is greater than or equal to 2000, then there is no internal
! surface.
! ---------------------------------------------------------------->>>

      DO L = 1, DIMENSION_IS
! Make sure an IS has been specified
         IF (IS_DEFINED(L)) THEN
            IF (IS_TYPE(L)=='IMPERMEABLE' .OR. &
                IS_TYPE(L)(3:13)=='IMPERMEABLE') THEN
               FLAGX = 0
            ELSEIF (IS_TYPE(L)=='SEMIPERMEABLE' .OR. &
                    IS_TYPE(L)(3:15)=='SEMIPERMEABLE') THEN
               FLAGX = 1000 + L
            ELSE
               IF(DMP_LOG)WRITE (UNIT_LOG, 1100) L
               call mfix_exit(myPE)
            ENDIF

            IF (IS_X_W(L)==IS_X_E(L) .AND. DO_I) THEN
               IS_PLANE(L) = 'E'
               I = IS_I_W(L)
               DO K = IS_K_B(L), IS_K_T(L)
                  DO J = IS_J_S(L), IS_J_N(L)
                     IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                     IJK = FUNIJK(I,J,K)
                     FLAG_E(IJK) = FLAGX
                  ENDDO
               ENDDO
            ELSEIF (IS_TYPE(L)(1:1) == 'X') THEN
               IS_PLANE(L) = 'E'
               DO I = IS_I_W(L), IS_I_E(L)
                  DO K = IS_K_B(L), IS_K_T(L)
                     DO J = IS_J_S(L), IS_J_N(L)
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        FLAG_E(IJK) = FLAGX
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

            IF (IS_Y_S(L)==IS_Y_N(L) .AND. DO_J) THEN
               IS_PLANE(L) = 'N'
               J = IS_J_S(L)
               DO K = IS_K_B(L), IS_K_T(L)
                  DO I = IS_I_W(L), IS_I_E(L)
                     IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                     IJK = FUNIJK(I,J,K)
                     FLAG_N(IJK) = FLAGX
                  ENDDO
               ENDDO
            ELSEIF (IS_TYPE(L)(1:1) == 'Y') THEN
               IS_PLANE(L) = 'N'
               DO J = IS_J_S(L), IS_J_N(L)
                  DO K = IS_K_B(L), IS_K_T(L)
                     DO I = IS_I_W(L), IS_I_E(L)
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        FLAG_N(IJK) = FLAGX
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

            IF (IS_Z_B(L)==IS_Z_T(L) .AND. DO_K) THEN
               IS_PLANE(L) = 'T'
               K = IS_K_B(L)
               DO J = IS_J_S(L), IS_J_N(L)
                  DO I = IS_I_W(L), IS_I_E(L)
                     IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                     IJK = FUNIJK(I,J,K)
                     FLAG_T(IJK) = FLAGX
                  ENDDO
               ENDDO
            ELSEIF (IS_TYPE(L)(1:1) == 'Z') THEN
               IS_PLANE(L) = 'T'
               DO K = IS_K_B(L), IS_K_T(L)
                  DO J = IS_J_S(L), IS_J_N(L)
                     DO I = IS_I_W(L), IS_I_E(L)
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        FLAG_T(IJK) = FLAGX
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ENDIF
      call send_recv(flag,2)
      call send_recv(flag_t,2)
      call send_recv(flag_n,2)
      call send_recv(flag_e,2)
      ENDDO    ! end do loop (l = 1, dimension_is)
! ----------------------------------------------------------------<<<

      IF (MYPE.EQ.PE_IO) THEN
         ALLOCATE (ARR1(IJKMAX3))
      ELSE
         ALLOCATE (ARR1(1))
      ENDIF

      CALL GATHER(FLAG,ARR1,ROOT)
      CALL SCATTER(FLAG,ARR1,ROOT)

      DEALLOCATE (ARR1)


      RETURN
 1000 FORMAT(/1X,70('*')//' From: SET_FLAGS',/&
         ' Message: ICBC_FLAG(',I3,') = ',&
         A3,' is illegal',/1X,70('*')/)
 1100 FORMAT(/1X,70('*')//' From: SET_FLAGS',/&
         ' Message: Unknown IS_TYPE(',I3,&
         ')',/1X,70('*')/)

      END SUBROUTINE SET_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_FLAGS1                                             C
!  Purpose: Assign IP flag to the faces of wall cells                  C
!                                                                      C
!  Notes: This routine may still leave flag_e, flag_n and flag_t       C
!         undefined in some boundary cells                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified: FLAG_E, FLAG_N, FLAG_T                          C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLAGS1

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE bc
      USE is
      USE indices
      USE physprop
      USE funits
      USE compar
      USE sendrecv
      USE mpi_utility
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP
      INTEGER :: I, J, K
!
      INTEGER, DIMENSION(:), allocatable :: FLAG_TEMP
      INTEGER :: flag_size
!-----------------------------------------------


! Allocate storage for temporary flag arrays
      flag_size = ijkmax3
      if (myPE.eq.root) then
          flag_size = ijkmax3
      endif
      allocate( flag_temp(flag_size) )


      DO IJK = ijkstart3,ijkend3
         IMJK = IM_OF(IJK)
         IJMK = JM_OF(IJK)
         IJKM = KM_OF(IJK)
         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

! If the flag is greater than or equal to 2000, there is no
! internal surface.
         IF (WALL_AT(IJK)) THEN
! ---------------------------------------------------------------->>>
! the default is equivalent to an impermeable surface and these cells
! will be treated as such in the momentum routines
            FLAG_E(IJK) = 0
            FLAG_N(IJK) = 0
            FLAG_T(IJK) = 0
            FLAG_E(IMJK) = 0
            FLAG_N(IJMK) = 0
            FLAG_T(IJKM) = 0

            IF (CYCLIC_AT(IJK)) THEN
! make the upper (E, N, T) boundary permeable
               IF (I == IMAX2) THEN
                  IF ((J/=1.AND.J/=0.) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF (NO_K) THEN
                        IF(.NOT.WALL_AT(IMJK)) FLAG_E(IMJK) = 2000
                     ELSEIF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(.NOT.WALL_AT(IMJK)) FLAG_E(IMJK) = 2000
                     ENDIF
                  ENDIF
               ENDIF
               IF (J == JMAX2) THEN
                  IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3)) THEN
                     IF (NO_K) THEN
                        IF(.NOT.WALL_AT(IJMK)) FLAG_N(IJMK) = 2000
                     ELSE IF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(.NOT.WALL_AT(IJMK)) FLAG_N(IJMK) = 2000
                     ENDIF
                  ENDIF
                ENDIF
               IF (K == KMAX2) THEN
                  IF ((J/=1.AND.J/=0.) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3) .AND. &
                       .NOT.WALL_AT(IJKM)) FLAG_T(IJKM) = 2000
                  ENDIF
               ENDIF

            ENDIF   ! end if cyclic_at(ijk)

! ----------------------------------------------------------------<<<
         ELSEIF (FLUID_AT(IJK)) THEN
! ---------------------------------------------------------------->>>

            IF ( .NOT.WALL_AT(IMJK) .AND. FLAG_E(IMJK)==UNDEFINED_I) &
               FLAG_E(IMJK) = 2000 + FLAG(IMJK)
            IF ( .NOT.WALL_AT(IJMK) .AND. FLAG_N(IJMK)==UNDEFINED_I) &
               FLAG_N(IJMK) = 2000 + FLAG(IJMK)
            IF ( .NOT.WALL_AT(IJKM) .AND. FLAG_T(IJKM)==UNDEFINED_I) &
               FLAG_T(IJKM) = 2000 + FLAG(IJKM)
            IF ( .NOT.WALL_AT(IPJK) .AND. FLAG_E(IJK)==UNDEFINED_I) &
               FLAG_E(IJK) = 2000 + FLAG(IPJK)
            IF ( .NOT.WALL_AT(IJPK) .AND. FLAG_N(IJK)==UNDEFINED_I) &
               FLAG_N(IJK) = 2000 + FLAG(IJPK)
            IF ( .NOT.WALL_AT(IJKP) .AND. FLAG_T(IJK)==UNDEFINED_I) &
               FLAG_T(IJK) = 2000 + FLAG(IJKP)

         ENDIF   ! end if/else (wall_at(ijk)/fluid_at(ijk))

      ENDDO    ! end do loop (ijk = ijkstart3,ijkend3)
! ----------------------------------------------------------------<<<

! Fill the ghost layers using gather and scatter
      call gather( flag_e, flag_temp )
      call scatter( flag_e, flag_temp )
      call gather( flag_n, flag_temp )
      call scatter( flag_n, flag_temp )
      call gather( flag_t, flag_temp )
      call scatter( flag_t, flag_temp )

! deallocate storage of temporary flag arrays
      deallocate( flag_temp )
      call send_recv(flag_t,2)
      call send_recv(flag_n,2)
      call send_recv(flag_e,2)

      RETURN
      END SUBROUTINE SET_FLAGS1
