!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: PARTICLES_IN_CELL                                       !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PARTICLES_IN_CELL

! The number of particles on the current process.
      use discretelement, only: PIP, MAX_PIP
! The I/J/K, IJK, and phase index of each particle
      use discretelement, only: PIJK
! The number and list of particles in each fluid cell IJK.
      use derived_types, only: PIC
      use discretelement, only: PINC
! The East/North/Top face location of a given I/J/K index.
      use discretelement, only: XE, YN, ZT
! Flag for 2D simulations.
      use geometry, only: NO_K
! The start and end indices of IJK loops
      use compar, only: IJKStart3, IJKEnd3
! The Upper and Loper indices covered by the current process.
      use compar, only: ISTART3, IEND3
      use compar, only: JSTART3, JEND3
      use compar, only: KSTART3, KEND3
! Fluid grid cell dimensions and mesh size
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: JMIN2, JMAX2
      USE geometry, only: KMIN2, KMAX2
! Fixed array sizes in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      use param, only: DIMENSION_3
! Function to conpute IJK from I/J/K
      use functions, only: FUNIJK

      !     The accumulated number of particles in each IJK.
      use generate_particles, only: particle_count

      use discretelement, only: DES_POS_NEW
      use functions, only: IS_NONEXISTENT, IS_GHOST, IS_ENTERING_GHOST, IS_EXITING_GHOST

      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! ijk indices
      INTEGER I, J, K, IJK
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos
!......................................................................!

! following quantities are reset every call to particles_in_cell
      PINC(:) = 0

!      allocate(PARTICLE_COUNT(DIMENSION_3))
! Use an incremental approach to determine the new particle location.
!-----------------------------------------------------------------------
!!$omp parallel default(shared) private(L, I, J, K, IJK)
!!$omp do reduction(+:PINC) schedule (guided,50)

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(IS_NONEXISTENT(L)) CYCLE

         I = PIJK(L,1)
         IF(I <= ISTART3 .OR. I >= IEND3) THEN
            CALL PIC_SEARCH(I, DES_POS_NEW(L,1), XE,                   &
               DIMENSION_I, IMIN2, IMAX2)
         ELSE
            IF((DES_POS_NEW(L,1) >= XE(I-1)) .AND.                     &
               (DES_POS_NEW(L,1) <  XE(I))) THEN
               I = I
            ELSEIF((DES_POS_NEW(L,1) >= XE(I)) .AND.                   &
               (DES_POS_NEW(L,1) < XE(I+1))) THEN
              I = I+1
            ELSEIF((DES_POS_NEW(L,1) >= XE(I-2)) .AND.                 &
               (DES_POS_NEW(L,1) < XE(I-1))) THEN
               I = I-1
            ELSE
               CALL PIC_SEARCH(I, DES_POS_NEW(L,1), XE,                &
                  DIMENSION_I, IMIN2, IMAX2)
            ENDIF
         ENDIF

         J = PIJK(L,2)
         IF(J <= JSTART3 .OR. J >= JEND3) THEN
            CALL PIC_SEARCH(J, DES_POS_NEW(L,2), YN,                   &
               DIMENSION_J, JMIN2, JMAX2)
         ELSE
            IF((DES_POS_NEW(L,2) >= YN(J-1)) .AND.                     &
               (DES_POS_NEW(L,2) < YN(J))) THEN
               J = J
            ELSEIF((DES_POS_NEW(L,2) >= YN(J)) .AND.                   &
               (DES_POS_NEW(L,2) < YN(J+1))) THEN
               J = J+1
            ELSEIF((DES_POS_NEW(L,2) >= YN(J-2)) .AND.                 &
               (DES_POS_NEW(L,2) < YN(J-1)))THEN
               J = J-1
            ELSE
               CALL PIC_SEARCH(J, DES_POS_NEW(L,2), YN,                &
                  DIMENSION_J, JMIN2, JMAX2)
            ENDIF
         ENDIF


         IF(NO_K) THEN
            K = 1
         ELSE
            K = PIJK(L,3)
            IF(K <= KSTART3 .OR. K >= KEND3) THEN
               CALL PIC_SEARCH(K, DES_POS_NEW(L,3), ZT,                &
                  DIMENSION_K, KMIN2, KMAX2)
            ELSE
               IF((DES_POS_NEW(L,3) >= ZT(K-1)) .AND.                  &
                  (DES_POS_NEW(L,3) < ZT(K))) THEN
                  K = K
                ELSEIF((DES_POS_NEW(L,3) >= ZT(K)) .AND.               &
                  (DES_POS_NEW(L,3) < ZT(K+1))) THEN
                  K = K+1
               ELSEIF((DES_POS_NEW(L,3) >= ZT(K-2)) .AND.              &
                  (DES_POS_NEW(L,3) < ZT(K-1))) THEN
                  K = K-1
               ELSE
                  CALL PIC_SEARCH(K, DES_POS_NEW(L,3), ZT,             &
                     DIMENSION_K, KMIN2, KMAX2)
               ENDIF
            ENDIF
         ENDIF

! Calculate the fluid cell index.
         IJK = FUNIJK(I,J,K)

! Assign PIJK(L,1:4)
         PIJK(L,1) = I
         PIJK(L,2) = J
         PIJK(L,3) = K
         PIJK(L,4) = IJK

! Increment the number of particles in cell IJK
         IF(.NOT.IS_GHOST(L) .AND. .NOT.IS_ENTERING_GHOST(L) .AND. &
            .NOT.IS_EXITING_GHOST(L)) PINC(IJK) = PINC(IJK) + 1

      ENDDO
!!$omp end parallel

      CALL CHECK_CELL_MOVEMENT

! Assigning the variable PIC(IJK)%p(:). For each computational fluid
! cell compare the number of current particles in the cell to what was
! in the cell previously. If different reallocate. Store the particle
! ids
! ---------------------------------------------------------------->>>
!!$omp parallel do if(ijkend3 .ge. 2000) default(shared)           &
!!$omp private(ijk,npic) !schedule (guided,50)
      DO IJK = IJKSTART3, IJKEND3

! checking all cells (including ghost cells); updating entering/exiting
! particle regions
         NPIC =  PINC(IJK)
         IF (ASSOCIATED(PIC(IJK)%p)) THEN
            IF (NPIC.NE.SIZE(PIC(IJK)%p)) THEN
               DEALLOCATE(PIC(IJK)%p)
               IF (NPIC.GT.0) ALLOCATE(PIC(IJK)%p(NPIC))
            ENDIF
         ELSE
            IF (NPIC.GT.0) ALLOCATE(PIC(IJK)%p(NPIC))
         ENDIF
      ENDDO
!!$omp end parallel do

      PARTICLE_COUNT(:) = 1
      PC = 1
      DO L = 1, MAX_PIP
! exiting loop if reached max number of particles in processor
         IF(PC.GT.PIP) exit
! skipping indices with no particles (non-existent particles)
         IF(IS_NONEXISTENT(L)) CYCLE
! incrementing particle account when particle exists
         PC = PC+1
! skipping ghost particles
         IF(IS_GHOST(L) .OR. IS_ENTERING_GHOST(L) .OR. IS_EXITING_GHOST(L)) CYCLE
         IJK = PIJK(L,4)
         POS = PARTICLE_COUNT(IJK)
         PIC(IJK)%P(POS) = L
         PARTICLE_COUNT(IJK) = PARTICLE_COUNT(IJK) + 1
      ENDDO

!      deallocate(PARTICLE_COUNT)

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_PARTICLES_IN_CELL                                  !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_PARTICLES_IN_CELL

      use discretelement, only: PIJK, PINC
      USE discretelement, only: DES_POS_NEW
      USE discretelement, only: MAX_PIP
      USE discretelement, only: XE, YN, ZT
      USE functions, only: IS_NONEXISTENT, IS_GHOST, IS_ENTERING_GHOST, IS_EXITING_GHOST
      use mpi_funs_des, only: des_par_exchange

! Number of particles in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K

      use mpi_utility
      use sendrecv

      USE error_manager
      USE desgrid, only: desgrid_pic

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER :: L
! ijk indices
      INTEGER :: I, J, K, IJK

      CALL INIT_ERR_MSG("INIT_PARTICLES_IN_CELL")

! following quantities are reset every call to particles_in_cell
      PINC(:) = 0

! Bin the particles to the DES grid.
      CALL DESGRID_PIC(.TRUE.)
! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
      CALL DES_PAR_EXCHANGE

! Assigning PIJK(L,1), PIJK(L,2) and PIJK(L,3) the i, j, k indices
! of particle L (locating particle on fluid grid). Also determine
! composite ijk index. If first_pass, also assigning PIJK(L,5) the
! solids phase index of particle.
! ---------------------------------------------------------------->>>
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(IS_NONEXISTENT(L)) CYCLE

! Use a brute force technique to determine the particle locations in
! the Eulerian fluid grid.

         CALL PIC_SEARCH(I, DES_POS_NEW(L,1), XE,                      &
            DIMENSION_I, IMIN2, IMAX2)
         PIJK(L,1) = I

         CALL PIC_SEARCH(J, DES_POS_NEW(L,2), YN,                      &
            DIMENSION_J, JMIN2, JMAX2)
         PIJK(L,2) = J

         IF(NO_K) THEN
            K=1
            PIJK(L,3) = 1
         ELSE
            CALL PIC_SEARCH(K, DES_POS_NEW(L,3), ZT,                   &
               DIMENSION_K, KMIN2, KMAX2)
            PIJK(L,3) = K
         ENDIF

! Assigning PIJK(L,4) now that particles have been located on the fluid
         IJK = FUNIJK(I,J,K)
         PIJK(L,4) = IJK

! Enumerate the number of 'real' particles in the ghost cell.
         IF(.NOT.IS_GHOST(L) .AND. .NOT.IS_ENTERING_GHOST(L) .AND. &
            .NOT.IS_EXITING_GHOST(L)) PINC(IJK) = PINC(IJK) + 1
      ENDDO

! Bin the particles to the DES grid.
      CALL DESGRID_PIC(.TRUE.)
! Calling exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
! unclear why this needs to be called again.
      CALL DES_PAR_EXCHANGE

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE INIT_PARTICLES_IN_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PIC_SEARCH                                              !
!                                                                      !
!  Purpose: Identify the I (or J or K) index of the fluid cell that    !
!  contains the particle centroid.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_SEARCH(IDX, lPOS, ENT_POS, lDIMN, lSTART, lEND)

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index being searched for (I, J, or K)
      INTEGER, INTENT(OUT) :: IDX
! Particle x,y,z position
      DOUBLE PRECISION, INTENT(IN) :: lPOS
! Dimension of ENT_POS array
      INTEGER, INTENT(IN) :: lDIMN
! East, North, or Top cell face location
      DOUBLE PRECISION, INTENT(IN) :: ENT_POS(0:lDIMN)
! Search bounds (by rank)
      INTEGER, INTENT(IN) :: lSTART, lEND

      DO IDX = lSTART,lEND
         IF(lPOS >= ENT_POS(IDX-1) .AND. lPOS < ENT_POS(IDX)) EXIT
      ENDDO

      RETURN
      END SUBROUTINE PIC_SEARCH
