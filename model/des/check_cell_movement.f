!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT                                     !
!                                                                      !
!  Purpose: Check to see if particles have moved into ghost cells.     !
!                                                                      !
!  Note: This routine needs a global communicator to identify errors.  !
!  The collection could get expensive so the call frequency of this    !
!  routine should probably be reduced.                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT

! Global Variables:
!---------------------------------------------------------------------//
! Max number of particles in process
      use discretelement, only: MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! Number of particles in fluid cells
      use discretelement, only: PINC
! Run time flag indicating DEM or PIC solids.
      use run, only: DEM_SOLIDS, PIC_SOLIDS
! Flag for fluid cells
      use functions, only: FLUID_AT

      use functions, only: IS_NORMAL
      use mpi_utility
      use error_manager

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:
      INTEGER :: L, I, J, K, IJK
! Integer error flag.
      INTEGER :: IER
! Flag to remove rogue DEM particles
      LOGICAL, PARAMETER :: REMOVE_ROGUE_PARTICLES = .FALSE.

! Initialize local variables.
      IER = 0

! Set an error flag if any errors are found. Preform a global collection
! to sync error flags. If needed, reort errors.
!.......................................................................
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(IS_NORMAL(L)) THEN

            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
            IJK = PIJK(L,4)

            IF((I > IEND1 .OR. I < ISTART1) .OR. &
               (J > JEND1 .OR. J < JSTART1) .OR. &
               (K > KEND1 .OR. K < KSTART1)) THEN

! This particle cannot say in the current cell.
               PINC(IJK) = PINC(IJK) - 1
! PIC parcels can be relocated.
               IF(PIC_SOLIDS) THEN
                  CALL RECOVER_PARCEL(L)
! DEM particles can be removed.
               ELSEIF(REMOVE_ROGUE_PARTICLES) THEN
                  CALL DELETE_PARTICLE(L)
! Otherwise, this is a fatal error.
               ELSE
                  IER = 1
               ENDIF

            ENDIF
         ENDIF
      ENDDO

      IF(DEM_SOLIDS) THEN
         IF(REMOVE_ROGUE_PARTICLES) RETURN
         CALL GLOBAL_ALL_SUM(IER)
         IF(IER > 0) CALL CHECK_CELL_MOVEMENT_DEM
!     ELSE
!        CALL CHECK_CELL_MOVEMENT_PIC
      ENDIF


      RETURN
      END SUBROUTINE CHECK_CELL_MOVEMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT_DEM                                 !
!                                                                      !
!  Purpose: Report which DEM particles have moved into ghost cells.    !
!  This is a dead-end routine. Once called, the simulation will exit.  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_DEM

! Global Variables:
!---------------------------------------------------------------------//
! The global ID of a particle.
      use discretelement, only: iGlobal_ID
! Max number of particles in process
      use discretelement, only: MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! Particle positions and velocities
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW

      use functions, only: IS_NORMAL
      use mpi_utility
      USE error_manager

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:.
      INTEGER :: L, I, J, K
! Integer error flag
      INTEGER :: IER


      CALL INIT_ERR_MSG("CHECK_CELL_MOVEMENT_DEM")
      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Particles detected in a ghost cell:',/' ')

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'I',        &
               trim(iVal(I)),'X',DES_POS_NEW(L,1),'X',DES_VEL_NEW(L,1)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_id(L))),'J',        &
               trim(iVal(J)),'Y',DES_POS_NEW(L,2),'Y',DES_VEL_NEW(L,2)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF (DO_K .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'K',        &
               trim(iVal(K)),'Z',DES_POS_NEW(L,3),'Z',DES_VEL_NEW(L,3)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

 1101 FORMAT('Particle ',A,' moved into cell with ',A,' index ',A,/    &
         3x,A,'-Position: ',g11.4,6x,A,'-Velocity:',g11.4,/' ')

      WRITE(ERR_MSG, 1102)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
 1102 FORMAT('This is a fatal error. A particle output file (vtp) ',   &
         'will be written',/'to aid debugging.')


      CALL WRITE_DES_DATA
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE CHECK_CELL_MOVEMENT_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: RECOVER_PARCEL                                          !
!                                                                      !
!  Purpose: Try to recover the parcel. First move it back to its       !
!  previous position, then do a full search to bin the parcel in a     !
!  fluid cell. Delete the parcel if it still is located outside of a   !
!  fluid cell.                                                         !
!                                                                      !
!  Possible Future Work:                                               !
!                                                                      !
!  1. Redistribute particle's weight among other particles in the      !
!     domain to conserve mass.                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE RECOVER_PARCEL(NP)

! Global Variables:
!---------------------------------------------------------------------//
! The (max) number of particles in process
      use discretelement, only: PIP, MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! The number of particles in cell IJK
      use discretelement, only: PINC
! Particle position and velocity
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW
! The East/North/Top face loctions of fluid cells
      use discretelement, only: XE, YN, ZT
! Particle velocities
      use discretelement, only: DTSOLID
! Fluid grid cell dimensions and mesh size
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: JMIN2, JMAX2
      USE geometry, only: KMIN2, KMAX2
! The East/North/Top face location of a given I/J/K index.
      use discretelement, only: XE, YN, ZT
! Fixed array sizes in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K


      use cutcell
      use error_manager
      use functions, only: funijk, fluid_at
      use mpi_utility
      use discretelement, only: iglobal_id

      IMPLICIT NONE


! Dummy Arguments
!----------------------------------------------------------------------!
! Parcel ID
      INTEGER, INTENT(IN) :: NP


! Local Variables:
!----------------------------------------------------------------------!
! particle no.
      INTEGER :: I, J, K, IJK
! Integer error flag.
      INTEGER :: IER
! Local parameter to print verbose messages about particles.
      LOGICAL, PARAMETER :: lDEBUG = .FALSE.

      DOUBLE PRECISION :: oPOS(3)
!.......................................................................

      CALL INIT_ERR_MSG("RECOVER_PARCEL")
      IF(lDEBUG) CALL OPEN_PE_LOG(IER)

      oPOS = DES_POS_NEW(NP,:)

! Reflect the parcle.
      DES_VEL_NEW(NP,:) = -DES_VEL_NEW(NP,:)

! Move the particle back to the previous position.
      DES_POS_NEW(NP,:) = DES_POS_NEW(NP,:) + &
         DES_VEL_NEW(NP,:) * DTSOLID

! Rebin the particle to the fluid grid.
      CALL PIC_SEARCH(I,DES_POS_NEW(NP,1),XE,DIMENSION_I,IMIN2,IMAX2)
      CALL PIC_SEARCH(J,DES_POS_NEW(NP,2),YN,DIMENSION_J,JMIN2,JMAX2)
      CALL PIC_SEARCH(K,DES_POS_NEW(NP,3),ZT,DIMENSION_K,KMIN2,KMAX2)

! Calculate the fluid cell index.
      IJK = FUNIJK(I,J,K)
      IF(FLUID_AT(IJK)) THEN

! Assign PIJK(L,1:4)
         PIJK(NP,1) = I
         PIJK(NP,2) = J
         PIJK(NP,3) = K
         PIJK(NP,4) = IJK

         PINC(IJK) = PINC(IJK) + 1
      ELSE

         write(*,*) 'Still not cool -->', iGLOBAL_ID(NP)

!         write(*,*) 'POS:',oPOS
!         call write_des_data
!         stop 'killer stop'
         CALL DELETE_PARTICLE(NP)

      ENDIF


 1100 FORMAT('Warninge 1100: Particle ',A,' was recovered from a ',    &
         'ghost cell.',2/2x,'Moved into cell with ',A1,' index: ',A,   &
         /2x,A1,'-Position OLD:',g11.4,/2x,A1,'-Position NEW:',g11.4,  &
         /2x,A1,'-Velocity:',g11.4)

 1110 FORMAT('Warninge 1110: Particle ',A,' was deleted from a ',      &
         'ghost cell.',2/2x,'Moved into cell with ',A1,' index: ',A,   &
         /2x,A1,'-Position OLD:',g11.4,/2x,A1,'-Position NEW:',g11.4,  &
         /2X,A1,'-Velocity:',g11.4,/2x,'Fluid Cell: ',A,/2x,           &
         'Cut cell? ',L1,/2x,'Fluid at? ',L1)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE RECOVER_PARCEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT_PIC                                 !
!                                                                      !
!  Purpose: Report which PIC particles have moved into ghost cells.    !
!  Unlike DEM particles, this routine either deletes particles that    !
!  are out of the domain, or reassigns the index to try and recover    !
!  the particle.                                                       !
!                                                                      !
!  Notes:                                                              !
!                                                                      !
!  PIC particles may end up in ghost cells if the cell adjacent to the !
!  ghost cell is a cut-cell the particle was not detected outside the  !
!  system because of tolerances.                                       !
!                                                                      !
!  Future Work:                                                        !
!                                                                      !
!  1. Redistribute particle's weight among other particles in the      !
!     domain to conserve mass.                                         !
!                                                                      !
!  2. Rather than deactivating the particle, reflect the particle      !
!     inside the domain using the ghost cell bc's                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_PIC

! Global Variables:
!---------------------------------------------------------------------//
! The (max) number of particles in process
      use discretelement, only: PIP, MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! The number of particles in cell IJK
      use discretelement, only: PINC
! Particle positions New/Previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! The East/North/Top face loctions of fluid cells
      use discretelement, only: XE, YN, ZT
! Particle velocities
      use discretelement, only: DES_VEL_NEW
! Flag: identifies a fluid cell as a cut cell.
      use cutcell, only: CUT_CELL_AT

      use mpi_utility
      use error_manager
      use functions

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! particle no.
      INTEGER :: L, I, J, K, IJK
! Integer error flag.
      INTEGER :: IER
! Position
      DOUBLE PRECISION :: lPOS
! Number of deleted particles found on the local process
      INTEGER :: lDELETED, gDELETED
! Number of recovered particles found on the local process
      INTEGER :: lRECOVERED, gRECOVERED
! Local parameter to print verbose messages about particles.
      LOGICAL, PARAMETER :: lDEBUG = .FALSE.

!.......................................................................


      CALL INIT_ERR_MSG("CHECK_CELL_MOVEMENT_PIC")
      IF(lDEBUG) CALL OPEN_PE_LOG(IER)


! Initialize the counters for deleted/recovered parcels.
      lDELETED = 0
      lRECOVERED = 0

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

! this ijk is still an old value as it has not been updated
         IJK=PIJK(L,4)

! in MPPIC a particle can lie on the surface of the wall as only the
! centers are tracked.
         IF (I > IEND1 .OR. I < ISTART1) THEN

            lPOS = DES_POS_NEW(L,1)
            IF(I.EQ.IEND1+1 .AND. &
               (lPOS >= XE(IEND1-1) .AND. lPOS <= XE(IEND1)) )THEN

               lRECOVERED = lRECOVERED + 1
               PIJK(L,1) = IEND1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1100) trim(iVal(L)),'I',trim(iVal(I)), &
                  'X',DES_POS_OLD(L,1),'X',lPOS,'X',DES_VEL_NEW(L,1)
                  CALL FLUSH_ERR_MSG
               ENDIF
            ELSE

               lDELETED = lDELETED + 1
               CALL SET_NONEXISTENT(L)
               PINC(IJK) = PINC(IJK) - 1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1110) trim(iVal(L)),'I',trim(iVal(I)), &
                  'X',DES_POS_OLD(L,1),'X',lPOS,'X',DES_VEL_NEW(L,1),  &
                  trim(iVal(IJK)), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  CALL FLUSH_ERR_MSG
               ENDIF
               CYCLE
            ENDIF
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            lPOS = DES_POS_NEW(L,2)
            IF(J.EQ.JEND1+1.AND.&
              (lPOS >= YN(JEND1-1) .AND. lPOS <= YN(JEND1)) ) THEN

               lRECOVERED = lRECOVERED + 1
               PIJK(L,2) = JEND1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1100) trim(iVal(L)),'J',trim(iVal(J)), &
                  'Y',DES_POS_OLD(L,2),'Y',lPOS,'Y',DES_VEL_NEW(L,2)
                  CALL FLUSH_ERR_MSG
               ENDIF

            ELSE

               lDELETED = lDELETED + 1
               CALL SET_NONEXISTENT(L)
               PINC(IJK) = PINC(IJK) - 1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1110) trim(iVal(L)),'J',trim(iVal(J)), &
                  'Y',DES_POS_OLD(L,2),'Y',lPOS,'Y',DES_VEL_NEW(L,2),  &
                  trim(iVal(IJK)), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  CALL FLUSH_ERR_MSG
               ENDIF
               CYCLE
            ENDIF
         ENDIF

         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) THEN
            lPOS = DES_POS_NEW(L,3)
            IF(K == KEND1+1 .AND. &
              (lPOS >= ZT(KEND1-1) .AND. lPOS <= ZT(KEND1)) ) THEN

               lRECOVERED = lRECOVERED + 1
               PIJK(L,3) = KEND1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1100) trim(iVal(L)),'K',trim(iVal(K)), &
                  'Z',DES_POS_OLD(L,3),'Z',lPOS,'Z',DES_VEL_NEW(L,3)
                  CALL FLUSH_ERR_MSG
               ENDIF
            ELSE

               lDELETED = lDELETED + 1
               CALL SET_NONEXISTENT(L)
               PINC(IJK) = PINC(IJK) - 1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1110) trim(iVal(L)),'K',trim(iVal(K)), &
                  'Z',DES_POS_OLD(L,3),'Z',lPOS,'Z',DES_VEL_NEW(L,3),  &
                  trim(iVal(IJK)), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  CALL FLUSH_ERR_MSG
               ENDIF
               CYCLE
            ENDIF
         ENDIF
      ENDDO

 1100 FORMAT('Warninge 1100: Particle ',A,' was recovered from a ',    &
         'ghost cell.',2/2x,'Moved into cell with ',A1,' index: ',A,   &
         /2x,A1,'-Position OLD:',g11.4,/2x,A1,'-Position NEW:',g11.4,  &
         /2x,A1,'-Velocity:',g11.4)

 1110 FORMAT('Warninge 1110: Particle ',A,' was deleted from a ',      &
         'ghost cell.',2/2x,'Moved into cell with ',A1,' index: ',A,   &
         /2x,A1,'-Position OLD:',g11.4,/2x,A1,'-Position NEW:',g11.4,  &
         /2X,A1,'-Velocity:',g11.4,/2x,'Fluid Cell: ',A,/2x,           &
         'Cut cell? ',L1,/2x,'Fluid at? ',L1)

! Update the number of particles
      PIP = PIP - lDELETED

! Send the numbers to the IO process.
      CALL GLOBAL_SUM(lRECOVERED, gRECOVERED, PE_IO)
      CALL GLOBAL_SUM(lDELETED, gDELETED, PE_IO)

      IF(gRECOVERED + gDELETED > 0) THEN
         WRITE(ERR_MSG,1115) trim(iVal(gDELETED + gRECOVERED)),        &
            trim(iVal(gDELETED)), trim(iVal(gRECOVERED))
         CALL FLUSH_ERR_MSG
      ENDIF

 1115 FORMAT('Warning 1115: ',A,' particles detected outside the ',    &
         'domain.',/2x,A,' particles were deleted.',/2x,A,' particles',&
         ' were recovered.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_CELL_MOVEMENT_PIC
