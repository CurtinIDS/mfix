!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_MARCH                                       !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_MARCH

! Modules
!---------------------------------------------------------------------//
      use calc_collision_wall, only: calc_dem_thermo_with_wall_stl
      use des_bc, only: DEM_BCMI, DEM_BCMO
      use des_thermo, only: CALC_RADT_DES
      use desgrid, only: desgrid_pic
      use discretelement
      use error_manager
      use functions
      use machine
      use mpi_funs_des, only: DESMPI_SEND_RECV_FIELD_VARS
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      use mpi_utility
      use output_man, only: OUTPUT_MANAGER
      use run, only: ANY_SPECIES_EQ
      use run, only: CALL_USR
      use run, only: ENERGY_EQ
      use run, only: NSTEP
      use run, only: TIME, TSTOP, DT
      use sendrecv
! Include mfix.dat input vars that set printing frequency
      use particle_filter, only: DT_DRAG_PRINT, &
                & DT_CONTACT_PRINT, DT_VEL_PRINT
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Total number of particles
      INTEGER, SAVE :: NP=0

      !type(sap_t) :: sap

! time step loop counter index
      INTEGER :: NN,ii,nnn
! loop counter index for any initial particle settling incoupled cases
      INTEGER :: FACTOR
! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      DOUBLE PRECISION :: DTSOLID_TMP
! Numbers to calculate wall time spent in DEM calculations.
      DOUBLE PRECISION :: TMP_WALL

!......................................................................!

! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME
      DTSOLID_TMP = DTSOLID
      TMP_WALL = WALL_TIME()

! Initialize time stepping variables for coupled gas/solids simulations.
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID = DT
         ENDIF

! Initialize time stepping variable for pure granular simulations.
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
      ENDIF   ! end if/else (des_continuum_coupled)

      NP = PIP - IGHOST_CNT
      CALL GLOBAL_ALL_SUM(NP)

      IF(DES_CONTINUUM_COUPLED) THEN
         WRITE(ERR_MSG, 1000) trim(iVal(factor)), trim(iVAL(NP))
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ELSE
         WRITE(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(factor))
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ENDIF
 1000 FORMAT(/'DEM NITs: ',A,3x,'Total PIP: ', A)
 1100 FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

      IF(CALL_USR) CALL USR0_DES

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_EXPLICITLY_COUPLED) THEN
            CALL DRAG_GS_DES1
            IF(ENERGY_EQ) CALL CONV_GS_DES1
            IF(ANY_SPECIES_EQ) CALL DES_REACTION_MODEL
         ELSE
            IF(ANY_SPECIES_EQ) CALL ZERO_RRATE_DES
            IF(ENERGY_EQ) CALL ZERO_ENERGY_SOURCE
         ENDIF
         CALL CALC_PG_GRAD
      ENDIF

      IF(any(CALC_RADT_DES)) CALL CALC_avgTs


! Main DEM time loop
!----------------------------------------------------------------->>>
      DO NN = 1, FACTOR

         IF(DES_CONTINUUM_COUPLED) THEN
! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the discrete loop
            IF(S_TIME.GE.(TIME+DT)) EXIT
! If next time step in the discrete loop will exceed the current time
! in the continuum simulation, modify the discrete time step so final
! time will match
            IF((S_TIME+DTSOLID).GT.(TIME+DT)) &
               DTSOLID = TIME + DT - S_TIME
         ENDIF

! Calculate inter particle forces acting (collisional, cohesion)
         CALL CALC_FORCE_DEM
! Calculate or distribute fluid-particle drag force.
         CALL CALC_DRAG_DES
! Calculate heat conduction to/from wall
         IF(ENERGY_EQ)CALL CALC_DEM_THERMO_WITH_WALL_STL

! Update the old values of particle position and velocity with the new
! values computed
         IF (DO_OLD) CALL CFUPDATEOLD
! Calculate thermochemical sources (energy and  rates of formation).
         CALL CALC_THERMO_DES
! Call user functions.

         IF(CALL_USR) CALL USR1_DES

         ! sap = multisap%saps(0)
         ! ! CHECK SORT
         ! do ii=2, sap%x_endpoints_len
         !    if (sap%x_endpoints(ii)%value < sap%x_endpoints(ii-1)%value) then
         !       print *,"****************************************************************************************"
         !       print *,"ii:",ii,"  endpoints(ii):",sap%x_endpoints(ii)%box_id,sap%x_endpoints(ii)%value
         !       print *,"****************************************************************************************"
         !       stop __LINE__
         !    endif
         ! enddo

         ! do nnn=0, size(multisap%saps)-1
         !    !print *,"nnn = ",nnn
         !    if (.not.check_boxes(multisap%saps(nnn))) stop __LINE__
         !    if (.not.check_sort(multisap%saps(nnn))) stop __LINE__
         ! enddo

! Update position and velocities
         CALL CFNEWVALUES

         ! do nnn=0, size(multisap%saps)-1
         !    !print *,"nnn = ",nnn
         !    if (.not.check_boxes(multisap%saps(nnn))) stop __LINE__
         !    if (.not.check_sort(multisap%saps(nnn))) stop __LINE__
         ! enddo

         ! sap = multisap%saps(0)
         ! ! CHECK SORT
         ! do ii=2, sap%x_endpoints_len
         !    if (sap%x_endpoints(ii)%value < sap%x_endpoints(ii-1)%value) then
         !       print *,"****************************************************************************************"
         !       print *,"ii:",ii,"  endpoints(ii):",sap%x_endpoints(ii)%box_id,sap%x_endpoints(ii)%value
         !       print *,"****************************************************************************************"
         !       stop __LINE__
         !    endif
         ! enddo


! Update particle temperatures
         CALL DES_THERMO_NEWVALUES

! Set DO_NSEARCH before calling DES_PAR_EXCHANGE.
         DO_NSEARCH = (NN == 1 .OR. MOD(NN,NEIGHBOR_SEARCH_N) == 0)

! Add/Remove particles to the system via flow BCs.
         IF(DEM_BCMI > 0) CALL MASS_INFLOW_DEM
         IF(DEM_BCMO > 0) CALL MASS_OUTFLOW_DEM(DO_NSEARCH)

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
         IF (DO_NSEARCH .OR. (numPEs>1) .OR. DES_PERIODIC_WALLS) THEN
            CALL DESGRID_PIC(.TRUE.)
            CALL DES_PAR_EXCHANGE
         ENDIF

         IF(DO_NSEARCH) CALL NEIGHBOUR

! Explicitly coupled simulations do not need to rebin particles to
! the fluid grid every time step. However, this implies that the
! fluid cell information and interpolation weights become stale.
         IF(DES_CONTINUUM_COUPLED .AND. &
            .NOT.DES_EXPLICITLY_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS
         ENDIF

! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME and number of steps for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES and SPx data.
            CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
         ENDIF  ! end if (.not.des_continuum_coupled)

         IF(CALL_USR) CALL USR2_DES


! Write Drag force, Particle #, velocity, Contact force
! after every timestep. This routine only works for
! OpenMP and serial code
! Write only during the last iteration
        IF (NN .EQ. FACTOR) THEN
! NSTEP used instead of TIMESTEP_CUST
!                TIMESTEP_CUST=TIMESTEP_CUST+1
! Used for testing
!              PRINT *,SIZE(DRG_FC),DT_DRAG_PRINT,DT_VEL_PRINT
               CALL WRITE_CUST(NSTEP+1,DT_DRAG_PRINT,'DRAG',& 
                        &PART_INFO,DRG_FC,SIZE(DRG_FC),&
                        &TIMESTEP_CUST)
               CALL WRITE_CUST(NSTEP+1,DT_CONTACT_PRINT,&
                        &'CONTACT',PART_INFO,CONTACT_FC,&
                        &SIZE(CONTACT_FC),TIMESTEP_CUST)
               CALL WRITE_CUST(NSTEP+1,DT_VEL_PRINT,&
                        &'VELOCITY',PART_INFO,PART_VEL,&
                        & SIZE(PART_VEL),TIMESTEP_CUST)
               IF (TIMESTEP_CUST.EQ.0) THEN
                        TIMESTEP_CUST=1
               END IF
        END IF
      ENDDO ! end do NN = 1, FACTOR

! END DEM time loop
!-----------------------------------------------------------------<<<

      IF(CALL_USR) CALL USR3_DES

! Reset the discrete time step to original value.
      DTSOLID = DTSOLID_TMP

      IF(DES_CONTINUUM_COUPLED) CALL DESMPI_SEND_RECV_FIELD_VARS

      TMP_WALL = WALL_TIME() - TMP_WALL
      IF(TMP_WALL > 1.0d-10) THEN
         WRITE(ERR_MSG, 9000) trim(iVal(dble(FACTOR)/TMP_WALL))
      ELSE
         WRITE(ERR_MSG, 9000) '+Inf'
      ENDIF
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

 9000 FORMAT('    NITs/SEC = ',A)

      RETURN
      END SUBROUTINE DES_TIME_MARCH
