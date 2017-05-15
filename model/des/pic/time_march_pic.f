!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: PIC_TIME_MARCH                                          !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Main PIC driver routine.                                   !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PIC_TIME_MARCH

! Global variables
!---------------------------------------------------------------------//
! Fluid time, simulation end time, time step size, number of time steps
      use run, only: TIME, TSTOP, DT, NSTEP
! Discrete particle time, time step size
      use discretelement, only: S_TIME, DTSOLID
! MPPIC model step-size bounds
      use mfix_pic, only: DTPIC_MAX, DTPIC_CFL, DTPIC_TAUP
! Local particle count
      use discretelement, only: PIP
! Flag: Coupled fluid-solids simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag: Store _OLD arrays
      use discretelement, only: DO_OLD
! Flag: Call user defined subroutines
      use run, only: CALL_USR
! Flag: Explicitly coupled gas-solids drag
      use discretelement, only: DES_EXPLICITLY_COUPLED
! Number of mass outflows/inflows
      use pic_bc, only: PIC_BCMO, PIC_BCMI

! Module procedures
!---------------------------------------------------------------------//
      use desgrid, only: DESGRID_PIC
      use error_manager
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      use mpi_utility, only: GLOBAL_ALL_SUM
      use output_man, only: output_manager

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! time till which the PIC loop will be run
      double precision :: TEND_PIC_LOOP
! number of PIC time steps
      Integer :: PIC_ITERS
! Global number of parcels.
      INTEGER :: gPIP
!......................................................................!


! Set solids time to fluid time.
      S_TIME = TIME

      IF(DES_CONTINUUM_COUPLED) THEN
         TEND_PIC_LOOP = TIME+DT
         DTSOLID = min(DTPIC_MAX, DT)
      ELSE
         TEND_PIC_LOOP = TSTOP
         DTSOLID = DT
      ENDIF
      PIC_ITERS = 0


      IF(CALL_USR) CALL USR0_DES

! Compute the gas-phase pressure gradient
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_EXPLICITLY_COUPLED) CALL DRAG_GS_DES1
         CALL CALC_PG_GRAD
      ENDIF


! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the lagrangian loop
      DO WHILE(S_TIME.LT.TEND_PIC_LOOP)

         PIC_ITERS  = PIC_ITERS + 1

! Set the solids time step
!         DTSOLID = MERGE(MIN(DTPIC_MAX, DT), DTPIC_MAX,                &
!            DES_CONTINUUM_COUPLED)

! If next time step in the discrete loop will exceed the current time
! in the continuum simulation, modify the discrete time step so final
! time will match
         IF(S_TIME + DTSOLID > TEND_PIC_LOOP) &
            DTSOLID = TEND_PIC_LOOP - S_TIME

! Calculate the solids pressure
         CALL CALC_PS_PIC
         CALL CALC_PS_GRAD_PIC
         CALL INTERPOLATE_PIC

         IF(DES_CONTINUUM_COUPLED) CALL CALC_DRAG_DES

         IF (DO_OLD) CALL CFUPDATEOLD

         CALL INTEGRATE_TIME_PIC 

!         CALL WRITE_PARTICLE(6010)

! Apply mass outflow/inflow boundary conditions
         IF(PIC_BCMO > 0) CALL MASS_OUTFLOW_PIC
         IF(PIC_BCMI > 0) CALL MASS_INFLOW_PIC

! Impose the wall-particle boundary condition
         CALL APPLY_WALL_BC_PIC

! Exchange particle crossing processor boundaries
         CALL DESGRID_PIC(.TRUE.)
         CALL DES_PAR_EXCHANGE

         IF(S_TIME + DTSOLID < TEND_PIC_LOOP .OR. &
            .NOT.DES_EXPLICITLY_COUPLED ) THEN
! Bin particles to the fluid grid
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields
            CALL COMP_MEAN_FIELDS
         ENDIF

! This was moved from particles in cell and the passed variables should
! be added to particles in cell or made global.
         CALL REPORT_STATS_PIC
! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

! When coupled, all write calls are made in time_march (the continuum
! portion) according to user settings for spx_time and res_time.
! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES and SPx data.
            CALL OUTPUT_MANAGER(.FALSE., .FALSE.)
         ENDIF  ! end if (.not.des_continuum_coupled)

      ENDDO

      CALL GLOBAL_ALL_SUM(PIP, gPIP)
      WRITE(ERR_MSG, 3000) trim(iVal(PIC_ITERS)), trim(iVal(gPIP))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 3000 FORMAT(/'PIC NITs: ',A,3x,'Total PIP: ', A)

      RETURN
      END SUBROUTINE PIC_TIME_MARCH




!         !DTPIC_MAX = MIN( 1e-04, DTPIC_MAX)
!         IF(MOD(PIC_ITERS, 10).eq.0) then
!            IF(DES_CONTINUUM_COUPLED) then
!               WRITE(ERR_MSG, 2000) DTSOLID, DTPIC_CFL, DTPIC_TAUP, DT
!            ELSE
!               WRITE(ERR_MSG, 2001) S_TIME, DTSOLID, DTPIC_CFL, DTPIC_TAUP, DT
!            ENDIF
!            CALL FLUSH_ERR_MSG(HEADER = .FALSE., FOOTER = .FALSE.)
!         ENDIF
!
! 2000 FORMAT(/5x,'DTSOLID CURRENT  = ',g17.8,/5x,'DTPIC_CFL',8x,'= ',  &
!         g17.8, /5x,'DTPIC TAUP',7x,'= ',g17.8,/5x,'DT FLOW',10x,'= ', &
!         g17.8)
!
! 2001 FORMAT(/5x,'TIME',13X,'= ',g17.8,/5x,'DTSOLID CURRENT  = ',g17.8,&
!         /5x,'DTPIC_CFL',8X,'= ', g17.8,/5x,'DTPIC TAUP',7x,'= ',g17.8,&
!         /5x,'DT FLOW',10X,'= ', g17.8)




      SUBROUTINE WRITE_PARTICLE(NP)

      Use usr
      use compar
      use discretelement

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP

      INTEGER, SAVE :: CALLS = 0
      CHARACTER(len=128) :: FNAME

      FNAME=''; WRITE(FNAME, 2000) NP, myPE, CALLS
 2000 FORMAT('DBG/DBG_',I9.9,'_',I4.4,'_',I5.5,'.vtp')

      OPEN(UNIT=555, FILE=trim(FNAME), STATUS='UNKNOWN')

      write(*,"('Saving: ',A,' at ',F15.8)") trim(FNAME), S_TIME


      WRITE(555, 3000)
 3000 FORMAT('<?xml version="1.0"?>')

      WRITE(555, 3001)
 3001 FORMAT('<VTKFile type="PolyData" ' &
         'version="0.1" byte_order="LittleEndian">')

      WRITE(555,"('<PolyData>')")

      WRITE(555, 3002)
 3002 FORMAT('<Piece NumberOfPoints="1" ',         &
         'NumberOfVerts="0" NumberOfLines="0" ', &
         'NumberOfStrips="0" ',                    &
         'NumberOfPolys="0">')

      WRITE(555,"('<Points>')")

 3003 FORMAT('<DataArray type="Float32" Name="Position" ', &
         'NumberOfComponents="3" format="ascii">')
      WRITE(555, 3003)
      WRITE(555,"(3(3x,F15.8))") DES_POS_NEW(NP,:)
      WRITE(555,"('</DataArray>')")

      WRITE(555,"('</Points>')")

 3004 FORMAT('<PointData Scalars="Diameter" Vectors="Velocity">')
      WRITE(555, 3004)

 3005 FORMAT('<DataArray type="Float32" ', &
         'Name="Diameter" format="ascii">')
      WRITE(555, 3005)
      WRITE(555,"(3x,F15.8)") DES_RADIUS(NP)*2.0d0
      WRITE(555,"('</DataArray>')")

 3006 FORMAT('<DataArray type="Float32" Name="Velocity" ',&
         'NumberOfComponents="3" format="ascii">')
      WRITE(555, 3006)
      WRITE(555,"(3(3x,F15.8))") DES_VEL_NEW(NP,:)
      WRITE(555,"('</DataArray>')")

      WRITE(555,"('</PointData>')")
      WRITE(555,"('<CellData></CellData>')")
      WRITE(555,"('<Verts></Verts>')")
      WRITE(555,"('<Lines></Lines>')")
      WRITE(555,"('<Strips></Strips>')")
      WRITE(555,"('<Polys></Polys>')")
      WRITE(555,"('</Piece>')")
      WRITE(555,"('</PolyData>')")
      WRITE(555,"('</VTKFile>')")

      close(555)

      CALLS = CALLS+1

      RETURN
      END SUBROUTINE WRITE_PARTICLE

