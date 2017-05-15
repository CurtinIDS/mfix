! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  MODULE: MAIN                                                        !
!                                                                      !
!  Purpose: Main module for top level mfix subroutines.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE MAIN

      use exit, only: mfix_exit

!-----------------------------------------------
! Module variables
!-----------------------------------------------
! Final value of CPU time.
      DOUBLE PRECISION :: CPU1
! time used for computations.
      DOUBLE PRECISION :: CPUTIME_USED, WALLTIME_USED
! DISTIO variable for specifying the mfix version
      CHARACTER(LEN=512) :: version
! environment variable
!$ CHARACTER(LEN=512) :: omp_num_threads
!$ INTEGER :: length
!$ INTEGER :: status

! Number of iterations
      INTEGER :: NIT_TOTAL
! used for activating check_data_30
      INTEGER :: NCHECK, DNCHECK

! Flag to save results and cleanly exit.
      LOGICAL :: EXIT_SIGNAL = .FALSE.

      CHARACTER(LEN=80), DIMENSION(100) :: CMD_LINE_ARGS
      INTEGER :: CMD_LINE_ARGS_COUNT = 0

      CONTAINS

      SUBROUTINE INITIALIZE

#ifdef MPI
      USE mpi, only: mpi_comm_world, mpi_barrier  ! ignore-depcomp
#endif
      USE cdist, only: bdoing_postmfix
      USE cdist, only: bglobalnetcdf, bstart_with_one_res, bdist_io, bwrite_netcdf
      USE check, only: check_mass_balance
      USE check_data_cg, only: check_bc_flags, report_best_processor_size
      USE coeff, only: init_coeff
      USE compar, only: mpierr, mype, pe_io
      USE cont, only: do_cont
      USE cutcell, only: cartesian_grid, re_indexing, set_corner_cells
      USE discretelement, only: discrete_element
      USE drag, only: f_gs
      USE error_manager, only: err_msg, flush_err_msg
      USE error_manager, only: init_err_msg, finl_err_msg
      USE fldvar, only: rop_g, rop_s
      USE funits, only: dmp_log, unit_log, unit_res
      USE machine, only: start_log, end_log
      USE machine, only: wall_time, pc_quickwin, machine_cons, get_run_id, start_log, end_log
      USE mfix_netcdf, only: mfix_usingnetcdf
      USE output, only: dbgprn_layout
      USE output_man, only: init_output_vars, output_manager
      USE parallel_mpi, only: parallel_init, parallel_fin
      USE param1, only: n_spx, undefined, zero
      USE pgcor, only: d_e, d_n, d_t, phase_4_p_g, switch_4_p_g
      USE physprop, only: mmax
      USE pscor, only: e_e, e_n, e_t, do_p_s, phase_4_p_s, mcp, switch_4_p_s
      USE qmom_kinetic_equation, only: qmomk
      USE read_input, only: get_data
      USE run, only: id_version, ier
      USE run, only: automatic_restart, call_usr, dem_solids, dt_max, dt_min
      USE run, only: iter_restart, nstep, pic_solids, run_type, dt, shear, time, v_sh
      USE time_cpu, only: CPU00, wall0
      USE time_cpu, only: cpu_io, cpu_nlog, cpu0, cpuos, time_nlog
      USE vtk, only: write_vtk_files

      IMPLICIT NONE

!$    INTEGER num_threads, threads_specified, omp_id
!$    INTEGER omp_get_num_threads
!$    INTEGER omp_get_thread_num

      ! Temporary storage for DT
      DOUBLE PRECISION :: DT_tmp
      ! Save TIME in input file for RESTART_2
      DOUBLE PRECISION :: TIME_SAVE

      INTEGER :: LL, MM

! DISTIO
! If you change the value below in this subroutine, you must also
! change it in write_res0.f and the value should also be consistent
! with the check in read_res0
      version = 'RES = 01.6'

      bDoing_postmfix = .false.

! Invoke MPI initialization routines and get rank info.
      CALL PARALLEL_INIT
      CALL GEN_LOG_BASENAME

! we want only PE_IO to write out common error messages
      DMP_LOG = (myPE == PE_IO)

! set the version.release of the software
      ID_VERSION = '2016-1'

! set automatic restart flag to false
!      AUTOMATIC_RESTART = .FALSE.
!      ITER_RESTART      = 1

! specify the number of processors to be used
!$        call get_environment_variable("OMP_NUM_THREADS",omp_num_threads,length,status, .true.)
!$      if (status.eq.0 .and. length.ne.0) then
!$        read(omp_num_threads,*) threads_specified
!$      else
!$        WRITE(*,'(A,$)') 'Enter the number of threads to be used for SMP: '
!$        READ(*,*) threads_specified
!$      endif

!$      call omp_set_num_threads(threads_specified)

! Find the number of processors used
!$omp  parallel
!$      num_threads = omp_get_num_threads()
!$      omp_id = omp_get_thread_num()
!$      if(omp_id.eq.0) Write(*,*)' Number of threads used for SMP = ',  num_threads
!$omp  end parallel

! Set machine dependent constants
      CALL MACHINE_CONS

! Get the date and time. They give the unique run_id in binary output
! files
      CALL GET_RUN_ID

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates. timestep at the beginning of execution
      CALL CPU_TIME (CPU00)
      WALL0 = WALL_TIME()

! Read input data, check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      CALL GET_DATA

! Write the initial part of the standard output file
      CALL WRITE_OUT0
      IF(.NOT.CARTESIAN_GRID)  CALL WRITE_FLAGS

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

!$    CALL START_LOG
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' Number of processors used = ', threads_specified
!$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
!$    CALL END_LOG

!  setup for PC quickwin application
      CALL PC_QUICKWIN

      CALL INIT_ERR_MSG('MFIX')


! if not netcdf writes asked for ... globally turn off netcdf
      if(MFIX_usingNETCDF()) then
         bGlobalNetcdf = .false.
         do LL = 1,20
            if (bWrite_netcdf(LL)) bGlobalNetcdf = .true.
         enddo
      endif

      DT_TMP = DT
      SELECT CASE (TRIM(RUN_TYPE))

      CASE ('NEW')
! Write the initial part of the restart files
         CALL WRITE_RES0
         DO LL = 1, N_SPX
            CALL WRITE_SPX0 (LL, 0)
         ENDDO

      CASE ('RESTART_1')
! Read the time-dependent part of the restart file
         CALL READ_RES1
         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

      CASE ('RESTART_2')
         TIME_SAVE = TIME
! DISTIO
         if (myPE .ne. PE_IO .and. bDist_IO .and. bStart_with_one_res) then
            write (unit_res,rec=1) version
            write (unit_res,rec=2) 4
            write (unit_res,rec=3) 4
         endif

         CALL READ_RES1
         TIME = TIME_SAVE

1010     FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
              G12.5,/'Time step number (NSTEP) =',I7)

         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

         CALL WRITE_RES0

! Writing the RES1 and SPX1 can only be done here when re-indexing is turned off
! This will be done after the cell re-indexing is done later in this file.
! This allows restarting independently of the re-indexing setting between
! the previous and current run.
         IF(.NOT.RE_INDEXING) THEN
            CALL WRITE_RES1
            DO LL = 1, N_SPX
               CALL WRITE_SPX0 (LL, 0)
               CALL WRITE_SPX1 (LL, 0)
            END DO
            call write_netcdf(0,0,time)
         ENDIF

      CASE DEFAULT
         CALL START_LOG
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
              ' MFIX: Do not know how to process'
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file'
         CALL END_LOG
         call mfix_exit(myPE)

      END SELECT

#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

      IF (DT_TMP /= UNDEFINED) THEN
         DT = MAX(DT_MIN,MIN(DT_MAX,DT))
      ELSE
         DT = DT_TMP
      ENDIF

! Set arrays for computing indices. A secondary call is made
! after cut cell-preprocessing to update array indices.
      IF(CARTESIAN_GRID) THEN
         CALL SET_INCREMENTS
         CALL SET_INCREMENTS3
      ENDIF

!      IF(.NOT.RE_INDEXING) CALL WRITE_IJK_VALUES

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1

!  Update flags for Cartesian_GRID.
      IF(CARTESIAN_GRID) CALL CHECK_BC_FLAGS

! Calculate cell volumes and face areas
      IF(.NOT.CARTESIAN_GRID) CALL SET_GEOMETRY1

! Find corner cells and set their face areas to zero
      IF(.NOT.CARTESIAN_GRID)  THEN
         CALL GET_CORNER_CELLS()
      ELSE
         IF (SET_CORNER_CELLS)  CALL GET_CORNER_CELLS ()
      ENDIF

! Set constant physical properties
      CALL SET_CONSTPROP

! Set initial conditions
      CALL SET_IC

! Set point sources.
      CALL SET_PS

! Set boundary conditions
      CALL ZERO_NORM_VEL
      CALL SET_BC0

! Cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_BC0

! Set gas mixture molecular weight
      CALL SET_MW_MIX_G

! Set the pressure field for a fluidized bed
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P

! Initialize densities.
      IF (RUN_TYPE == 'NEW') CALL SET_RO_G
      IF (RUN_TYPE == 'NEW') CALL SET_RO_S

! Initialize time dependent boundary conditions
      CALL SET_BC1

! Check the field variable data and report errors.
      IF(.NOT.CARTESIAN_GRID)  CALL CHECK_DATA_20

!=======================================================================
! JFD: START MODIFICATION FOR RE-INDEXING CELLS
!=======================================================================
      IF(CARTESIAN_GRID.AND.RE_INDEXING) THEN

         IF(myPE == PE_IO) THEN
            WRITE(*,"(72('='))")
            WRITE(*,*)' RE-INDEXING CELLS FOR CARTESIAN GRID...'
         ENDIF
         CALL RE_INDEX_ARRAYS


         !IF(myPE == PE_IO)print*,'Calling REPORT_BEST_IJK_SIZE:'
         !CALL REPORT_BEST_IJK_SIZE
         CALL REPORT_BEST_PROCESSOR_SIZE
         !IF(myPE == PE_IO)print*,'Exiting MFIX after REPORT_BEST_IJK_SIZE.'

         IF(myPE == PE_IO) WRITE(*,"(72('='))")

! In case of a RESTART_2, write the RES1 and SPX1 files here
! This was commented out earlier in this file.
         IF(RUN_TYPE == 'RESTART_2') THEN
            CALL WRITE_RES1
            DO LL = 1, N_SPX
               CALL WRITE_SPX0 (LL, 0)
               CALL WRITE_SPX1 (LL, 0)
            END DO
            call write_netcdf(0,0,time)
         ENDIF
      ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR RE-INDEXING CELLS
!=======================================================================

! Setup VTK data for regular (no cut cells) grid
      IF(.NOT.CARTESIAN_GRID.AND.WRITE_VTK_FILES) CALL SETUP_VTK_NO_CUTCELL

      IF(DISCRETE_ELEMENT) CALL MAKE_ARRAYS_DES
      IF(QMOMK) CALL QMOMK_MAKE_ARRAYS

! Set the inflow/outflow BCs for DEM solids
      IF(DEM_SOLIDS) CALL SET_BC_DEM
! Set the inflow/outflow BC for PIC solids
      IF(PIC_SOLIDS) CALL SET_BC_PIC

! Set the inital properties of each particle.
      IF(DEM_SOLIDS) CALL SET_IC_DEM

! AEOLUS: debug prints
      if (DBGPRN_LAYOUT .or. bdist_io) then
         !write (*,*) myPE , ' E.4 ... version = ' , version(1:33)
         call debug_write_layout()
         call write_parallel_info()
      endif

! Initializations for CPU time calculations in iterate
      CPUOS = 0.
      CALL CPU_TIME (CPU1)
      CPU_NLOG = CPU1
      TIME_NLOG = TIME - DT

! Get the initial value of CPU time
      CALL CPU_TIME (CPU0)

! Find the solution of the equations from TIME to TSTOP at
! intervals of DT

!-----------------------------------------------

      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0

      CALL INIT_OUTPUT_VARS

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

      CALL RRATES_INIT()

! Calculate all the coefficients once before entering the time loop
      CALL INIT_COEFF(IER)

      DO MM=1, MMAX
         F_gs(1,MM) = ZERO
      ENDDO

! Remove undefined values at wall cells for scalars
      CALL UNDEF_2_0 (ROP_G)
      DO MM = 1, MMAX
         CALL UNDEF_2_0 (ROP_S(1,MM))
      ENDDO

! Initialize d's and e's to zero
      DO MM = 0, MMAX
         D_E(1,MM) = ZERO
         D_N(1,MM) = ZERO
         D_T(1,MM) = ZERO
      ENDDO
      E_E(:) = ZERO
      E_N(:) = ZERO
      E_T(:) = ZERO

! calculate shear velocities if periodic shear BCs are used
      IF(SHEAR) CALL CAL_D(V_sh)

! Initialize check_mass_balance.  This routine is not active by default.
! Specify a reporting interval (hard-wired in the routine) to activate
! the routine.
      Call check_mass_balance (0)

! sof modification: now it's only needed to do this once before time-loop
! Mark the phase whose continuity will be solved and used to correct
! void/volume fraction in calc_vol_fr (see subroutine for details)
      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
           DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S)

      END SUBROUTINE INITIALIZE

      SUBROUTINE FINALIZE

      USE cutcell, only: cartesian_grid
      USE dashboard
      USE cut_cell_preproc, only: close_cut_cell_files
      USE error_manager, only: finl_err_msg
      USE machine, only: wall_time
      USE parallel_mpi, only: parallel_fin
      USE run, only: dt, call_usr, dt_min, get_tunit, tunit
      USE time_cpu
      IMPLICIT NONE

! Call user-defined subroutine after time-loop.
      IF (CALL_USR) CALL USR3

! Get the final value of CPU time.  The difference gives the
! CPU time used for the computations.
      CALL CPU_TIME (CPU1)

! Compute the CPU time and write it out in the .OUT file.
      CPUTIME_USED = CPU1 - CPU0 - CPU_IO
      WALLTIME_USED = WALL_TIME() - WALL0
      CALL WRITE_OUT3 (CPUTIME_USED, WALLTIME_USED, CPU_IO)

! JFD: cartesian grid implementation
      IF(WRITE_DASHBOARD) THEN
         IF(DT>=DT_MIN) THEN
            RUN_STATUS = 'Complete.'
         ELSE
            RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
         ENDIF
         CALL GET_TUNIT(CPUTIME_USED,TUNIT)
         CALL UPDATE_DASHBOARD(0,CPUTIME_USED,TUNIT)
      ENDIF
      IF(CARTESIAN_GRID)  CALL CLOSE_CUT_CELL_FILES

! Finalize and terminate MPI
      call parallel_fin

      CALL FINL_ERR_MSG

      END SUBROUTINE FINALIZE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GEN_LOG_BASENAME                                        !
!  Author: Aytekin Gel                                Date: 19-SEP-03  !
!                                                                      !
!  Purpose: Generate the file base for DMP logs.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GEN_LOG_BASENAME

      use compar, only: myPE
      use compar, only: fbname

      implicit none

! Variables for generating file basename with processor id
      INTEGER :: i1, i10, i100, i1000, i10000

! PAR_I/O Generate file basename for LOG files
      i10000 = int(myPE/10000)
      i1000  = int((myPE-i10000*10000)/1000)
      i100   = int((myPE-i10000*10000-i1000*1000)/100)
      i10    = int((myPE-i10000*10000-i1000*1000-i100*100)/10)
      i1     = int((myPE-i10000*10000-i1000*1000-i100*100-i10*10)/1)

      i10000 = i10000 + 48
      i1000  = i1000  + 48
      i100   = i100   + 48
      i10    = i10    + 48
      i1     = i1     + 48

      fbname=char(i10000)//char(i1000)//char(i100)//char(i10)//char(i1)

      RETURN
      END SUBROUTINE GEN_LOG_BASENAME

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: debug_write()                                          C
!  Purpose: Write out full geometry index setup information for the
!  case
!                                                                      C
!  Author: Aytekin Gel                                Date: 19-SEP-03  C
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

   SUBROUTINE debug_write_layout()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE cdist
      USE compar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE leqsol
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE run
      USE sendrecv
      USE sendrecv3
      USE time_cpu
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! phase index
      INTEGER :: M
! indices
      INTEGER :: i, j, k, ijk, ijk_GL, ijk_PROC, ijk_IO
!
      integer :: indxA, indxA_gl, indxB, indxB_gl, indxC, indxC_gl
      integer :: indxD, indxD_gl, indxE, indxE_gl, indxF, indxF_gl
      integer :: indxG, indxG_gl, indxH, indxH_gl
!
      logical :: amgdbg = .TRUE.

      character(LEN=80) :: fname

!DISTIO
!      fname = "layout_xxxx.txt"
!      write (fname(8:11),'(i4.4)') myPE
      fname = "layout_xxxxx.txt"
      write (fname(8:12),'(i5.5)') myPE
      open (unit=11,file=fname,status='unknown')

      write (11,*) ' ********************************************'
      write (11,*) ' ********************************************'
      write (11,*) ' ********************************************'
      write (11,*) ' ********************************************'
      write (11,*) ' '
      write (11,*) ' '
      write (11,*) ' myPE =           ' , myPE
      write (11,*) ' '
      write (11,*) ' '


      IF (AMGDBG .OR. bDist_IO) THEN
         write(11,"('BLK1: Running from istart3,iend3 .AND. jstart3, jend3 .AND. kstart3, kend3')")
         write(11,"(' (   i ,    j,     k) =>    ijk      ijk_GL     ijk_PROC    ijk_IO')")
         write(11,"(' ====================      =====     =======    ========    ======')")
         DO k = kstart3, kend3
            DO i = istart3,iend3
               DO j = jstart3, jend3
                  ijk = FUNIJK(i,j,k)
                  ijk_GL = FUNIJK_GL(i,j,k)
                  ijk_PROC = FUNIJK_PROC(i,j,k,myPE)
                  ijk_IO = FUNIJK_IO(i,j,k)
                  write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',4(I8,' , '))") &
                       i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
               ENDDO
            ENDDO
         ENDDO

         write(11,"(/,/,'BLK2: Print out Bottom, South, West, East, North, Top neighbors')")
         write(11,"(' (   i ,    j,     k) =>    ijk    ijk_GL    B_of    S_of    W_of    E_of    N_of    T_of')")
         write(11,"(' ====================      =====   =======  ======  ======  ======  ======  ======  ======')")
         DO k = kstart3, kend3
            DO i = istart3,iend3
               DO j = jstart3, jend3
                  ijk = FUNIJK(i,j,k)
                  ijk_GL = FUNIJK_GL(i,j,k)
                  write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
                       i,j,k,ijk,ijk_GL,bottom_of(ijk),south_of(ijk),west_of(ijk),&
                       east_of(ijk),north_of(ijk),top_of(ijk)
               ENDDO
            ENDDO
         ENDDO

         write(11,"(/,/,'BLK3: Print out km, jm, im, ip, jp, kp neighbors')")
         write(11,"(' (   i ,    j,     k) =>    ijk    ijk_GL    km_of   jm_of   im_of   ip_of   jp_of   kp_of')")
         write(11,"(' ====================      =====   =======  ======  ======  ======  ======  ======  ======')")
         DO k = kstart3, kend3
            DO i = istart3,iend3
               DO j = jstart3, jend3
                  ijk = FUNIJK(i,j,k)
                  ijk_GL = FUNIJK_GL(i,j,k)
                  write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
                       i,j,k,ijk,ijk_GL,km_of(ijk),jm_of(ijk),im_of(ijk),&
                       ip_of(ijk),jp_of(ijk),kp_of(ijk)
               ENDDO
            ENDDO
         ENDDO

         write(11,"(/,'BLK4a: Active Fluid Cells:FLUID_AT(ijk)=.T.',/,&
              &           ' (   i ,    j,     k) =>    ijk  [   x ,     ,     z]')")
         write(11,"(' ====================      =====  ====================')")
         DO ijk = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            !         IF (FLOW_AT_E(IJK)) THEN
            IF (FLUID_AT(IJK)) THEN
               !          write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
               write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8,' [',E12.5,',',E12.5,' ]')") I,J,K,ijk,X(i),Z(k)
            ENDIF
         ENDDO

         write(11,"(/,'BLK4b: Cells that are (.NOT.WALL_AT(IJK)) = .T.',/,&
              &           ' (   i ,    j,     k) =>    ijk  [   x ,     ,     z]')")
         write(11,"(' ====================      =====  ====================')")
         DO ijk = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IF (.NOT.WALL_AT(IJK)) THEN
               !          write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
               write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8,' [',E12.5,',',E12.5,' ]')") I,J,K,ijk,X(i),Z(k)
            ENDIF
         ENDDO

         DO k = kstart3, kend3
            DO i = istart3,iend3
               DO j = jstart3, jend3
                  ijk = FUNIJK(i,j,k)
                  ijk_GL = FUNIJK_GL(i,j,k)

                  if (i == istart2 .AND. j == jstart2) then
                     indxA = ijk
                     indxA_gl = ijk_GL
                  endif
                  if (i == istart1 .AND. j == jstart1) then
                     indxE = ijk
                     indxE_gl = ijk_GL
                  endif
                  if (i == istart2 .AND. j == jend2) then
                     indxB = ijk
                     indxB_gl = ijk_GL
                  endif
                  if (i == istart1 .AND. j == jend1) then
                     indxF = ijk
                     indxF_gl = ijk_GL
                  endif
                  if (i == iend1 .AND. j == jstart1) then
                     indxH = ijk
                     indxH_gl = ijk_GL
                  endif
                  if (i == iend2 .AND. j == jstart2) then
                     indxD = ijk
                     indxD_gl = ijk_GL
                  endif
                  if (i == iend1 .AND. j == jend1) then
                     indxG = ijk
                     indxG_gl = ijk_GL
                  endif
                  if (i == iend2 .AND. j == jend2) then
                     indxC = ijk
                     indxC_gl = ijk_GL
                  endif
               ENDDO
            ENDDO
            write(11,"('BLK5:')")
            write(11,"(57('='))")
            write(11,"('k= ',I5,/,57('='))") k
            write(11,"('B= ',I5,' (',I7,')',20X,'C= ',I5,' (',I7,')',/)") indxB, indxB_gl, &
                 indxC, indxC_gl
            !        write(UNIT_LOG,"(' \',34X,'/')")
            !        write(UNIT_LOG,"(2X,'\',32X,'/')")
            write(11,"(3X,'F= ',I5,' (',I7,')',12X,'G= ',I5,' (',I7,')')") indxF, indxF_gl, &
                 indxG, indxG_gl
            write(11,"(4(9X,'|',29X,'|',/))")
            write(11,"(3X,'E= ',I5,' (',I7,')',12X,'H= ',I5,' (',I7,')',/)") indxE, indxE_gl, &
                 indxH, indxH_gl
            !        write(UNIT_LOG,"(2X,'/',32X,'\')")
            !        write(UNIT_LOG,"('/',34X,'\')")
            write(11,"('A= ',I5,' (',I7,')',20X,'D= ',I5,' (',I7,')',/,/)") indxA, indxA_gl, &
                 indxD, indxD_gl

            !        write(UNIT_LOG,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
            !                                         i,j,k,ijk,ijk_GL,bottom_of(ijk),south_of(ijk),west_of(ijk),&
            !                                        east_of(ijk),north_of(ijk),top_of(ijk)

         ENDDO

         !      write(UNIT_LOG,"(/,' (   i ,    j,     k) =>    ijk (Active Fluid)')")
         !      write(UNIT_LOG,"(' ====================      =====')")
         !       DO ijk = ijkstart3, ijkend3
         !         I = I_OF(IJK)
         !         J = J_OF(IJK)
         !         K = K_OF(IJK)

         !         IF (FLOW_AT_E(IJK)) THEN
         !         IF (FLUID_AT(IJK)) THEN
         !           write(UNIT_LOG,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
         !         ENDIF
         !      END DO


      endif   ! end if(amgdbg .or. bdist_io)

      M = 0
      !      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER)

      IF (AMGDBG .OR. bDist_IO) THEN
         write(11,"(/,/,'BLK6: ========= ORIGINAL MFIX VARIABLES ===========')")
         write(11,"('PE ',I5,': imin1  = ',I6,3X,'imax1= ',I6,/,'PE ',I5,': jmin1  = ',I6,3X,'jmax1= ',I6)") &
              myPE,imin1,imax1,myPE,jmin1,jmax1
         write(11,"('PE ',I5,': kmin1  = ',I6,3X,'kmax1= ',I6)") myPE,kmin1,kmax1
         write(11,"('-----')")
         write(11,"('PE ',I5,': imin2  = ',I6,3X,'imax2= ',I6,/,'PE ',I5,': jmin2  = ',I6,3X,'jmax2= ',I6)") &
              myPE,imin2,imax2,myPE,jmin2,jmax2
         write(11,"('PE ',I5,': kmin2  = ',I6,3X,'kmax2= ',I6)") myPE,kmin2,kmax2
         write(11,"('----- Below xxx3 set is DMP extension ------------')")
         write(11,"('PE ',I5,': imin3  = ',I6,3X,'imax3= ',I6,/,'PE ',I5,': jmin3  = ',I6,3X,'jmax3= ',I6)") &
              myPE,imin3,imax3,myPE,jmin3,jmax3
         write(11,"('PE ',I5,': kmin3  = ',I6,3X,'kmax3= ',I6)") myPE,kmin3,kmax3
         write(11,"('----- End of Below xxx3 set is DMP extension -----')")
         !      write(11,"('PE ',I5,': ijkmax2= ',I6)") myPE,ijkmax2
         write(11,"('PE ',I5,': ijmax2 = ',I6)") myPE,ijmax2
         write(11,"('PE ',I5,': ijkmin1= ',I6,' ijkmax1= ',I12)") myPE,ijkmin1, ijkmax1
         write(11,"('PE ',I5,':          ',6X,' ijkmax2= ',I12)") myPE,ijkmax2
         write(11,"('PE ',I5,':          ',6X,' ijkmax3= ',I12)") myPE,ijkmax3
         write(11,"('PE ',I5,': ijkmin4= ',I6,' ijkmax4= ',I12)") myPE,ijkmin4, ijkmax4


         write(11,"(/,/,' ========= DMP EXTENSION VARIABLES ===========')")
         !      write(UNIT_LOG,"('PE ',I5,': ijksize  = ',I6)") myPE,ijksize
         write(11,"('PE ',I5,': ijksize3 = ',I6,3X,'ijksize3_all = ',I6)") myPE,ijksize3,ijksize3_all(myPE)
         write(11,"('PE ',I5,': ijksize4 = ',I6,3X,'ijksize4_all = ',I6)") myPE,ijksize4,ijksize4_all(myPE)
         write(11,"('PE ',I5,': ijkstart3  = ',I6,3X,'ijkend3  = ',I6)") myPE,ijkstart3, ijkend3
         write(11,"('PE ',I5,': ijkstart3_all = ',I6,3X,'ijkstart4_all = ',I6)") myPE,ijkstart3_all(myPE),ijkstart4_all(myPE)
         write(11,"('PE ',I5,': istart_all = ',I6,3X,'iend_all = ',I6,/,'PE ',I5,': jstart_all = ',I6,3X,'jend_all = ',I6)") &
              myPE,istart_all(myPE),iend_all(myPE),myPE,jstart_all(myPE),jend_all(myPE)
         write(11,"('PE ',I5,': kstart_all = ',I6,3X,'kend_all = ',I6,/,'----------------------')") &
              myPE,kstart_all(myPE),kend_all(myPE)

         write(11,"('PE ',I5,': istart1_all= ',I6,3X,'iend1_all= ',I6,/,'PE ',I5,': jstart1_all= ',I6,3X,'jend3_all= ',I6)") &
              myPE,istart1_all(myPE),iend1_all(myPE),myPE,jstart1_all(myPE),jend1_all(myPE)
         write(11,"('PE ',I5,': kstart1_all= ',I6,3X,'kend1_all= ',I6,/,'----------------------')") &
              myPE,kstart1_all(myPE),kend1_all(myPE)

         write(11,"('PE ',I5,': istart2_all= ',I6,3X,'iend2_all= ',I6,/,'PE ',I5,': jstart2_all= ',I6,3X,'jend3_all= ',I6)") &
              myPE,istart2_all(myPE),iend2_all(myPE),myPE,jstart2_all(myPE),jend2_all(myPE)
         write(11,"('PE ',I5,': kstart2_all= ',I6,3X,'kend2_all= ',I6,/,'----------------------')") &
              myPE,kstart2_all(myPE),kend2_all(myPE)

         write(11,"('PE ',I5,': istart3_all= ',I6,3X,'iend3_all= ',I6,/,'PE ',I5,': jstart3_all= ',I6,3X,'jend3_all= ',I6)") &
              myPE,istart3_all(myPE),iend3_all(myPE),myPE,jstart3_all(myPE),jend3_all(myPE)
         write(11,"('PE ',I5,': kstart3_all= ',I6,3X,'kend3_all= ',I6,/,'----------------------')") &
              myPE,kstart3_all(myPE),kend3_all(myPE)

         write(11,"('PE ',I5,': istart1= ',I6,3X,'iend1= ',I6,/,'PE ',I5,': jstart1= ',I6,3X,'jend1= ',I6)") &
              myPE,istart1,iend1,myPE,jstart1,jend1
         write(11,"('PE ',I5,': kstart1= ',I6,3X,'kend1= ',I6,/,'----------------------')") &
              myPE,kstart1,kend1
         write(11,"('PE ',I5,': istart2= ',I6,3X,'iend2= ',I6,/,'PE ',I5,': jstart2= ',I6,3X,'jend2= ',I6)") &
              myPE,istart2,iend2,myPE,jstart2,jend2
         write(11,"('PE ',I5,': kstart2= ',I6,3X,'kend2= ',I6,/,'----------------------')") &
              myPE,kstart2,kend2
         write(11,"('PE ',I5,': istart3= ',I6,3X,'iend3= ',I6,/,'PE ',I5,': jstart3= ',I6,3X,'jend3= ',I6)") &
              myPE,istart3,iend3,myPE,jstart3,jend3
         write(11,"('PE ',I5,': kstart3= ',I6,3X,'kend3= ',I6,/,'----------------------')") &
              myPE,kstart3,kend3

      ENDIF   ! end if(amgdbg .or. bdist_io)

      close(unit=11)


      RETURN
   END SUBROUTINE DEBUG_WRITE_LAYOUT



   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   SUBROUTINE write_parallel_info()

      !-----------------------------------------------
      !   M o d u l e s
      !-----------------------------------------------
      USE compar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE leqsol
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE run
      USE sendrecv
      USE sendrecv3
      USE time_cpu
      IMPLICIT NONE
      !-----------------------------------------------
      ! Dummy arguments
      !-----------------------------------------------
      ! Local Variables
      !-----------------------------------------------
      ! phase index
      INTEGER :: M
      ! indices
      INTEGER :: i, j, k, ijk, ijk_GL, ijk_PROC, ijk_IO
      !
      character(LEN=80) :: fname
      !-----------------------------------------------

      !DISTIO
      !      fname = "p_info_xxxx.txt"
      !      write (fname(8:11),'(i4.4)') myPE
      fname = "p_info_xxxxx.txt"
      write (fname(8:12),'(i5.5)') myPE
      open (unit=11,file=fname,status='unknown')

      write (11,*) myPe , ' = myPE'

      write (11,*) myPE , istart3,iend3
      write (11,*) myPE , jstart3,jend3
      write (11,*) myPE , kstart3,kend3

      write(11,"('BLK1: Running from istart3,iend3 .AND. jstart3, jend3 .AND. kstart3, kend3')")
      write(11,"(' (   i ,    j,     k)       ijk      ijk_GL     ijk_PROC    ijk_IO')")
      write(11,"(' ====================      =====     =======    ========    ======')")
      DO k = kstart3, kend3
         DO i = istart3,iend3
            DO j = jstart3, jend3
               ijk = FUNIJK(i,j,k)
               ijk_GL = FUNIJK_GL(i,j,k)
               ijk_PROC = FUNIJK_PROC(i,j,k,myPE)
               ijk_IO = FUNIJK_IO(i,j,k)
               write(11,"('  ',I4,'   ',I4,'   ',I4,'     ',4(I8,'   '))" ) &
                    i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
            ENDDO
         ENDDO
      ENDDO

      M = 0
      !      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER)

      close(unit=11)

      RETURN
   END SUBROUTINE write_parallel_info

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  SUBROUTINE: DO_MPI_BCAST                                            !
   !                                                                      !
   !  Purpose: Used by pymfix for broadcasting commands.                  !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   function do_mpi_bcast(str)
#ifdef MPI
      use mpi ! ignore-depcomp
#endif
      implicit none

      ! TODO: return a dynamically allocated string, instead of fixed size of 100,000 bytes
      character(len=100000),intent(in) :: str
      character :: aa(100000)
      character :: do_mpi_bcast(100000)
      integer :: ii
      integer :: ierr

      do ii = 1,len(str)
         aa(ii) = str(ii:ii)
      end do

#ifdef MPI
      call mpi_bcast(aa,100000,mpi_character,0,mpi_comm_world,ierr)
#endif

      do ii = 1,100000
         do_mpi_bcast(ii:ii) = aa(ii)
      end do

   end function do_mpi_bcast

   subroutine do_write_dbg_vtu_and_vtp_files
      implicit none
      call write_dbg_vtu_and_vtp_files
   end subroutine do_write_dbg_vtu_and_vtp_files

   subroutine do_backupres
      use output_man, only: backup_res
      implicit none
      call backup_res
   end subroutine do_backupres

   subroutine do_reinit(filename)
      use reinit, only: reinitialize
      implicit none
      ! filename of uploaded mfix.dat file
      character(len=*), intent(in) :: filename
      call reinitialize(filename)
   end subroutine do_reinit

   subroutine do_abort
      use compar, only: mype
      implicit none
      call mfix_exit(mype)
   end subroutine do_abort

   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !  Subroutine: ADD_COMMAND_LINE_ARGUMENT                               !
   !  Author: M.Meredith                                 Date: 03-FEB-16  !
   !                                                                      !
   !  Purpose: Save command line arguments in CMD_LINE_ARGS array.        !
   !           Used by both mfix.f and pymfix.                            !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE ADD_COMMAND_LINE_ARGUMENT(ARG)
      implicit none
      CHARACTER(LEN=80), INTENT(IN) :: ARG

      CMD_LINE_ARGS_COUNT = CMD_LINE_ARGS_COUNT + 1

      if (CMD_LINE_ARGS_COUNT > 100) THEN
         print *,"TOO MANY COMMAND LINE ARGUMENTS"
         stop
      ENDIF

      CMD_LINE_ARGS(CMD_LINE_ARGS_COUNT) = arg

   END SUBROUTINE ADD_COMMAND_LINE_ARGUMENT

END MODULE MAIN
