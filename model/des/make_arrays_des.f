!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES

      USE calc_collision_wall
      USE compar
      use constant, only: PI
      USE cutcell
      USE des_rxns
      USE des_thermo
      USE desgrid
      USE discretelement
      USE error_manager
      USE functions
      USE funits
      use fldvar, only: ro_s
      USE GENERATE_PARTICLES, only: GENERATE_PARTICLE_CONFIG
      USE geometry
      use mfix_pic, only:  MPPIC
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      USE mpi_utility
      USE param1
      use physprop, only: mmax, ro_s0
      USE run
      USE stl
      USE stl_functions_des
      use stl_preproc_des, only: add_facet

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L, IJK
      INTEGER :: I1, I2, J1, J2, K1, K2, II, JJ, KK, IJK2
      INTEGER :: lcurpar, lpip_all(0:numpes-1), lglobal_id
      INTEGER :: CELL_ID, I_CELL, J_CELL, K_CELL, COUNT, NF
      INTEGER :: IMINUS1, IPLUS1, JMINUS1, JPLUS1, KMINUS1, KPLUS1

! MPPIC related quantities
      CALL INIT_ERR_MSG("MAKE_ARRAYS_DES")

! Check interpolation input.
      CALL SET_FILTER_DES

! cfassign and des_init_bc called before reading the particle info
      CALL CFASSIGN

      VOL_SURR(:) = ZERO

      ! initialize VOL_SURR array
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = funijk(I,J,K)
               I1 = I
               I2 = I+1
               J1 = J
               J2 = J+1
               K1 = K
               K2 = merge(K, K+1, NO_K)

! looping over stencil points (node values)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells
                        IJK2 = funijk_map_c(II, JJ, KK)
                        IF(FLUID_AT(IJK2)) VOL_SURR(IJK) = &
                           VOL_SURR(IJK)+VOL(IJK2)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO



! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN
         IF(PARTICLES /= 0) THEN
            IF(GENER_PART_CONFIG) THEN
               CALL GENERATE_PARTICLE_CONFIG
            ELSE
               CALL READ_PAR_INPUT
            ENDIF
         ENDIF

! Set the global ID for the particles and set the ghost cnt
         ighost_cnt = 0
         lpip_all = 0
         lpip_all(mype) = pip
         call global_all_sum(lpip_all)
         lglobal_id = sum(lpip_all(0:mype-1))
         imax_global_id = 0
         do lcurpar  = 1,pip
            lglobal_id = lglobal_id + 1
            iglobal_id(lcurpar) = lglobal_id
            imax_global_id = iglobal_id(pip)
         end do
         call global_all_max(imax_global_id)

! Initialize old values
         omega_new(:,:)   = zero

! Particle orientation
         IF(PARTICLE_ORIENTATION) THEN
            ORIENTATION(1,:) = INIT_ORIENTATION(1)
            ORIENTATION(2,:) = INIT_ORIENTATION(2)
            ORIENTATION(3,:) = INIT_ORIENTATION(3)
         ENDIF


         IF (DO_OLD) THEN
            omega_old(:,:)   = zero
            des_pos_old(:,:) = des_pos_new(:,:)
            des_vel_old(:,:) = des_vel_new(:,:)
         ENDIF

! Read the restart file.
      ELSEIF(RUN_TYPE == 'RESTART_1' .OR. RUN_TYPE == 'RESTART_2') THEN

         CALL READ_RES0_DES
         imax_global_id = maxval(iglobal_id(1:pip))
         call global_all_max(imax_global_id)

! Initizlie the old values.
         IF (DO_OLD) THEN
            omega_old(:,:)   = omega_new(:,:)
            des_pos_old(:,:) = des_pos_new(:,:)
            des_vel_old(:,:) = des_vel_new(:,:)
         ENDIF

      ELSE

         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1100 FORMAT('Error 1100: Unsupported RUN_TYPE for DES.')

      ENDIF

      IF(RUN_TYPE == 'RESTART_2') VTP_FINDEX=0

! setting additional particle properties now that the particles
! have been identified
      DO L = 1, MAX_PIP
! Skip 'empty' locations when populating the particle property arrays.
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_GHOST(L) .OR. IS_ENTERING_GHOST(L) .OR. IS_EXITING_GHOST(L)) CYCLE
         PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L)
         OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
      ENDDO

      CALL SET_PHASE_INDEX
      CALL INIT_PARTICLES_IN_CELL

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.
! Bin the particles to the DES grid.
      CALL DESGRID_PIC(PLOCATE=.TRUE.)
      CALL DES_PAR_EXCHANGE
      CALL PARTICLES_IN_CELL

      IF(DEM_SOLIDS) THEN
         CALL NEIGHBOUR
         CALL INIT_SETTLING_DEM
      ENDIF

      IF(RUN_TYPE == 'NEW') CALL SET_IC_DEM


! Calculate interpolation weights
      CALL CALC_INTERP_WEIGHTS
! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS

      IF(MPPIC) CALL CALC_DTPIC


      IF(RUN_TYPE /= 'RESTART_1' .AND. PRINT_DES_DATA) THEN
         S_TIME = TIME
         CALL WRITE_DES_DATA
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES
