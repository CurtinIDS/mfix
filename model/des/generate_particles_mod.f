      MODULE GENERATE_PARTICLES

        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PARTICLE_COUNT

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 19-Mar-14  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GENERATE_PARTICLE_CONFIG

      use mfix_pic, only: MPPIC
      use discretelement, only: PIP, PARTICLES
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! Parameter for detecting unspecified values, zero, and one
      use param1, only: ONE
! Maximum number of initial conditions
      use param, only: DIMENSION_IC
! IC Region gas volume fraction.
      use ic, only: IC_EP_G

      use mpi_utility

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

      INTEGER :: ICV

! Initialize the error manager.
      CALL INIT_ERR_MSG("Generate_Particle_Config")

      DO ICV = 1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE
         IF(IC_EP_G(ICV) == ONE) CYCLE

         IF(MPPIC) THEN
            CALL GENERATE_PARTICLE_CONFIG_MPPIC(ICV)
         ELSE
            CALL GENERATE_PARTICLE_CONFIG_DEM(ICV)
         ENDIF

      ENDDO

      CALL GLOBAL_SUM(PIP,PARTICLES)

      WRITE(ERR_MSG, 1004) PARTICLES
 1004 FORMAT(/,'Total number of particles in the system: ',I15)

      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GENERATE_PARTICLE_CONFIG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
!  Authors: Rahul Garg                               Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: Generate particle configuration for DEM solids for each IC !
!           region. Now using the particle linked lists for initial    !
!           build                                                      !
!           This routine will ultimately supersede the older rouine    !
!           that has not been deleted yet                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM(ICV)

! Global Variables:
!---------------------------------------------------------------------//
! particle radius and density
      use discretelement, only: DES_RADIUS, RO_Sol
! particle position new and old
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! particle velocity new and old
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Simulation dimension (2D/3D)
      use discretelement, only: DIMN
! Number of particles in the system (current)
      use discretelement, only: PIP
! Number of DEM solids phases.
      use discretelement, only: DES_MMAX
! Flag to use _OLD variables
      use discretelement, only: DO_OLD
! Angular velocity
      use discretelement, only: OMEGA_OLD, OMEGA_NEW, PIJK
! solid phase diameters and densities.
      use physprop, only: D_p0, RO_s0, MMAX
! IC Region solids volume fraction.
      use ic, only: IC_EP_S

! Constant: 3.14159...
      use constant, only: PI
! min and max physical co-ordinates of IC regions in each direction
      use ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initally specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
! Flag to extend the lattice distribution in a given IC to available area
      use ic, only: IC_DES_FIT_TO_REGION
! Parameter for detecting unspecified values, zero, and one
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Parameter for small and large numbers
      use param1, only: SMALL_NUMBER, LARGE_NUMBER

! to access random number generator subroutines
      use randomno
      use mpi_utility
      use functions, only: SET_NORMAL

      use desgrid, only: dg_xstart, dg_ystart, dg_zstart
      use desgrid, only: dg_xend, dg_yend, dg_zend

! direction wise spans of the domain and grid spacing in each direction
      use geometry, only: xlength, ylength, zlength

      use cutcell, only : CARTESIAN_GRID
      use stl_functions_des, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW

      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_OWNs
      use toleranc, only: compare

      use discretelement, only: max_pip, max_radius, xe, yn, zt
      use error_manager
      use functions
      use param, only: dim_m
      use param, only: dimension_i, dimension_j, dimension_k

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Starting positions in the axial directions
      DOUBLE PRECISION :: xINIT, yINIT, zINIT
! Fractor used to scale particle diameter
      DOUBLE PRECISION :: lFAC
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Number of particles in the lattice
      INTEGER :: SEED_X, SEED_Y, SEED_Z
! Loop indices phase/fluid cell
      INTEGER :: M, MM, I, J, K, IJK, LB, UB
! Loop indicies for seeding
      INTEGER :: II, JJ, KK
! Start and end bound for IC region.
      DOUBLE PRECISION :: IC_START(3), IC_END(3)
! Volume and lengths of the IC Region
      DOUBLE PRECISION :: DOM_VOL, DOML(3)
! Flag to skip the particle
      LOGICAL :: SKIP
! Diameter adjusted for space padding
      DOUBLE PRECISION :: ADJ_DIA
! Number of particles calculated from volume fracton
      INTEGER :: rPARTS(DIM_M), tPARTS
! Spacing between particles.
      DOUBLE PRECISION :: lDEL, lDX, lDY, lDZ
! Flag that the setup failed to fit the particles to the IC region
      LOGICAL :: FIT_FAILED
! Number of seeded particles
      INTEGER :: pCOUNT(DIM_M), tCOUNT

      DOUBLE PRECISION :: SOLIDS_DATA(0:DIM_M)

      LOGICAL :: VEL_FLUCT
      DOUBLE PRECISION :: VEL_SIG
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)

      logical :: report = .true.
      logical :: found

!......................................................................!

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_DEM")

      WRITE(ERR_MSG,"(2/,'Generating initial particle configuration:')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))

! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lFAC = 1.05D0

! Setup local arrays with IC region bounds.
      IC_START(1)=IC_X_W(ICV);   IC_END(1)=IC_X_E(ICV)
      IC_START(2)=IC_Y_S(ICV);   IC_END(2)=IC_Y_N(ICV)
      IC_START(3)=IC_Z_B(ICV);   IC_END(3)=IC_Z_T(ICV)

      DOML = IC_END-IC_START
      IF(NO_K) DOML(3)=DZ(1)

! Volume of the IC region
      DOM_VOL = DOML(1)*DOML(2)*DOML(3)

      rPARTS=0
      DO M=MMAX+1,MMAX+DES_MMAX
         IF(SOLIDS_MODEL(M) == 'DEM') THEN
! Number of particles for phase M
            rPARTS(M) = &
               floor((6.0d0*IC_EP_S(ICV,M)*DOM_VOL)/(PI*(D_P0(M)**3)))
         ENDIF
      ENDDO

! Total number of particles in this IC region.
      tPARTS = sum(rPARTS)
      IF(tPARTS == 0) RETURN

      ADJ_DIA = 2.0d0*MAX_RADIUS*lFAC

! Attempt to seed particle throughout the IC region
      FIT_FAILED=.FALSE.
      IF(IC_DES_FIT_TO_REGION(ICV)) THEN
         IF(NO_K) THEN
            lDEL = (DOML(1)-ADJ_DIA)*(DOML(2)-ADJ_DIA)
            lDEL = (lDEL/dble(tPARTS))**(1.0/2.0)
            SEED_X = max(1,ceiling((DOML(1)-ADJ_DIA)/lDEL))
            SEED_Y = max(1,ceiling((DOML(2)-ADJ_DIA)/lDEL))
            SEED_Z = 1
         ELSE
            lDEL = (DOML(1)-ADJ_DIA)*(DOML(2)-ADJ_DIA)*(DOML(3)-ADJ_DIA)
            lDEL = (lDEL/dble(tPARTS))**(1.0/3.0)
            SEED_X = max(1,ceiling((DOML(1)-ADJ_DIA)/lDEL))
            SEED_Y = max(1,ceiling((DOML(2)-ADJ_DIA)/lDEL))
            SEED_Z = max(1,ceiling((DOML(3)-ADJ_DIA)/lDEL))
         ENDIF
         FIT_FAILED=(dble(SEED_X*SEED_Y*SEED_Z) < tPARTS)
      ENDIF

! Generic filling
      IF(.NOT.IC_DES_FIT_TO_REGION(ICV) .OR. FIT_FAILED) THEN
         SEED_X = max(1,floor((IC_END(1)-IC_START(1)-ADJ_DIA)/ADJ_DIA))
         SEED_Y = max(1,floor((IC_END(2)-IC_START(2)-ADJ_DIA)/ADJ_DIA))
         SEED_Z = max(1,floor((IC_END(3)-IC_START(3)-ADJ_DIA)/ADJ_DIA))
      ENDIF

      lDX = DOML(1)/dble(SEED_X)
      lDY = DOML(2)/dble(SEED_Y)
      IF(DO_K) THEN
         lDZ = DOML(3)/dble(SEED_Z)
      ELSE
         lDZ = 0.0d0
      ENDIF

      xINIT = IC_START(1)+HALF*lDX
      yINIT = IC_START(2)+HALF*lDY
      zINIT = IC_START(3)+HALF*lDZ

      M=1
      pCOUNT = 0
      tCOUNT = 0

      VEL_FLUCT = SET_VEL_FLUCT(ICV,M)

      JJ_LP: DO JJ=1, SEED_Y
         POS(2) = YINIT + (JJ-1)*lDY
         IF(compare(POS(2),dg_ystart) .OR. compare(POS(2),dg_yend))    &
            POS(2) = POS(2) + SMALL_NUMBER

      KK_LP: DO KK=1, SEED_Z
         POS(3) = ZINIT + (KK-1)*lDZ
         IF(DO_K) THEN
            IF(compare(POS(3),dg_zstart) .OR. compare(POS(3),dg_zend)) &
               POS(3) = POS(3) + SMALL_NUMBER
         ENDIF

      II_LP: DO II=1, SEED_X
         POS(1) = xINIT + (II-1)*lDX
         IF(compare(POS(1),dg_xstart) .OR. compare(POS(1),dg_xend))    &
            POS(1) = POS(1) + SMALL_NUMBER

! Exit if all particles were seeded.
         IF(tCOUNT > int(tPARTS)) THEN
            EXIT JJ_LP
! Find the next phase that needs to be seeded
         ELSEIF(pCOUNT(M) > int(rPARTS(M))) THEN
            MM_LP: DO MM=M+1,MMAX+DES_MMAX
               IF(rPARTS(MM) > 0.0) THEN
                  M=MM
                  EXIT MM_LP
               ENDIF
            ENDDO MM_LP
            IF(M > MMAX+DES_MMAX) EXIT JJ_LP
            VEL_FLUCT = SET_VEL_FLUCT(ICV,M)
         ENDIF

         pCOUNT(M) = pCOUNT(M) + 1
         tCOUNT = tCOUNT + 1

! Keep only particles that belong to this process.
         IF(.NOT.dg_is_ON_myPE_OWNs(IofPOS(POS(1)), &
            JofPOS(POS(2)),KofPOS(POS(3)))) CYCLE

! Bin the parcel to the fuild grid.
         K=1
         IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)
         CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
         CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)

! Skip cells that return invalid IJKs.
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

! Skip cells that are not part of the local fuild domain.
         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF(CARTESIAN_GRID) THEN
            CALL CHECK_IF_PARTICLE_OVERLAPS_STL(POS, I, J, K, SKIP)
            IF(SKIP) CYCLE
         ENDIF

         PIP = PIP + 1
         CALL PARTICLE_GROW(PIP)
         MAX_PIP = max(PIP,MAX_PIP)

         CALL SET_NORMAL(PIP)

         IF(VEL_FLUCT) THEN
            VEL(1) = randVEL(pCOUNT(M),1)
            VEL(2) = randVEL(pCOUNT(M),2)
            VEL(3) = randVEL(pCOUNT(M),3)
         ELSE
            VEL(1) = IC_U_s(ICV,M)
            VEL(2) = IC_V_s(ICV,M)
            VEL(3) = IC_W_s(ICV,M)
         ENDIF
         IF(NO_K) VEL(3) = 0.0d0


         DES_POS_NEW(PIP,:) = POS(:)
         DES_VEL_NEW(PIP,:) = VEL(:)
         OMEGA_NEW(PIP,:) = 0.0d0

         DES_RADIUS(PIP) = D_P0(M)*HALF
         RO_SOL(PIP) =  RO_S0(M)

         PIJK(PIP,1) = I
         PIJK(PIP,2) = J
         PIJK(PIP,3) = K
         PIJK(PIP,4) = IJK
         PIJK(PIP,5) = M

         IF(DO_OLD) THEN
            DES_VEL_OLD(PIP,:) = DES_VEL_NEW(PIP,:)
            DES_POS_OLD(PIP,:) = DES_POS_NEW(PIP,:)
            OMEGA_OLD(PIP,:) = ZERO
         ENDIF

         SOLIDS_DATA(M) = SOLIDS_DATA(M) + 1.0

      ENDDO II_LP
      ENDDO KK_LP
      ENDDO JJ_LP

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

1000 FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=',&
         ES15.4,/'Please correct the mfix.dat file.')

      WRITE(ERR_MSG,2000) ICV
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      DO M=MMAX+1, MMAX+DES_MMAX
         IF(SOLIDS_DATA(M) < SMALL_NUMBER) CYCLE
         WRITE(ERR_MSG,2010) M, int(SOLIDS_DATA(M)), IC_EP_S(ICV,M),   &
            (dble(SOLIDS_DATA(M))*(Pi/6.0d0)*D_P0(M)**3)/SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      IF(allocated(randVEL)) deallocate(randVEL)

      CALL FINL_ERR_MSG

      RETURN

 2000 FORMAT(/2x,'|',43('-'),'|',/2x,'| IC Region: ',I3,28x,'|',/2x,   &
         '|',43('-'),'|',/2x,'| Phase | Number of |    EPs    |    EP',&
         's    |',/2x,'|   ID  | Particles | Specified |   Actual  |', &
         /2x,'|-------|',3(11('-'),'|'))

 2010 FORMAT(2x,'|  ',I3,'  |',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
         '|-------|',3(11('-'),'|'))


      CONTAINS

!......................................................................!
! Function: SET_VEL_FLUCT                                              !
! Purpose: Set the flag for random velocity fluctuations. If needed    !
! the random velocities are calculated.                                !
!......................................................................!
      LOGICAL FUNCTION SET_VEL_FLUCT(lICV, lM)
      INTEGER, INTENT(IN) :: lICV, lM
      DOUBLE PRECISION :: VEL_SIG
      SET_VEL_FLUCT=(IC_Theta_M(lICV,lM) /= 0.0d0)
      IF(SET_VEL_FLUCT) THEN
         if(allocated(randVEL)) deallocate(randVEL)
         allocate(randVEL(100+int(rPARTS(lM)),3))
         VEL_SIG = sqrt(IC_Theta_M(lICV,lM))
         CALL NOR_RNO(randVEL(:,1), IC_U_s(lICV,lM),VEL_SIG)
         CALL NOR_RNO(randVEL(:,2), IC_V_s(lICV,lM),VEL_SIG)
         IF(DO_K) CALL NOR_RNO(randVEL(:,3),IC_W_s(lICV,lM),VEL_SIG)
      ENDIF

      END FUNCTION SET_VEL_FLUCT

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         !
!  Author: Rahul Garg                                 Date: 3-May-2011 !
!                                                                      !
!  Purpose: Generates particle position distribution for MPPIC.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC(ICV)

! Global variables
!---------------------------------------------------------------------//
! Number of DES solids phases.
      use discretelement, only: DES_MMAX
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_ROP_s
! IC Region gas volume fraction.
      use ic, only: IC_EP_G
! MPPIC specific IC region specification.
      use ic, only: IC_PIC_CONST_NPC, IC_PIC_CONST_STATWT

      use param1, only: UNDEFINED, UNDEFINED_I
      use param1, only: ZERO, ONE, HALF
! Maximum number of IC regions and solids phases
      use param, only: DIMENSION_IC
      use param, only: DIM_M
      use physprop, only: mmax

! The accumulated number of particles in each IJK.
      use mpi_utility, only: GLOBAL_ALL_SUM

      use error_manager

      use run, only: solids_model
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Generic loop counters
      INTEGER :: M
! Actual volume of IC region
      DOUBLE PRECISION :: IC_VOL
! Solids data in IC Region by phase:
      DOUBLE PRECISION :: SOLIDS_DATA(0:4*DIM_M)
!......................................................................!

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_MPPIC")

      WRITE(ERR_MSG,"(2/,'Generating initial parcel configuration:')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)


      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))

      WRITE(ERR_MSG,2000) ICV
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

! Set up the individual solids phases.
      DO M=MMAX+1, DES_MMAX+MMAX
         IF(SOLIDS_MODEL(M) == 'PIC') THEN
            IF(IC_ROP_S(ICV,M) == ZERO) CYCLE
! Seed parcels with constant stastical weight.
            CALL GPC_MPPIC_CONST_NPC(ICV, M, SOLIDS_DATA(0), &
            SOLIDS_DATA((4*M-3):(4*M)))
         ENDIF
      ENDDO

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=',&
         ES15.4,/'Please correct the mfix.dat file.')

! Report solids information for the IC region.
      DO M=MMAX+1, DES_MMAX+MMAX
         IF(SOLIDS_MODEL(M) == 'PIC') THEN
            WRITE(ERR_MSG,2010) M, int(SOLIDS_DATA(4*M-3)),&
               int(SOLIDS_DATA(4*M-3)*SOLIDS_DATA(4*M-2)), &
               SOLIDS_DATA(4*M-2), SOLIDS_DATA(4*M-1),     &
               SOLIDS_DATA(4*M)/SOLIDS_DATA(0)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 2000 FORMAT(/2x,'|',67('-'),'|',/2x,'| IC Region: ',I3,52x,'|',/2x,   &
         '|',67('-'),'|',/2x,'| Phase | Num. Comp | Num. Real ',       &
         '| Stastical |    EPs    |    EPs    |',/2x,'|   ID  |  ',    &
         'Parcels  | Particles |   Weight  | Specified |   Actual  |', &
         /2x,'|-------|',5(11('-'),'|'))

 2010 FORMAT(2x,'|  ',I3,'  |',2(1x,I9,1x,'|'),3(1x,ES9.2,1x,'|'),/2x,&
         '|-------|',5(11('-'),'|'))

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         !
!  Author: Rahul Garg                                 Date: 3-May-2011 !
!                                                                      !
!  Purpose: generates particle position distribution for MPPIC         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GPC_MPPIC_CONST_NPC(ICV, M, IC_VOL, sDATA)

! Global variables
!---------------------------------------------------------------------//
! Constant: 3.14159...
      use constant, only: PI
! Cut_cell identifier array
      use cutcell, only: cut_cell_at
      use cutcell, only : CARTESIAN_GRID
      use discretelement, only: XE, YN, ZT

      use des_allocate, only: PARTICLE_GROW
      use discretelement

! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_EP_s
      use ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

! initally specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
      use ic, only: IC_PIC_CONST_NPC
      use ic, only: IC_PIC_CONST_STATWT

      use mfix_pic, only: des_stat_wt
      use mpi_utility

! Maximum number of IC regions and solids phases
      use param, only: DIMENSION_IC
      use param, only: DIM_M
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, HALF

! solid phase diameters and densities.
      use physprop, only: D_p0, RO_s0
      use run, only: solids_model

      use randomno
      use error_manager
      use functions
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW

      IMPLICIT NONE

! Dummy Arguments
!----------------------------------------------------------------------//
! Index of IC region and solids phase
      INTEGER, INTENT(IN) :: ICV, M
! Specific volume of IC region (accounts for blocked cells)
      DOUBLE PRECISION, INTENT(IN) :: IC_VOL
! Data about solids in the IC region.
      DOUBLE PRECISION, INTENT(OUT) :: sDATA(4)

! Local variables
!----------------------------------------------------------------------//
      DOUBLE PRECISION :: EP_SM

! Number of real and comp. particles in a cell.
      DOUBLE PRECISION ::  rPARTS
      INTEGER :: maxPARTS
      DOUBLE PRECISION :: DOML(3), IC_START(3)
! Parcel position with random
      DOUBLE PRECISION :: POS(3)
! Average velocity and standard deivation
      DOUBLE PRECISION :: IC_VEL(3), VEL_SIG
! Flag to not keep parcel.
      LOGICAL :: SKIP
! Arrasy for assigning random position and velocities
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)
      DOUBLE PRECISION :: RAND(3)
! Statistical weights
      DOUBLE PRECISION :: STAT_WT, SUM_STAT_WT
! Volume of a parcel and total solids volume
      DOUBLE PRECISION :: sVOL, sVOL_TOT
! Counter for seeded parcles.
      INTEGER :: SEEDED
! Generic loop indices and loop counters
      INTEGER :: I, J, K, IJK, LC, LC_MAX
!......................................................................!

      CALL INIT_ERR_MSG("GPC_MPPIC_CONST_NPC")

      maxPARTS=25
      allocate(randVEL(maxPARTS,3))

      IC_VEL(1) = IC_U_s(ICV,M)
      IC_VEL(2) = IC_V_s(ICV,M)
      IC_VEL(3) = merge(IC_W_s(ICV,M),0.0d0,DO_K)

      VEL_SIG = sqrt(IC_Theta_M(ICV,M))

! Volume occupied by one particle
      sVOL = (Pi/6.0d0)*(D_P0(M)**3.d0)

      SEEDED = 0
      sVOL_TOT = 0.0d0
      SUM_STAT_WT = 0.0d0

      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle
         IF(DEAD_CELL_AT(I,J,K)) cycle

         IJK = FUNIJK(I,J,K)
         IF(.not.FLUID_AT(IJK)) cycle

         rPARTS = IC_EP_s(ICV,M)*VOL(IJK)/sVOL

! Seed parcels with a constant stastical weight
         IF(IC_PIC_CONST_STATWT(ICV,M) /= ZERO) THEN
            STAT_WT = IC_PIC_CONST_STATWT(ICV,M)
            LC_MAX = int(rPARTS/STAT_WT)

! Seed parcels with a constant number per cell
         ELSEIF(IC_PIC_CONST_NPC(ICV,M) /= 0) THEN
            LC_MAX = IC_PIC_CONST_NPC(ICV,M)
            STAT_WT = rPARTS/dble(LC_MAX)
            IF(CUT_CELL_AT(IJK)) LC_MAX = max(1,int(VOL(IJK)*dble(&
               IC_PIC_CONST_NPC(ICV,M))/(DX(I)*DY(J)*DZ(K))))
         ENDIF

! Increase particle buffer
         IF(LC_MAX > maxPARTS) THEN
            maxPARTS = 2*LC_MAX
            if(allocated(randVEL)) deallocate(randVEL)
            allocate(randVEL(maxPARTS,3))
         ENDIF

         DO LC=1, merge(2,3,NO_K)
            IF(VEL_SIG > ZERO) THEN
               CALL NOR_RNO(randVEL(1:LC_MAX,LC), IC_VEL(LC), VEL_SIG)
            ELSE
               randVEL(1:LC_MAX,LC) = IC_VEL(LC)
            ENDIF
         ENDDO
         IF(NO_K) randVEL(1:LC_MAX,3) = 0.0d0

         IC_START(1) = XE(I-1)
         IC_START(2) = YN(J-1)
         IC_START(3) = ZERO;  IF(DO_K) IC_START(3) = ZT(K-1)

         DOML(1) = DX(I)
         DOML(2) = DY(J)
         DOML(3) = ZERO;  IF(DO_K) DOML(3) = DZ(K)

         DO LC=1,LC_MAX

            CALL RANDOM_NUMBER(RAND)
            POS(:) = IC_START + DOML*RAND

            IF(CARTESIAN_GRID) THEN
               CALL CHECK_IF_PARCEL_OVERLAPS_STL(POS, SKIP)
               DO WHILE(SKIP)
                  CALL RANDOM_NUMBER(RAND)
                  POS(:) = IC_START + DOML*RAND
                  CALL CHECK_IF_PARCEL_OVERLAPS_STL(POS, SKIP)
               ENDDO
            ENDIF

            PIP = PIP + 1
            CALL PARTICLE_GROW(PIP)
            MAX_PIP = max(PIP,MAX_PIP)

            DES_POS_NEW(PIP,:) = POS(:)
            DES_VEL_NEW(PIP,:) = randVEL(LC,:)

            DES_RADIUS(PIP) = D_P0(M)*HALF
            RO_SOL(PIP) =  RO_S0(M)

            PIJK(PIP,1) = I
            PIJK(PIP,2) = J
            PIJK(PIP,3) = K
            PIJK(PIP,4) = IJK
            PIJK(PIP,5) = M

            SUM_STAT_WT = SUM_STAT_WT + STAT_WT

            DES_STAT_WT(PIP) = STAT_WT
            sVOL_TOT = sVOL_TOT + sVOL*STAT_WT

            CALL SET_NORMAL(PIP)

            SEEDED = SEEDED + 1

         ENDDO

      ENDDO
      ENDDO
      ENDDO

      IF(allocated(randVEL)) deallocate(randVEL)

      sDATA(1) = dble(SEEDED)
      sDATA(2) = SUM_STAT_WT/dble(SEEDED)
      sDATA(3) = IC_EP_S(ICV,M)
      sDATA(4) = sVOL_TOT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GPC_MPPIC_CONST_NPC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GET_IC_VOLUME                                           !
!  Author: J.Musser                                 Date: 26-Aug-2015  !
!                                                                      !
!  Purpose: Calculate the actual volume of the IC region.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_IC_VOLUME(ICV, IC_VOL)

! IC region index bounds
      use ic, only: IC_I_w, IC_I_e
      use ic, only: IC_J_s, IC_J_n
      use ic, only: IC_K_b, IC_K_t

! Volume of computational cells.
      use geometry, only: VOL

      use functions
      use compar, only: dead_cell_at

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Index of IC region
      INTEGER, INTENT(IN) :: ICV
! Total calculated volume of IC region
      DOUBLE PRECISION, INTENT(OUT) :: IC_VOL

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK
!......................................................................!


      IC_VOL = 0.0d0
      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I,J,K)
         IF(FLUID_AT(IJK)) IC_VOL = IC_VOL + VOL(IJK)

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_IC_VOLUME

      END MODULE GENERATE_PARTICLES
