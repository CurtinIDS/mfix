!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: MFIX_PIC                                                    !
!  Purpose: MP-PIC related data                                        !
!  Author: R. Garg                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE MFIX_PIC

      use param, only: dim_m
!......................................................................!

! flag for turning on MP-PIC method
      LOGICAL :: MPPIC

      DOUBLE PRECISION :: PSFAC_FRIC_PIC
      DOUBLE PRECISION :: FRIC_EXP_PIC
      DOUBLE PRECISION :: FRIC_NON_SING_FAC

      LOGICAL :: MPPIC_SOLID_STRESS_SNIDER
      LOGICAL :: MPPIC_CORR_VOLFRAC
      INTEGER :: ITER_VOL_FRAC_CORR

! force to solids pressure gradient
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PS_FORCE_PIC

! avg solid velocity at particle position (used for MP-PIC)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AVGSOLVEL_P

! EP_g interpolated at particle position (used for MP-PIC)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EPG_P

! to initially seed the particles based on constant number of
! particles per cell or constant statistical weight of each particle
! number of particles per cell for the case of CONSTANTNPC
      LOGICAL :: MPPIC_CONSTANTNPC, MPPIC_CONSTANTWT

      INTEGER NPC_PIC(DIM_M)

! coefficeient of restituion used in MPPIC case in the
! frictional regime
      DOUBLE PRECISION :: MPPIC_COEFF_EN1, MPPIC_COEFF_EN2
      DOUBLE PRECISION :: MPPIC_COEFF_EN_WALL, MPPIC_COEFF_ET_WALL

! statistical weight or number of real particles per computational
! particle for the case of CONSTANTWT
      DOUBLE PRECISION STATWT_PIC(DIM_M)

! # of computational particles for the entire grid
! keep it double precision for inlow BC
      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE :: CNP_ARRAY

! Statistical weight of each particle
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_STAT_WT

!     logical variable to decide if special treatment is needed or not
!     in the direction of gravity in the frictional stress tensor
      LOGICAL :: MPPIC_GRAV_TREATMENT

! Particle response time scale for each phase
      DOUBLE PRECISION :: DES_TAU_P(DIM_M)

! The maximum dt for point particles based on particle response time
! (taup) and cfl. See cfassign for its computation
      DOUBLE PRECISION :: DTPIC_MAX

! CFL value for point-particles that is used to control DTSOLID
! DTPP_CFL, dt for point-particles based on user specified CFL_PP
! and maximum velocity of particles
      DOUBLE PRECISION CFL_PIC, DTPIC_CFL, DTPIC_TAUP

      DOUBLE PRECISION :: DTSOLID_ORIG
! solid pressure gradient
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PS_GRAD

! flag to turn on implicit treatment of drag force term in particle
! trajectory evolution equation
      LOGICAL MPPIC_PDRAG_IMPLICIT

! Solids pressure as a result of granular motion
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: PIC_P_s
!Face centered u-velocity required by PIC model.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PIC_U_S
! Face centered v-velocity required by PIC model. Not using the U_s arrays anymore
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PIC_V_S
! Face centered z-velocity required by PIC model. Not using the U_s arrays anymore
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PIC_W_S

! A run time flag to report minimum value of gas voidage
      LOGICAL :: PIC_REPORT_MIN_EPG

      END MODULE MFIX_PIC

