!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO                                             !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM heat transfer.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_THERMO

      USE param, only: dim_m

! Heat transfer correlation specified in mfix.dat
! Default [RANZ_1952]
      CHARACTER(LEN=24) :: DES_CONV_CORR

      INTEGER :: DES_CONV_CORR_ENUM
      INTEGER, PARAMETER :: RANZ_1952 = 0

! Run time flags for calculating the various modes of heat transfer
      LOGICAL :: CALC_CONV_DES = .FALSE.
      LOGICAL :: CALC_COND_DES(DIM_M) = .FALSE.
      LOGICAL :: CALC_RADT_DES(DIM_M) = .FALSE.

! Particle properties
!-----------------------------------------------------------------------
! Particle temperature
      DOUBLE PRECISION, ALLOCATABLE :: DES_T_s(:)

! DES specific heat of particles by particle
      DOUBLE PRECISION, ALLOCATABLE :: DES_C_ps(:)

! Emissivity of particles
      DOUBLE PRECISION :: DES_Em(DIM_M)
! Stefan-Boltzmann Constant
      DOUBLE PRECISION :: SB_CONST
! Bulk solids temperature for radiative heat transfer
      DOUBLE PRECISION, ALLOCATABLE :: avgDES_T_s(:)

! Convective heat transfer coefficient TIMES particle surface area
      DOUBLE PRECISION, ALLOCATABLE :: GAMMAxSA(:)


! Thermodynamic Neighborhood
!-----------------------------------------------------------------------
! Fluid Lens Proportion Constant used to calculate the radius of the
! fluid lens that surrounds the particle for particle-fluid-particle
! conduction.  Default [ 0.2 ]
      DOUBLE PRECISION :: FLPC

! Mininum separation distance between the surface of two contacting
! particles. This value is used to remove the singluarity that the
! particle-fluid-particle conduciton model develops at the contact
! interface. [4.0x10^(-10) meters]
      DOUBLE PRECISION :: DES_MIN_COND_DIST



! Rates of heat transfer
!-----------------------------------------------------------------------
! Generic heat transfer source to the particle.
      DOUBLE PRECISION, ALLOCATABLE :: Q_Source(:), Q_Source0(:)
! Convective heat transfer source
      DOUBLE PRECISION, ALLOCATABLE :: CONV_Qs(:)
! Heat source resulting from chemical reactions
      DOUBLE PRECISION, ALLOCATABLE :: RXNS_Qs(:)


! Fluid/Particle coupling
!---------------------------------------------------------------------//
! Gas Phase Energy Eq source terms for gas-particle convection
      DOUBLE PRECISION, ALLOCATABLE :: CONV_Sc(:), CONV_Sp(:)

      contains


!----------------------------------------------------------------------!
! Function: CALC_Cp_DES                                                !
! Author: J.Musser                                    Date: 11-DEC-15  !
!                                                                      !
! Purpose: Calculate the specific heat of a particle.                  !
!----------------------------------------------------------------------!
      PURE DOUBLE PRECISION FUNCTION CALC_Cp_DES(pNP)

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Calculate the specific heat from polynomical data obtained from the
! thermodynamic databases.
      use read_thermochemical, only: calc_CpoR
! Number of species comprising species
      use physprop, only: NMAX_s
      use physprop, only: C_PS0
      use physprop, only: MW_s
      use discretelement, only: PIJK
      use des_rxns, only: DES_X_s
      use run, only: UNITS

      use param1, only: ZERO, UNDEFINED

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: pNP
!     DOUBLE PRECISION, INTENT(IN) :: pTs
! Local Variables
!---------------------------------------------------------------------//
! loop counter and error indicator
      INTEGER :: NN, IER, MM
      DOUBLE PRECISION :: lTs
!......................................................................!

! Temperature of particle
      lTs = DES_T_s(pNP)
! Phase of particle.
      MM=PIJK(pNP,5)

! Calculate the specific heat based on the species composition of the
! particle and the data from the thermodynamic databases.
      IF(C_PS0(MM) == UNDEFINED) THEN
         CALC_Cp_DES = ZERO
         DO NN = 1, NMAX_s(MM)
            CALC_Cp_DES = CALC_Cp_DES + calc_CpoR(lTs, MM, NN) *  &
               DES_X_s(pNP,NN)*RGAS / MW_s(MM,NN)
         ENDDO
! Convert to SI units if needed.
         IF (UNITS == 'SI') CALC_Cp_DES = 4183.925D0*CALC_Cp_DES
      ELSE
! If a constant value specific heat has been assigned to the particle
! in the mfix.dat file, use this value.
         CALC_CP_DES = C_PS0(MM)
      ENDIF
      END FUNCTION CALC_Cp_DES

      END MODULE DES_THERMO
