      MODULE rxns

        use param, only: dim_m, dim_n_all, dim_n_g, dim_n_s
        Use rxn_com, only: reaction_block

! reaction rates
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ReactionRates

! number of ReactionRates
      INTEGER nRR
! total number of species
      INTEGER N_all

      LOGICAL rDatabase(0:DIM_M, DIM_N_g)

!-----------------------------------------------------------------------

! Indicates that reaction rates are to be calculated.
      LOGICAL :: RRATE
! Indicates if the legacy reaction rates file (rrates.f) is used.
      LOGICAL :: USE_RRATES

! Rate of production of gas species
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R_gp
! Rate of consumption of gas species/X_g
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RoX_gc
! Net production of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SUM_R_g

! Rate of production of solids species
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: R_sp
! Rate of consumption of solids species/X_s
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: RoX_sc
! Net production of solids
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  SUM_R_s

! Rate of mass transfer from phase M to Phase L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  R_phase

! Actual number of Reactions
      INTEGER NO_OF_RXNS

! Species names -- Legacy variable
      CHARACTER(len=18) SPECIES_NAME(DIM_N_ALL)

! Gas phase species names (database) and aliases
      CHARACTER(len=18) SPECIES_g(DIM_N_g) ! database name
      CHARACTER(len=32)  SPECIES_ALIAS_g(DIM_N_g) ! alias

! Solids phase species names (database) and aliases
      CHARACTER(len=18) SPECIES_s(DIM_M, DIM_N_s) ! database name
      CHARACTER(len=32)  SPECIES_ALIAS_s(DIM_M, DIM_N_s) ! alias

! Array linking all of the reaction data.
      TYPE(REACTION_BLOCK), DIMENSION(:), TARGET, ALLOCATABLE :: Reaction

      END MODULE rxns
