  MODULE ghdtheory


      Use param
      Use param1

!
!     Zeroth order dissipation term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zeta0
!
!     cooling rate transport coefficient (1st order)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZetaU
!
!     Thermal diffusivity DiT
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DiT
!
!     Mass mobility coefficient
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: DijF
!
!     Thermal mobility
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Lij
!
!     Ordinary diffusion
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Dij
!
!     Dufour coefficient
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: DijQ
!
!     Species mass flux in X-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiX
!
!     Species mass flux in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiY
!
!     Species mass flux in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiZ
!
!     external force in X-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiX
!
!     external force in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiY
!
!     external force in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiZ

!     external force in X-direction--flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiXvel
!
!     external force in Y-direction--flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiYvel
!
!     external force in Z-direction--flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiZvel
!
!     external force Minus Drag term in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiMinusDragX
!
!     Species mass flux Without Drag term in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiMinusDragX
!
!     external force Minus Drag term in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiMinusDragY
!
!     Species mass flux Without Drag term in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiMinusDragY
!
!     external force Minus Drag term in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: FiMinusDragZ
!
!     Species mass flux Without Drag term in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiMinusDragZ


!     Everything that does not depend on the velocity of ith particle
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DELTAU
!
!     Everything that does not depend on the velocity of the ith particle
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DELTAV
!
!     Everything that does not depned on the velocity of the ith particle
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DELTAW

!     Everything that does not depend on the velocity of ith particle
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DEL_DOT_J

!     drag force in X-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dragFx
!
!     drag force in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dragFy
!
!     drag force in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dragFz
!     drag force in X-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dragFxflux
!
!     drag force in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dragFyflux
!
!     drag force in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dragFzflux
!     drag force in X-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: beta_cell_X
!
!     drag force in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: beta_cell_Y
!
!     drag force in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: beta_cell_Z
!     drag force in X-direction
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: beta_ij_cell_X
!
!     drag force in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: beta_ij_cell_Y
!
!     drag force in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: beta_ij_cell_Z

!       decide whether to do harmonic or arithmetic thermal diffusivity east
       LOGICAL, DIMENSION(:), ALLOCATABLE :: DiT_HarmE
!       decide whether to do harmonic or arithmetic thermal diffusivity north
       LOGICAL, DIMENSION(:), ALLOCATABLE :: DiT_HarmN
!       decide whether to do harmonic or arithmetic thermal diffusivity top
       LOGICAL, DIMENSION(:), ALLOCATABLE :: DiT_HarmT

!       decide whether to do harmonic or arithmetic mass mobility east
       LOGICAL, DIMENSION(:,:), ALLOCATABLE :: DijF_HarmE
!       decide whether to do harmonic or arithmetic mass mobility north
       LOGICAL, DIMENSION(:,:), ALLOCATABLE :: DijF_HarmN
!       decide whether to do harmonic or arithmetic mass mobility top
       LOGICAL, DIMENSION(:,:), ALLOCATABLE :: DijF_HarmT

!       decide whether to do harmonic or arithmetic ordinary diffusivity east
       LOGICAL, DIMENSION(:,:), ALLOCATABLE :: Dij_HarmE
!       decide whether to do harmonic or arithmetic ordinary diffusivity north
       LOGICAL, DIMENSION(:,:), ALLOCATABLE :: Dij_HarmN
!       decide whether to do harmonic or arithmetic ordinary diffusivity top
       LOGICAL, DIMENSION(:,:), ALLOCATABLE :: Dij_HarmT

      END MODULE ghdtheory
