
      MODULE scalars

      use param, only: dim_scalar

! Number of scalar equations solved
      INTEGER :: NScalar

! Index of phase associated with scalar n
      INTEGER, DIMENSION(1:DIM_Scalar) :: Phase4Scalar

! Source term for User-defined Scalars is linearized as
! S = Scalar_c - Scalar_p * Scalar
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Scalar_c
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Scalar_p

! Diffusion coefficient for User-defined Scalars
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Dif_Scalar

! New scalars for population balance equation
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  Source_a
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE ::  S_bar
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE ::  omega
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Matrix_a
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Matrix_b
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Matrix_c
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Inv_a
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ystart
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: A
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: beta_a
      INTEGER :: IJK_INDEX


      END MODULE scalars
