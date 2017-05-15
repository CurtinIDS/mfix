      MODULE visc_s


! Granular first coefficient of (shear) viscosity(i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s

! 'Solids' first coefficient of (shear) viscosity(i,j,k) multiplied
! by the volume fraction of that phase if ishii otherwise by 1
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  EPMU_s

! Granular first coefficient of (shear) viscosity(i,j,k)
! Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_v

! Granular first coefficient of (shear) viscosity(i,j,k)
! Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_p

! Granular first coefficient of (shear) viscosity(i,j,k)
! Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_f

! Bulk viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mu_b_v

! Granular first coefficient of (shear) viscosity(i,j,k)
! Collisional contribution of viscous component
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s_c

! Granular second coefficient of viscosity(i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s

! 'Solids' second coefficient of viscosity(i,j,k) multiplied by
! the volume fraction of that phase if ishii otherwise by 1
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  EPLAMBDA_s

! Granular second coefficient of (shear) viscosity(i,j,k)
! Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_v

! Granular second coefficient of (shear) viscosity(i,j,k)
! Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_p

! Granular second coefficient of (shear) viscosity(i,j,k)
! Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_f

! Granular second coefficient of viscosity
! Collisional contribution of viscous component
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s_c

! Packed bed (close packed) void fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_star_array

! Start of Blending
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_g_blend_start

! End of Blending
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_g_blend_end

! trace of D_s (rate of strain tensor) at i, j, k
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s

! Second invariant of the deviator of D_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: I2_devD_s

! Boyle-Massoudi stress tensor coefficient (i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ALPHA_s

! For Boyle-Massoudi stress tensor
! Trace of M_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trM_s

! Trace of (D_s)(M_s)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trDM_s


            END MODULE visc_s
