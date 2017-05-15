      MODULE tau_g

! cross terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_W_g

! divergence of complete gas phase stress tensor in
! the x, y, z momentum eqs
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cTAU_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cTAU_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cTAU_W_g

      END MODULE tau_g
