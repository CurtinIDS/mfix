
      MODULE turb

! Source term for K & Epsilon eq is linearized as
!    Source = K_Turb_G_c - K_Turb_G_p * K_Turb_G
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_Turb_G_c
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_Turb_G_p
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  E_Turb_G_c
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  E_Turb_G_p

! Diffusion coefficient for K & Epsilon
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Dif_K_Turb_G
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Dif_E_Turb_G

! Cross-corelation K_12 and time-scales in Simonin and Ahmadi
! models
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_12
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Tau_12
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Tau_1

      END MODULE turb
