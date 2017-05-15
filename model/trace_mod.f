      MODULE trace


! Trace of D_s
! Note same quantity in visc_s_mod: trd_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s_C

! Trace of the square of D_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s2

! Trace of D_s at previous timestep
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s_Co

! Trace of D_s at previous time
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s_Co2

! Trace of the dot of D_s (Mth solids phase) and
! D_sl (Lth solids phase)
      DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: trD_s2_ip

      END MODULE trace
