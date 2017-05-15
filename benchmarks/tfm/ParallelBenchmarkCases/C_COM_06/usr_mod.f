      MODULE usr


      Use param
      Use param1


! Dummy variable.
      DOUBLE PRECISION :: DUMMY_DP

! Sherwood number
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: N_Sh

! Function of ash-layer void fraction
      DOUBLE PRECISION :: f_EP_A

! Prox. Analysis of char and ash.
      DOUBLE PRECISION :: PAFC, PAA



      END MODULE usr
