!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: indices                                                     !
!  Author: M. Syamlal                                 Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Global arrays for index computations.                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE indices

      Use param1, only: MAX_CLASS

! Increments used for index computation
      INTEGER :: INCREMENT_FOR_n (MAX_CLASS)
      INTEGER :: INCREMENT_FOR_s (MAX_CLASS)
      INTEGER :: INCREMENT_FOR_e (MAX_CLASS)
      INTEGER :: INCREMENT_FOR_w (MAX_CLASS)
      INTEGER :: INCREMENT_FOR_t (MAX_CLASS)
      INTEGER :: INCREMENT_FOR_b (MAX_CLASS)
      INTEGER :: INCREMENT_FOR_im(MAX_CLASS)
      INTEGER :: INCREMENT_FOR_ip(MAX_CLASS)
      INTEGER :: INCREMENT_FOR_jm(MAX_CLASS)
      INTEGER :: INCREMENT_FOR_jp(MAX_CLASS)
      INTEGER :: INCREMENT_FOR_km(MAX_CLASS)
      INTEGER :: INCREMENT_FOR_kp(MAX_CLASS)

      INTEGER :: INCREMENT_FOR_MP(6,MAX_CLASS) = 0
      INTEGER :: INCREMENT_FOR_NB(6,MAX_CLASS) = 0

! Increments used for index computation of 3rd layer
      INTEGER :: INCREMENT3_FOR_im(MAX_CLASS)
      INTEGER :: INCREMENT3_FOR_ip(MAX_CLASS)
      INTEGER :: INCREMENT3_FOR_jm(MAX_CLASS)
      INTEGER :: INCREMENT3_FOR_jp(MAX_CLASS)
      INTEGER :: INCREMENT3_FOR_km(MAX_CLASS)
      INTEGER :: INCREMENT3_FOR_kp(MAX_CLASS)

! Store LM index values
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: STORE_LM

! Identification of the cell class
      INTEGER, DIMENSION(:), ALLOCATABLE :: CELL_CLASS

! Maps between IJK and the base I/J/K indices
      INTEGER, DIMENSION(:), ALLOCATABLE :: I_OF ! IJK --> I
      INTEGER, DIMENSION(:), ALLOCATABLE :: J_OF ! IJK --> J
      INTEGER, DIMENSION(:), ALLOCATABLE :: K_OF ! IJK --> K

! +/- increments, shifted for cyclic BCs and DMP partitions.
      INTEGER, DIMENSION(:), ALLOCATABLE :: Im1, Ip1 ! I-1, I+1
      INTEGER, DIMENSION(:), ALLOCATABLE :: Jm1, Jp1 ! J-1, J+1
      INTEGER, DIMENSION(:), ALLOCATABLE :: Km1, Kp1 ! K-1, K+1

! Identification of the cell class for higher order scheme
      INTEGER, DIMENSION(:), ALLOCATABLE :: CELL_CLASS3

      INTEGER, DIMENSION(:), ALLOCATABLE :: I3_OF
      INTEGER, DIMENSION(:), ALLOCATABLE :: J3_OF
      INTEGER, DIMENSION(:), ALLOCATABLE :: K3_OF

      INTEGER, DIMENSION(:), ALLOCATABLE :: Im1_3, Ip1_3
      INTEGER, DIMENSION(:), ALLOCATABLE :: Jm1_3, Jp1_3
      INTEGER, DIMENSION(:), ALLOCATABLE :: Km1_3, Kp1_3

! Save original IJK value of Background grid (new to old mapping)
      INTEGER, DIMENSION(:), ALLOCATABLE :: BACKGROUND_IJK_OF
! Save new IJK value of Background grid (old to new mapping)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJK_OF_BACKGROUND
! Save original IJKEND3 value of Background grid
      INTEGER :: BACKGROUND_IJKEND3

      END MODULE indices
