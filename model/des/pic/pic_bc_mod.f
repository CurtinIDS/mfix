!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: pic_bc_mod.f                                           !
!                                                                      !
!  Purpose: Common elements needed for the pic mass inflow boundary    !
!  condition.                                                          !
!                                                                      !
!  Author: R. Garg                                   Date: 11-Jun-14   !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE PIC_BC

      use param, only: dimension_bc
!......................................................................!

      INTEGER :: PIC_BCMI
      INTEGER :: PIC_BCMO

      LOGICAL :: PIC_MIO  ! either inlet or outlet exists

! Map between PIC MI/MO IDs and the user input BC index.
      INTEGER :: PIC_BCMI_MAP(DIMENSION_BC)
      INTEGER :: PIC_BCMO_MAP(DIMENSION_BC)


      INTEGER, ALLOCATABLE :: PIC_BCMO_IJKSTART(:)
      INTEGER, ALLOCATABLE :: PIC_BCMO_IJKEND(:)

      INTEGER, ALLOCATABLE :: PIC_BCMO_IJK(:)


      INTEGER, ALLOCATABLE :: PIC_BCMI_IJKSTART(:)
      INTEGER, ALLOCATABLE :: PIC_BCMI_IJKEND(:)

      INTEGER, ALLOCATABLE :: PIC_BCMI_IJK(:)

! Direction of the MI BC plane
      INTEGER, ALLOCATABLE :: PIC_BCMI_NORMDIR(:,:)
! BC planes will be calculated from the fluid node as
! XE-dx, YN-dy, and ZT-dz.
! For BC-plane 'S', YN-dy will provide the wrong plane, as
! it shud be just yn. The array below is used to provide this correction.
! BC_plane coords will be calculated as x_i - offset dx_i and offset will
! have default value of 1. It will be set to 0 for a case described above.
      INTEGER, ALLOCATABLE :: PIC_BCMI_OFFSET(:,:)


! This will store cumulative number of computational parcels
      Double precision, Allocatable :: PIC_BCMI_CNP(:, :)
! This will store cumulative number of implied real particles
      Double precision, Allocatable :: PIC_BCMI_RNP(:, :)

! To include or not to include cut-cells in the MI BC

      LOGICAL, Allocatable :: PIC_BCMI_INCL_CUTCELL(:)

! Logicals to print seeding and deletion of particles based on user needs
      LOGICAL :: PIC_REPORT_SEEDING_STATS, PIC_REPORT_DELETION_STATS

! Mininum velocity for parcles at rest of offset gravitational forces
      DOUBLE PRECISION :: minVEL(3), minVEL_MAG, OoMinVEL_MAG

      END MODULE PIC_BC

