module resize

  public :: byte_grow
  public :: integer_grow
  public :: integer_grow2_reverse
  public :: integer_grow2
  public :: logical_grow
  public :: logical_grow2
  public :: real_grow
  public :: real_grow2
  public :: real_grow2_reverse
  public :: real_grow3
  public :: logical_grow2_reverse

  contains

      SUBROUTINE BYTE_GROW(byte_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: byte_array
        INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: byte_tmp
        INTEGER lSIZE

        lSIZE = size(byte_array,1)
        allocate(byte_tmp(new_size))
        byte_tmp(1:lSIZE) = byte_array(1:lSIZE)
        call move_alloc(byte_tmp,byte_array)

      END SUBROUTINE BYTE_GROW

      SUBROUTINE INTEGER_GROW(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE

        lSIZE = size(integer_array,1)
        allocate(integer_tmp(new_size))
        integer_tmp(1:lSIZE) = integer_array(1:lSIZE)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW

      SUBROUTINE INTEGER_GROW2_reverse(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(new_size,lSIZE2))
        integer_tmp(1:lSIZE,:) = integer_array(1:lSIZE,:)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2_reverse

      SUBROUTINE INTEGER_GROW2(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        allocate(integer_tmp(lSIZE,new_size))
        integer_tmp(:,1:lSIZE2) = integer_array(:,1:lSIZE2)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2

      SUBROUTINE LOGICAL_GROW(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE

        lSIZE = size(logical_array,1)
        allocate(logical_tmp(new_size))
        logical_tmp(1:lSIZE) = logical_array(1:lSIZE)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW

      SUBROUTINE LOGICAL_GROW2(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(logical_array,1)
        lSIZE2 = size(logical_array,2)
        allocate(logical_tmp(lSIZE,new_size))
        logical_tmp(:,1:lSIZE2) = logical_array(:,1:lSIZE2)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW2

      SUBROUTINE REAL_GROW(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE

        lSIZE = size(real_array,1)
        allocate(real_tmp(new_size))
        real_tmp(1:lSIZE) = real_array(1:lSIZE)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW

      SUBROUTINE REAL_GROW2(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(lSIZE,new_size))
        real_tmp(:,1:lSIZE2) = real_array(:,1:lSIZE2)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2

      SUBROUTINE REAL_GROW2_reverse(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(new_size,lSIZE2))
        real_tmp(1:lSIZE,:) = real_array(1:lSIZE,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2_REVERSE

      SUBROUTINE REAL_GROW3(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2, lSIZE3

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lSIZE3 = size(real_array,3)
        allocate(real_tmp(lSIZE,lSIZE2,new_size))
        real_tmp(:,:,1:lSIZE3) = real_array(:,:,1:lSIZE3)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW3

      SUBROUTINE LOGICAL_GROW2_REVERSE(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        allocate(real_tmp(new_size,lSIZE2))
        real_tmp(1:lSIZE,:) = real_array(1:lSIZE,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE LOGICAL_GROW2_REVERSE

end module resize
