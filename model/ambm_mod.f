! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ambm                                                        !
!  Purpose:                                                            !
!     IMPORTANT:  For using these arrays in a subroutine               !
!     -lock the module in the beginning of the subroutine              !
!      call lock_ambm                                                  !
!     -and unlock the module at the end of the subroutine              !
!      call unlock_ambm                                                !
! Contains the following subroutines:                                  !
!      lock_ambm, unlock_ambm                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE ambm

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE exit, only: mfix_exit
      USE funits
!-----------------------------------------------

! linear equation matrix and vector
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: A_m
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: B_m

      LOGICAL :: ambm_locked = .false.

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE lock_ambm
         IMPLICIT NONE
      IF(ambm_locked) THEN
         IF (DMP_LOG) WRITE(*,*) &
            'Error:  Multiple use of ambm (ambm_mod.f)'
         CALL MFIX_EXIT(myPE)
      ELSE
         ambm_locked = .true.
      ENDIF
      END SUBROUTINE lock_ambm

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE unlock_ambm
      ambm_locked = .false.
      END SUBROUTINE unlock_ambm

      END MODULE ambm
