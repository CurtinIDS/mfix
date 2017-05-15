!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS                                                    !
!  Author: J.Musser                                   Date: 17-Nov-14  !
!                                                                      !
!  REF: Di Renzo, A. and Di Maio F.P. "Comparison of contact-force     !
!       models for the simulation of collisions in DEM-based granular  !
!       flow codes," Chemical Engineering Science, 59(3), pg 525-541.  !
!                                                                      !
!  REF: Kharaz, A.H., Gorham, D.A., and Salman, A.D. "An experimental  !
!       study of the elastic rebound of spheres," Powder Technology,   !
!       120(3), pg 281-291.                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE usr

      Use param
      Use param1

! a dummy variable listed in usrnlst.inc
      DOUBLE PRECISION DUMMY_DP

! Initial collision angles
      DOUBLE PRECISION :: INIT_ANGLE(50)
! Initial Tangential velocity
      DOUBLE PRECISION :: INIT_Vel_T(50)

! The start/end indices for particles in rebound angle test
      INTEGER, PARAMETER :: ANGLE_START=1, ANGLE_END=9
! Experimentally measured rebound angles.
      DOUBLE PRECISION, PARAMETER :: EXP_ANGLE(ANGLE_START:ANGLE_END) =&
         (/  0.1194,  1.7630, 2.9577, 6.6917, 12.3645, 21.7809,        &
            32.5455, 45.1051, 57.3652 /)

! The start/end indices for particles in angular velocity tests
      INTEGER, PARAMETER :: AVEL_START=10, AVEL_END=16
! Experimentally measured post-collision angluar velocities
      DOUBLE PRECISION, PARAMETER :: EXP_OMEGA(AVEL_START:AVEL_END) =  &
         (/ 135.5703, 264.5493, 577.0096, 612.6481, 595.4843,          &
            452.9444, 381.3185 /)

! The start/end indices for particles in angular velocity tests
      INTEGER, PARAMETER :: COEFF_START=17, COEFF_END=23
! Experimentally measured post-collision
      DOUBLE PRECISION, PARAMETER :: EXP_COEFF(COEFF_START:COEFF_END) =&
         (/ 0.4902, 0.3439, 0.1915, 0.3756, 0.5427, 0.6927, 0.8024 /)

      END MODULE usr
