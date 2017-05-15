!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS3_DES                                               !
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
      SUBROUTINE USR3_DES

      use discretelement, only: PIP
      use discretelement, only: DES_VEL_NEW
      use discretelement, only: OMEGA_NEW
      use constant, only: PI

      use usr

      IMPLICIT NONE


! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER(LEN=64) :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: FEXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: UDF_UNIT = 2030
! Calculated tangential restitution coefficient.
      DOUBLE PRECISION :: RST_COEFF(2)
! Rebound Angle (degrees)
      DOUBLE PRECISION :: RBND_ANGLE(2)
! Particle loop counter
      INTEGER :: LC, NP
! XZ Velocity Mag
      DOUBLE PRECISION :: VEL_XZ(2), ROT_XZ(2)

! Particle-Wall Rebound Angle (degrees)
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_ALPHA.dat'
      OPEN(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

      DO LC=1,31

! Calculate the particle-wall rebound angle.
         NP = LC+31
         VEL_XZ(1) = sqrt(DES_VEL_NEW(NP,1)**2 + DES_VEL_NEW(NP,3)**2)
         RBND_ANGLE(1) = atan(VEL_XZ(1)/DES_VEL_NEW(NP,2))*180.0/PI

! Calculate the particle-wall rebound angle.
         NP = LC
         VEL_XZ(2) = sqrt(DES_VEL_NEW(NP,1)**2 + DES_VEL_NEW(NP,3)**2)
         RBND_ANGLE(2) = atan(VEL_XZ(2)/DES_VEL_NEW(NP,2))*180.0/PI

! Write the results to a file.
         WRITE(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(NP), RBND_ANGLE(1:2)
      ENDDO
      CLOSE(UDF_UNIT)



! Particle-wall Angular velocity (rad/sec)
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_OMEGA.dat'
      OPEN(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

      DO LC=1,31

! Calculate particle-particle angular velocity
         NP=LC+31
         ROT_XZ(1) = sqrt(OMEGA_NEW(NP,1)**2 + OMEGA_NEW(NP,3)**2)

! Calculate particle-wall angular velocity
         NP=LC
         ROT_XZ(2) = sqrt(OMEGA_NEW(NP,1)**2 + OMEGA_NEW(NP,3)**2)

! Write the results to a file.
         WRITE(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(NP), ROT_XZ(1:2)
      ENDDO
      CLOSE(UDF_UNIT)


! Particle-wall Tangential Restitution Coefficient
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_COEFF.dat'
      OPEN(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

      DO LC=2, 31
! Calculate the particle-wall restitution coefficient.
         NP=LC+31
         VEL_XZ(1) = sqrt(DES_VEL_NEW(NP,1)**2 + DES_VEL_NEW(NP,3)**2)
         RST_COEFF(1) =  VEL_XZ(1)/INIT_VEL_T(NP)
! Calculate the particle-wall restitution coefficient.
         NP=LC
         VEL_XZ(2) = sqrt(DES_VEL_NEW(NP,1)**2 + DES_VEL_NEW(NP,3)**2)
         RST_COEFF(2) =  VEL_XZ(2)/INIT_VEL_T(NP)
! Write the results to a file.
         WRITE(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(NP), RST_COEFF(1:2)
      ENDDO
      CLOSE(UDF_UNIT)


      RETURN
      END SUBROUTINE USR3_DES
