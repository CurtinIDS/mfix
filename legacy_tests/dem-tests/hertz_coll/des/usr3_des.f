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
      DOUBLE PRECISION :: RST_COEFF
! Rebound Angle (degrees)
      DOUBLE PRECISION :: RBND_ANGLE
! Calculated absolute percent relative error
      DOUBLE PRECISION :: rErr
! Particle loop counter
      INTEGER :: NP


! Rebound Angle (degrees)
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_ALPHA.dat'
      INQUIRE(FILE=FNAME,EXIST=FEXISTS)
      IF (.NOT.FEXISTS) THEN
         OPEN(UNIT=UDF_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(UDF_UNIT,1000)
      ELSE
         OPEN(UNIT=UDF_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      DO NP=ANGLE_START, ANGLE_END
! Calculate the rebound angle.
         RBND_ANGLE = atan(DES_VEL_NEW(NP,1)/DES_VEL_NEW(NP,2))*180.0/PI
! Calculate the absolute relative error.
         rErr = (ABS(EXP_ANGLE(NP) - RBND_ANGLE)/                      &
            ABS(EXP_ANGLE(NP)))*100
! Write the results to a file.
         WRITE(UDF_UNIT,"(4(3x,F11.4))") INIT_ANGLE(NP),               &
             EXP_ANGLE(NP), RBND_ANGLE, rErr
      ENDDO
      CLOSE(UDF_UNIT)


! Angular velocity
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_OMEGA.dat'
      INQUIRE(FILE=FNAME,EXIST=FEXISTS)
      IF (.NOT.FEXISTS) THEN
         OPEN(UNIT=UDF_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(UDF_UNIT,1000)
      ELSE
         OPEN(UNIT=UDF_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      DO NP=AVEL_START, AVEL_END
! Calculate the absolute relative error.
         rErr = (ABS(EXP_OMEGA(NP) - OMEGA_NEW(NP,3)) /                &
            ABS(EXP_OMEGA(NP)))*100
! Write the results to a file.
         WRITE(UDF_UNIT,"(4(3x,F11.4))") INIT_ANGLE(NP),               &
            EXP_OMEGA(NP), OMEGA_NEW(NP,3), rErr
      ENDDO
      CLOSE(UDF_UNIT)


! Tangential Restitution Coefficient
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_COEFF.dat'
      INQUIRE(FILE=FNAME,EXIST=FEXISTS)
      IF (.NOT.FEXISTS) THEN
         OPEN(UNIT=UDF_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(UDF_UNIT,1000)
      ELSE
         OPEN(UNIT=UDF_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      DO NP=COEFF_START, COEFF_END
! Calculate the restitution coefficient.
         RST_COEFF = DES_VEL_new(NP,1)/INIT_VEL_T(NP)
! Calculate the absolute relative error.
         rErr = abs(EXP_COEFF(NP)-RST_COEFF)/abs(EXP_COEFF(NP))*100
! Write the results to a file.
         WRITE(UDF_UNIT,"(4(3x,F11.4))") INIT_ANGLE(NP),               &
            EXP_COEFF(NP), RST_COEFF, rErr
      ENDDO
      CLOSE(UDF_UNIT)


 1000 FORMAT(/4x,'Init Angle',4x,'Experiment',7x,'MFIX',8x,'|%Rel Err|')

      RETURN
      END SUBROUTINE USR3_DES
