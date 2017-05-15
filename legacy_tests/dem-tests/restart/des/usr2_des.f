!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE


!-----------------------------------------------
! Local variables
!-----------------------------------------------
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: TP_UNIT = 2030

      DOUBLE PRECISION Tg, Tp0, Coeff, Tp,REL_ERR


      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
      DOUBLE PRECISION, SAVE :: OUT_TIME
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 0.0001d0

      if ( 0.335 < time .and. time < 0.336 ) then
         CALL WRITE_DES_Tp(S_TIME)
      endif

      RETURN

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_URS2                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_DES_Tp(lTime)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION lTime

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: TP_UNIT = 2030

      DOUBLE PRECISION Tg, Tp0, Coeff, Tp,REL_ERR


      FNAME = 'POST_posvel.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=TP_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(TP_UNIT,"(7X,'Time',12X,'x1',10X,'y1',8X,'u1',10x,'v1',10x,'x2',10x,'y2',10x,'u2',10x,'v2')")
      ELSE
         OPEN(UNIT=TP_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Calculate the value for the analytic solution.

! Calculate the relative error.
!      REL_ERR = (ABS(Tp - DES_T_s(1))/ABS(Tp))*100

! Write the data to a file.
      WRITE(TP_UNIT,"(9(3X,F12.5))")lTime, &
           des_pos_new(1,1),des_pos_new(1,2), &
           des_vel_new(1,1),des_vel_new(1,2), &
           des_pos_new(2,1),des_pos_new(2,2), &
           des_vel_new(2,1),des_vel_new(2,2)
      CLOSE(TP_UNIT)

      END SUBROUTINE WRITE_DES_Tp

      END SUBROUTINE USR2_DES
