!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR0                                                   !
!  Author: J.Musser                                   Date: dd-mmm-yy  !
!  Purpose: This routine is called before the time loop starts and is  !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting constants and checking errors in   !
!           data.  This routine is not called from an IJK loop, hence  !
!           all indices are undefined.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use constant
      use discretelement
      use exit, only: mfix_exit
      use run, only: TSTOP
      use usr

      IMPLICIT NONE

      if(particles /= 1) then
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      endif

! Store the initial particle velocity.
      u0 = des_vel_new(1,1) 

! Calculate the time at which slipping ends.
      tsA = (-2.0d0*u0)/(7.0d0*MEW_W*grav(2))

      if(TSTOP < tsA) then
         write(*,"(3x,'simulation not long enough.')")
         call mfix_exit(0)
      endif

      return
      END SUBROUTINE USR0
