!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use discretelement, only: DTSOLID
      use constant, only: C

      IMPLICIT NONE

! Set the solids time step to a larger value to a millisecond. Note
! that this likely results in greater numerical error, but at least
! the case runs quickly. This simulation can get away with this large
! DEM time step size because the particles' positions are fixed.
      DTSOLID = C(1)


      RETURN
      END SUBROUTINE USR0
