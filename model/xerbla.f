      SUBROUTINE XERBLA ( SRNAME, INFO )
!     ..    Modules
      USE compar
      USE exit, only: mfix_exit

!     ..    Scalar Arguments ..
      INTEGER            INFO
      CHARACTER(LEN=6)        SRNAME
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the Level 2 BLAS routines.
!
!  It is called by the Level 2 BLAS routines if an input parameter is
!  invalid.
!
!  Installers should consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!  Parameters
!  ==========
!
!  SRNAME - CHARACTER*6.
!           On entry, SRNAME specifies the name of the routine which
!           called XERBLA.
!
!  INFO   - INTEGER.
!           On entry, INFO specifies the position of the invalid
!           parameter in the parameter-list of the calling routine.
!
!
!  Auxiliary routine for Level 2 Blas.
!
!  Written on 20-July-1986.
!
!     .. Executable Statements ..
!
      WRITE (*,99999) myPE,SRNAME, INFO

!     STOP
      call mfix_exit(myPE)
99999 FORMAT ( '(PE ',I6,'): ** On entry to ', A6, ' parameter number ', I2, &
         ' had an illegal value' )
!     End of XERBLA.
      END SUBROUTINE XERBLA

