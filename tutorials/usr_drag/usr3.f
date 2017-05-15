!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: USR3                                                     !
!                                                                      !
! Purpose: Collect some information at the end of the run.             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR3

      use discretelement, only: DES_MMAX
      use discretelement, only: DES_VEL_NEW
      use discretelement, only: DES_POS_NEW
      use discretelement, only: PIJK, PIP
      use functions, only: IS_NONEXISTENT

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//

      INTEGER :: lUNIT
      INTEGER :: NP, M

      INTEGER, DIMENSION(5) :: rID = (/1,25,50,75,100/)

      DOUBLE PRECISION, DIMENSION(5) :: lEXP = &
         (/39.45d0, 22.42d0, 29.88d0, 13.07d0, 21.83d0/)

      DOUBLE PRECISION, DIMENSION(5) :: lTHL = &
         (/35.40d0, 24.27d0, 34.43d0, 13.62d0, 23.52d0/)

      DO NP=1, PIP
         IF(IS_NONEXISTENT(NP)) CYCLE
         M =  PIJK(NP,5)
         lUNIT = 750 + M
         WRITE(lUNIT,"(2/,4x,'Run ID: ',I3)") rID(M)
         WRITE(lUNIT,"(   4x,'Theoretical:  ',F5.2)") lTHL(M)
         WRITE(lUNIT,"(   4x,'Experimental: ',F5.2)") lEXP(M)
         WRITE(lUNIT,"(   4x,'Simulated:    ',F5.2)") &
            sqrt(DOT_PRODUCT(DES_VEL_NEW(NP,:),DES_VEL_NEW(NP,:)))

         close(lUNIT)
      ENDDO

      RETURN
      END SUBROUTINE USR3
