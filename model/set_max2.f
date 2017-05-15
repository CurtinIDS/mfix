!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_MAX2                                               C
!  Purpose: calculate IMAX1,IMAX2,JMAX1,JMAX2,KMAX1,KMAX2,IJMAX2       C
!                     IJKMAX2                                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX, JMAX, KMAX, NO_I, NO_J, NO_K            C
!  Variables modified: IMAX1, IMAX2, JMAX1, JMAX2, KMAX1, KMAX2        C
!                      IJMAX2, IJKMAX2, IMIN1, JMIN1, KMIN1            C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_MAX2

! X-Axix partition specifications.
      USE compar, only: NODESI
      USE geometry, only: IMAX,  NO_I,  DO_I
      USE geometry, only: IMIN1, IMIN2, IMIN3, IMIN4
      USE geometry, only: IMAX1, IMAX2, IMAX3, IMAX4

! Y-Axix partition specifications.
      USE compar, only: NODESJ
      USE geometry, only: JMAX,  NO_J,  DO_J
      USE geometry, only: JMIN1, JMIN2, JMIN3, JMIN4
      USE geometry, only: JMAX1, JMAX2, JMAX3, JMAX4

! Z-Axix partition specifications.
      USE compar, only: NODESK
      USE geometry, only: KMAX,  NO_K,  DO_K
      USE geometry, only: KMIN1, KMIN2, KMIN3, KMIN4
      USE geometry, only: KMAX1, KMAX2, KMAX3, KMAX4

! Calculated array sizes.
      USE geometry, only: IJMAX2
      USE geometry, only: IJKMAX1, IJKMIN1
      USE geometry, only: IJKMAX2, IJKMAX3, IJKMAX4

      IMPLICIT NONE


! Initialize I's
      IMIN1=1;  IMIN2=1;  IMIN3=1;  IMIN4=1
      IMAX1=1;  IMAX2=1;  IMAX3=1;  IMAX4=1

      DO_I=.NOT.NO_I

! Set the domain specific values.
      IF(DO_I) THEN
         IMIN1 = 2
         IMAX1 = IMAX + 1
         IMAX2 = IMAX + 2
         IMIN2 = 1
         IF(NODESI.NE.1) THEN
            IMIN3 = 0
            IMAX3 = IMAX + 3
            IMIN4 = -1
            IMAX4 = IMAX + 4
         ELSE
            IMIN3 = IMIN2
            IMAX3 = IMAX2
            IMIN4 = IMIN3
            IMAX4 = IMAX3
         ENDIF
      ENDIF


! Initialize J's
      JMIN1=1;  JMIN2=1;  JMIN3=1;  JMIN4=1
      JMAX1=1;  JMAX2=1;  JMAX3=1;  JMAX4=1

      DO_J=.NOT.NO_J

! Set the domain specific values.
      IF(DO_J) THEN
         JMIN1 = 2
         JMAX1 = JMAX + 1
         JMAX2 = JMAX + 2
         JMIN2 = 1
         IF(NODESJ.NE.1) THEN
            JMIN3 = 0
            JMAX3 = JMAX + 3
            JMIN4 = -1
            JMAX4 = JMAX + 4
         ELSE
            JMIN3 = JMIN2
            JMAX3 = JMAX2
            JMIN4 = JMIN3
            JMAX4 = JMAX3
         ENDIF
      ENDIF


! Initialize J's
      KMIN1=1;  KMIN2=1;  KMIN3=1;  KMIN4=1
      KMAX1=1;  KMAX2=1;  KMAX3=1;  KMAX4=1

      DO_K=.NOT.NO_K

! Set the domain specific values.
      IF (DO_K) THEN
         KMIN1 = 2
         KMAX1 = KMAX + 1
         KMAX2 = KMAX + 2
         KMIN2 = 1
         IF(NODESK.NE.1) THEN
            KMIN3 = 0
            KMAX3 = KMAX + 3
            KMIN4 = -1
            KMAX4 = KMAX + 4
         ELSE
            KMIN3 = KMIN2
            KMAX3 = KMAX2
            KMIN4 = KMIN3
            KMAX4 = KMAX3
         ENDIF
      ENDIF

! Number of cells in I/J plane.
      IJMAX2 = IMAX2*JMAX2
! Totoal number of possible fluid cells.
      IJKMAX2 = IMAX2*JMAX2*KMAX2

      IF (DO_K) THEN
         IJKMIN1 = IJMAX2 + 1
         IJKMAX1 = IJKMAX2 - IJMAX2
      ELSE
         IJKMIN1 = IMAX2 + 1
         IJKMAX1 = IJKMAX2 - IMAX2
      ENDIF

! Max cell count with one layer of ghost cells.
      IJKMAX3 = (IMAX3-IMIN3+1)*(JMAX3-JMIN3+1)*(KMAX3-KMIN3+1)
! Max cell count with two layers of ghost cells.
      IJKMAX4 = (IMAX4-IMIN4+1)*(JMAX4-JMIN4+1)*(KMAX4-KMIN4+1)

      RETURN
      END SUBROUTINE SET_MAX2
