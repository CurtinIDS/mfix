!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NSQUARE                                                C
!>  Purpose: DES - N-Square neighbor search
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NSQUARE

      use des_allocate
      use des_bc
      use des_thermo
      use discretelement
      use functions
      use geometry, only: DO_K, xlength, ylength, zlength
      use param1, only: zero

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I, J, K
      INTEGER L, LL, CC
      DOUBLE PRECISION DISTVEC(3), DIST, R_LM
! Temporary variables to adjust particle position in event of periodic
! boundaries
      DOUBLE PRECISION XPOS(2), YPOS(2), ZPOS(2), TMPPOS(3)
! Max loop limit in each coordinate direction
      INTEGER II, JJ, KK
! Index to track accounted for particles
      INTEGER PC
! Index to track unaccounted for particles
      INTEGER PNPC
!-----------------------------------------------

      PC=1
      DO L=1, MAX_PIP

         NEIGHBOR_INDEX(L) = 1
         if (L .gt. 1) NEIGHBOR_INDEX(L) = NEIGHBOR_INDEX(L-1)

         IF(PC .GE. PIP ) EXIT
         IF(IS_NONEXISTENT(L)) CYCLE

         PNPC = PIP - PC
         DO LL = L+1, MAX_PIP

            IF(PNPC .LE. 0) EXIT
            IF(IS_NONEXISTENT(LL)) CYCLE

            R_LM = DES_RADIUS(L) + DES_RADIUS(LL)
            R_LM = FACTOR_RLM*R_LM

! the following section adjusts the neighbor check routine when
! any boundary is periodic
! ------------------------------
            IF (DES_PERIODIC_WALLS) THEN
               XPOS(:) = DES_POS_NEW(LL,1)
               YPOS(:) = DES_POS_NEW(LL,2)
               II = 1
               JJ = 1
               KK = 1

               IF(DES_PERIODIC_WALLS_X) THEN
                  IF (DES_POS_NEW(L,1) + R_LM > XLENGTH) THEN
                     II = 2
                     XPOS(II) = DES_POS_NEW(LL,1) + XLENGTH
                  ELSEIF (DES_POS_NEW(L,1) - R_LM < ZERO) THEN
                     II = 2
                     XPOS(II) = DES_POS_NEW(LL,1) - XLENGTH
                  ENDIF
               ENDIF
               IF(DES_PERIODIC_WALLS_Y) THEN
                  IF (DES_POS_NEW(L,2) + R_LM > YLENGTH) THEN
                     JJ = 2
                     YPOS(JJ) = DES_POS_NEW(LL,2) + YLENGTH
                  ELSEIF (DES_POS_NEW(L,2) - R_LM < YLENGTH) THEN
                     JJ = 2
                     YPOS(JJ) = DES_POS_NEW(LL,2) - YLENGTH
                  ENDIF
               ENDIF
               IF(DO_K) THEN
                  ZPOS(:) = DES_POS_NEW(LL,3)
                  IF(DES_PERIODIC_WALLS_Z) THEN
                     IF (DES_POS_NEW(L,3) + R_LM > ZLENGTH) THEN
                        KK = 2
                        ZPOS(KK) = DES_POS_NEW(LL,3) + ZLENGTH
                     ELSEIF (DES_POS_NEW(L,3) - R_LM < ZERO) THEN
                        KK = 2
                        ZPOS(KK) = DES_POS_NEW(LL,3) - ZLENGTH
                     ENDIF
                  ENDIF
               ENDIF

! if particle L is within R_LM of a periodic boundary then check
! particles LL current position and its position shifted to the
! opposite boundary for neighbor contact
               OUTER: DO I = 1,II
                  DO J = 1,JJ
                     TMPPOS(1) = XPOS(I)
                     TMPPOS(2) = YPOS(J)
                     IF (DO_K) THEN
                        DO K = 1,KK
                           TMPPOS(3) = ZPOS(K)
                           DISTVEC(:) = TMPPOS(:) - DES_POS_NEW(L,:)
                           DIST = dot_product(DISTVEC,DISTVEC)
                           IF (DIST.LE.R_LM) EXIT OUTER
                        ENDDO
                     ELSE
                        DISTVEC(:) = TMPPOS(:) - DES_POS_NEW(L,:)
                        DIST = dot_product(DISTVEC,DISTVEC)
                        IF (DIST.LE.R_LM) EXIT OUTER
                     ENDIF
                  ENDDO
               ENDDO OUTER

            ELSE   ! if .not.des_periodic_walls
               DISTVEC(:) = DES_POS_NEW(LL,:) - DES_POS_NEW(L,:)
               DIST = dot_product(DISTVEC,DISTVEC)
            ENDIF    ! endif des_periodic_walls
! ------------------------------

            IF (DIST < R_LM**2) THEN
               cc = add_pair(L, LL)
            ENDIF
            PNPC = PNPC - 1
         ENDDO   ! end loop over LL

         PC = PC + 1
      ENDDO   ! end loop over L

      RETURN
      END SUBROUTINE NSQUARE
