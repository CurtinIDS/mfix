!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFWALLCONTACT(WALL, L, WALLCONTACT)
!  Purpose: Check if particle L is in contact with WALL.  If so, set
!           WALLCONTACT to 1, else WALLCONTACT is 0
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFWALLCONTACT(WALL, L, WALLCONTACT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE constant
      USE parallel
      USE compar
      Use discretelement
      USE des_bc
      use geometry, only: DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Given wall ID number (1=west, 2=east, 3=south, 4=north, 5=bottom,
! 6=top)
      INTEGER, INTENT (IN) :: WALL
! Given particle ID number
      INTEGER, INTENT (IN) :: L
! Flag to indicate whether given particle is in contact with given wall
! (1=contact, 0 = no contact)
      INTEGER, INTENT (INOUT) :: WALLCONTACT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variables for x, y, z position of the particle
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS
! local variables to define system dimensions
      DOUBLE PRECISION :: LXE, LXW, LYN, LYS, LZT, LZB
! local variables: distance between particle surface and wall
      DOUBLE PRECISION :: DistApart
!-----------------------------------------------

! assign temporary local variables for quick reference
      LXE = EX2
      LXW = WX1
      LYN = TY2
      LYS = BY1
      LZT = NZ2
      LZB = SZ1

! assign temporary local variables for manipulation/use
      XPOS = DES_POS_NEW(L,1)
      YPOS = DES_POS_NEW(L,2)
      IF(DO_K) ZPOS = DES_POS_NEW(L,3)


! initialize
      WALLCONTACT = 0

      IF (DES_LE_BC) THEN
! Current implementation of Lees & Edwards boundaries implies all other
! boundaries are periodic (i.e. no walls in system)
         RETURN
      ELSEIF (DES_PERIODIC_WALLS) THEN
! Check if current wall corresponds to a periodic boundary (i.e. no wall)
         IF( (DES_PERIODIC_WALLS_X .AND. (WALL.EQ.1.OR.WALL.EQ.2)).OR.&
             (DES_PERIODIC_WALLS_Y .AND. (WALL.EQ.3.OR.WALL.EQ.4)).OR.&
             (DO_K.AND.DES_PERIODIC_WALLS_Z .AND. &
             (WALL.EQ.5.OR.WALL.EQ.6)) ) THEN
            RETURN
         ENDIF
      ENDIF

! Note that if no cohesion is used WALL_VDW_OUTER_CUTOFF = zero.
! west wall (X)
      IF(WALL.EQ.1) THEN
         DistApart = XPOS-LXW-DES_RADIUS(L)
! consider wall contact calculations for particle wall distances less
! than the cutoff
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) WALLCONTACT = 1

! east wall (X)
      ELSEIF(WALL.EQ.2) THEN
         DistApart = LXE-XPOS-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) WALLCONTACT = 1

! south wall (Y)
      ELSEIF(WALL.EQ.3) THEN
         DistApart = YPOS-(LYS)-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) WALLCONTACT = 1

! north wall (Y)
      ELSEIF(WALL.EQ.4) THEN
         DistApart = LYN-YPOS-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) WALLCONTACT = 1

! bottom wall (Z)
      ELSEIF(WALL.EQ.5) THEN
         DistApart = ZPOS-LZB-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) WALLCONTACT = 1

! top wall (Z)
      ELSEIF(WALL.EQ.6) THEN
         DistApart = LZT-ZPOS-DES_RADIUS(L)
         IF( DistApart <= WALL_VDW_OUTER_CUTOFF ) WALLCONTACT = 1
      ENDIF

      RETURN
      END SUBROUTINE CFWALLCONTACT
