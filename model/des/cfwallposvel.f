!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFWALLPOSVEL
!  Purpose: Calculate the position and velocity of wall particle
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFWALLPOSVEL(L, IW, WALL_POS, WALL_VEL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      USE des_bc
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle no. index
      INTEGER, INTENT(IN) :: L
! index of wall no. (1-6)
      INTEGER, INTENT(IN) :: IW
! wall position is returned, set based on particle position and radius
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_POS(3)
! wall velocity is returned, set based on any wall bc settings
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_VEL(3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! solids phase index of particle L
      INTEGER M
! local value for particle radius
      DOUBLE PRECISION DES_R
!-----------------------------------------------

! initialize wall velocity
      WALL_VEL = ZERO
      WALL_POS = ZERO

      DES_R = DES_RADIUS(L)

! Assigning wall position and velocity
! Find which solids phase the particle belongs to
      M = PIJK(L,5)

! west wall; in 3D on yz plane (x=wx1->0)
      IF(IW.EQ.1) THEN
         WALL_POS(1) = WX1 - DES_R
         WALL_POS(2) = DES_POS_NEW(L,2)
         IF(DO_K) WALL_POS(3) = DES_POS_NEW(L,3)

         WALL_VEL(1) = ZERO
         WALL_VEL(2) = DES_BC_Vw_s(IW,M)
         IF(DO_K) WALL_VEL(3) = DES_BC_Ww_s(IW,M)

         WALL_NORMAL(1,1) = -ONE
         WALL_NORMAL(1,2) = ZERO
         WALL_NORMAL(1,3) = ZERO

! east wall; in 3D on yz plane (x=ex2->xlength)
      ELSEIF(IW.EQ.2) THEN
         WALL_POS(1) = EX2 + DES_R
         WALL_POS(2) = DES_POS_NEW(L,2)
         IF(DO_K) WALL_POS(3) = DES_POS_NEW(L,3)

         WALL_VEL(1) = ZERO
         WALL_VEL(2) = DES_BC_Vw_s(IW,M)
         IF(DO_K) WALL_VEL(3) = DES_BC_Ww_s(IW,M)

         WALL_NORMAL(2,1) = ONE
         WALL_NORMAL(2,2) = ZERO
         WALL_NORMAL(2,3) = ZERO

! bottom wall; in 3D on xz plane (y=by1->0)
      ELSEIF(IW.EQ.3) THEN
         WALL_POS(1) = DES_POS_NEW(L,1)
         WALL_POS(2) = BY1 - DES_R
         IF(DO_K) WALL_POS(3) = DES_POS_NEW(L,3)

         WALL_VEL(1) = DES_BC_Uw_s(IW,M)
         WALL_VEL(2) = ZERO
         IF(DO_K) WALL_VEL(3) = DES_BC_Ww_s(IW,M)

         WALL_NORMAL(3,1) = ZERO
         WALL_NORMAL(3,2) = -ONE
         WALL_NORMAL(3,3) = ZERO

! top wall; in 3D on xz plane (y=ty2->ylength)
      ELSEIF(IW.EQ.4) THEN
         WALL_POS(1) = DES_POS_NEW(L,1)
         WALL_POS(2) = TY2 + DES_R
         IF(DO_K) WALL_POS(3) = DES_POS_NEW(L,3)

         WALL_VEL(1) = DES_BC_Uw_s(IW,M)
         WALL_VEL(2) = ZERO
         IF(DO_K) WALL_VEL(3) = DES_BC_Ww_s(IW,M)

         WALL_NORMAL(4,1) = ZERO
         WALL_NORMAL(4,2) = ONE
         WALL_NORMAL(4,3) = ZERO

! south wall; in 3D on xy plane (z=sz1->0)
      ELSEIF(IW.EQ.5) THEN
         WALL_POS(1) = DES_POS_NEW(L,1)
         WALL_POS(2) = DES_POS_NEW(L,2)
         WALL_POS(3) = SZ1 - DES_R

         WALL_VEL(1) = DES_BC_Uw_s(IW,M)
         WALL_VEL(2) = DES_BC_Vw_s(IW,M)
         WALL_VEL(3) = ZERO

         WALL_NORMAL(5,1) = ZERO
         WALL_NORMAL(5,2) = ZERO
         WALL_NORMAL(5,3) = -ONE

! north wall; in 3D on xy plane (z=nz2->zlength)
      ELSEIF(IW.EQ.6) THEN
         WALL_POS(1) = DES_POS_NEW(L,1)
         WALL_POS(2) = DES_POS_NEW(L,2)
         WALL_POS(3) = NZ2 + DES_R

         WALL_VEL(1) = DES_BC_Uw_s(IW,M)
         WALL_VEL(2) = DES_BC_Vw_s(IW,M)
         WALL_VEL(3) = ZERO

         WALL_NORMAL(6,1) = ZERO
         WALL_NORMAL(6,2) = ZERO
         WALL_NORMAL(6,3) = ONE

      ENDIF

      RETURN
      END SUBROUTINE CFWALLPOSVEL



