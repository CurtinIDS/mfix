!----------------------------------------------------------------------!
! ******************************************************************** !
! *  There are minor differences between this file and the file of   * !
! *  the same name in the model directory:                           * !
! *  1) FUNIJK_IO is calcluated differently                          * !
! *  2) FUNIJK is set to FUNIJK_IO                                   * !
! *  3) FUNIJK_0 is set to FUNIJK_IO                                 * !
! *  4) FUNIJK_GL is set to FUNIJK_IO                                * !
! ******************************************************************** !
!----------------------------------------------------------------------!
      MODULE functions

      CONTAINS

!//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
!     FUNIJK (LI, LJ, LK) = c0 + LI + (LJ-jstart3_all(myPE))*c1 + (LK-kstart3_all(myPE))* c2
!      funijk(li,lj,lk) = lj + c0 + li*c1 + lk*c2
      INTEGER FUNCTION funijk_0(li,lj,lk)
      USE compar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      funijk_0 = funijk_io(li,lj,lk)
      END FUNCTION funijk_0

      INTEGER FUNCTION funijk(li,lj,lk)
      USE compar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      funijk = funijk_io(li,lj,lk)
      END FUNCTION funijk

! Function for generating the LOCAL 3-D array index IJK from the
! the 1-D indices I, J, K and IPROC.
!     FUNIJK_PROC(LI, LJ, LK, LIPROC) = 1 + (LI - istart3_all(LIPROC))+ &
!     (LJ-jstart3_all(LIPROC))*(iend3_all(LIPROC)-istart3_all(LIPROC)+1) &
!     + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
!     (iend3_all(LIPROC)-istart3_all(LIPROC)+1)
      INTEGER FUNCTION FUNIJK_PROC(LI, LJ, LK, LIPROC)
      USE compar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK, LIPROC
      FUNIJK_PROC = 1 + (LJ - jstart3_all(LIPROC))+ &
         (LI-Istart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1) &
         + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
         (iend3_all(LIPROC)-istart3_all(LIPROC)+1)
      END FUNCTION FUNIJK_PROC

      INTEGER FUNCTION FUNIJK_GL (LI, LJ, LK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      FUNIJK_GL = FUNIJK_IO(LI, LJ, LK)
      END FUNCTION FUNIJK_GL

      INTEGER FUNCTION FUNIJK_IO(LI, LJ, LK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      FUNIJK_IO = LI + (LJ-1)*IMAX2 + (LK-1)*IJMAX2
      END FUNCTION FUNIJK_IO


!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_OWNS                                           !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process.                          !
!                                                                      !
!  o Ownership is defined as belonging to the current PE's domain but  !
!    as a cell in any of the PE's ghost layers.                        !
!                                                                      !
!  o Each computational cell is owned by one -and only one- PE.        !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_OWNS(LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_OWNS = &
         LI >= ISTART .AND. LI <= IEND .AND. &
         LJ >= JSTART .AND. LJ <= JEND .AND. &
         LK >= KSTART .AND. LK <= KEND

      RETURN
      END FUNCTION IS_ON_MYPE_OWNS


!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_WOBND                                          !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process and not a exterior ghost  !
!  cell.                                                               !
!                                                                      !
!  o This is a subset of IS_ON_myPE_OWNS.                              !
!                                                                      !
!  o Exterior ghost cells are those in cells surrounding the domain.   !
!    These are cells created to fully define boundary conditions       !
!    (e.g., I == 1 where X_E(1) == ZERO).                              !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_wobnd (LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_WOBND = &
         LI >= ISTART1 .AND. LI <= IEND1 .AND. &
         LJ >= JSTART1 .AND. LJ <= JEND1 .AND. &
         LK >= KSTART1 .AND. LK <= KEND1   !.AND. &
!        (.NOT.DEAD_CELL_AT(LI,LJ,LK))

      RETURN
      END FUNCTION IS_ON_myPE_wobnd

!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_Plus1Layer                                     !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process or contained in the fisrt !
!  layer of ghost cells seen by the current PE.                        !
!                                                                      !
!  o This is a superset of IS_ON_myPE_OWNS.                            !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_plus1layer (LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_PLUS1LAYER = &
         LI >= ISTART2 .AND. LI <= IEND2 .AND. &
         LJ >= JSTART2 .AND. LJ <= JEND2 .AND. &
         LK >= KSTART2 .AND. LK <= KEND2

      RETURN
      END FUNCTION IS_ON_myPE_plus1layer


!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_Plus2Layer                                     !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process or contained in the fisrt !
!  two layers of ghost cells seen by the current PE.                   !
!                                                                      !
!  o This is a superset of IS_ON_Plus1Layer.                           !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_plus2layers (LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_PLUS2LAYERS = &
         LI >= ISTART3 .AND. LI <= IEND3 .AND. &
         LJ >= JSTART3 .AND. LJ <= JEND3 .AND. &
         LK >= KSTART3 .AND. LK <= KEND3  !.AND. &
!        (.NOT.DEAD_CELL_AT(LI,LJ,LK))

      RETURN
      END FUNCTION IS_ON_myPE_plus2layers

!---------------------------------------------------------------------//
! WEST_OF  (IJK)   = IJK + INCREMENT_FOR_w (CELL_CLASS(IJK))
! EAST_OF  (IJK)   = IJK + INCREMENT_FOR_e (CELL_CLASS(IJK))
! SOUTH_OF (IJK)   = IJK + INCREMENT_FOR_s (CELL_CLASS(IJK))
! NORTH_OF (IJK)   = IJK + INCREMENT_FOR_n (CELL_CLASS(IJK))
! BOTTOM_OF(IJK)   = IJK + INCREMENT_FOR_b (CELL_CLASS(IJK))
! TOP_OF   (IJK)   = IJK + INCREMENT_FOR_t (CELL_CLASS(IJK))

!      WEST_OF  (IJK)   = WEST_ARRAY_OF(IJK)
!      EAST_OF  (IJK)   = EAST_ARRAY_OF(IJK)
!      SOUTH_OF  (IJK)   = SOUTH_ARRAY_OF(IJK)
!      NORTH_OF  (IJK)   = NORTH_ARRAY_OF(IJK)
!      BOTTOM_OF  (IJK)   = BOTTOM_ARRAY_OF(IJK)
!      TOP_OF  (IJK)   = TOP_ARRAY_OF(IJK)

! Function for calculating IJKE:  EAST_OF, EAST_OF_0
! Returns IPJK if IPJK is not a wall cell else IJK
      INTEGER FUNCTION EAST_OF  (IJK)
      USE indices, only: increment_for_nb, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      EAST_OF = IJK + INCREMENT_FOR_NB (1,CELL_CLASS(IJK))
      END FUNCTION EAST_OF

! Function for calculating IJKW:  WEST_OF, WEST_OF_0
! Returns IMJK if IMJK is not a wall cell else IJK
      INTEGER FUNCTION WEST_OF  (IJK)
      USE indices, only: increment_for_nb, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      WEST_OF = IJK + INCREMENT_FOR_NB (2,CELL_CLASS(IJK))
      END FUNCTION WEST_OF

! Function for calculating IJKN:  NORTH_OF, NORTH_OF_0
! Returns IJPK if IJPK is not a wall cell else IJK
      INTEGER FUNCTION NORTH_OF (IJK)
      USE indices, only: increment_for_nb, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      NORTH_OF = IJK + INCREMENT_FOR_NB (4,CELL_CLASS(IJK))
      END FUNCTION NORTH_OF

! Function for calculating IJKS:  SOUTH_OF, SOUTH_OF_0
! Returns IJMK if IJMK is not a wall cell else IJK
      INTEGER FUNCTION SOUTH_OF (IJK)
      USE indices, only: increment_for_nb, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      SOUTH_OF = IJK + INCREMENT_FOR_NB (3,CELL_CLASS(IJK))
      END FUNCTION SOUTH_OF

! Function for calculating IJKT:  TOP_OF, TOP_OF_0
! Returns IJKP if IJKP is not a wall cell else IJK
      INTEGER FUNCTION TOP_OF   (IJK)
      USE indices, only: increment_for_nb, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      TOP_OF = IJK + INCREMENT_FOR_NB (6,CELL_CLASS(IJK))
      END FUNCTION TOP_OF

! Function for calculating IJKB:  BOTTOM_OF, BOTTOM_OF_0
! Returns IJKM if IJKM is not a wall cell else IJK
      INTEGER FUNCTION BOTTOM_OF(IJK)
      USE indices, only: increment_for_nb, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      BOTTOM_OF = IJK + INCREMENT_FOR_NB (5,CELL_CLASS(IJK))
      END FUNCTION BOTTOM_OF


      INTEGER FUNCTION WEST_OF_0  (IJK)
      USE indices, only: increment_for_w, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      WEST_OF_0 = IJK + INCREMENT_FOR_w (CELL_CLASS(IJK))
      END FUNCTION WEST_OF_0

      INTEGER FUNCTION EAST_OF_0  (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      EAST_OF_0 = IJK + INCREMENT_FOR_e (CELL_CLASS(IJK))
      END FUNCTION EAST_OF_0

      INTEGER FUNCTION SOUTH_OF_0 (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      SOUTH_OF_0 = IJK + INCREMENT_FOR_s (CELL_CLASS(IJK))
      END FUNCTION SOUTH_OF_0

      INTEGER FUNCTION NORTH_OF_0 (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      NORTH_OF_0 = IJK + INCREMENT_FOR_n (CELL_CLASS(IJK))
      END FUNCTION NORTH_OF_0

      INTEGER FUNCTION BOTTOM_OF_0(IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      BOTTOM_OF_0 = IJK + INCREMENT_FOR_b (CELL_CLASS(IJK))
      END FUNCTION BOTTOM_OF_0

      INTEGER FUNCTION TOP_OF_0   (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      TOP_OF_0 = IJK + INCREMENT_FOR_t (CELL_CLASS(IJK))
      END FUNCTION TOP_OF_0

!---------------------------------------------------------------------//
! IM_OF  (IJK)     = IJK + INCREMENT_FOR_im(CELL_CLASS(IJK))
! IP_OF  (IJK)     = IJK + INCREMENT_FOR_ip(CELL_CLASS(IJK))
! JM_OF (IJK)      = IJK + INCREMENT_FOR_jm(CELL_CLASS(IJK))
! JP_OF (IJK)      = IJK + INCREMENT_FOR_jp(CELL_CLASS(IJK))
! KM_OF(IJK)       = IJK + INCREMENT_FOR_km(CELL_CLASS(IJK))
! KP_OF   (IJK)    = IJK + INCREMENT_FOR_kp(CELL_CLASS(IJK))
!      IM_OF  (IJK)   = IM_ARRAY_OF(IJK)
!      IP_OF  (IJK)   = IP_ARRAY_OF(IJK)
!      JM_OF  (IJK)   = JM_ARRAY_OF(IJK)
!      JP_OF  (IJK)   = JP_ARRAY_OF(IJK)
!      KM_OF  (IJK)   = KM_ARRAY_OF(IJK)
!      KP_OF  (IJK)   = KP_ARRAY_OF(IJK)

! Function for calculating IMJK:  IM_OF, IM_OF_0
! Returns composite ijk index for i-1, j, k
      INTEGER FUNCTION IM_OF  (IJK)
      USE indices, only: increment_for_mp, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      IM_OF = IJK + INCREMENT_FOR_MP(1,CELL_CLASS(IJK))
      END FUNCTION IM_OF

! Function for calculating IPJK:  IP_OF, IP_OF_0
! Returns composite ijk index for i+1, j, k
      INTEGER FUNCTION IP_OF  (IJK)
      USE indices, only: increment_for_mp, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      IP_OF = IJK + INCREMENT_FOR_MP(2,CELL_CLASS(IJK))
      END FUNCTION IP_OF

! Function for calculating IJMK:  JM_OF, JM_OF_0
! Returns composite ijk index for i, j-1, k
      INTEGER FUNCTION JM_OF (IJK)
      USE indices, only: increment_for_mp, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      JM_OF = IJK + INCREMENT_FOR_MP(3,CELL_CLASS(IJK))
      END FUNCTION JM_OF

! Function for calculating IJPK:  JP_OF, JP_OF_0
! Returns composite ijk index for i, j+1, k
      INTEGER FUNCTION JP_OF (IJK)
      USE indices, only: increment_for_mp, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      JP_OF = IJK + INCREMENT_FOR_MP(4,CELL_CLASS(IJK))
      END FUNCTION JP_OF

! Function for calculating IJKM:  KM_OF, KM_OF_0
! Returns composite ijk index for i, j, k-1
      INTEGER FUNCTION KM_OF(IJK)
      USE indices, only: increment_for_mp, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      KM_OF = IJK + INCREMENT_FOR_MP(5,CELL_CLASS(IJK))
      END FUNCTION KM_OF

! Function for calculating IJKP:  KP_OF, KP_OF_0
! Returns composite ijk index for i, j, k+1
      INTEGER FUNCTION KP_OF   (IJK)
      USE indices, only: increment_for_mp, cell_class
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      KP_OF = IJK + INCREMENT_FOR_MP(6,CELL_CLASS(IJK))
      END FUNCTION KP_OF

      INTEGER FUNCTION IM_OF_0  (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      IM_OF_0 = IJK + INCREMENT_FOR_im(CELL_CLASS(IJK))
      END FUNCTION IM_OF_0

      INTEGER FUNCTION IP_OF_0  (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      IP_OF_0 = IJK + INCREMENT_FOR_ip(CELL_CLASS(IJK))
      END FUNCTION IP_OF_0

      INTEGER FUNCTION JM_OF_0 (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      JM_OF_0 = IJK + INCREMENT_FOR_jm(CELL_CLASS(IJK))
      END FUNCTION JM_OF_0

      INTEGER FUNCTION JP_OF_0 (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      JP_OF_0 = IJK + INCREMENT_FOR_jp(CELL_CLASS(IJK))
      END FUNCTION JP_OF_0

      INTEGER FUNCTION KM_OF_0(IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      KM_OF_0 = IJK + INCREMENT_FOR_km(CELL_CLASS(IJK))
      END FUNCTION KM_OF_0

      INTEGER FUNCTION KP_OF_0   (IJK)
      USE indices
      IMPLICIT NONE
      INTEGER IJK
      KP_OF_0 = IJK + INCREMENT_FOR_kp(CELL_CLASS(IJK))
      END FUNCTION KP_OF_0


! logical function to identify various fluid/flow cells
!---------------------------------------------------------------------//
! logical function to identify a fluid cell
      LOGICAL FUNCTION FLUID_AT(IJK)
      USE geometry, only: flag
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      FLUID_AT    = FLAG(IJK) .EQ. 1
      END FUNCTION FLUID_AT

! logical function to identify a specified pressure inflow cell
      LOGICAL FUNCTION P_FLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      P_FLOW_AT = FLAG(IJK) .EQ. 10 .OR. &
         FLAG(IJK) .EQ. 11
      END FUNCTION P_FLOW_AT

! logical function to identify a specified pressure outflow cell
      LOGICAL FUNCTION P_OUTFLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      P_OUTFLOW_AT= FLAG(IJK) .EQ. 11
      END FUNCTION P_OUTFLOW_AT

! logical function to identify either a specified pressure inflow
! or outflow cell or a fluid cell (simplified check)
! FLUID_AT or P_FLOW_AT (simplified check)
      LOGICAL FUNCTION FLUIDorP_FLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FLUIDorP_FLOW_AT = FLAG(IJK) .LE. 11
      END FUNCTION FLUIDorP_FLOW_AT

! logical function to identify a specified mass outflow cell
      LOGICAL FUNCTION MASS_OUTFLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MASS_OUTFLOW_AT= FLAG(IJK) .EQ. 21
      END FUNCTION MASS_OUTFLOW_AT

! logical function to identify a specified outflow cell
      LOGICAL FUNCTION OUTFLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      OUTFLOW_AT  = FLAG(IJK) .EQ. 31
      END FUNCTION OUTFLOW_AT

! logical function to identify any type of flow in/out at cell
! pressure inflow/outflow, mass inflow/outflow or outflow
      LOGICAL FUNCTION FLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FLOW_AT     = FLAG(IJK) .GE. 10 .AND. FLAG(IJK) .LE. 31
      END FUNCTION FLOW_AT

! Logical function to identify default walls
      LOGICAL FUNCTION WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      WALL_AT     = FLAG(IJK) .GE. 100
      END FUNCTION WALL_AT

! Logical function to identify a No-slip wall cell
      LOGICAL FUNCTION NS_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NS_WALL_AT  = FLAG(IJK) .EQ. 100
      END FUNCTION NS_WALL_AT

! Logical function to identify a Free-slip wall cell
      LOGICAL FUNCTION FS_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FS_WALL_AT  = FLAG(IJK) .EQ. 101
      END FUNCTION FS_WALL_AT

! Logical function to identify a Partial-slip wall cell
      LOGICAL FUNCTION PS_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PS_WALL_AT  = FLAG(IJK) .EQ. 102
      END FUNCTION PS_WALL_AT

! Logical function to identify wall ICBC_FLAG
      LOGICAL FUNCTION WALL_ICBC_FLAG(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      WALL_ICBC_FLAG = ICBC_FLAG(IJK)(1:1) .EQ. 'W' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'S' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 's' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'c' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'C'
      END FUNCTION WALL_ICBC_FLAG

      LOGICAL FUNCTION DEFAULT_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      DEFAULT_WALL_AT = ICBC_FLAG(IJK)(2:3) .EQ. '--' .AND. &
         (ICBC_FLAG(IJK)(1:1) .NE. 'c'  .AND. &
         ICBC_FLAG(IJK)(1:1) .NE. 'C')
      END FUNCTION DEFAULT_WALL_AT


! Cyclic
!---------------------------------------------------------------------//
      LOGICAL FUNCTION CYCLIC_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      CYCLIC_AT = FLAG(IJK) .EQ. 106 .OR. &
         FLAG(IJK) .EQ. 107
      END FUNCTION CYCLIC_AT

! logical function to identify cyclic condition at east boundary
      LOGICAL FUNCTION CYCLIC_AT_E(IJK)
      USE geometry, only: flag_e
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      CYCLIC_AT_E   = FLAG_E(IJK) .EQ. 2000
      END FUNCTION CYCLIC_AT_E

! logical function to identify cyclic condition at north boundary
      LOGICAL FUNCTION CYCLIC_AT_N(IJK)
      USE geometry, only: flag_n
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      CYCLIC_AT_N   = FLAG_N(IJK) .EQ. 2000
      END FUNCTION CYCLIC_AT_N

! logical function to identify cyclic condition at top boundary
      LOGICAL FUNCTION CYCLIC_AT_T(IJK)
      USE geometry, only: flag_t
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      CYCLIC_AT_T   = FLAG_T(IJK) .EQ. 2000
      END FUNCTION CYCLIC_AT_T


! Flow boundaries
!---------------------------------------------------------------------//
! identify flow at east boundary
      LOGICAL FUNCTION FLOW_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FLOW_AT_E   = FLAG_E(IJK) .GE. 2000 .AND.&
         FLAG_E(IJK) .LE. 2011
      END FUNCTION FLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION FLOW_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FLOW_AT_N   = FLAG_N(IJK) .GE. 2000 .AND.&
         FLAG_N(IJK) .LE. 2011
      END FUNCTION FLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION FLOW_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FLOW_AT_T   = FLAG_T(IJK) .GE. 2000 .AND.&
         FLAG_T(IJK) .LE. 2011
      END FUNCTION FLOW_AT_T

! identify const. pressure flow top boundary
      LOGICAL FUNCTION PFLOW_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PFLOW_AT_E  = FLAG_E(IJK) .EQ. 2010 .OR.&
         FLAG_E(IJK) .EQ. 2011
      END FUNCTION PFLOW_AT_E

! identify const. pressure flow north boundary
      LOGICAL FUNCTION PFLOW_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PFLOW_AT_N  = FLAG_N(IJK) .EQ. 2010 .OR.&
         FLAG_N(IJK) .EQ. 2011
      END FUNCTION PFLOW_AT_N

! identify const. pressure flow east boundary
      LOGICAL FUNCTION PFLOW_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PFLOW_AT_T  = FLAG_T(IJK) .EQ. 2010 .OR.&
         FLAG_T(IJK) .EQ. 2011
      END FUNCTION PFLOW_AT_T

! identify specified flow east boundary
      LOGICAL FUNCTION MFLOW_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MFLOW_AT_E  = FLAG_E(IJK) .EQ. 2020 .OR. &
         FLAG_E(IJK) .EQ. 2021 .OR. &
         FLAG_E(IJK) .EQ. 2031
      END FUNCTION MFLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION MFLOW_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MFLOW_AT_N  = FLAG_N(IJK) .EQ. 2020 .OR. &
         FLAG_N(IJK) .EQ. 2021 .OR. &
         FLAG_N(IJK) .EQ. 2031
      END FUNCTION MFLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION MFLOW_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MFLOW_AT_T  = FLAG_T(IJK) .EQ. 2020 .OR. &
         FLAG_T(IJK) .EQ. 2021 .OR. &
         FLAG_T(IJK) .EQ. 2031
      END FUNCTION MFLOW_AT_T


! Functions to identify a impermeable and/or semi-permeable surface at
! indicated boundary (specific type of internal surface)
!---------------------------------------------------------------------//
! Logical function to identify IP (impermeable surface) at East
! of the cell
      LOGICAL FUNCTION IP_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IP_AT_E     = FLAG_E(IJK) .LT. 1000
      END FUNCTION IP_AT_E

! Logical function to identify IP (impermeable surface) at North
! of the cell
      LOGICAL FUNCTION IP_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IP_AT_N     = FLAG_N(IJK) .LT. 1000
      END FUNCTION IP_AT_N

! Logical function to identify IP (impermeable surface) at Top
! of the cell
      LOGICAL FUNCTION IP_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IP_AT_T     = FLAG_T(IJK) .LT. 1000
      END FUNCTION IP_AT_T

! Logical function to identify SP or IP (semi or impermeable surface)
! at east of the cell
      LOGICAL FUNCTION SIP_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SIP_AT_E    = (FLAG_E(IJK) .LT. 2000)
      END FUNCTION SIP_AT_E

! Logical function to identify SP or IP (semi or impermeable surface)
! at north of the cell
      LOGICAL FUNCTION SIP_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SIP_AT_N    = (FLAG_N(IJK) .LT. 2000)
      END FUNCTION SIP_AT_N

! Logical function to identify SP or IP (semi or impermeable surface)
! at top of the cell
      LOGICAL FUNCTION SIP_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SIP_AT_T    = (FLAG_T(IJK) .LT. 2000)
      END FUNCTION SIP_AT_T

! Logical function to identify SP (semi-permeable surface) at east
! of cell
      LOGICAL FUNCTION SP_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SP_AT_E     = (FLAG_E(IJK) .LT. 2000) .AND. &
         (FLAG_E(IJK) .GE. 1000)
      END FUNCTION SP_AT_E

! Logical function to identify SP (semi-permeable surface) at north
! of cell
      LOGICAL FUNCTION SP_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SP_AT_N     = (FLAG_N(IJK) .LT. 2000) .AND. &
         (FLAG_N(IJK) .GE. 1000)
      END FUNCTION SP_AT_N

! Logical function to identify SP (semi-permeable surface) at top
! of cell
      LOGICAL FUNCTION SP_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SP_AT_T     = (FLAG_T(IJK) .LT. 2000) .AND. &
         (FLAG_T(IJK) .GE. 1000)
      END FUNCTION SP_AT_T


! Logical functions concerning general internal surfaces
! Integer functions to return internal surface ID
!---------------------------------------------------------------------//
! Internal surface ID for east face
      INTEGER FUNCTION IS_ID_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_ID_AT_E = FLAG_E(IJK) - 1000
      END FUNCTION IS_ID_AT_E

! Internal surface ID for north face
      INTEGER FUNCTION IS_ID_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_ID_AT_N  = FLAG_N(IJK) - 1000
      END FUNCTION IS_ID_AT_N

! Internal surface ID for top face
      INTEGER FUNCTION IS_ID_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_ID_AT_T  = FLAG_T(IJK) - 1000
      END FUNCTION IS_ID_AT_T

! Logical function to identify IS at East of the cell
      LOGICAL FUNCTION IS_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_AT_E     = FLAG_E(IJK) .LT. 2000
      END FUNCTION IS_AT_E

! Logical function to identify IS at North of the cell
      LOGICAL FUNCTION IS_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_AT_N     = FLAG_N(IJK) .LT. 2000
      END FUNCTION IS_AT_N

! Logical function to identify IS at Top of the cell
      LOGICAL FUNCTION IS_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_AT_T     = FLAG_T(IJK) .LT. 2000
      END FUNCTION IS_AT_T

! Logical function to identify No IS at East of the cell
      LOGICAL FUNCTION NO_IS_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NO_IS_AT_E  = FLAG_E(IJK) .GE. 2000
      END FUNCTION NO_IS_AT_E

! Logical function to identify No IS at North of the cell
      LOGICAL FUNCTION NO_IS_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NO_IS_AT_N  = FLAG_N(IJK) .GE. 2000
      END FUNCTION NO_IS_AT_N

! Logical function to identify No IS at Top of the cell
      LOGICAL FUNCTION NO_IS_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NO_IS_AT_T  = FLAG_T(IJK) .GE. 2000
      END FUNCTION NO_IS_AT_T

! Misc
!---------------------------------------------------------------------//
! Function for generating the index for the entries to the upper
! triangle (excluding the diagonal) of an (L,M) matrix.
      INTEGER FUNCTION FUNLM (L1, L2)
      USE indices, only: store_lm
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L1, L2
      FUNLM = STORE_LM (L1, L2)
      END FUNCTION FUNLM

! Function that returns the maximum of zero or input
      DOUBLE PRECISION FUNCTION ZMAX(XXX)
      USE param1, only: zero
      IMPLICIT NONE
      DOUBLE PRECISION XXX
      ZMAX       = MAX(XXX, ZERO)
      END FUNCTION ZMAX

      LOGICAL FUNCTION IS_NONEXISTENT(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NONEXISTENT
        INTEGER, INTENT(IN) :: PP
        IS_NONEXISTENT = (PARTICLE_STATE(PP)==NONEXISTENT)
      END FUNCTION IS_NONEXISTENT

      LOGICAL FUNCTION IS_NORMAL(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_PARTICLE
        INTEGER, INTENT(IN) :: PP
        IS_NORMAL = (PARTICLE_STATE(PP)==NORMAL_PARTICLE)
      END FUNCTION IS_NORMAL

      LOGICAL FUNCTION IS_ENTERING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        IS_ENTERING = (PARTICLE_STATE(PP)==ENTERING_PARTICLE)
      END FUNCTION IS_ENTERING

      LOGICAL FUNCTION IS_EXITING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        IS_EXITING = (PARTICLE_STATE(PP)==EXITING_PARTICLE)
      END FUNCTION IS_EXITING

      LOGICAL FUNCTION IS_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_GHOST = (PARTICLE_STATE(PP)==NORMAL_GHOST)
      END FUNCTION IS_GHOST

      LOGICAL FUNCTION IS_ENTERING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_ENTERING_GHOST = (PARTICLE_STATE(PP)==ENTERING_GHOST)
      END FUNCTION IS_ENTERING_GHOST

      LOGICAL FUNCTION IS_EXITING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_EXITING_GHOST = (PARTICLE_STATE(PP)==EXITING_GHOST)
      END FUNCTION IS_EXITING_GHOST

      LOGICAL FUNCTION IS_ANY_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_GHOST
        USE discretelement, ONLY: ENTERING_GHOST, EXITING_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_ANY_GHOST = ((PARTICLE_STATE(PP)==NORMAL_GHOST) .OR.        &
             (PARTICLE_STATE(PP)==ENTERING_GHOST) .OR.                   &
             (PARTICLE_STATE(PP)==EXITING_GHOST))
      END FUNCTION IS_ANY_GHOST

      SUBROUTINE SET_NONEXISTENT(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NONEXISTENT
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=NONEXISTENT
      END SUBROUTINE SET_NONEXISTENT

      SUBROUTINE SET_NORMAL(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_PARTICLE
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=NORMAL_PARTICLE
      END SUBROUTINE SET_NORMAL

      SUBROUTINE SET_ENTERING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=ENTERING_PARTICLE
      END SUBROUTINE SET_ENTERING

      SUBROUTINE SET_EXITING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=EXITING_PARTICLE
      END SUBROUTINE SET_EXITING

      SUBROUTINE SET_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_GHOST
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=NORMAL_GHOST
      END SUBROUTINE SET_GHOST

      SUBROUTINE SET_ENTERING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_GHOST
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=ENTERING_GHOST
      END SUBROUTINE SET_ENTERING_GHOST

      SUBROUTINE SET_EXITING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_GHOST
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=EXITING_GHOST
      END SUBROUTINE SET_EXITING_GHOST

      INTEGER FUNCTION BOUND_FUNIJK(pLI, pLJ, pLK)
        USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
        IMPLICIT NONE
        INTEGER          pLI, pLJ, pLK
        BOUND_FUNIJK  = FUNIJK ( MIN( IEND3, MAX (ISTART3, pLI) ),&
             MIN( JEND3, MAX (JSTART3, pLJ) ),&
             MIN( KEND3, MAX (KSTART3, pLK) ) )
      END FUNCTION BOUND_FUNIJK

      END MODULE functions
