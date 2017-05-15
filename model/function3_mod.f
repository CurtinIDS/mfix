MODULE function3

  USE functions

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  !
  !                      Function for generating the LOCAL 3-D array index IJK
  !                      from the 1-D indices I, J, and K.
  !//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
  !      INTEGER          FUNIJK3
  !                      Function for generating the LOCAL 3-D array index IJK
  !                      from the 1-D indices I, J, K and IPROC.
  !      INTEGER          FUNIJK3_PROC
  !                      Function for generating the GLOBAL 3-D array index IJK
  !                      from the 1-D indices I, J, and K.
  !      INTEGER          FUNIJK3_GL
  !
  !                      Function for calculating IMJK
  !      INTEGER          IM3_OF
  !
  !                      Function for calculating IPJK
  !      INTEGER          IP3_OF
  !
  !                      Function for calculating IJMK
  !      INTEGER          JM3_OF
  !
  !                      Function for calculating IJPK
  !      INTEGER          JP3_OF
  !
  !                      Function for calculating IJKM
  !      INTEGER          KM3_OF
  !
  !                      Function for calculating IJKP
  !      INTEGER          KP3_OF
  !
  !                      Logical function to identify a fluid cell
  !      LOGICAL          FLUID3_AT
  !
  !                      Logical function to identify a cyclic cell
  !      LOGICAL          CYCLIC3_AT
  !
  !                      DUMMY INDICES
  !      INTEGER          LI3, LJ3, LK3, LIPROC3, IJK3
  !

contains

  INTEGER FUNCTION FUNIJK3(LI3,LJ3,LK3)
    USE compar
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI3,LJ3,LK3
    FUNIJK3 = LJ3 + c0_3 + LI3*c1_3 + LK3*c2_3
  END FUNCTION FUNIJK3

  INTEGER FUNCTION FUNIJK3_PROC(LI3, LJ3, LK3, LIPROC3)
    USE compar
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI3,LJ3,LK3,LIPROC3
    FUNIJK3_PROC = 1 + (LJ3 - jstart4_all(LIPROC3))+ &
         (LI3-Istart4_all(LIPROC3))*(jend4_all(LIPROC3)-jstart4_all(LIPROC3)+1) &
         + (LK3-kstart4_all(LIPROC3))*(jend4_all(LIPROC3)-jstart4_all(LIPROC3)+1)* &
         (iend4_all(LIPROC3)-istart4_all(LIPROC3)+1)
  END FUNCTION FUNIJK3_PROC

  INTEGER FUNCTION FUNIJK3_GL (LI3, LJ3, LK3)
    USE geometry
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LI3,LJ3,LK3
    FUNIJK3_GL = 1 + (LJ3 - jmin4) + (LI3-imin4)*(jmax4-jmin4+1) &
         + (LK3-kmin4)*(jmax4-jmin4+1)*(imax4-imin4+1)
  END FUNCTION FUNIJK3_GL

  INTEGER FUNCTION IM3_OF  (IJK3)
    USE indices
    IMPLICIT NONE
    INTEGER IJK3
    IM3_OF       = IJK3 + INCREMENT3_FOR_im(CELL_CLASS3(IJK3))
  END FUNCTION IM3_OF

  INTEGER FUNCTION IP3_OF  (IJK3)
    USE indices
    IMPLICIT NONE
    INTEGER IJK3
    IP3_OF       = IJK3 + INCREMENT3_FOR_ip(CELL_CLASS3(IJK3))
  END FUNCTION IP3_OF

  INTEGER FUNCTION JM3_OF (IJK3)
    USE indices
    IMPLICIT NONE
    INTEGER IJK3
    JM3_OF       = IJK3 + INCREMENT3_FOR_jm(CELL_CLASS3(IJK3))
  END FUNCTION JM3_OF

  INTEGER FUNCTION JP3_OF (IJK3)
    USE indices
    IMPLICIT NONE
    INTEGER IJK3
    JP3_OF       = IJK3 + INCREMENT3_FOR_jp(CELL_CLASS3(IJK3))
  END FUNCTION JP3_OF

  INTEGER FUNCTION KM3_OF(IJK3)
    USE indices
    IMPLICIT NONE
    INTEGER IJK3
    KM3_OF       = IJK3 + INCREMENT3_FOR_km(CELL_CLASS3(IJK3))
  END FUNCTION KM3_OF

  INTEGER FUNCTION KP3_OF   (IJK3)
    USE indices
    IMPLICIT NONE
    INTEGER IJK3
    KP3_OF       = IJK3 + INCREMENT3_FOR_kp(CELL_CLASS3(IJK3))
  END FUNCTION KP3_OF

  LOGICAL FUNCTION FLUID3_AT(IJK3)
    USE geometry
    IMPLICIT NONE
    INTEGER IJK3
    FLUID3_AT    = FLAG3(IJK3) .EQ. 1
  END FUNCTION FLUID3_AT

  LOGICAL FUNCTION CYCLIC3_AT(IJK3)
    USE geometry
    IMPLICIT NONE
    INTEGER IJK3
    CYCLIC3_AT   = FLAG3(IJK3) .EQ. 106 .OR. &
         FLAG3(IJK3) .EQ. 107
  END FUNCTION CYCLIC3_AT

  INTEGER FUNCTION BOUND_FUNIJK3(pLI, pLJ, pLK)
    USE compar, only: istart4, iend4, jstart4, jend4, kstart4, kend4
    IMPLICIT NONE
    INTEGER          pLI, pLJ, pLK
    BOUND_FUNIJK3  = FUNIJK3 ( MIN( IEND4, MAX (ISTART4, pLI) ),&
         MIN( JEND4, MAX (JSTART4, pLJ) ),&
         MIN( KEND4, MAX (KSTART4, pLK) ) )
  END FUNCTION BOUND_FUNIJK3

END MODULE function3
