!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RO_g(L)                                           C
!  Purpose: Calculate gas density                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 22-NOV-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      REAL FUNCTION CALC_RO_g(L)
!
      USE compar
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop

      IMPLICIT NONE
!
!              Passed value of IJK index
      INTEGER  L, IJK
!
      DOUBLE PRECISION MW
!
        IF(RO_g0 .EQ. UNDEFINED .AND. .NOT.WALL_AT(L)) THEN
          IF(MW_AVG .EQ. UNDEFINED) THEN
            MW = CALC_MW(X_g, DIMENSION_3, L, NMAX(0), MW_g)
          ELSE
            MW = MW_AVG
          ENDIF
          CALC_RO_g = EOSG(MW, P_g(L), T_g(L))
        ELSE
          CALC_RO_g = RO_g0
        ENDIF
!
      RETURN
      END
