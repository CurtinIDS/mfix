!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VFLOW_gy(I, J, K, IJK)                                 C
!  Purpose: Calculate volumetric flow of gas in y direction            C
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
      REAL FUNCTION VFLOW_gy(I, J, K, IJK)
!
      Use param
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE
!
      INTEGER I, J, K, IJK, IJPK

      IF(V_g(IJK) .GT. ZERO) THEN
        VFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJK)
      ELSE
        IJPK = JP_OF(IJK)
        VFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJPK)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFLOW_gy(I, J, K, IJK)                                 C
!  Purpose: Calculate mass flow of gas in y direction                  C
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
      REAL FUNCTION MFLOW_gy(I, J, K, IJK)
!
      Use param
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER I, J, K, IJK, IJPK
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(V_g(IJK) .GT. ZERO) THEN
        MFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK)
      ELSE
        IJPK = JP_OF(IJK)
        MFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJPK)&
                   * CALC_RO_g(IJPK)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLUX_gy(IJK)                                           C
!  Purpose: Calculate mass flux of gas in y direction                  C
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
      REAL FUNCTION FLUX_gy(IJK)
!
      Use param
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER IJK, IJPK
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(V_g(IJK) .GT. ZERO) THEN
        FLUX_gy = V_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK)
      ELSE
        IJPK = JP_OF(IJK)
        FLUX_gy = V_g(IJK) * EP_g(IJPK)&
                   * CALC_RO_g(IJPK)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XFLOW_gy(I, J, K, IJK, N)                              C
!  Purpose: Calculate gas species mass flow in y direction             C
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
      REAL FUNCTION XFLOW_gy(I, J, K, IJK, N)
!
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER I, J, K, IJK, N, IJPK
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(V_g(IJK) .GT. ZERO) THEN
        XFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK) * X_g(IJK, N)
      ELSE
        IJPK = JP_OF(IJK)
        XFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJPK)&
                   * CALC_RO_g(IJPK) * X_g(IJPK, N)
      ENDIF
!
      RETURN
      END
