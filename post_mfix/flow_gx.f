!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VFLOW_gx(I, J, K, IJK)                                 C
!  Purpose: Calculate volumetric flow of gas in x direction            C
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
      REAL FUNCTION VFLOW_gx(I, J, K, IJK)
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
      INTEGER I, J, K, IJK, IPJK

      IF(U_g(IJK) .GT. ZERO) THEN
        VFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IJK)
      ELSE
        IPJK = IP_OF(IJK)
        VFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IPJK)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFLOW_gx(I, J, K, IJK)                                 C
!  Purpose: Calculate mass flow of gas in x direction                  C
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
      REAL FUNCTION MFLOW_gx(I, J, K, IJK)
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

      INTEGER I, J, K, IJK, IPJK
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(U_g(IJK) .GT. ZERO) THEN
        MFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IJK)&
                   * CALC_RO_g(IJK)
      ELSE
        IPJK = IP_OF(IJK)
        MFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IPJK)&
                   * CALC_RO_g(IPJK)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLUX_gx(IJK)                                           C
!  Purpose: Calculate mass flux of gas in x direction                  C
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
      REAL FUNCTION FLUX_gx(IJK)
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

      INTEGER IJK, IPJK
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(U_g(IJK) .GT. ZERO) THEN
        FLUX_gx = U_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK)
      ELSE
        IPJK = IP_OF(IJK)
        FLUX_gx = U_g(IJK) * EP_g(IPJK)&
                   * CALC_RO_g(IPJK)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XFLOW_gx(I, J, K, IJK, N)                              C
!  Purpose: Calculate gas species mass flow in x direction             C
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
      REAL FUNCTION XFLOW_gx(I, J, K, IJK, N)
!
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER I, J, K, IJK, N, IPJK
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(U_g(IJK) .GT. ZERO) THEN
        XFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK) * X_g(IJK, N)
      ELSE
        IPJK = IP_OF(IJK)
        XFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IPJK)&
                   * CALC_RO_g(IPJK) * X_g(IPJK, N)
      ENDIF
!
      RETURN
      END
