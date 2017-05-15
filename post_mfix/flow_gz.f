!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VFLOW_gz(I, J, K, IJK)                                 C
!  Purpose: Calculate volumetric flow of gas in z direction            C
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
      REAL FUNCTION VFLOW_gz(I, J, K, IJK)
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
      INTEGER I, J, K, IJK, IJKP
!
      IF(W_g(IJK) .GT. ZERO) THEN
        VFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJK)
      ELSE
        IJKP = KP_OF(IJK)
        VFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJKP)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFLOW_gz(I, J, K, IJK)                                 C
!  Purpose: Calculate mass flow of gas in z direction                  C
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
      REAL FUNCTION MFLOW_gz(I, J, K, IJK)
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

      INTEGER I, J, K, IJK, IJKP
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(W_g(IJK) .GT. ZERO) THEN
        MFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK)
      ELSE
        IJKP = KP_OF(IJK)
        MFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJKP)&
                   * CALC_RO_g(IJKP)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLUX_gz(IJK)                                           C
!  Purpose: Calculate mass flux of gas in z direction                  C
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
      REAL FUNCTION FLUX_gz(IJK)
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

      INTEGER IJK, IJKP
!
!  Function subroutines
!
      REAL CALC_RO_g

      IF(W_g(IJK) .GT. ZERO) THEN
        FLUX_gz = W_g(IJK) * EP_g(IJK) &
                   * CALC_RO_g(IJK)
      ELSE
        IJKP = KP_OF(IJK)
        FLUX_gz = W_g(IJK) * EP_g(IJKP)&
                   * CALC_RO_g(IJKP)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XFLOW_gz(I, J, K, IJK, N)                              C
!  Purpose: Calculate gas species mass flow in z direction             C
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
      REAL FUNCTION XFLOW_gz(I, J, K, IJK, N)
!
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE
!
      INTEGER I, J, K, IJK, N, IJKP
!
!  Function subroutines
!
      REAL CALC_RO_g
!
      IF(W_g(IJK) .GT. ZERO) THEN
        XFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJK)&
                   * CALC_RO_g(IJK) * X_g(IJK, N)
      ELSE
        IJKP = KP_OF(IJK)
        XFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJKP)&
                   * CALC_RO_g(IJKP) * X_g(IJKP, N)
      ENDIF
!
      RETURN
      END
