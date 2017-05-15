!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VFLOW_sy(I, J, K, IJK, M)                              C
!  Purpose: Calculate volumetric flow of solids in x direction         C
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
      REAL FUNCTION VFLOW_sy(I, J, K, IJK, M)
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

      INTEGER I, J, K, IJK, M, IJPK

      IF(V_s(IJK, M) .GT. ZERO) THEN
        VFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * EP_s(IJK, M)
      ELSE
        IJPK = JP_OF(IJK)
        VFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * EP_s(IJPK, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFLOW_sy(I, J, K, IJK, M)                              C
!  Purpose: Calculate mass flow of solids in y direction               C
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
      REAL FUNCTION MFLOW_sy(I, J, K, IJK, M)
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

      INTEGER I, J, K, IJK, M, IJPK

      IF(V_s(IJK, M) .GT. ZERO) THEN
        MFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * ROP_s(IJK, M)
      ELSE
        IJPK = JP_OF(IJK)
        MFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * ROP_s(IJPK, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLUX_sy(IJK, M)                                        C
!  Purpose: Calculate mass flux of solids in y direction               C
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
      REAL FUNCTION FLUX_sy(IJK, M)
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

      INTEGER IJK, M, IJPK

      IF(V_s(IJK, M) .GT. ZERO) THEN
        FLUX_sy = V_s(IJK, M) * ROP_s(IJK, M)
      ELSE
        IJPK = JP_OF(IJK)
        FLUX_sy = V_s(IJK, M) * ROP_s(IJPK, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XFLOW_sy(I, J, K, IJK, M, N)                           C
!  Purpose: Calculate solids species mass flow in x direction          C
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
      REAL FUNCTION XFLOW_sy(I, J, K, IJK, M, N)
!
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER I, J, K, IJK, M, N, IJPK

      IF(V_s(IJK, M) .GT. ZERO) THEN
        XFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M)&
                   * ROP_s(IJK, M) * X_s(IJK, M, N)
      ELSE
        IJPK = JP_OF(IJK)
        XFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M)&
                   * ROP_s(IJPK, M) * X_s(IJPK, M, N)
      ENDIF
!
      RETURN
      END
