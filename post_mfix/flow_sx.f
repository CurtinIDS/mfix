!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VFLOW_sx(I, J, K, IJK, M)                              C
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
      REAL FUNCTION VFLOW_sx(I, J, K, IJK, M)
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

      INTEGER I, J, K, IJK, M, IPJK

      IF(U_s(IJK, M) .GT. ZERO) THEN
        VFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * EP_s(IJK, M)
      ELSE
        IPJK = IP_OF(IJK)
        VFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * EP_s(IPJK, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFLOW_sx(I, J, K, IJK, M)                              C
!  Purpose: Calculate mass flow of solids in x direction               C
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
      REAL FUNCTION MFLOW_sx(I, J, K, IJK, M)
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

      INTEGER I, J, K, IJK, M, IPJK

      IF(U_s(IJK, M) .GT. ZERO) THEN
        MFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * ROP_s(IJK, M)
      ELSE
        IPJK = IP_OF(IJK)
        MFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * ROP_s(IPJK, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLUX_sx(IJK, M)                                        C
!  Purpose: Calculate mass flux of solids in x direction               C
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
      REAL FUNCTION FLUX_sx(IJK, M)
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

      INTEGER IJK, M, IPJK

      IF(U_s(IJK, M) .GT. ZERO) THEN
        FLUX_sx = U_s(IJK, M) * ROP_s(IJK, M)
      ELSE
        IPJK = IP_OF(IJK)
        FLUX_sx = U_s(IJK, M) * ROP_s(IPJK, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XFLOW_sx(I, J, K, IJK, M, N)                           C
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
      REAL FUNCTION XFLOW_sx(I, J, K, IJK, M, N)
!
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER I, J, K, IJK, M, N, IPJK

      IF(U_s(IJK, M) .GT. ZERO) THEN
        XFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M)&
                   * ROP_s(IJK, M) * X_s(IJK, M, N)
      ELSE
        IPJK = IP_OF(IJK)
        XFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M)&
                   * ROP_s(IPJK, M) * X_s(IPJK, M, N)
      ENDIF
!
      RETURN
      END
