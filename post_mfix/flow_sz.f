!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VFLOW_sz(I, J, K, IJK, M)                              C
!  Purpose: Calculate volumetric flow of solids in z direction         C
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
      REAL FUNCTION VFLOW_sz(I, J, K, IJK, M)
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

      INTEGER I, J, K, IJK, M, IJKP

      IF(W_s(IJK, M) .GT. ZERO) THEN
        VFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * EP_s(IJK, M)
      ELSE
        IJKP = KP_OF(IJK)
        VFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * EP_s(IJKP, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MFLOW_sz(I, J, K, IJK, M)                              C
!  Purpose: Calculate mass flow of solids in z direction               C
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
      REAL FUNCTION MFLOW_sz(I, J, K, IJK, M)
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

      INTEGER I, J, K, IJK, M, IJKP

      IF(W_s(IJK, M) .GT. ZERO) THEN
        MFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * ROP_s(IJK, M)
      ELSE
        IJKP = KP_OF(IJK)
        MFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * ROP_s(IJKP, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FLUX_sz(IJK, M)                                        C
!  Purpose: Calculate mass flux of solids in z direction               C
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
      REAL FUNCTION FLUX_sz(IJK, M)
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

      INTEGER IJK, M, IJKP

      IF(W_s(IJK, M) .GT. ZERO) THEN
        FLUX_sz = W_s(IJK, M) * ROP_s(IJK, M)
      ELSE
        IJKP = KP_OF(IJK)
        FLUX_sz = W_s(IJK, M) * ROP_s(IJKP, M)
      ENDIF
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: XFLOW_sz(I, J, K, IJK, M, N)                           C
!  Purpose: Calculate solids species mass flow in z direction          C
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
      REAL FUNCTION XFLOW_sz(I, J, K, IJK, M, N)
!
      Use param1
      Use fldvar
      Use indices
      Use physprop
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE

      INTEGER I, J, K, IJK, M, N, IJKP

      IF(W_s(IJK, M) .GT. ZERO) THEN
        XFLOW_sz = DX(I) * DY(J) * W_s(IJK, M)  &
                   * ROP_s(IJK, M) * X_s(IJK, M, N)
      ELSE
        IJKP = KP_OF(IJK)
        XFLOW_sz = DX(I) * DY(J) * W_s(IJK, M)&
                   * ROP_s(IJKP, M) * X_s(IJKP, M, N)
      ENDIF
!
      RETURN
      END
