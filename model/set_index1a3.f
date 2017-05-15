!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_INDEX1A3(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK,      C
!                            IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,       C
!                            IJKB, IJKT)                               C
!  Purpose: Set the indices of the neighbors of cell ijk (brute force) C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modify index computations for K for setting periodic       C
!           boundary conditions in a cylindrical geometry where z goes C
!           from 0 to 2 pi                                             C
!  Author: M. Syamlal                                 Date: 10-MAR-92  C
!  Revision Number: 2                                                  C
!  Purpose:  Calculate only the nearest neighbor indices.( for code    C
!            optimization)                                             C
!  Author: M. Syamlal                                 Date: 23-SEP-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: I, J, K, IJK                                  C
!                                                                      C
!  Variables modified: IJKM, IJMK, IMJK, IPJK, IJPK, IJKP, IJKW, IJKE, C
!                      IJKS, IJKN, IJKB, IJKT                          C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_INDEX1A3(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
         IJKW, IJKE, IJKS, IJKN, IJKB, IJKT)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE functions
      USE function3
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, IJKW, IJKE, &
         IJKS, IJKN, IJKB, IJKT
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------

      IMJK = UNDEFINED_I
      IPJK = UNDEFINED_I
      IJMK = UNDEFINED_I
      IJPK = UNDEFINED_I
      IJKM = UNDEFINED_I
      IJKP = UNDEFINED_I


      IF(IM1_3(I).NE.UNDEFINED_I) THEN
        IMJK = BOUND_FUNIJK3(IM1_3(I),J,K)
      ENDIF

      IF(IP1_3(I).NE.UNDEFINED_I) THEN
        IPJK = BOUND_FUNIJK3(IP1_3(I),J,K)
      ENDIF

      IF(JM1_3(J).NE.UNDEFINED_I) THEN
        IJMK = BOUND_FUNIJK3(I,JM1_3(J),K)
      ENDIF

      IF(JP1_3(J).NE.UNDEFINED_I) THEN
        IJPK = BOUND_FUNIJK3(I,JP1_3(J),K)
      ENDIF

      IF(KM1_3(K).NE.UNDEFINED_I) THEN
        IJKM = BOUND_FUNIJK3(I,J,KM1_3(K))
      ENDIF

      IF(KP1_3(K).NE.UNDEFINED_I) THEN
        IJKP = BOUND_FUNIJK3(I,J,KP1_3(K))
      ENDIF
!
      RETURN
      END SUBROUTINE SET_INDEX1A3

!// Comments on the modifications for DMP version implementation
!// Modified calls to BOUND_FUNIJK to have a self consistent formulation
