!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK,       C
!                            IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,       C
!                           IJKB, IJKT, IM, JM, KM)                    C
!  Purpose: Set the indices of the first neighbors of cell ijk         C
!           This version adds 'increments' stored in STORE_INCREMENTS  C
!           to IJK to find indices of the neighbors.                   C
!                                                                      C
!  Author: M. Syamlal, W. A. Rogers                   Date: 17-Dec-91  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Second index of store_increments changed to PARAMETER's.   C
!           Use IM_OF etc. instead of BOUND_FUNIJK for computing       C
!           IMJK etc.                                                  C
!  Author: M. Syamlal                                 Date: 18-FEB-92  C
!  Revision Number:2                                                   C
!  Purpose: change STORE_INCREMENTS to INCREMENT_FOR_xx.  Do only      C
!           calculation of the nearest neighbors.  Remove MIN and MAX  C
!           from IM, IP etc. calculations.                             C
!  Author: M. Syamlal                                 Date: 18-SEP-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: CELL_CLASS, INCREMENTS_FOR_xx, IJK, I, J, K   C
!                                                                      C
!  Variables modified: IJKN,IJKS,IJKE,IJKW,IJKT,IJKB, IJKM, IJMK, IMJK,C
!                      IPJK, IJPK, IJKP                                C
!                                                                      C
!  Local variables:                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
         IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, IM, JM, KM)
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
      USE fldvar
      USE geometry
      USE constant
      USE indices
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                      IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, &
                      IM, JM, KM
!
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
!
      IM = IM1(I)
      JM = JM1(J)
      KM = KM1(K)
!
!
!     Determine the true indices of neighboring cells
!
!
      IJKW = WEST_OF(IJK)
      IJKE = EAST_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKB = BOTTOM_OF(IJK)
      IJKT = TOP_OF(IJK)
      IMJK = IM_OF(IJK)
      IPJK = IP_OF(IJK)
      IJMK = JM_OF(IJK)
      IJPK = JP_OF(IJK)
      IJKM = KM_OF(IJK)
      IJKP = KP_OF(IJK)
!
      RETURN
      END SUBROUTINE SET_INDEX1

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
