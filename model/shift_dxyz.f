!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SHIFT_DXYZ                                             C
!  Purpose:  shift the data in the dx,dy,dz arrays from 1:IMAX to      C
!            IMIN1:IMAX1,  1:JMAX to JMIN1:JMAX1 ,                     C
!            1:KMAX to KMIN1:KMAX1                                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX, IMAX1, IMAX2, JMAX, JMAX1, JMAX2, KMAX  C
!                        KMAX1 , KMAX2, IMIN1, JMIN1, KMIN1, NO_I,     C
!                        NO_J, NO_K                                    C
!  Variables modified:  DX, DY, DZ                                     C
!                                                                      C
!  Local variables: LC                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SHIFT_DXYZ
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!              loop counter
      INTEGER :: LC
!-----------------------------------------------
!
!
      IF (DO_I) THEN
         DX(IMAX3) = DX(IMAX-1)
         DX(IMAX2) = DX(IMAX-1)
         DO LC = IMAX1, IMIN1, -1
            DX(LC) = DX(LC-2)
         ENDDO
         DX(IMIN2) = DX(IMIN1)
         DX(IMIN3) =DX(IMIN2)
      ENDIF
!

      IF (DO_J) THEN
         DY(JMAX3) = DY(JMAX-1)
         DY(JMAX2) = DY(JMAX-1)
         DO LC = JMAX1, JMIN1, -1
            DY(LC) = DY(LC-2)
         ENDDO
         DY(JMIN2) = DY(JMIN1)
         DY(JMIN3) =DY(JMIN2)

      ENDIF
!
      IF (DO_K) THEN

         DZ(KMAX3) = DZ(KMAX-1)
         DZ(KMAX2) = DZ(KMAX-1)
         DO LC = KMAX1, KMIN1, -1
            DZ(LC) = DZ(LC-2)
         ENDDO
         DZ(KMIN2) = DZ(KMIN1)
         DZ(KMIN3) =DZ(KMIN2)
      ENDIF
!
      RETURN
      END SUBROUTINE SHIFT_DXYZ

!// Comments on the modifications for DMP version implementation
!// 120 Added new initializations at IMAX2, JMAX2, KMAX2 etc.
