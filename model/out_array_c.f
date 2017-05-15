!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_C (ARRAY,MESSAGE)                            C
!  Purpose: print out a 3D array to standard output (character)        C
!                                                                      C
!  Author: P.Nicoletti                                Date: 10-JAN-92  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: KMAX2                                         C
!  Variables modified: K                                               C
!                                                                      C
!  Local variables: IJK                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_ARRAY_C(ARRAY, MESSAGE)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE indices
      USE funits
      USE compar
      USE functions
      USE in_binary_512i
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                       array to print out
      CHARACTER(LEN=4) :: ARRAY(*)
!
!                       message to print out
      CHARACTER(LEN=*) :: MESSAGE
!
! local variables
!
!                       pointer into array (points to start of a k-plane)
      INTEGER           IJK
!
!                       K
      INTEGER           K

      character(LEN=4),  allocatable :: array1c(:)

!-----------------------------------------------
!
!//d      call lock_tmp_array

      allocate (array1c(ijkmax2))
      call convert_to_io_c(array,array1c,ijkmax2)
!
      DO K = 1, KMAX2
         IJK = FUNIJK_IO(1,1,K)
!
         WRITE (UNIT_OUT, 1100) MESSAGE, K
         CALL OUT_ARRAY_KC (ARRAY1C(IJK))
      END DO
 1100 FORMAT(/,1X,A,' at K = ',I4,/)
!
!//d      call unlock_tmp_array
      deallocate (array1c)

      RETURN
      END SUBROUTINE OUT_ARRAY_C

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization, array1c(ijkmax2)
