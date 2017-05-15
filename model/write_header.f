!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_HEADER                                           C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-APR-97  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_HEADER
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE funits
      USE machine
      USE output
      USE param
      USE param1
      USE run
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
!
!
!                      Memory required for the run
      DOUBLE PRECISION :: MEMORY
!-----------------------------------------------
!
!
      CALL START_LOG
!
      IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
      IF(DMP_LOG)WRITE (UNIT_LOG, 1005) ID_VERSION, ID_NODE
      IF(DMP_LOG)WRITE (UNIT_LOG,1010)RUN_NAME,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR
!
      IF (FULL_LOG .and. myPE.eq.PE_IO) THEN    !//d
         WRITE (*, *) ' '
         WRITE (*, 1005) ID_VERSION, ID_NODE
         WRITE(*,1010)RUN_NAME,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR
      ENDIF
!
!   Calculate the memory requirement for the present run
!
      MEMORY = 9. + (8.*DIMENSION_3/ONEMEG)*(95. + 32.*DIMENSION_M + 4.*&
         DIMENSION_N_G + 4.*DIMENSION_M*DIMENSION_N_S)
      IF(DMP_LOG)WRITE (UNIT_LOG, '(1X,A,F7.2,A)') 'Memory required: ', MEMORY, ' Mb'
      IF (FULL_LOG .and. myPE.eq.PE_IO) THEN     !//d
         WRITE (*, '(1X,A,F7.2,A)') 'Memory required: ', MEMORY, ' Mb'
         WRITE (*, 1015)
      ENDIF
!
      IF(DMP_LOG)WRITE (UNIT_LOG, 1015)
      CALL END_LOG
!
      RETURN
 1005 FORMAT(1X,'MFIX (',A10,') simulation on computer: ',A20)
 1010 FORMAT(1X,'Run name: ',A20,2X,'Time: ',I2,':',I2.0,20X,'Date: ',I2,'-',I2&
         ,'-',I4)
 1015 FORMAT(72('_'))
      END SUBROUTINE WRITE_HEADER
