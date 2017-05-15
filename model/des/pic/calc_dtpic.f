!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DTPIC

      USE compar
      USE cutcell
      USE des_rxns
      USE stl_functions_des
      USE des_thermo
      USE discretelement
      USE functions
      USE funits
      USE geometry
      USE param1
      USE run
      use desmpi
      use error_manager
      use mfix_pic, only: CFL_PIC, DTPIC_CFL, DTPIC_TAUP
      use mfix_pic, only: DTPIC_MAX
      use mpi_utility

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: L, PC


! MPPIC related quantities
      DOUBLE PRECISION :: DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------

      DTPIC_CFL = LARGE_NUMBER

      PC = 1
      DO L = 1, MAX_PIP
      IF(PC.GT.PIP) EXIT
         IF(IS_NONEXISTENT(L)) CYCLE
         PC = PC+1
         IF(IS_GHOST(L) .or. IS_ENTERING_GHOST(L) .or. IS_EXITING_GHOST(L)) CYCLE

         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/&
            (ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/&
            (ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DO_K) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/&
            (ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)

         DTPIC_CFL = MIN(DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
      ENDDO

      CALL global_all_max(DTPIC_CFL)

      DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)
      DTSOLID = DTPIC_MAX

      WRITE(ERR_MSG,2000) dtpic_cfl, dtpic_taup, DTSOLID
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 2000 FORMAT('DTPIC BASED ON CFL AND TAUP:', 2x, 2(2x,g11.4),          &
         /'DTSOLID set to ', g11.4)

      RETURN
      END SUBROUTINE CALC_DTPIC
