!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_PIC_MO                                           !
!                                                                      !
!                                                                      !
!  Author: R.Garg                                     Date: 23-June-14 !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_PIC_MO

      use bc, only: BC_PLANE
      use bc, only: BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t

      use pic_bc, only: PIC_BCMO, PIC_BCMO_MAP, PIC_BCMO_IJK
      use pic_bc, only: PIC_BCMO_IJKSTART, PIC_BCMO_IJKEND

      use funits, only: DMP_LOG

      use mpi_utility
      use error_manager
      use functions

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: BCV, BCV_I      ! BC loop counter

      INTEGER :: LC

      LOGICAL, parameter :: setDBG = .false.
      LOGICAL :: dFlag

      INTEGER :: MAX_CELLS, BND1, BND2

      INTEGER, ALLOCATABLE :: LOC_PIC_BCMO_IJK(:)

      INTEGER :: I,J,K,IJK
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t

      CALL INIT_ERR_MSG("SET_BC_PIC_MO")

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'PIC outlet count: ',I4)") PIC_BCMO

! Loop over the outflow BCs to get an approximate count of the number
! of fluid cells that are adjacent to the outlet.
      MAX_CELLS = 0
      DO BCV_I = 1, PIC_BCMO
         BCV = PIC_BCMO_MAP(BCV_I)

! Set the search area to the dimensions of the inlet.
         if(dFlag) WRITE(*,"(/2x,'Adding cells for BCV: ',I3)") BCV
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S')
            BND1 = BC_I_e(BCV) - BC_I_w(BCV)
            BND2 = BC_K_t(BCV) - BC_K_b(BCV)

         CASE('E','W')
            BND1 = BC_J_n(BCV) - BC_J_s(BCV)
            BND2 = BC_K_t(BCV) - BC_K_b(BCV)

         CASE('T','B')
            BND1 = BC_I_e(BCV) - BC_I_w(BCV)
            BND2 = BC_J_n(BCV) - BC_J_s(BCV)
         END SELECT

         MAX_CELLS = MAX_CELLS +                                      &
            2*(BND1+1)*(BND2+1) + 2*(BND1+2) + 2*(BND2+2)

         if(dFlag) WRITE(*,"(4x,'Plane:   ',A)") BC_PLANE(BCV)
         if(dFlag) WRITE(*,"(4x,'Cells: ', I8)") (BND1+1)*(BND2+1)
      ENDDO

      if(dFlag) write(*,"(2x,'Max Cells: ',I8)") MAX_CELLS

! Allocate an array to hold the IJK values. This array should be
! more than enough to store all the IJKs.
      allocate( LOC_PIC_BCMO_IJK(MAX_CELLS) )

! Loop over the IJKs for each BC and store only the IJKs that you
! own as well as the start/end array positions for each BC.
      LC = 1
      DO BCV_I = 1, PIC_BCMO

         PIC_BCMO_IJKSTART(BCV_I) = LC
         BCV = PIC_BCMO_MAP(BCV_I)

         if(dFlag) write(*,"(/2x,'Searching for fluid cells:',I3)") BCV

         I_w = BC_I_w(BCV); I_e = BC_I_e(BCV)
         J_s = BC_J_s(BCV); J_n = BC_J_n(BCV)
         K_b = BC_K_b(BCV); K_t = BC_K_t(BCV)

! Depending on the flow plane, the 'common' index needs shifted to
! reference the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); J_s = J_s+1;  J_n = J_s
         CASE('S'); J_s = J_s-1;  J_n = J_s
         CASE('E'); I_w = I_w+1;  I_e = I_w
         CASE('W'); I_w = I_w-1;  I_e = I_w
         CASE('T'); K_b = K_b+1;  K_t = K_b
         CASE('B'); K_b = K_b-1;  K_t = K_b
        END SELECT

         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
! Limit only to fluid cells
            IF(.NOT.IS_ON_myPE_wobnd(I,J,K)) CYCLE
            IJK = FUNIJK(I,J,K)

            IF(.NOT.FLUID_AT(IJK)) CYCLE

            LOC_PIC_BCMO_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         PIC_BCMO_IJKEND(BCV_I) = LC-1

         if(dFLAG) write(*,1111) BCV, BCV_I,                           &
            PIC_BCMO_IJKSTART(BCV_I),PIC_BCMO_IJKEND(BCV_I)

      ENDDO

 1111 FORMAT(/2x,'PIC Mass Outflow:',/4x,'BC:',I4,3x,'MAP:',I4,&
         /4x,'IJKSTART:',I6,/4x,'IJKEND:  ',I6)

! Allocate the global store arrary array. This changes across MPI ranks.
      IF(LC > 1) THEN
         allocate( PIC_BCMO_IJK(LC-1) )
         PIC_BCMO_IJK(1:LC-1) = LOC_PIC_BCMO_IJK(1:LC-1)
      ELSE
         allocate( PIC_BCMO_IJK(1) )
         PIC_BCMO_IJK(1) = LOC_PIC_BCMO_IJK(1)
      ENDIF

      deallocate(LOC_PIC_BCMO_IJK)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_PIC_MO


