!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_PIC_MI                                           !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_PIC_MI

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE bc
      USE pic_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run

      USE error_manager
      USE toleranc

      IMPLICIT NONE

      INTEGER :: BCV

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV_I      ! BC loop counter
      INTEGER M           ! Mass phase loop counter
      INTEGER PHASE_CNT        ! Number of solid phases at bc
      INTEGER PHASE_LIST(DIM_M) ! List of phases used in current bc

! the number of particles injected in a solids time step
      DOUBLE PRECISION MAX_DIA ! Max diameter of incoming particles at bc

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

!-----------------------------------------------

      CALL INIT_ERR_MSG("SET_BC_PIC_MI")

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'PIC inlet count: ',I4)") PIC_BCMI

      PIC_BCMI_INCL_CUTCELL(:) = .false.

! Loop over BCs that are flagged for PIC mass inflow.
      DO BCV_I = 1, PIC_BCMI

! Get the user defined BC ID.
         BCV = PIC_BCMI_MAP(BCV_I)

         PIC_BCMI_OFFSET (BCV_I,:) = 1
         PIC_BCMI_NORMDIR(BCV_I,:) = 0

         SELECT CASE(BC_PLANE(BCV))
         CASE('E'); PIC_BCMI_NORMDIR(BCV_I,1) = 1
         CASE('W')
            PIC_BCMI_NORMDIR(BCV_I,1) = -1
            PIC_BCMI_OFFSET (BCV_I,1) = 0

         CASE('N'); PIC_BCMI_NORMDIR(BCV_I,2) = 1
         CASE('S')
            PIC_BCMI_NORMDIR(BCV_I,2) = -1
            PIC_BCMI_OFFSET (BCV_I,2) = 0

         CASE('T'); PIC_BCMI_NORMDIR(BCV_I,3) = 1
         CASE('B')
            PIC_BCMI_NORMDIR(BCV_I,3) = -1
            PIC_BCMI_OFFSET (BCV_I,3) =  0
         END SELECT

         if(dFlag) write(*,"(2/,'Setting PIC_MI:',I3)") BCV_I

! The number of mass phases at this inlet.  While a system may be
! polydisperse, the inlet could consist of a single mass phase
         PHASE_CNT = 0
! The mass phase indices of incoming particles at this inlet
         PHASE_LIST(:) = -1
! The max diameter of incoming particles at this inlet
         MAX_DIA = ZERO

! Determine if the inlet is mono or polydisperse
         DO M=1, DES_MMAX+MMAX
            IF(SOLIDS_MODEL(M) /= 'PIC') CYCLE
            IF(BC_ROP_s(BCV,M) == UNDEFINED) CYCLE
            IF(COMPARE(BC_ROP_s(BCV,M),ZERO)) CYCLE
            PHASE_CNT = PHASE_CNT + 1
            PHASE_LIST(PHASE_CNT) = M
            MAX_DIA = MAX(MAX_DIA,D_P0(M))
         ENDDO

      ENDDO

      CALL SET_PIC_BCMI_IJK

      CALL FINL_ERR_MSG


      RETURN
      END SUBROUTINE SET_BC_PIC_MI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_PIC                                              !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the provided information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_PIC_BCMI_IJK

! Modules
!---------------------------------------------------------------------//
      use bc, only: BC_PLANE
      use bc, only: BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t
      USE cutcell, only: CUT_CELL_AT
      USE discretelement, only : DES_MMAX
      use error_manager
      use functions
      use funits, only: DMP_LOG
      use mpi_utility
      use param, only: dim_m
      use pic_bc, only: PIC_BCMI, PIC_BCMI_MAP, PIC_BCMI_IJK
      use pic_bc, only: PIC_BCMI_IJKSTART, PIC_BCMI_IJKEND
      use pic_bc, only: pic_bcmi_cnp
      use pic_bc, only: pic_bcmi_rnp
      use pic_bc, only: PIC_BCMI_INCL_CUTCELL
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER, ALLOCATABLE :: LOC_PIC_BCMI_IJK(:)
      INTEGER :: BCV, BCV_I
      INTEGER :: LC
      INTEGER :: MAX_CELLS
      INTEGER :: BND1, BND2
      LOGICAL, parameter :: setDBG = .false.
      LOGICAL :: dFlag
      INTEGER :: I,J,K,IJK
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t
!......................................................................!


      CALL INIT_ERR_MSG("SET_PIC_BCMI_IJK")

      dFlag = (DMP_LOG .AND. setDBG)

      if(dFlag) write(*,"(2/,2x,'From: SET_PIC_BCMI_IJK')")

! Loop over all inflow BCs to get an approximate count of the number
! of fluid cells that are adjacent to them.
      MAX_CELLS = 0
      DO BCV_I = 1, PIC_BCMI
         BCV = PIC_BCMI_MAP(BCV_I)

! Set the search area a little bigger than the inlet area.
         if(dFlag) WRITE(*,"(/2x,'Adding cells for BCV: ',I3)") BCV
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S')
            BND1 = min(BC_I_e(BCV)+1,IMAX1) - max(BC_I_w(BCV)-1,IMIN1)
            BND2 = min(BC_K_t(BCV)+1,KMAX1) - max(BC_K_b(BCV)-1,KMIN1)

         CASE('E','W')
            BND1 = min(BC_J_n(BCV)+1,JMAX1) - max(BC_J_s(BCV)-1,JMIN1)
            BND2 = min(BC_K_t(BCV)+1,KMAX1) - max(BC_K_b(BCV)-1,KMIN1)

         CASE('T','B')
            BND1 = min(BC_I_e(BCV)+1,IMAX1) - max(BC_I_w(BCV)-1,IMIN1)
            BND2 = min(BC_J_n(BCV)+1,JMAX1) - max(BC_J_s(BCV)-1,JMIN1)
         END SELECT

         MAX_CELLS = MAX_CELLS + (BND1 + 1)*(BND2 + 1)
         if(dFlag) WRITE(*,"(4x,'Plane:   ',A)") BC_PLANE(BCV)
         if(dFlag) WRITE(*,"(4x,'Cells: ', I8)") (BND1+1)*(BND2+1)
      ENDDO

      if(dFlag) write(*,"(2x,'Max Cells: ',I8)") MAX_CELLS

! Allocate an array to hold the IJK values. This array should be
! more than enough to store all the IJKs.
      allocate( LOC_PIC_BCMI_IJK(MAX_CELLS) )


! Loop over the IJKs for each BC and store only the IJKs that you
! own as well as the start/end array positions for each BC.
      LC = 1
      DO BCV_I = 1, PIC_BCMI

         PIC_BCMI_IJKSTART(BCV_I) = LC
         BCV = PIC_BCMI_MAP(BCV_I)

         if(dFlag) write(*,"(/2x,'Searching for fluid cells:',I3)") BCV

         I_w = max(BC_I_w(BCV),IMIN1); I_e = min(BC_I_e(BCV),IMAX1)
         J_s = max(BC_J_s(BCV),JMIN1); J_n = min(BC_J_n(BCV),JMAX1)
         K_b = max(BC_K_b(BCV),KMIN1); K_t = min(BC_K_t(BCV),KMAX1)

! Depending on the flow plane, the 'common' index needs set to reference
! the fluid cell.
         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); J_s = BC_J_s(BCV)+1;   J_n = BC_J_n(BCV)+1
         CASE('S'); J_s = BC_J_s(BCV)-1;   J_n = BC_J_n(BCV)-1
         CASE('E'); I_w = BC_I_w(BCV)+1;   I_e = BC_I_e(BCV)+1
         CASE('W'); I_w = BC_I_w(BCV)-1;   I_e = BC_I_e(BCV)-1
         CASE('T'); K_b = BC_K_b(BCV)+1;   K_t = BC_K_t(BCV)+1
         CASE('B'); K_b = BC_K_b(BCV)-1;   K_t = BC_K_t(BCV)-1
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

            !do not include this cell if the IO BC has been set to
            !not include cutcells
            IF(CUT_CELL_AT(IJK).AND.(.NOT.PIC_BCMI_INCL_CUTCELL(BCV_I))) CYCLE
            LOC_PIC_BCMI_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO

         PIC_BCMI_IJKEND(BCV_I) = LC-1

         IF(dFLAG) write(*,1111) BCV, BCV_I,                           &
            PIC_BCMI_IJKSTART(BCV_I), PIC_BCMI_IJKEND(BCV_I)

      ENDDO

 1111 FORMAT(/2x,'PIC Mass Inflow:',/4x,'BC:',I4,3x,'MAP:',I4,         &
         /4x,'IJKSTART:',I6,/4x,'IJKEND:  ',I6)


! Allocate the global store arrary array. This changes across MPI ranks.
      IF(LC > 1) THEN
         allocate( PIC_BCMI_IJK(LC-1) )
         allocate(pic_bcmi_cnp(LC-1, DIM_M))
         allocate(pic_bcmi_rnp(LC-1, DIM_M))

         PIC_BCMI_IJK(1:LC-1) = LOC_PIC_BCMI_IJK(1:LC-1)

         pic_bcmi_cnp(1:LC-1,:) = 0.d0
         pic_bcmi_rnp(1:LC-1,:) = 0.d0
      ELSE
         allocate( PIC_BCMI_IJK(1) )
         allocate(pic_bcmi_cnp(1,1))
         allocate(pic_bcmi_rnp(1,1))

         PIC_BCMI_IJK(1) = LOC_PIC_BCMI_IJK(1)
      ENDIF

      deallocate(LOC_PIC_BCMI_IJK)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_PIC_BCMI_IJK
