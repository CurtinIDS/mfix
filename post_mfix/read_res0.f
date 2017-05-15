!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: READ_RES0                                               C
!  Purpose: read the initial restart records (namelist data)           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE READ_RES0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE fldvar
      USE funits
      USE geometry
      USE ic
      USE in_binary_512
      USE in_binary_512i
      USE is
      USE leqsol
      USE machine
      USE mpi_utility
      USE output
      USE param, only: dimension_3, dimension_i, dimension_j, dimension_k, dimension_m, dimension_n_g, dimension_scalar
      USE param1
      USE physprop
      USE run
      USE rxns
      USE scalars
      USE scales
      USE stiff_chem
      USE toleranc
      USE ur_facs

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counters
      INTEGER    LC, L , N, M
! Pointer to the next record
      INTEGER    NEXT_RECA
! file version id
      CHARACTER(LEN=512) :: VERSION
! version number
      REAL       VERSION_NUMBER
! Temporary arrays
      DOUBLE PRECISION IC_Tmp(DIMENSION_IC), BC_Tmp(DIMENSION_BC)
!
      INTEGER    DIM_IC , DIM_BC , DIM_C , DIM_IS
! Read an array dimension
      INTEGER :: DIM_tmp
!
      logical :: doingPost
      integer :: n_spx_res

! Place holder for deprecated keywords:
      LOGICAL :: CALL_ISAT

! declare global scratch arrays:
! declare integer Global SCRatch array
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IGTEMP1, iGTEMP2
! declare real*4 Global SCRatch array
      REAL, ALLOCATABLE, DIMENSION(:) :: rGTEMP
! declare real*8 Global SCRatch array
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dGTEMP
! declaring following arrays to pack scalar variables when
! BCASTing to reduce the number of BCAST calls:
! packing array for integers
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INTPACK
! packing array for doubles
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DBLPACK
!-----------------------------------------------

      doingPost = .false.

!  1) Check to ensure that this subroutine was updated.
!  2) Initialize missing constants from earlier versions.
!  3) Add new read statements at the end of the file.

! only PE_IO reads the restart file
    if (myPE == PE_IO ) then
      READ (UNIT_RES, REC=1) VERSION
      READ (VERSION(6:512), *) VERSION_NUMBER

      IF (VERSION_NUMBER > 1.8) THEN
         WRITE (*, *) ' Update Subroutine read_res0'
         CALL SLUMBER
!         STOP
         call exitMPI(myPE)  ! Abort all PEs, not only the current one
      ENDIF
    endif

! Initialize required constants missing from earlier versions
      P_REF = ZERO
      P_SCALE = ONE
      DIM_IC = 5
      DIM_BC = 5
      DIM_C = 5
      DIM_IS = 5
      C_E = 1.0
      C_F = 0.0
      PHI = 0.0
      PHI_W = 0.0

! only PE_IO reads the restart
    if (myPE == PE_IO ) then
      READ (UNIT_RES, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, &
         ID_MINUTE, ID_SECOND
      READ (UNIT_RES, REC=3) NEXT_RECA

      IF (VERSION == 'RES = 01.00') THEN
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DT, &
            XLENGTH, YLENGTH, ZLENGTH
      ELSEIF (VERSION=='RES = 01.01' .OR. VERSION=='RES = 01.02') THEN
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DT, XLENGTH, YLENGTH, ZLENGTH
      ELSEIF (VERSION == 'RES = 01.03') THEN
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH
      ELSEIF (VERSION == 'RES = 01.04') THEN
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DIM_C, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH
      ELSEIF (VERSION == 'RES = 01.05') THEN
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DIM_C, DIM_IS, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH
      ELSE
         READ (UNIT_RES, REC=4) IMIN1, JMIN1, KMIN1, IMAX, JMAX, KMAX, IMAX1, &
            JMAX1, KMAX1, IMAX2, JMAX2, KMAX2, IJMAX2, IJKMAX2, MMAX, DIM_IC, &
            DIM_BC, DIM_C, DIM_IS, DT, XMIN, XLENGTH, YLENGTH, ZLENGTH, C_E, &
            C_F, PHI, PHI_W
      ENDIF
    endif

    Allocate( INTPACK(30))  ! ALLOCate packing array
    Allocate( DBLPACK(20))  ! ALLOCate packing array
    INTPACK = 0
    DBLPACK = 0.d0

! Root PE (PE_IO) : Pack and broadcast ; others unpack
! WARNING: This implementation assumes VERSION > 1.05, filtering needed for
!          earlier versions
    if (myPE == PE_IO ) then
       call bcast(VERSION_NUMBER, PE_IO)  ! BCAST0r
       INTPACK(1) = ID_MONTH
       INTPACK(2) = ID_DAY
       INTPACK(3) = ID_YEAR
       INTPACK(4) = ID_HOUR
       INTPACK(5) = ID_MINUTE
       INTPACK(6) = ID_SECOND
       INTPACK(7) = IMIN1
       INTPACK(8) = JMIN1
       INTPACK(9) = KMIN1
       INTPACK(10) = IMAX
       INTPACK(11) = JMAX
       INTPACK(12) = KMAX
       INTPACK(13) = IMAX1
       INTPACK(14) = JMAX1
       INTPACK(15) = KMAX1
       INTPACK(16) = IMAX2
       INTPACK(17) = JMAX2
       INTPACK(18) = KMAX2
       INTPACK(19) = IJMAX2
       INTPACK(20) = IJKMAX2
       INTPACK(21) = MMAX
       INTPACK(22) = DIM_IC
       INTPACK(23) = DIM_BC
       INTPACK(24) = DIM_C
       INTPACK(25) = DIM_IS
! In spite of the overhead, using MPI_Pack/MPI_Unpack may be worth using here
       call bcast(INTPACK(1:25),PE_IO)  ! BCAST1i
       DBLPACK(1) = DT
       DBLPACK(2) = XMIN
       DBLPACK(3) = XLENGTH
       DBLPACK(4) = YLENGTH
       DBLPACK(5) = ZLENGTH
       DBLPACK(6) = C_E
       DBLPACK(7) = C_F
       DBLPACK(8) = PHI
       DBLPACK(9) = PHI_W
       call bcast(DBLPACK,PE_IO) ! BCAST1d
       call bcast(VERSION,PE_IO) ! BCAST0c
       call bcast(RUN_NAME,PE_IO) ! BCAST0c
    else
       call bcast(VERSION_NUMBER, PE_IO)  ! BCAST0r
       call bcast(INTPACK(1:25),PE_IO)    ! BCAST1i (receive)
       ID_MONTH = INTPACK(1)
       ID_DAY = INTPACK(2)
       ID_YEAR = INTPACK(3)
       ID_HOUR = INTPACK(4)
       ID_MINUTE = INTPACK(5)
       ID_SECOND = INTPACK(6)
       IMIN1 = INTPACK(7)
       JMIN1 = INTPACK(8)
       KMIN1 = INTPACK(9)
       IMAX = INTPACK(10)
       JMAX = INTPACK(11)
       KMAX = INTPACK(12)
       IMAX1 = INTPACK(13)
       JMAX1 = INTPACK(14)
       KMAX1 = INTPACK(15)
       IMAX2 = INTPACK(16)
       JMAX2 = INTPACK(17)
       KMAX2 = INTPACK(18)
       IJMAX2 = INTPACK(19)
       IJKMAX2 = INTPACK(20)
       MMAX = INTPACK(21)
       DIM_IC = INTPACK(22)
       DIM_BC = INTPACK(23)
       DIM_C = INTPACK(24)
       DIM_IS = INTPACK(25)
       call bcast(DBLPACK,PE_IO)  ! BCAST1d (recv)
       DT = DBLPACK(1)
       XMIN = DBLPACK(2)
       XLENGTH = DBLPACK(3)
       YLENGTH = DBLPACK(4)
       ZLENGTH = DBLPACK(5)
       C_E = DBLPACK(6)
       C_F = DBLPACK(7)
       PHI = DBLPACK(8)
       PHI_W = DBLPACK(9)
       call bcast(VERSION,PE_IO)  ! BCAST0c (recv)
       call bcast(RUN_NAME,PE_IO) ! BCAST0c (recv)
    endif

    DeAllocate (INTPACK)



! CHECK DIMENSIONS

      IF (.NOT. ( (DIM_IC <= DIMENSION_IC)  .AND. &
                  (DIM_BC <= DIMENSION_BC)  .AND. &
                  (DIM_C <= DIMENSION_C)  .AND. &
                  (DIM_IS <= DIMENSION_IS) )) GOTO 900

      IF (MMAX + 1 > 0) THEN
         NMAX(:MMAX) = 1
      ENDIF

      NEXT_RECA = 5


      IF (VERSION_NUMBER >= 1.04) THEN

         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, C, DIM_C, NEXT_RECA)
! work around for -O3 compiler bug
             NEXT_RECA = 1 + NEXT_RECA
             NEXT_RECA = NEXT_RECA - 1
         endif
         call bcast(C,PE_IO)   ! BCAST1d user defined constants,C

         if (myPE == PE_IO) then
            DO LC = 1, DIM_C
               READ (UNIT_RES, REC=NEXT_RECA) C_NAME(LC)
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
         endif
         call bcast(C_NAME,PE_IO)  ! BCAST1c user defined constant names

         IF (myPE == PE_IO) THEN
            IF (VERSION_NUMBER < 1.12) THEN
               CALL IN_BIN_512I (UNIT_RES, NMAX, MMAX + 1, NEXT_RECA)
            ELSE
               READ (UNIT_RES, REC=NEXT_RECA) (NMAX(L),L=0, MMAX)
               NEXT_RECA = NEXT_RECA + 1
            ENDIF
         ENDIF
         call bcast(NMAX,PE_IO)   !BCAST1i total # of gas OR solid species
      ENDIF

! The following occurs when mfix.dat has not been read and the
! dimensions are from the .RES file.  Note that check on NMAX for
! solids is not done.
      IF ( .NOT. ( (IMAX2 <= DIMENSION_I)    .AND. &
                   (JMAX2 <= DIMENSION_J)    .AND. &
                   (KMAX2 <= DIMENSION_K)    .AND. &
                   (IJKMAX2 <= DIMENSION_3)  .AND. &
                   (MMAX <= DIMENSION_M)     .AND. &
                   (NMAX(0) <= DIMENSION_N_G) ) ) then

         IF(IMAX2 == 1)NO_I=.TRUE.
         IF(JMAX2 == 1)NO_J=.TRUE.
         IF(KMAX2 == 1)NO_K=.TRUE.

! This is a "fix" for post_mfix ... assumes only 1 processor
         write (*,*) ' '
         write (*,*) ' ********************************************************'
         write (*,*) ' '
         write (*,*) ' read_res0 : code valid for running on 1 processor only'
         write (*,*) ' '
         write (*,*) ' ********************************************************'
         write (*,*) ' '
         IMAX3 = imax2
         JMAX3 = jmax2
         KMAX3 = kmax2
         kend3 = kmax2
         kstart3 = 1
         jend3 = jmax2
         jstart3 = 1
         iend3 = imax2
         istart3 = 1
         ijkmax3 = imax3*jmax3*kmax3
         allocate (ijksize3_all(0:1))
         ijksize3_all(:) = ijkmax3

         nScalar = 0                ! since NScalar
         nRR     = 0                ! and NRR not read
         doingPost = .true.         ! until later

         call set_parameters
         call allocate_arrays_geometry
         call allocate_arrays_increments
         call allocate_arrays       ! do for mfix/post_mfix
         deallocate(ijksize3_all)   ! post_mfix "fix"

      ENDIF

! this read statement needs to be consistent with the write_res from
! mfix in terms of the index start. however, care must then be used
! when calc_distance is called to adjust how the array is read
      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, DX(1), IMAX2, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, DY(1), JMAX2, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, DZ(1), KMAX2, NEXT_RECA)
      endif


! the RES file for version <= 1.4 write out : dx(1) to dx(imax2)
! but reads starting at dx(0). Need to shift the arrays for
! proper pos-processing etc.
      if (version_number < 1.41) then
         do L = imax2,1,-1
            dx(L) = dx(L-1)
         end do
         do L = jmax2,1,-1
            dy(L) = dy(L-1)
         end do
         do L = kmax2,1,-1
            dz(L) = dz(L-1)
         end do
      endif
      call bcast(dx, PE_IO)   !//PAR_I/O BCAST1d
      call bcast(dy, PE_IO)   !//PAR_I/O BCAST1d
      call bcast(dz, PE_IO)   !//PAR_I/O BCAST1d


      IF (myPE == PE_IO) THEN
         READ (UNIT_RES, REC=NEXT_RECA) RUN_NAME, &
               DESCRIPTION, UNITS, RUN_TYPE, &
               COORDINATES
         NEXT_RECA = NEXT_RECA + 1
      ENDIF
      call bcast(RUN_NAME,PE_IO)    ! BCAST0c
      call bcast(DESCRIPTION,PE_IO) ! BCAST0c
      call bcast(UNITS,PE_IO)       ! BCAST0c
      call bcast(RUN_TYPE,PE_IO)    ! BCAST0c
      call bcast(COORDINATES,PE_IO) ! BCAST0c

      IF (myPE == PE_IO) THEN
         IF (VERSION=='RES = 01.00' .OR. VERSION=='RES = 01.01') THEN
            READ (UNIT_RES, REC=NEXT_RECA) (D_P0(L),L=1,&
               MMAX), (RO_S0(L),L=1,MMAX), EP_STAR, &
               MU_G0, MW_AVG
         ELSEIF (VERSION == 'RES = 01.02') THEN
            READ (UNIT_RES, REC=NEXT_RECA) (D_P0(L),L=1,&
               MMAX), (RO_S0(L),L=1,MMAX), EP_STAR, &
               RO_G0, MU_G0, MW_AVG
         ELSEIF (VERSION == 'RES = 01.03') THEN
            READ (UNIT_RES, REC=NEXT_RECA) (D_P0(L),L=1,&
                MMAX), (RO_S0(L),L=1,MMAX), EP_STAR, &
                RO_G0, MU_G0, MW_AVG
         ELSEIF (VERSION_NUMBER >= 1.04) THEN
            READ (UNIT_RES, REC=NEXT_RECA) (D_P0(L),L=1,&
                MMAX), (RO_S0(L),L=1,MMAX), EP_STAR, &
                RO_G0, MU_G0, MW_AVG
         ENDIF
         NEXT_RECA = NEXT_RECA + 1
      ENDIF
      call bcast(D_P0, PE_IO)      ! BCAST1d
      call bcast(RO_S0, PE_IO)     ! BCAST1d
      call bcast(EP_STAR, PE_IO)   ! BCAST0d
      call bcast(RO_G0, PE_IO)     ! BCAST0d
      call bcast(MU_G0, PE_IO)     ! BCAST0d
      call bcast(MW_AVG, PE_IO)    ! BCAST0d

      IF (VERSION_NUMBER >= 1.04) THEN
         IF (myPE == PE_IO) THEN
            CALL IN_BIN_512 (UNIT_RES, MW_G, NMAX(0), NEXT_RECA)
            DO LC = 1, MMAX
               READ (UNIT_RES, REC=NEXT_RECA) (MW_S(LC,N),&
                  N=1,NMAX(LC))
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
         ENDIF
         call bcast(MW_G, PE_IO)     ! BCAST1d
         call bcast(MW_S, PE_IO)     ! BCAST2d
      ENDIF

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, IC_X_W, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_X_E, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_Y_S, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_Y_N, DIM_IC, NEXT_RECA)
      endif
      call bcast(IC_X_W, PE_IO)     ! BCAST1d
      call bcast(IC_X_E, PE_IO)     ! BCAST1d
      call bcast(IC_Y_S, PE_IO)     ! BCAST1d
      call bcast(IC_Y_N, PE_IO)     ! BCAST1d

      if (myPE == PE_IO)then
         CALL IN_BIN_512 (UNIT_RES, IC_Z_B, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_Z_T, DIM_IC, NEXT_RECA)
      endif
      call bcast(IC_Z_B, PE_IO)     ! BCAST1d
      call bcast(IC_Z_T, PE_IO)     ! BCAST1d

      if (myPE == PE_IO)then
         CALL IN_BIN_512I (UNIT_RES, IC_I_W, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512I (UNIT_RES, IC_I_E, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512I (UNIT_RES, IC_J_S, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512I (UNIT_RES, IC_J_N, DIM_IC, NEXT_RECA)
      endif
      call bcast(IC_I_W, PE_IO)     ! BCAST1i
      call bcast(IC_I_E, PE_IO)     ! BCAST1i
      call bcast(IC_J_S, PE_IO)     ! BCAST1i
      call bcast(IC_J_N, PE_IO)     ! BCAST1i

      if (myPE == PE_IO)then
         CALL IN_BIN_512I (UNIT_RES, IC_K_B, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512I (UNIT_RES, IC_K_T, DIM_IC, NEXT_RECA)
      endif
      call bcast(IC_K_B, PE_IO)     ! BCAST1i
      call bcast(IC_K_T, PE_IO)     ! BCAST1i

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, IC_EP_G, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_P_G, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_T_G, DIM_IC, NEXT_RECA)
      endif
      call bcast(IC_EP_G, PE_IO)    ! BCAST1d
      call bcast(IC_P_G, PE_IO)     ! BCAST1d
      call bcast(IC_T_G, PE_IO)     ! BCAST1d

      IF (VERSION_NUMBER < 1.15) THEN
         IF (myPE == PE_IO) THEN
            CALL IN_BIN_512 (UNIT_RES, IC_T_S(1,1), DIM_IC, NEXT_RECA)
            IF (MMAX >= 2) THEN
               CALL IN_BIN_512 (UNIT_RES, IC_T_S(1,2), DIM_IC, NEXT_RECA)
            ELSE
               CALL IN_BIN_512 (UNIT_RES, IC_TMP, DIM_IC, NEXT_RECA)
            ENDIF
         ENDIF
            call bcast(IC_T_S, PE_IO)    ! BCAST2d
      ENDIF

      IF (VERSION_NUMBER >= 1.04) THEN
         IF (myPE == PE_IO) THEN
            DO N = 1, NMAX(0)
               CALL IN_BIN_512 (UNIT_RES, IC_X_G(1,N), DIM_IC, &
                  NEXT_RECA)
            ENDDO
         ENDIF
         call bcast(IC_X_G, PE_IO)    ! BCAST2d
      ENDIF

      if (myPE == PE_IO)then
         CALL IN_BIN_512 (UNIT_RES, IC_U_G, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_V_G, DIM_IC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, IC_W_G, DIM_IC, NEXT_RECA)
      ENDIF
      call bcast(IC_U_G, PE_IO)    ! BCAST1d
      call bcast(IC_V_G, PE_IO)    ! BCAST1d
      call bcast(IC_W_G, PE_IO)    ! BCAST1d

      IF (myPE == PE_IO) THEN
         DO LC = 1, MMAX
            CALL IN_BIN_512 (UNIT_RES, IC_ROP_S(1,LC), &
                 DIM_IC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IC_U_S(1,LC), &
                 DIM_IC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IC_V_S(1,LC), &
                 DIM_IC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IC_W_S(1,LC), &
                 DIM_IC, NEXT_RECA)
            IF (VERSION_NUMBER >= 1.15) THEN
               CALL IN_BIN_512(UNIT_RES, IC_T_S(1,LC), DIM_IC, &
                 NEXT_RECA)
            ENDIF

            IF (VERSION_NUMBER >= 1.04) THEN
               DO N = 1, NMAX(LC)
                 CALL IN_BIN_512 (UNIT_RES, IC_X_S(1,LC,N), &
                       DIM_IC, NEXT_RECA)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      call bcast(IC_ROP_S, PE_IO)  ! BCAST2d
      call bcast(IC_U_S, PE_IO)    ! BCAST2d
      call bcast(IC_V_S, PE_IO)    ! BCAST2d
      call bcast(IC_W_S, PE_IO)    ! BCAST2d

      if (VERSION_NUMBER >= 1.15) call bcast(IC_T_S, PE_IO)  ! BCAST2d
      if (VERSION_NUMBER >= 1.04) call bcast(IC_X_S, PE_IO)  ! BCAST3d

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, BC_X_W, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_X_E, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_Y_S, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_Y_N, DIM_BC, NEXT_RECA)
      endif
      call bcast(BC_X_W, PE_IO)    ! BCAST1d
      call bcast(BC_X_E, PE_IO)    ! BCAST1d
      call bcast(BC_Y_S, PE_IO)    ! BCAST1d
      call bcast(BC_Y_N, PE_IO)    ! BCAST1d

      if (myPE == PE_IO) then
         CALL IN_BIN_512 (UNIT_RES, BC_Z_B, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_Z_T, DIM_BC, NEXT_RECA)
      endif
      call bcast(BC_Z_B, PE_IO)    ! BCAST1d
      call bcast(BC_Z_T, PE_IO)    ! BCAST1d

      IF (myPE == PE_IO) THEN
         CALL IN_BIN_512I (UNIT_RES, BC_I_W, DIM_BC, NEXT_RECA)
          CALL IN_BIN_512I (UNIT_RES, BC_I_E, DIM_BC, NEXT_RECA)
          CALL IN_BIN_512I (UNIT_RES, BC_J_S, DIM_BC, NEXT_RECA)
          CALL IN_BIN_512I (UNIT_RES, BC_J_N, DIM_BC, NEXT_RECA)
      ENDIF
      call bcast(BC_I_W, PE_IO)    ! BCAST1i
      call bcast(BC_I_E, PE_IO)    ! BCAST1i
      call bcast(BC_J_S, PE_IO)    ! BCAST1i
      call bcast(BC_J_N, PE_IO)    ! BCAST1i

      IF (myPE == PE_IO) THEN
           CALL IN_BIN_512I (UNIT_RES, BC_K_B, DIM_BC, NEXT_RECA)
           CALL IN_BIN_512I (UNIT_RES, BC_K_T, DIM_BC, NEXT_RECA)
      ENDIF
      call bcast(BC_K_B, PE_IO)    ! BCAST1i
      call bcast(BC_K_T, PE_IO)    ! BCAST1i

      IF (myPE == PE_IO) THEN
         CALL IN_BIN_512 (UNIT_RES, BC_EP_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_P_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_T_G, DIM_BC, NEXT_RECA)
      ENDIF
      call bcast(BC_EP_G, PE_IO)   ! BCAST1d
      call bcast(BC_P_G, PE_IO)    ! BCAST1d
      call bcast(BC_T_G, PE_IO)    ! BCAST1d

      IF (VERSION_NUMBER < 1.15) THEN
         IF (myPE == PE_IO) THEN
            CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,1), DIM_BC, NEXT_RECA)
            IF (MMAX >= 2) THEN
               CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,2), DIM_BC, NEXT_RECA)
            ELSE
! dummy read, no need to broadcast
               CALL IN_BIN_512 (UNIT_RES, BC_TMP, DIM_BC, NEXT_RECA)
            ENDIF
         ENDIF
         call bcast(BC_T_S, PE_IO)   ! BCAST2d
      ENDIF

      IF (VERSION_NUMBER >= 1.04) THEN
         IF (myPE == PE_IO) THEN
            DO N = 1, NMAX(0)
               CALL IN_BIN_512 (UNIT_RES, BC_X_G(1,N), &
                    DIM_BC, NEXT_RECA)
            ENDDO
         ENDIF
         call bcast(BC_X_G, PE_IO)   ! BCAST2d
      ENDIF

      IF (myPE == PE_IO) THEN
         CALL IN_BIN_512 (UNIT_RES, BC_U_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_V_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_W_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_RO_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_ROP_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_VOLFLOW_G, DIM_BC, NEXT_RECA)
         CALL IN_BIN_512 (UNIT_RES, BC_MASSFLOW_G, DIM_BC, NEXT_RECA)
      ENDIF
      call bcast(BC_U_G, PE_IO)       ! BCAST1d
      call bcast(BC_V_G, PE_IO)       ! BCAST1d
      call bcast(BC_W_G, PE_IO)       ! BCAST1d
      call bcast(BC_RO_G, PE_IO)      ! BCAST1d
      call bcast(BC_ROP_G, PE_IO)     ! BCAST1d
      call bcast(BC_VOLFLOW_G, PE_IO) ! BCAST1d
      call bcast(BC_MASSFLOW_G, PE_IO)! BCAST1d

      IF (myPE == PE_IO) THEN
         DO LC = 1, MMAX
            CALL IN_BIN_512 (UNIT_RES, BC_ROP_S(1,LC), DIM_BC, &
               NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_U_S(1,LC), DIM_BC, &
               NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_V_S(1,LC), DIM_BC, &
               NEXT_RECA)

! Note : previous versions did not write out BC_W_s, X_S
            IF (VERSION_NUMBER >= 1.04) THEN
               CALL IN_BIN_512 (UNIT_RES, BC_W_S(1,LC), DIM_BC, &
                  NEXT_RECA)
! Note : previous versions did not write out BC_T_s
               IF (VERSION_NUMBER >= 1.15) &
                  CALL IN_BIN_512 (UNIT_RES, BC_T_S(1,LC), DIM_BC, &
                     NEXT_RECA)
               DO N = 1, NMAX(LC)
                  CALL IN_BIN_512 (UNIT_RES, BC_X_S(1,LC,N), &
                     DIM_BC, NEXT_RECA)
               ENDDO
            ENDIF
            CALL IN_BIN_512 (UNIT_RES, BC_VOLFLOW_S(1,LC), DIM_BC,&
               NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_MASSFLOW_S(1,LC), DIM_BC,&
                NEXT_RECA)
         ENDDO
      ENDIF
      call bcast(BC_ROP_S, PE_IO) ! BCAST2d
      call bcast(BC_U_S, PE_IO)   ! BCAST2d
      call bcast(BC_V_S, PE_IO)   ! BCAST2d
      if (VERSION_NUMBER >= 1.04) then
         call bcast(BC_W_S, PE_IO)   ! BCAST2d
         if (VERSION_NUMBER >= 1.15)  call bcast(BC_T_S, PE_IO)   ! BCAST2d
         call bcast(BC_X_S, PE_IO)   ! BCAST2d
      endif
      call bcast(BC_VOLFLOW_S, PE_IO)   ! BCAST2d
      call bcast(BC_MASSFLOW_S, PE_IO)  ! BCAST2d

      IF (myPE == PE_IO) THEN
         IF (VERSION == 'RES = 01.00') THEN
            L = 10
         ELSE
            L = DIM_BC
         ENDIF
         DO LC = 1, L
            READ (UNIT_RES, REC=NEXT_RECA) BC_TYPE(LC)
            NEXT_RECA = NEXT_RECA + 1
         ENDDO
      ENDIF
      call bcast(BC_TYPE, PE_IO)   ! BCAST1c

      if (myPE == PE_IO) then
         Allocate(IGTEMP1(IJKMAX2))   ! ALLOCate INT Global scratch
         Allocate(iGTEMP2(IJKMAX3))   ! ALLOCate INT Global scratch
         CALL IN_BIN_512I (UNIT_RES, IGTEMP1, IJKMAX2, NEXT_RECA)
         call convert_from_io_i(IGTEMP1,iGTEMP2,ijkmax2)
      else
         Allocate(IGTEMP1(1))   ! ALLOCate INT Global scratch
         Allocate(iGTEMP2(1))   ! ALLOCate INT Global scratch
      endif
      flag = iGTEMP2
      DeAllocate (IGTEMP1, iGTEMP2)


! ------------------------------------------------------------------------
      IF (VERSION_NUMBER >= 1.04) THEN
         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, IS_X_W, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_X_E, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_Y_S, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_Y_N, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_Z_B, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_Z_T, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512I (UNIT_RES, IS_I_W, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512I (UNIT_RES, IS_I_E, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512I (UNIT_RES, IS_J_S, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512I (UNIT_RES, IS_J_N, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512I (UNIT_RES, IS_K_B, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512I (UNIT_RES, IS_K_T, DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_PC(1,1), DIM_IS, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IS_PC(1,2), DIM_IS, NEXT_RECA)
         endif
         call bcast(IS_X_W, PE_IO)   ! BCAST1d
         call bcast(IS_X_E, PE_IO)   ! BCAST1d
         call bcast(IS_Y_S, PE_IO)   ! BCAST1d
         call bcast(IS_Y_N, PE_IO)   ! BCAST1d
         call bcast(IS_Z_B, PE_IO)   ! BCAST1d
         call bcast(IS_Z_T, PE_IO)   ! BCAST1d
         call bcast(IS_I_W, PE_IO)   ! BCAST1i
         call bcast(IS_I_E, PE_IO)   ! BCAST1i
         call bcast(IS_J_S, PE_IO)   ! BCAST1i
         call bcast(IS_J_N, PE_IO)   ! BCAST1i
         call bcast(IS_K_B, PE_IO)   ! BCAST1i
         call bcast(IS_K_T, PE_IO)   ! BCAST1i
         call bcast(IS_PC, PE_IO)    ! BCAST1i

         IF (VERSION_NUMBER >= 1.07) THEN
            if (myPE == PE_IO) then
               DO LC = 1, MMAX
                 CALL IN_BIN_512 (UNIT_RES, IS_VEL_S(1,LC), &
                       DIM_IS, NEXT_RECA)
               ENDDO
            endif
            call bcast(IS_VEL_S, PE_IO)   ! BCAST2d
         ENDIF

         if (myPE == PE_IO) then
            DO LC = 1, DIM_IS
               READ (UNIT_RES, REC=NEXT_RECA) IS_TYPE(LC)
               NEXT_RECA = NEXT_RECA + 1
             ENDDO
         endif
         call bcast(IS_TYPE, PE_IO)   ! BCAST1c
      ENDIF

! ------------------------------------------------------------------------
! Additions from new versions of .RES file
      IF (VERSION_NUMBER >= 1.08) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) CYCLIC_X, &
               CYCLIC_Y, CYCLIC_Z, CYCLIC_X_PD, &
               CYCLIC_Y_PD, CYCLIC_Z_PD, DELP_X, DELP_Y, &
               DELP_Z, U_G0, U_S0, V_G0, V_S0, W_G0, W_S0
            NEXT_RECA = NEXT_RECA + 1
         ENDIF
         call bcast(CYCLIC_X,PE_IO)     ! BCAST0l
         call bcast(CYCLIC_Y,PE_IO)     ! BCAST0l
         call bcast(CYCLIC_Z,PE_IO)     ! BCAST0l
         call bcast(CYCLIC_X_PD,PE_IO)  ! BCAST0l
         call bcast(CYCLIC_Y_PD,PE_IO)  ! BCAST0l
         call bcast(CYCLIC_Z_PD,PE_IO)  ! BCAST0l
         call bcast(DELP_X,PE_IO)       ! BCAST1d
         call bcast(DELP_Y,PE_IO)       ! BCAST1d
         call bcast(DELP_Z,PE_IO)       ! BCAST1d
         call bcast(U_G0,PE_IO)         ! BCAST1d
         call bcast(U_S0,PE_IO)         ! BCAST1d
         call bcast(V_G0,PE_IO)         ! BCAST1d
         call bcast(V_S0,PE_IO)         ! BCAST1d
         call bcast(W_G0,PE_IO)         ! BCAST1d
         call bcast(W_S0,PE_IO)         ! BCAST1d
      ENDIF

! ------------------------------------------------------------------------
      IF (VERSION_NUMBER >= 1.09) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) TIME, TSTOP, &
                ENERGY_EQ, RES_DT, OUT_DT, NLOG, &
                L_SCALE0, NO_I, NO_J, NO_K, CALL_USR
             NEXT_RECA = NEXT_RECA + 1
         ENDIF
         call bcast(TIME,PE_IO)      ! BCAST0d
         call bcast(TSTOP,PE_IO)     ! BCAST0d
         call bcast(ENERGY_EQ,PE_IO) ! BCAST0l
         call bcast(RES_DT,PE_IO)    ! BCAST0d
         call bcast(OUT_DT,PE_IO)    ! BCAST0d
         call bcast(NLOG,PE_IO)      ! BCAST0i
         call bcast(L_SCALE0,PE_IO)  ! BCAST0d
         call bcast(NO_I,PE_IO)      ! BCAST0l
         call bcast(NO_J,PE_IO)      ! BCAST0l
         call bcast(NO_K,PE_IO)      ! BCAST0l
         call bcast(CALL_USR,PE_IO)  ! BCAST0l

         IF (myPE == PE_IO) THEN
            if (version_number >= 1.50) then
               read (unit_res,rec=next_reca) n_spx_res
               next_reca = next_reca + 1
               if (n_spx_res > n_spx) then
                  write (*,*) ' n_spx too small '
                  write (*,*) ' n_spx = ' , n_spx
                  write (*,*) ' n_spx must equal ' , n_spx_res
                  call exitMPI(myPE)  ! Abort all PEs, not only the current one
               endif
            else
               n_spx_res = 9
            endif

            DO LC = 1, N_SPX_RES
               READ (UNIT_RES, REC=NEXT_RECA) SPX_DT(LC)
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
            DO LC = 0, MMAX
               READ (UNIT_RES, REC=NEXT_RECA) SPECIES_EQ(LC)
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
         ENDIF
         call bcast(SPX_DT,PE_IO)     ! BCAST1d
         call bcast(SPECIES_EQ,PE_IO) ! BCAST1l (recv)

         IF (myPE == PE_IO) THEN
            CALL IN_BIN_512 (UNIT_RES, USR_DT, DIMENSION_USR, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, USR_X_W, DIMENSION_USR, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, USR_X_E, DIMENSION_USR, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, USR_Y_S, DIMENSION_USR, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, USR_Y_N, DIMENSION_USR, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, USR_Z_B, DIMENSION_USR, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, USR_Z_T, DIMENSION_USR, NEXT_RECA)
         ENDIF
         call bcast(USR_DT,PE_IO)  ! BCAST1d
         call bcast(USR_X_W,PE_IO) ! BCAST1d
         call bcast(USR_X_E,PE_IO) ! BCAST1d
         call bcast(USR_Y_S,PE_IO) ! BCAST1d
         call bcast(USR_Y_N,PE_IO) ! BCAST1d
         call bcast(USR_Z_B,PE_IO) ! BCAST1d
         call bcast(USR_Z_T,PE_IO) ! BCAST1d

         IF (myPE == PE_IO) THEN
            DO LC = 1, DIMENSION_USR
               READ (UNIT_RES, REC=NEXT_RECA) USR_FORMAT(LC), &
                  USR_EXT(LC), USR_TYPE(LC), USR_VAR(LC)
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
         ENDIF
         call bcast(USR_FORMAT,PE_IO) ! BCAST1c
         call bcast(USR_EXT,PE_IO)    ! BCAST1c
         call bcast(USR_TYPE,PE_IO)   ! BCAST1c
         call bcast(USR_VAR,PE_IO)    ! BCAST1c

         if (myPE == PE_IO) then
            CALL IN_BIN_512 (UNIT_RES, IC_P_STAR, DIM_IC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, IC_L_SCALE, DIM_IC, NEXT_RECA)
         endif
         call bcast(IC_P_STAR,PE_IO)  ! BCAST1d
         call bcast(IC_L_SCALE,PE_IO) ! BCAST1d

         IF (myPE == PE_IO) THEN
            DO LC = 1, DIM_IC
               READ (UNIT_RES, REC=NEXT_RECA) IC_TYPE(LC)
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
         ENDIF
         call bcast(IC_TYPE,PE_IO) ! BCAST1c

         IF (myPE == PE_IO) THEN
            CALL IN_BIN_512 (UNIT_RES, BC_DT_0, DIM_BC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_JET_G0, DIM_BC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_DT_H, DIM_BC , NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_JET_GH, DIM_BC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_DT_L, DIM_BC , NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_JET_GL, DIM_BC, NEXT_RECA)
         ENDIF
         call bcast(BC_DT_0,PE_IO) !//PAR_I/O BCAST1d
         call bcast(BC_JET_G0,PE_IO) !//PAR_I/O BCAST1d
         call bcast(BC_DT_H,PE_IO) !//PAR_I/O BCAST1d
         call bcast(BC_JET_GH,PE_IO) !//PAR_I/O BCAST1d
         call bcast(BC_DT_L,PE_IO) !//PAR_I/O BCAST1d
         call bcast(BC_JET_GL,PE_IO) !//PAR_I/O BCAST1d
      ENDIF   ! endif (version_number >=1.09)

! ------------------------------------------------------------------------
      IF (VERSION_NUMBER >= 1.10) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) MU_GMAX
            NEXT_RECA = NEXT_RECA + 1
         ENDIF
         call bcast(MU_GMAX,PE_IO) ! BCAST0d
      ENDIF

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.11) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) V_EX, MODEL_B
           NEXT_RECA = NEXT_RECA + 1
        ENDIF
        call bcast(V_EX,PE_IO)    ! BCAST0d
        call bcast(MODEL_B,PE_IO) ! BCAST0l
     ENDIF

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.12) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) P_REF, &
              P_SCALE, UR_FAC, TOL_RESID, DT_MAX, &
              DT_MIN, DT_FAC, CLOSE_PACKED, GRAVITY, &
              MU_S0(1)
            NEXT_RECA = NEXT_RECA + 1
        ENDIF
        call bcast(P_REF,PE_IO)        ! BCAST0d
        call bcast(P_SCALE,PE_IO)      ! BCAST0d
        call bcast(UR_FAC,PE_IO)       ! BCAST0d
        call bcast(TOL_RESID,PE_IO)    ! BCAST0d
        call bcast(DT_MAX,PE_IO)       ! BCAST0d
        call bcast(DT_MIN,PE_IO)       ! BCAST0d
        call bcast(DT_FAC,PE_IO)       ! BCAST0d
        call bcast(CLOSE_PACKED,PE_IO) ! BCAST0l
        call bcast(GRAVITY,PE_IO)      ! BCAST0d
        call bcast(MU_S0(1),PE_IO)        ! BCAST0d

        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) LEQ_IT, LEQ_METHOD
           NEXT_RECA = NEXT_RECA + 1
        ENDIF
        call bcast(LEQ_IT,PE_IO)     ! BCAST1i
        call bcast(LEQ_METHOD,PE_IO) ! BCAST1i

        IF (myPE == PE_IO) THEN
           CALL IN_BIN_512 (UNIT_RES, BC_HW_G, DIM_BC, NEXT_RECA)
           CALL IN_BIN_512 (UNIT_RES, BC_UW_G, DIM_BC, NEXT_RECA)
           CALL IN_BIN_512 (UNIT_RES, BC_VW_G, DIM_BC, NEXT_RECA)
            CALL IN_BIN_512 (UNIT_RES, BC_WW_G, DIM_BC, NEXT_RECA)
        ENDIF
        call bcast(BC_HW_G,PE_IO) ! BCAST1d
        call bcast(BC_UW_G,PE_IO) ! BCAST1d
        call bcast(BC_VW_G,PE_IO) ! BCAST1d
        call bcast(BC_WW_G,PE_IO) ! BCAST1d

        IF (myPE == PE_IO) THEN
           DO LC = 1, MMAX
              CALL IN_BIN_512 (UNIT_RES, BC_HW_S(1,LC), DIM_BC, NEXT_RECA)
              CALL IN_BIN_512 (UNIT_RES, BC_UW_S(1,LC), DIM_BC, NEXT_RECA)
              CALL IN_BIN_512 (UNIT_RES, BC_VW_S(1,LC), DIM_BC, NEXT_RECA)
              CALL IN_BIN_512 (UNIT_RES, BC_WW_S(1,LC), DIM_BC, NEXT_RECA)
            ENDDO
        ENDIF
        call bcast(BC_HW_S,PE_IO) ! BCAST2d
        call bcast(BC_UW_S,PE_IO) ! BCAST2d
        call bcast(BC_VW_S,PE_IO) ! BCAST2d
        call bcast(BC_WW_S,PE_IO) ! BCAST2d
     ENDIF   ! endif (version_number >=1.12)

     LC = 0
     IF (MMAX + 1 > 0) THEN
        MOMENTUM_X_EQ(:MMAX) = .TRUE.
        MOMENTUM_Y_EQ(:MMAX) = .TRUE.
        MOMENTUM_Z_EQ(:MMAX) = .TRUE.
        LC = MMAX + 1
     ENDIF
     TOL_DIVERGE = 1.E+4

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.13) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) MOMENTUM_X_EQ, &
              MOMENTUM_Y_EQ, MOMENTUM_Z_EQ, TOL_DIVERGE, &
              DISCRETIZE, FULL_LOG
           NEXT_RECA = NEXT_RECA + 1
        ENDIF
        call bcast(MOMENTUM_X_EQ,PE_IO) ! BCAST1l
        call bcast(MOMENTUM_Y_EQ,PE_IO) ! BCAST1l
        call bcast(MOMENTUM_Z_EQ,PE_IO) ! BCAST1l
        call bcast(TOL_DIVERGE,PE_IO)   ! BCAST0d
        call bcast(DISCRETIZE,PE_IO)    ! BCAST1i
        call bcast(FULL_LOG,PE_IO)      ! BCAST0l
     ENDIF

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.14) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) DETECT_STALL
           NEXT_RECA = NEXT_RECA + 1
        ENDIF
        call bcast(DETECT_STALL,PE_IO) ! BCAST0l
     ENDIF

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.15) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) K_G0, K_S0(1), &
              C_PG0, C_PS0(1), TOL_RESID_T, TOL_RESID_X
           NEXT_RECA = NEXT_RECA + 1
        CALL IN_BIN_512 (UNIT_RES, IC_GAMA_RG, DIM_IC, NEXT_RECA)
        CALL IN_BIN_512 (UNIT_RES, IC_T_RG, DIM_IC, NEXT_RECA)
        ENDIF
        call bcast(K_G0,PE_IO)        ! BCAST0d
        call bcast(K_S0(1),PE_IO)        ! BCAST0d
        call bcast(C_PG0,PE_IO)       ! BCAST0d
        call bcast(C_PS0(1),PE_IO)    ! BCAST0d
        call bcast(TOL_RESID_T,PE_IO) ! BCAST0d
        call bcast(TOL_RESID_X,PE_IO) ! BCAST0d
        call bcast(IC_GAMA_RG,PE_IO)  ! BCAST1d
        call bcast(IC_T_RG,PE_IO)     ! BCAST1d

        IF (myPE == PE_IO) THEN
           DO LC = 1, MMAX
              CALL IN_BIN_512 (UNIT_RES, IC_GAMA_RS(1,LC), DIM_IC, &
                 NEXT_RECA)
              CALL IN_BIN_512 (UNIT_RES, IC_T_RS(1,LC), DIM_IC, &
                 NEXT_RECA)
           ENDDO
        ENDIF
        call bcast(IC_GAMA_RS,PE_IO) ! BCAST2d
        call bcast(IC_T_RS,PE_IO)    ! BCAST2d
     ENDIF

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.2) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) NORM_G, NORM_S
           NEXT_RECA = NEXT_RECA + 1
        ENDIF
        call bcast(NORM_G,PE_IO) ! BCAST0d
        call bcast(NORM_S,PE_IO) ! BCAST0d
     ENDIF

! ------------------------------------------------------------------------
     IF (VERSION_NUMBER >= 1.3) THEN
        IF (myPE == PE_IO) THEN
           READ (UNIT_RES, REC=NEXT_RECA) NScalar, TOL_RESID_Scalar, DIM_tmp
           NEXT_RECA = NEXT_RECA + 1
           CALL IN_BIN_512I (UNIT_RES, Phase4Scalar, DIM_tmp, NEXT_RECA)

! post mfix fix ...
           if (doingPost .and. nscalar.gt.0) then
              DIMENSION_Scalar = NScalar
              Allocate(  Scalar (DIMENSION_3,  DIMENSION_Scalar) )
              Allocate(  Scalaro (DIMENSION_3, DIMENSION_Scalar) )
           endif

        ENDIF
        call bcast(NScalar,PE_IO)          ! BCAST0d
        call bcast(TOL_RESID_Scalar,PE_IO) ! BCAST0d
        call bcast(Phase4Scalar,PE_IO)     ! BCAST0d
     ELSE
        NScalar = 0
     ENDIF

! ------------------------------------------------------------------------
! Version 1.4 -- read radiation variables in read_res1
      IF (VERSION_NUMBER >= 1.499) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) nRR
            NEXT_RECA = NEXT_RECA + 1
            if (doingPost .and. nRR.gt.0) then
               Allocate( ReactionRates(DIMENSION_3, nRR) )
            endif
         ENDIF
         call bcast(nRR,PE_IO) ! BCAST0d
      ELSE
         nRR = 0
      ENDIF

! ------------------------------------------------------------------------
! Version 1.6 -- read K_Epsilon and dqmom
      IF (VERSION_NUMBER >= 1.599) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) K_Epsilon, Call_DQMOM
            NEXT_RECA = NEXT_RECA + 1
            if (doingPost .and. K_epsilon) then
               Allocate( K_Turb_G(DIMENSION_3) )
               Allocate( E_Turb_G(DIMENSION_3) )
            end if
         ENDIF
         call bcast(K_Epsilon,PE_IO) !//PAR_I/O BCAST0d
         call bcast(Call_DQMOM,PE_IO) !//PAR_I/O BCAST0d
      ELSE
         K_Epsilon = .FALSE.
         Call_DQMOM =.FALSE.
      ENDIF

! ------------------------------------------------------------------------
! Version 1.7 -- Stiff Chemistry
      IF (VERSION_NUMBER >= 1.699) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) STIFF_CHEMISTRY, CALL_ISAT
            NEXT_RECA = NEXT_RECA + 1
         ENDIF
         call bcast(STIFF_CHEMISTRY,PE_IO)
         !call bcast(CALL_ISAT,PE_IO)
      ELSE
         STIFF_CHEMISTRY = .FALSE.
      ENDIF

! ------------------------------------------------------------------------
! Version 1.8 -- Variable solid density and  each solids species
      IF (VERSION_NUMBER >= 1.799) THEN
         IF (myPE == PE_IO) THEN
            READ (UNIT_RES, REC=NEXT_RECA) (SOLVE_ROs(LC),LC=1,MMAX)
            NEXT_RECA = NEXT_RECA + 1
            DO LC = 1, MMAX
               READ (UNIT_RES, REC=NEXT_RECA) INERT_SPECIES(LC), &
                  RO_s0(LC), (X_s0(LC,N),N=1,NMAX(LC))
               NEXT_RECA = NEXT_RECA + 1
            ENDDO
         ENDIF
         call bcast(SOLVE_ROs,PE_IO)
         call bcast(INERT_SPECIES,PE_IO)
         call bcast(RO_s0,PE_IO)
         call bcast(X_s0,PE_IO)

         ANY_SOLVE_ROs = ANY(SOLVE_ROs)
         IF (doingPost .AND. (.NOT.ANY_SOLVE_ROs)) THEN
            DO LC = 1, MMAX
               RO_S(:,LC) = RO_S0(LC)
            END DO
         ENDIF
      ELSE
         SOLVE_ROs = .FALSE.
         ANY_SOLVE_ROs = .FALSE.
         DO LC = 1, MMAX
            RO_S(:,LC) = RO_S0(LC)
         END DO
      ENDIF

! Add new read statements above this line.  Remember to update NEXT_RECA.
! Remember to update the version number check near beginning of this subroutine.
!------------------------------------------------------------------------------

      READ (UNIT_RES, REC=3) NEXT_RECA

! Since the value of UNDEFINED was changed ...
      IF (RO_G0 >= 1E30) RO_G0 = UNDEFINED
      IF (MU_G0 >= 1E30) MU_G0 = UNDEFINED
      IF (MW_AVG >= 1E30) MW_AVG = UNDEFINED
      IF (C_E >= 1E30) C_E = UNDEFINED

      RETURN

! HERE IF DIMENSION ERROR

  900 CONTINUE
      WRITE (*, *) ' '
      WRITE (*, *) ' **************************************'
      WRITE (*, "('(PE ',I6,'): From: READ_RES0')") myPE
      WRITE (*, *) ' DIMENSION ERROR ---'
      WRITE (*, *) ' '
      WRITE (*, *) ' DIMENSION_IC = ', DIMENSION_IC, ' DIM_IC       = ', DIM_IC
      WRITE (*, *) ' DIMENSION_BC = ', DIMENSION_BC, ' DIM_BC       = ', DIM_BC
      WRITE (*, *) ' DIMENSION_IS = ', DIMENSION_IS, ' DIM_IS       = ', DIM_IS
      WRITE (*, *) ' DIMENSION_C  = ', DIMENSION_C, ' DIM_C        = ', DIM_C
      WRITE (*, *) ' '

      call exitMPI(myPE)

      END SUBROUTINE READ_RES0

