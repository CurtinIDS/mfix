!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DBG_VTU_AND_VTP_FILES                            C
!  Purpose: Writes the cell and particle data in VTK format            C
!           for debug regions only.                                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 22-Jul-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DBG_VTU_AND_VTP_FILES

      use discretelement, only: DISCRETE_ELEMENT
      use vtp, only: write_vtp_file
      use vtk, only: DIMENSION_VTK

      IMPLICIT NONE
      INTEGER :: LC

      DO LC = 1, DIMENSION_VTK
         CALL WRITE_VTU_FILE(LC,1)
         IF(DISCRETE_ELEMENT) CALL WRITE_VTP_FILE(LC,1)
      ENDDO

      END SUBROUTINE WRITE_DBG_VTU_AND_VTP_FILES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VTU_FILE                                         C
!  Purpose: Writes the cell data grid in VTK format (Unstructured VTU) C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VTU_FILE(LCV,MODE)

      USE compar
      USE constant
      USE cutcell
      USE discretelement, Only :  DISCRETE_ELEMENT
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE mfix_pic
      USE mpi_utility
      USE output
      USE parallel
      USE parallel_mpi
      USE param
      USE param1
      USE pgcor
      USE pgcor
      USE physprop
      USE pscor
      USE quadric
      USE run
      USE rxns
      USE scalars
      USE sendrecv
      USE stl
      USE toleranc
      USE visc_s
      USE vtk
      USE vtp

      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,NN,R,IJK,LCV

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FACET_COUNT_DES, NEIGHBORING_FACET

      INTEGER :: SPECIES_COUNTER,LT

      CHARACTER (LEN=32) :: SUBM,SUBN,SUBR
      CHARACTER (LEN=64) :: VAR_NAME

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  DP_BC_ID, IJK_ARRAY

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)

      ! There is nothing to write if we are not in adefined vtk region
      VTK_REGION = LCV
      IF(.NOT.VTK_DEFINED(VTK_REGION)) RETURN
      IF(VTK_DATA(LCV)/='C') RETURN
      IF(MODE==0.AND.(VTK_DBG_FILE(LCV))) RETURN
      IF(MODE==1.AND.(.NOT.VTK_DBG_FILE(LCV))) RETURN

!     Location of U-momentum cells for original (uncut grid)
      IF (DO_I) THEN
        XG_E(1) = ZERO
        DO I = IMIN1, IMAX2
           XG_E(I) = XG_E(I-1) + DX(I)
        END DO
      ENDIF

!     Location of V-momentum cells for original (uncut grid)
      IF (DO_J) THEN
        YG_N(1) = ZERO
        DO J = JMIN1, JMAX2
           YG_N(J) = YG_N(J-1) + DY(J)
        END DO
      ENDIF

!     Location of W-momentum cells for original (uncut grid)
      IF (DO_K) THEN
        ZG_T(1) = ZERO
        DO K = KMIN1, KMAX2
           ZG_T(K) = ZG_T(K-1) + DZ(K)
        END DO
      ELSE
         ZG_T = ZERO
      ENDIF


      CALL SETUP_VTK_REGION

      CALL OPEN_VTU_FILE_BIN(MODE)

      IF(MODE==0) CALL OPEN_PVD_FILE

      CALL CLEAN_GEOMETRY

      DO PASS=WRITE_HEADER,WRITE_DATA


         CALL WRITE_GEOMETRY_IN_VTU_BIN(PASS)

         IF(VTK_EP_g(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('EP_G',EP_G,PASS)

         IF(VTK_P_g(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('P_G',P_G,PASS)

         IF(VTK_P_star(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('P_S',P_S,PASS)

         IF(VTK_VEL_g(VTK_REGION)) &
            CALL WRITE_VECTOR_IN_VTU_BIN('Gas_Velocity',U_G,V_G,W_G,PASS)

         IF(VTK_U_g(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('U_G',P_G,PASS)

         IF(VTK_V_g(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('V_G',P_G,PASS)

         IF(VTK_W_g(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('W_G',P_G,PASS)

         DO M = 1,MMAX
            IF(VTK_VEL_s(VTK_REGION,M)) THEN
               WRITE(SUBM,*)M
               CALL WRITE_VECTOR_IN_VTU_BIN('Solids_Velocity_'//ADJUSTL(SUBM),U_S(:,M),V_S(:,M),W_S(:,M),PASS)
            ENDIF
            IF(VTK_U_s(VTK_REGION,M)) THEN
               WRITE(SUBM,*)M
               CALL WRITE_SCALAR_IN_VTU_BIN('U_s_'//ADJUSTL(SUBM),U_S(:,M),PASS)
            ENDIF
            IF(VTK_V_s(VTK_REGION,M)) THEN
               WRITE(SUBM,*)M
               CALL WRITE_SCALAR_IN_VTU_BIN('V_s_'//ADJUSTL(SUBM),V_S(:,M),PASS)
            ENDIF
            IF(VTK_W_s(VTK_REGION,M)) THEN
               WRITE(SUBM,*)M
               CALL WRITE_SCALAR_IN_VTU_BIN('W_s_'//ADJUSTL(SUBM),W_S(:,M),PASS)
            ENDIF
         END DO

         DO M = 1,MMAX
            IF(VTK_ROP_s(VTK_REGION,M)) THEN
               WRITE(SUBM,*)M
               CALL WRITE_SCALAR_IN_VTU_BIN('Solids_density_'//ADJUSTL(SUBM),ROP_S(:,M),PASS)
            ENDIF
         END DO

         IF(VTK_T_g(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTU_BIN('Gas_temperature',T_g,PASS)

         DO M = 1,MMAX
            IF(VTK_T_s(VTK_REGION,M)) THEN
               WRITE(SUBM,*)M
               CALL WRITE_SCALAR_IN_VTU_BIN('Solids_temperature_'//ADJUSTL(SUBM),T_S(:,M),PASS)
            ENDIF
         END DO


         SPECIES_COUNTER = 0
         DO NN = 1,NMAX(0)
            IF(VTK_X_g(VTK_REGION,NN)) THEN
               WRITE(SUBN,*)NN
               IF(USE_RRATES) THEN
                  SPECIES_COUNTER = SPECIES_COUNTER + 1
                  VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                  LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
               ELSE
                  VAR_NAME = ADJUSTL(SPECIES_ALIAS_g(NN))
                  LT = LEN_TRIM(ADJUSTL(SPECIES_ALIAS_g(NN)))
               ENDIF
               VAR_NAME = VAR_NAME(1:LT)//'_Gas_mass_fractions_'//ADJUSTL(SUBN)
               CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,X_g(:,NN),PASS)
            ENDIF
         END DO

        DO M = 1, MMAX
           WRITE(SUBM,*)M
           DO NN = 1,NMAX(M)
              IF(VTK_X_s(VTK_REGION,M,NN)) THEN
                 WRITE(SUBN,*)NN
                 IF(USE_RRATES) THEN
                    SPECIES_COUNTER = SPECIES_COUNTER + 1
                    VAR_NAME = ADJUSTL(SPECIES_NAME(SPECIES_COUNTER))
                    LT = LEN_TRIM(ADJUSTL(SPECIES_NAME(SPECIES_COUNTER)))
                 ELSE
                    VAR_NAME = ADJUSTL(SPECIES_ALIAS_s(M,NN))
                    LT = LEN_TRIM(ADJUSTL(SPECIES_ALIAS_s(M,NN)))
                 ENDIF
                 VAR_NAME = VAR_NAME(1:LT)//'_Solids_mass_fractions_'//TRIM(ADJUSTL(SUBM))//'_'//ADJUSTL(SUBN)
                 CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,X_s(:,M,NN),PASS)
              ENDIF
           END DO
        END DO

        DO M = 1,MMAX
           IF(VTK_Theta_m(VTK_REGION,M)) THEN
              WRITE(SUBM,*)M
              CALL WRITE_SCALAR_IN_VTU_BIN('Granular_temperature_'//ADJUSTL(SUBM),Theta_m(:,M),PASS)
           ENDIF
        END DO

        DO NN = 1,NSCALAR
           IF(VTK_Scalar(VTK_REGION,NN)) THEN
              WRITE(SUBN,*)NN
              VAR_NAME = 'Scalar_'//ADJUSTL(SUBN)
              CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,Scalar(:,NN),PASS)
           ENDIF
        END DO

        DO R = 1,nRR
           IF(VTK_RRate(VTK_REGION,R)) THEN
              WRITE(SUBR,*)R
              VAR_NAME = 'RRates_'//ADJUSTL(SUBR)
              CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,ReactionRates(:, R),PASS)
           ENDIF
       END DO

       IF(K_EPSILON) THEN
          IF(VTK_K_Turb_G(VTK_REGION)) &
             CALL WRITE_SCALAR_IN_VTU_BIN('K_Turb_G',K_Turb_G,PASS)
          IF(VTK_E_Turb_G(VTK_REGION)) &
             CALL WRITE_SCALAR_IN_VTU_BIN('E_Turb_G',E_Turb_G,PASS)
       ENDIF


       IF(VTK_VORTICITY(VTK_REGION).OR.VTK_LAMBDA_2(VTK_REGION)) THEN
          CALL CALC_VORTICITY
       ENDIF

       IF(VTK_VORTICITY(VTK_REGION)) &
          CALL WRITE_SCALAR_IN_VTU_BIN('VORTICITY_MAG',VORTICITY,PASS)
       IF(VTK_LAMBDA_2(VTK_REGION)) &
          CALL WRITE_SCALAR_IN_VTU_BIN('LAMBDA_2',LAMBDA2,PASS)


       IF(VTK_PARTITION(VTK_REGION)) &
          CALL WRITE_SCALAR_IN_VTU_BIN('PARTITION',PARTITION,PASS)


       IF(VTK_BC_ID(VTK_REGION)) THEN
          Allocate(DP_BC_ID(DIMENSION_3))
          DP_BC_ID = DBLE(BC_ID)
          CALL WRITE_SCALAR_IN_VTU_BIN('BC_ID',DP_BC_ID,PASS)
          DeAllocate(DP_BC_ID)
       ENDIF


       IF(VTK_DWALL(VTK_REGION)) &
          CALL WRITE_SCALAR_IN_VTU_BIN('DISTANCE_TO_WALL',DWALL,PASS)

       IF(VTK_IJK(VTK_REGION)) THEN
         Allocate(IJK_ARRAY(DIMENSION_3))
         DO IJK = IJKSTART3, IJKEND3
            IJK_ARRAY(IJK) = DBLE(IJK)
         ENDDO
         CALL WRITE_SCALAR_IN_VTU_BIN('IJK',IJK_ARRAY,PASS)
         DeAllocate(IJK_ARRAY)
      ENDIF

       IF(VTK_IJK(VTK_REGION)) &
          CALL WRITE_VECTOR_IN_VTU_BIN('Scalar normal',NORMAL_S(:,1),NORMAL_S(:,2),NORMAL_S(:,3),PASS)

       DO NN=1,15
          IF(VTK_DEBUG(VTK_REGION,NN)) THEN
             WRITE(SUBN,*)NN
             VAR_NAME = 'DEBUG_'//ADJUSTL(SUBN)
             CALL WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,DEBUG_CG(:,NN),PASS)
          ENDIF
       ENDDO



      ENDDO ! PASS LOOP, EITHER HEADER OR DATA


      CALL CLOSE_VTU_FILE_BIN(MODE)
      IF(MODE==0) CALL UPDATE_AND_CLOSE_PVD_FILE

#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

! Update Frames
      IF (myPE == PE_IO.AND.TIME_DEPENDENT_FILENAME) THEN
         OPEN(CONVERT='BIG_ENDIAN',UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
         DO L = 1, DIMENSION_VTK
            IF(VTK_DEFINED(L)) WRITE(VTU_FRAME_UNIT,*) L,FRAME(L)
         ENDDO
         CLOSE(VTU_FRAME_UNIT)
      ENDIF

     IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,20)' DONE.'

20    FORMAT(A,1X/)
      RETURN

      END SUBROUTINE WRITE_VTU_FILE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_VTU_FILE                                          C
!  Purpose: Open a vtu file and writes the header                      C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_VTU_FILE_BIN(MODE)

      USE cdist
      USE compar
      USE constant
      USE cutcell
      USE exit, only: mfix_exit
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE mpi_utility
      USE output
      USE parallel
      USE param
      USE param1
      USE quadric
      USE run
      USE sendrecv
      USE toleranc
      USE vtk

      IMPLICIT NONE
      LOGICAL :: VTU_FRAME_FILE_EXISTS
      INTEGER :: ISTAT,BUFF1,BUFF2,L
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)

#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

! Only open the file from head node when not using distributed I/O
      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN

      IF(TRIM(VTU_DIR)/='.') CALL CREATE_DIR(trim(VTU_DIR))

      IF(TIME_DEPENDENT_FILENAME) THEN
         INQUIRE(FILE=VTU_FRAME_FILENAME,EXIST=VTU_FRAME_FILE_EXISTS)
         IF(VTU_FRAME_FILE_EXISTS) THEN
            OPEN(CONVERT='BIG_ENDIAN',UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
            DO L = 1, DIMENSION_VTK
               IF(VTK_DEFINED(L)) THEN
                  READ(VTU_FRAME_UNIT,*)BUFF1,BUFF2
                  FRAME(L)=BUFF2
               ENDIF
            ENDDO
            CLOSE(VTU_FRAME_UNIT)
         ENDIF
         IF(RESET_FRAME_AT_TIME_ZERO.AND.TIME==ZERO) THEN
            DO L = 1, DIMENSION_VTK
               IF(L==VTK_REGION) FRAME(L)=-1
            ENDDO
         ENDIF
         DO L = 1, DIMENSION_VTK
            IF(L==VTK_REGION) FRAME(L) = FRAME(L) + 1
         ENDDO
      ENDIF

      IF (BDIST_IO.AND.NUMBER_OF_VTK_CELLS>0) THEN


! For distributed I/O, define the file name for each processor
         IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
            WRITE(VTU_FILENAME,20) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION),MYPE
         ELSE
            WRITE(VTU_FILENAME,25) TRIM(VTK_FILEBASE(VTK_REGION)),MYPE
         ENDIF
      ELSE
         IF(MYPE.EQ.PE_IO) THEN
            IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
               WRITE(VTU_FILENAME,30) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION)
            ELSE
               WRITE(VTU_FILENAME,35) TRIM(VTK_FILEBASE(VTK_REGION))
            ENDIF
         END IF
      END IF

! Add the VTU directory path if necessary

      IF(TRIM(VTU_DIR)/='.') VTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//VTU_FILENAME

! Echo
      IF (FULL_LOG) THEN
         IF (.NOT.BDIST_IO) THEN
            WRITE(*,10,ADVANCE='NO')' WRITING VTU FILE : ', TRIM(VTU_FILENAME),' .'
         ELSE
            IF(myPE==PE_IO) WRITE(*,15,ADVANCE='NO')' EACH PROCESOR IS WRITING ITS OWN VTU FILE.'
         ENDIF
      ENDIF

! Open File
!      OPEN(CONVERT='BIG_ENDIAN',UNIT = VTU_UNIT, FILE = TRIM(VTU_FILENAME),FORM='BINARY',IOSTAT=ISTAT)


      IF(NUMBER_OF_VTK_CELLS>0) THEN

         OPEN(CONVERT='BIG_ENDIAN',UNIT     = VTU_UNIT,           &
              FILE     = TRIM(VTU_FILENAME), &
              FORM     = 'UNFORMATTED',      &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                                ! use 'BINARY' if 'UNFORMATTED' is not supported
              ACCESS   = 'STREAM',           &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                                ! use 'SEQUENTIAL' if 'STREAM' is not supported
              ACTION   = 'WRITE',            &
              IOSTAT=ISTAT)


         IF (ISTAT /= 0) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
            IF(FULL_LOG.AND.myPE == PE_IO) WRITE(*, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
            CALL MFIX_EXIT(myPE)
         ENDIF


    1001 FORMAT(/1X,70('*')//, ' From: OPEN_VTU_FILE',/,' Message: ',          &
            'Error opening vtu file. Terminating run.',/10X,          &
            'File name:  ',A,/10X,                                         &
            'DES_UNIT :  ',i4, /10X,                                       &
            'PLEASE VERIFY THAT VTU_DIR EXISTS: ', A, &
            /1X,70('*')/)


   ! Write file Header
         BUFFER='<?xml version="1.0"?>'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


         WRITE(BUFFER,110)'<!-- Time =',TIME,' sec. -->'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

         BUFFER='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

         BUFFER='  <UnstructuredGrid>'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

      ENDIF
! For distributed I/O, open .pvtu file that combines all *.vtu files for a given FRAME
! this is a simple ASCII file

      IF (myPE == PE_IO.AND.BDIST_IO) THEN

         IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
            WRITE(PVTU_FILENAME,40) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION)
         ELSE
            WRITE(PVTU_FILENAME,45) TRIM(VTK_FILEBASE(VTK_REGION))
         ENDIF

         IF(TRIM(VTU_DIR)/='.') PVTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//PVTU_FILENAME

         OPEN(CONVERT='BIG_ENDIAN',UNIT = PVTU_UNIT, FILE = TRIM(PVTU_FILENAME))

         WRITE(PVTU_UNIT,100) '<?xml version="1.0"?>'
         WRITE(PVTU_UNIT,110) '<!-- Time =',TIME,' sec. -->'
         WRITE(PVTU_UNIT,120) '<VTKFile type="PUnstructuredGrid"',&
                  ' version="0.1" byte_order="BigEndian">'

         WRITE(PVTU_UNIT,100) '  <PUnstructuredGrid GhostLevel="0">'
         WRITE(PVTU_UNIT,100) '      <PPoints>'
         WRITE(PVTU_UNIT,100) '        <PDataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
              &format="appended" offset=" 0" />'
         WRITE(PVTU_UNIT,100) '      </PPoints>'
         WRITE(PVTU_UNIT,100) ''
         WRITE(PVTU_UNIT,100) '      <PCellData Scalars="scalars">'

      ENDIF

100   FORMAT(A)
110   FORMAT(A,E14.7,A)
120   FORMAT(A,A)
10    FORMAT(/1X,3A)
15    FORMAT(/1X,A)
20    FORMAT(A,"_",I4.4,"_",I5.5,".vtu")
25    FORMAT(A,"_",I5.5,".vtu")
30    FORMAT(A,"_",I4.4,".vtu")
35    FORMAT(A,".vtu")
40    FORMAT(A,"_",I4.4,".pvtu")
45    FORMAT(A,".pvtu")

      RETURN

      END SUBROUTINE OPEN_VTU_FILE_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_GEOMETRY_IN_VTU_BIN                              C
!  Purpose: Write Geometry and connectivity in a vtu file              C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_GEOMETRY_IN_VTU_BIN(PASS)

      USE, INTRINSIC :: iso_c_binding
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      USE functions

      IMPLICIT NONE

      INTEGER :: IJK,L
      INTEGER :: OFFSET

      INTEGER :: CELL_TYPE

      REAL(c_float) :: float
      INTEGER(c_int) :: int

      INTEGER ::     nbytes_xyz,nbytes_connectivity,nbytes_offset,nbytes_type
      INTEGER ::     offset_xyz,offset_connectivity,offset_offset,offset_type

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2


! First a series of tags is written for the geometry (PASS=WRITE_HEADER)
!  - Coordinates
!  - Connectivity
!  - Connectivity offset
!  - cell types
!

! Since the data is appended (i.e., written after all tags), the
! offset, in number of bytes must be specified.  The offset includes
! the size of the data for each field, plus the size of the integer
! that stores the number of bytes.  this is why the offset of a field
! equals the offset of the previous field plus sizeof(int) plus the
! number of bytes of the field.

! Next, the actual data is written for the geometry (PASS=WRITE_DATA)
! The DATA is converted to single precision to save memory.

      IF (myPE == PE_IO.AND.(.NOT.BDIST_IO)) THEN
! The number of points and number of VTK cells is now computed in
! SETUP_VTK_REGION

! Number of bytes of each field
         nbytes_xyz          = NUMBER_OF_POINTS * 3 * sizeof(float)

         nbytes_connectivity = 0
         DO IJK = 1,IJKMAX3
            IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  nbytes_connectivity = nbytes_connectivity + GLOBAL_NUMBER_OF_NODES(IJK)
            ENDIF
         END DO
         nbytes_connectivity = nbytes_connectivity * sizeof(int)


         nbytes_offset       = NUMBER_OF_VTK_CELLS * sizeof(int)

         nbytes_type         = NUMBER_OF_VTK_CELLS * sizeof(int)


! Offset of each field
         offset_xyz = 0
         offset_connectivity = offset_xyz          + sizeof(int) + nbytes_xyz
         offset_offset       = offset_connectivity + sizeof(int) + nbytes_connectivity
         offset_type         = offset_offset       + sizeof(int) + nbytes_offset


         IF(PASS==WRITE_HEADER) THEN

            WRITE(BUFFER,100)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS,'" NumberOfCells="',NUMBER_OF_VTK_CELLS,'" >'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
                 &format="appended" offset="',offset_xyz,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="connectivity" format="appended" offset="', &
                 offset_connectivity,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="offsets" format="appended" offset="',offset_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="types" format="appended" offset="',offset_type,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            VTU_offset =  offset_type       + sizeof(int) + nbytes_type  ! Store offset for first variable to be written

            WRITE(BUFFER,110)'      <CellData>'                          ! Preparing CellData
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC




         ELSEIF(PASS==WRITE_DATA) THEN

            WRITE(BUFFER,110)'      </CellData>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'    </Piece>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  </UnstructuredGrid>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  <AppendedData encoding="raw">'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


! Starting raw binary data with an underscore

            WRITE(BUFFER,110)'_'
            WRITE(VTU_UNIT)TRIM(BUFFER)

! X,Y,Z coordinates
            WRITE(VTU_UNIT) nbytes_xyz, (GLOBAL_COORDS_OF_POINTS(1:3,L), L = 1,NUMBER_OF_POINTS)

! Connectivity
            WRITE(VTU_UNIT) nbytes_connectivity

            DO IJK = 1,IJKMAX3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  WRITE(VTU_UNIT) (GLOBAL_CLEANED_CONNECTIVITY(IJK,L)-1,L=1,GLOBAL_NUMBER_OF_NODES(IJK))
               ENDIF
            END DO

! Offsets
            WRITE(VTU_UNIT) nbytes_offset

            OFFSET = 0
            DO IJK = 1,IJKMAX3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  OFFSET = OFFSET + GLOBAL_NUMBER_OF_NODES(IJK)
                  WRITE(VTU_UNIT) OFFSET
               ENDIF
            END DO

! Types
            WRITE(VTU_UNIT)nbytes_type

            IF(NO_K) THEN
               CELL_TYPE = 7
            ELSE
               CELL_TYPE = 41
            ENDIF

            DO IJK = 1,IJKMAX3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK))  WRITE(VTU_UNIT) CELL_TYPE
            END DO


         ENDIF


      ELSEIF(BDIST_IO) THEN

! For distributed I/O, it works the same as above, except, the data is local to each processor
! First compute local number of cells and points

! The number of points and number of VTK cells is now computed in
! SETUP_VTK_REGION

! Number of bytes of each field
         nbytes_xyz          = NUMBER_OF_POINTS * 3 * sizeof(float)

         nbytes_connectivity = 0
         DO IJK = 1,IJKEND3
            IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  nbytes_connectivity = nbytes_connectivity + NUMBER_OF_NODES(IJK)
            ENDIF
         END DO
         nbytes_connectivity = nbytes_connectivity * sizeof(int)


         nbytes_offset       = NUMBER_OF_VTK_CELLS * sizeof(int)

         nbytes_type         = NUMBER_OF_VTK_CELLS * sizeof(int)


! Offset of each field
         offset_xyz = 0
         offset_connectivity = offset_xyz          + sizeof(int) + nbytes_xyz
         offset_offset       = offset_connectivity + sizeof(int) + nbytes_connectivity
         offset_type         = offset_offset       + sizeof(int) + nbytes_offset


         IF(PASS==WRITE_HEADER) THEN

            WRITE(BUFFER,100)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS,'" NumberOfCells="',NUMBER_OF_VTK_CELLS,'" >'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
                 &format="appended" offset="',offset_xyz,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      <Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="connectivity" format="appended" offset="', &
                 offset_connectivity,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="offsets" format="appended" offset="',offset_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,100)'        <DataArray type="Int32" Name="types" format="appended" offset="',offset_type,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'      </Cells>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            VTU_offset =  offset_type       + sizeof(int) + nbytes_type  ! Store offset for first variable to be written

            WRITE(BUFFER,110)'      <CellData>'                          ! Preparing CellData
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC




         ELSEIF(PASS==WRITE_DATA) THEN

            WRITE(BUFFER,110)'      </CellData>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,110)'    </Piece>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  </UnstructuredGrid>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            WRITE(BUFFER,110)'  <AppendedData encoding="raw">'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


! Starting raw binary data with an underscore

            WRITE(BUFFER,110)'_'
            WRITE(VTU_UNIT)TRIM(BUFFER)

! X,Y,Z coordinates
            WRITE(VTU_UNIT) nbytes_xyz, (COORDS_OF_POINTS(L,1:3), L = 1,NUMBER_OF_POINTS)

! Connectivity
            WRITE(VTU_UNIT) nbytes_connectivity

            DO IJK = 1,IJKEND3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  WRITE(VTU_UNIT) (CLEANED_CONNECTIVITY(IJK,L)-1,L=1,NUMBER_OF_NODES(IJK))
               ENDIF
            END DO

! Offsets
            WRITE(VTU_UNIT) nbytes_offset

            OFFSET = 0
            DO IJK = 1,IJKEND3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  OFFSET = OFFSET + NUMBER_OF_NODES(IJK)
                  WRITE(VTU_UNIT) OFFSET
               ENDIF
            END DO

! Types
            WRITE(VTU_UNIT)nbytes_type

            IF(NO_K) THEN
               CELL_TYPE = 7
            ELSE
               CELL_TYPE = 41
            ENDIF

            DO IJK = 1,IJKEND3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK))  WRITE(VTU_UNIT) CELL_TYPE
            END DO


         ENDIF


      ENDIF


100   FORMAT(A,I12,A,I12,A)
110   FORMAT(A)

      RETURN

      END SUBROUTINE WRITE_GEOMETRY_IN_VTU_BIN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SCALAR_IN_VTU_BIN                                C
!  Purpose: Write Scalar variable in a vtu file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_SCALAR_IN_VTU_BIN(VAR_NAME,VAR,PASS)

      USE, INTRINSIC :: iso_c_binding
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      USE output
      USE functions

      IMPLICIT NONE
      INTEGER :: I,IJK

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VAR
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VAR(:)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  TMP_VAR

      REAL(c_float) :: float

      INTEGER :: nbytes_scalar

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      IF (.NOT.BDIST_IO) THEN

! For each scalar, write a tag, with corresponding offset

         nbytes_scalar = NUMBER_OF_VTK_CELLS * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
!           For each scalar, write a tag, with corresponding offset

            DO I = 1,LEN_TRIM(VAR_NAME)
               IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
            ENDDO

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            VTU_offset = VTU_offset + sizeof(float) + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            IF (myPE == PE_IO) THEN
               allocate (GLOBAL_VAR(ijkmax3))
            ELSE
               allocate (GLOBAL_VAR(1))
            ENDIF

            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(VAR,TMP_VAR)
               CALL gather (TMP_VAR,GLOBAL_VAR,root)
            ELSE
               CALL gather (VAR,GLOBAL_VAR,root)
            ENDIF

            IF (myPE /= PE_IO) RETURN


            WRITE(VTU_UNIT) nbytes_scalar

            DO IJK = 1,IJKMAX3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) WRITE(VTU_UNIT) REAL(GLOBAL_VAR(IJK))
            ENDDO


            Deallocate (GLOBAL_VAR)

         ENDIF


      ELSE ! BDIST_IO=.TRUE.


         nbytes_scalar = NUMBER_OF_VTK_CELLS * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
!           For each scalar, write a tag, with corresponding offset

            DO I = 1,LEN_TRIM(VAR_NAME)
               IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
            ENDDO

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            VTU_offset = VTU_offset + sizeof(float) + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            WRITE(VTU_UNIT) nbytes_scalar

            DO IJK = 1,IJKEND3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) WRITE(VTU_UNIT) REAL(VAR(IJK))
            ENDDO

         ENDIF


         IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
            WRITE(PVTU_UNIT,90) '        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
         ENDIF


      ENDIF


      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I12,A)

      RETURN

      END SUBROUTINE WRITE_SCALAR_IN_VTU_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTU                                    C
!  Purpose: Write Vector variable in a vtu file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VECTOR_IN_VTU_BIN(VAR_NAME,VARX,VARY,VARZ,PASS)

      USE, INTRINSIC :: iso_c_binding
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      USE output
      USE functions

      IMPLICIT NONE
      INTEGER :: IJK

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  VARX,VARY,VARZ
      DOUBLE PRECISION, ALLOCATABLE :: GLOBAL_VARX(:),GLOBAL_VARY(:),GLOBAL_VARZ(:)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) ::  TMP_VAR

      REAL(c_float) :: float

      INTEGER :: nbytes_vector

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      IF (.NOT.BDIST_IO) THEN

         nbytes_vector = NUMBER_OF_VTK_CELLS * 3 * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
!           For each vector, write a tag, with corresponding offset

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            VTU_offset = VTU_offset + sizeof(float) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            IF (myPE == PE_IO) THEN
               allocate (GLOBAL_VARX(ijkmax3))
               allocate (GLOBAL_VARY(ijkmax3))
               allocate (GLOBAL_VARZ(ijkmax3))
            ELSE
               allocate (GLOBAL_VARX(1))
               allocate (GLOBAL_VARY(1))
               allocate (GLOBAL_VARZ(1))
            ENDIF

            IF(RE_INDEXING) THEN
               CALL UNSHIFT_DP_ARRAY(VARX,TMP_VAR)
               call gather (TMP_VAR,GLOBAL_VARX,root)

               CALL UNSHIFT_DP_ARRAY(VARY,TMP_VAR)
               call gather (TMP_VAR,GLOBAL_VARY,root)

               CALL UNSHIFT_DP_ARRAY(VARZ,TMP_VAR)
               call gather (TMP_VAR,GLOBAL_VARZ,root)

            ELSE
               call gather (VARX,GLOBAL_VARX,root)
               call gather (VARY,GLOBAL_VARY,root)
               call gather (VARZ,GLOBAL_VARZ,root)
            ENDIF


            IF (myPE /= PE_IO) RETURN


            WRITE(VTU_UNIT) nbytes_vector

            DO IJK = 1,IJKMAX3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  WRITE(VTU_UNIT) REAL(GLOBAL_VARX(IJK)),REAL(GLOBAL_VARY(IJK)),REAL(GLOBAL_VARZ(IJK))
               ENDIF
            ENDDO


            Deallocate (GLOBAL_VARX)
            Deallocate (GLOBAL_VARY)
            Deallocate (GLOBAL_VARZ)

         ENDIF


      ELSE ! BDIST_IO=.TRUE.


         nbytes_vector = NUMBER_OF_VTK_CELLS * 3 * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
!           For each vector, write a tag, with corresponding offset


            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


            VTU_offset = VTU_offset + sizeof(float) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN
!           and write the data, always preceded by its size in number of bytes

            WRITE(VTU_UNIT) nbytes_vector

            DO IJK = 1,IJKEND3
               IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                  WRITE(VTU_UNIT) REAL(VARX(IJK)),REAL(VARY(IJK)),REAL(VARZ(IJK))
               ENDIF
            ENDDO

         ENDIF


         IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
            WRITE(PVTU_UNIT,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
         ENDIF

      ENDIF


      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I12,A)

      RETURN

      END SUBROUTINE WRITE_VECTOR_IN_VTU_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_VTU_FILE_BIN                                     C
!  Purpose: Close a vtu file                                           C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_VTU_FILE_BIN(MODE)

      USE compar
      Use run
      USE vtk
      use cdist
      USE mpi_utility

      IMPLICIT NONE

      INTEGER:: N
      CHARACTER (LEN=32)  :: VTU_NAME
      INTEGER, DIMENSION(0:numPEs-1) :: ALL_VTK_CELL_COUNT
      INTEGER :: IERR
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)


      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN

      IF(NUMBER_OF_VTK_CELLS>0) THEN

! Write last tags and close the vtu file
          WRITE(BUFFER,110)'  </AppendedData>'
          WRITE(VTU_UNIT)END_REC//TRIM(BUFFER)//END_REC

          WRITE(BUFFER,110)'</VTKFile>'
          WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

         CLOSE(VTU_UNIT)

      ENDIF

! Update pvtu file and close

      IF(BDIST_IO)  CALL allgather_1i (NUMBER_OF_VTK_CELLS,ALL_VTK_CELL_COUNT,IERR)

      IF (myPE == PE_IO.AND.BDIST_IO) THEN
         WRITE(PVTU_UNIT,100) '      </PCellData>'

         DO N = 0,NumPEs-1
            IF(ALL_VTK_CELL_COUNT(N)>0) THEN
               IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
                  WRITE(VTU_NAME,20) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION),N
               ELSE
                  WRITE(VTU_NAME,25) TRIM(VTK_FILEBASE(VTK_REGION)),N
               ENDIF

               WRITE(PVTU_UNIT,110) '      <Piece Source="',TRIM(VTU_NAME),'"/>'
            ENDIF
         ENDDO


         WRITE(PVTU_UNIT,100) '  </PUnstructuredGrid>'
         WRITE(PVTU_UNIT,100) '</VTKFile>'
         CLOSE(PVTU_UNIT)
      ENDIF


20    FORMAT(A,"_",I4.4,"_",I5.5,".vtu")
25    FORMAT(A,"_",I5.5,".vtu")

100   FORMAT(A)
110   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE CLOSE_VTU_FILE_BIN






!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_PVD_FILE                                          C
!  Purpose: Open a PVD file and writes the header                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_PVD_FILE

      USE compar
      USE constant
      USE cutcell
      USE exit, only: mfix_exit
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE output
      USE parallel
      USE param
      USE param1
      USE quadric
      USE run
      USE sendrecv
      USE toleranc
      USE vtk

      IMPLICIT NONE
      LOGICAL :: PVD_EXISTS

      IF (myPE /= PE_IO) RETURN

      PVD_FILENAME = TRIM(VTK_FILEBASE(VTK_REGION)) // '.pvd'

! First, check if the file already exists.

      INQUIRE(FILE=PVD_FILENAME,EXIST=PVD_EXISTS)

! The first time this subroutine is executed, properly initialize the pvd file

      IF(.NOT.PVD_FILE_INITIALIZED(VTK_REGION)) THEN

         IF(RUN_TYPE == 'NEW'.OR.RUN_TYPE=='RESTART_2')THEN
            ! For a new or RESTART_2 run, the pvd file should not exist, and is created with appropriate header
            IF (.NOT.PVD_EXISTS) THEN
               OPEN(CONVERT='BIG_ENDIAN',UNIT = PVD_UNIT, FILE = TRIM(PVD_FILENAME))
               WRITE(PVD_UNIT,100) '<?xml version="1.0"?>'
               WRITE(PVD_UNIT,100) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
               WRITE(PVD_UNIT,100) '<Collection>'
!               CALL UPDATE_AND_CLOSE_PVD_FILE
               PVD_FILE_INITIALIZED(VTK_REGION)=.TRUE.
            ELSE ! If the pvd file exists, print error message and exits
               WRITE(*,1002) TRIM(PVD_FILENAME)
               WRITE(UNIT_LOG, 1002) TRIM(PVD_FILENAME)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
            ! For a restart_1 run, the pvd file must exist
            IF (.NOT.PVD_EXISTS) THEN
               ! If the pvd file does not exist, print error message and exits
               WRITE(*,1003) TRIM(PVD_FILENAME)
               WRITE(UNIT_LOG, 1003) TRIM(PVD_FILENAME)
               CALL MFIX_EXIT(myPE)
            ELSE
           ! If it already exists, go to the bottom of the file and prepare to append data (remove last two lines)
               OPEN(CONVERT='BIG_ENDIAN',UNIT=PVD_UNIT,FILE = TRIM(PVD_FILENAME),POSITION="APPEND",STATUS='OLD')
               BACKSPACE(PVD_UNIT)
               BACKSPACE(PVD_UNIT)
               PVD_FILE_INITIALIZED(VTK_REGION)=.TRUE.
            ENDIF
         ENDIF
      ELSE
         ! When properly initialized, open the file and go to the
         ! bottom of the file and prepare to append data (remove last two lines)
         OPEN(CONVERT='BIG_ENDIAN',UNIT=PVD_UNIT,FILE = TRIM(PVD_FILENAME),POSITION="APPEND",STATUS='OLD')
         BACKSPACE(PVD_UNIT)
         BACKSPACE(PVD_UNIT)
      ENDIF


100   FORMAT(A)

1002  FORMAT(/1X,70('*')/,' From: OPEN_PVD_FILE',/,' Message: ',       &
         A,' already exists in the run directory.',/10X,               &
           'This is not allowed for a new run.',/10X,                   &
         'Terminating run.',/1X,70('*')/)

1003  FORMAT(/1X,70('*')/,' From: OPEN_PVD_FILE',/,' Message: ',       &
         A,' is missing from the  the run directory,',/10X,            &
           ' and must be present for a restart run.',/10X,              &
         'Terminating run.',/1X,70('*')/)

      RETURN

      END SUBROUTINE OPEN_PVD_FILE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_AND_CLOSE_PVD_FILE                              C
!  Purpose: Updates and close a pvd file                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE UPDATE_AND_CLOSE_PVD_FILE

      USE compar
      use cdist
      USE run
      USE vtk

      IMPLICIT NONE

      CHARACTER (LEN=255)  :: FILENAME
      CHARACTER (LEN=5)   :: EXT

      IF (myPE /= PE_IO) RETURN

      IF(VTK_DATA(VTK_REGION)=='C') THEN
         EXT = '.pvtu'
      ELSEIF(VTK_DATA(VTK_REGION)=='P') THEN
         EXT = '.pvtp'
      ENDIF

      IF(.NOT.BDIST_IO) THEN
         FILENAME=VTU_FILENAME
      ELSE
         IF(TIME_DEPENDENT_FILENAME) THEN
            WRITE(FILENAME,40) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION),EXT
         ELSE
            WRITE(FILENAME,45) TRIM(VTK_FILEBASE(VTK_REGION)),EXT
         ENDIF
         IF(TRIM(VTU_DIR)/='.') FILENAME='./'//TRIM(VTU_DIR)//'/'//FILENAME
      ENDIF

! Write the data to the file
         WRITE(PVD_UNIT,100)&
         '<DataSet timestep="',TIME,'" ',& ! simulation time
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FILENAME),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(PVD_UNIT,110)'</Collection>'
         WRITE(PVD_UNIT,110)'</VTKFile>'

         CLOSE(PVD_UNIT)

! 40    FORMAT(A,"_",I4.4,".pvtu")
! 45    FORMAT(A,".pvtu")
40    FORMAT(A,"_",I4.4,A5)
45    FORMAT(A,A5)
100   FORMAT(6X,A,E14.7,5A)
110   FORMAT(A)

      RETURN

      END SUBROUTINE UPDATE_AND_CLOSE_PVD_FILE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_CUT_SURFACE_VTK                                  C
!  Purpose: Writes the cut cell surface in VTK format                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_CUT_SURFACE_VTK

      USE compar
      USE constant
      USE cut_cell_preproc, only: eval_f
      USE cutcell
      USE exit, only: mfix_exit
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE polygon
      USE quadric
      USE run
      USE sendrecv
      USE stl
      USE toleranc
      USE vtk

      IMPLICIT NONE

      INTEGER :: L,IJK,NODE
      INTEGER :: POINT_ID,POLY_COUNT,FACE_ID,Q_ID
      INTEGER :: N_CUT_FACE_NODES

      INTEGER NUMBER_OF_FACES
      INTEGER NUMBER_OF_SURFACE_POINTS

      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_CUT_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3)    :: NORMAL

      INTEGER, DIMENSION(DIMENSION_MAX_CUT_CELL,6) ::  FACE_CONNECTIVITY
      INTEGER, DIMENSION(DIMENSION_MAX_CUT_CELL)   ::  NUMBER_OF_CUT_FACE_POINTS

      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) ::  X_FACE_POINT
      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) ::  Y_FACE_POINT
      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) ::  Z_FACE_POINT

      DOUBLE PRECISION :: X_COPY,Y_COPY,Z_COPY,F_COPY

      LOGICAL :: CLIP_FLAG

      CHARACTER (LEN=255) :: FILENAME

      LOGICAL :: CORNER_POINT
      INTEGER :: NODE_OF_CORNER, IERROR

      IF(myPE/=0) RETURN

!======================================================================
!  Set-up connectivity for each cell, i.e., regular cells and cut cells
!======================================================================

      POLY_COUNT = 0

      NUMBER_OF_SURFACE_POINTS = 0

      NUMBER_OF_FACES = 0

      DO IJK = 1,IJKMAX3

         IF(GLOBAL_CUT_CELL_AT(IJK)) THEN

!======================================================================
!  Filter the connectivity to identify nodes belonging to cut face
!======================================================================


            NUMBER_OF_FACES = NUMBER_OF_FACES + 1

            N_CUT_FACE_NODES = 0

            CALL GET_GLOBAL_CELL_NODE_COORDINATES(IJK,'SCALAR')

            DO L = 1, GLOBAL_NUMBER_OF_NODES(IJK)
               IF(GLOBAL_CONNECTIVITY(IJK,L)>IJKMAX3) THEN   ! One of the new point
                  X_COPY = GLOBAL_X_NEW_POINT(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  Y_COPY = GLOBAL_Y_NEW_POINT(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  Z_COPY = GLOBAL_Z_NEW_POINT(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  CORNER_POINT = .FALSE.
               ELSE                                   ! An existing point
                  DO NODE = 1,8
                  CORNER_POINT = .TRUE.
                     IF(GLOBAL_CONNECTIVITY(IJK,L) == IJK_OF_NODE(NODE)) THEN
                        NODE_OF_CORNER = NODE
                        X_COPY = X_NODE(NODE)
                        Y_COPY = Y_NODE(NODE)
                        Z_COPY = Z_NODE(NODE)

                        IF (GLOBAL_SNAP(IJK_OF_NODE(NODE))) THEN ! One of the snapped corner point which now belongs to the cut face
                           N_CUT_FACE_NODES = N_CUT_FACE_NODES + 1
                           COORD_CUT_FACE_NODES(1,N_CUT_FACE_NODES) = X_COPY
                           COORD_CUT_FACE_NODES(2,N_CUT_FACE_NODES) = Y_COPY
                           COORD_CUT_FACE_NODES(3,N_CUT_FACE_NODES) = Z_COPY
                        ENDIF
                     ENDIF
                  END DO

               ENDIF

               IF(CORNER_POINT) THEN
                  Q_ID = 1

                  CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

                  CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

                  CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)

                  IF(USE_STL.OR.USE_MSH) F_COPY = GLOBAL_F_AT(IJK_OF_NODE(NODE_OF_CORNER))

!                  CALL EVAL_STL_FCT_AT('SCALAR',IJK,NODE_OF_CORNER,F_COPY,CLIP_FLAG,BCID2)
               ELSE
                  F_COPY = ZERO
               ENDIF

               IF (ABS(F_COPY) < TOL_F ) THEN ! belongs to cut face
                  N_CUT_FACE_NODES = N_CUT_FACE_NODES + 1
                  COORD_CUT_FACE_NODES(1,N_CUT_FACE_NODES) = X_COPY
                  COORD_CUT_FACE_NODES(2,N_CUT_FACE_NODES) = Y_COPY
                  COORD_CUT_FACE_NODES(3,N_CUT_FACE_NODES) = Z_COPY
               ENDIF

            END DO

            CALL REORDER_POLYGON(N_CUT_FACE_NODES,COORD_CUT_FACE_NODES,NORMAL,IERROR)

            NUMBER_OF_CUT_FACE_POINTS(NUMBER_OF_FACES) = N_CUT_FACE_NODES
            POLY_COUNT = POLY_COUNT + N_CUT_FACE_NODES + 1
            DO NODE = 1,N_CUT_FACE_NODES
               NUMBER_OF_SURFACE_POINTS = NUMBER_OF_SURFACE_POINTS + 1

               IF(NUMBER_OF_SURFACE_POINTS>=DIMENSION_MAX_CUT_CELL) THEN
                  WRITE(*,3000) 'ERROR IN SUBROUTINE WRITE_3DCUT_SURFACE_VTK:'
                  WRITE(*,3000) 'NUMBER_OF_SURFACE_POINTS>=DIMENSION_MAX_CUT_CELL:'
                  WRITE(*,3000) 'INCREASE VALUE OF FAC_DIM_MAX_CUT_CELL.'
                  WRITE(*,3010) 'CURRENT VALUE OF FAC_DIM_MAX_CUT_CELL =',FAC_DIM_MAX_CUT_CELL
                  WRITE(*,3020) 'CURRENT VALUE OF DIMENSION_MAX_CUT_CELL =',DIMENSION_MAX_CUT_CELL
                  WRITE(*,3000) 'MFiX will exit now.'
                  CALL MFIX_EXIT(myPE)
               ENDIF

               X_FACE_POINT(NUMBER_OF_SURFACE_POINTS) = COORD_CUT_FACE_NODES(1,NODE)
               Y_FACE_POINT(NUMBER_OF_SURFACE_POINTS) = COORD_CUT_FACE_NODES(2,NODE)
               Z_FACE_POINT(NUMBER_OF_SURFACE_POINTS) = COORD_CUT_FACE_NODES(3,NODE)
               FACE_CONNECTIVITY(NUMBER_OF_FACES,NODE) = NUMBER_OF_SURFACE_POINTS
            ENDDO

         ENDIF

      END DO

      FILENAME= TRIM(RUN_NAME) // '_boundary.vtk'
      FILENAME = TRIM(FILENAME)
      OPEN(CONVERT='BIG_ENDIAN',UNIT = 123, FILE = FILENAME)
      WRITE(123,1001)'# vtk DataFile Version 2.0'
      WRITE(123,1001)'3D CUT-CELL SURFACE'
      WRITE(123,1001)'ASCII'

      IF(NO_K) THEN   ! 2D GEOMETRY
         WRITE(123,1001)'DATASET UNSTRUCTURED_GRID'
      ELSE            ! 3D GEOMETRY
         WRITE(123,1001)'DATASET POLYDATA'
      ENDIF

      WRITE(123,1010)'POINTS ',NUMBER_OF_SURFACE_POINTS,' float'

      DO POINT_ID = 1,NUMBER_OF_SURFACE_POINTS
         WRITE(123,1020) X_FACE_POINT(POINT_ID),Y_FACE_POINT(POINT_ID),Z_FACE_POINT(POINT_ID)
      ENDDO

      IF(NO_K) THEN   ! 2D GEOMETRY

         WRITE(123,1030)'CELLS ',NUMBER_OF_FACES,POLY_COUNT
         DO FACE_ID = 1 , NUMBER_OF_FACES
            WRITE(123,1040) NUMBER_OF_CUT_FACE_POINTS(FACE_ID),(FACE_CONNECTIVITY(FACE_ID,L)-1,&
            L=1,NUMBER_OF_CUT_FACE_POINTS(FACE_ID))
         ENDDO
         WRITE(123,1030)'CELL_TYPES ',NUMBER_OF_FACES
         DO FACE_ID = 1 , NUMBER_OF_FACES
            WRITE(123,1040) 3
         ENDDO

      ELSE            ! 3D GEOMETRY

         WRITE(123,1030)'POLYGONS ',NUMBER_OF_FACES,POLY_COUNT
         DO FACE_ID = 1 , NUMBER_OF_FACES
            WRITE(123,1040) NUMBER_OF_CUT_FACE_POINTS(FACE_ID),(FACE_CONNECTIVITY(FACE_ID,L)-1,&
            L=1,NUMBER_OF_CUT_FACE_POINTS(FACE_ID))
         ENDDO

      ENDIF

1001  FORMAT(A)
1010  FORMAT(A,I8,A)
1020  FORMAT(3(E16.8,2X))
1030  FORMAT(A,2(I8,2X))
1040  FORMAT(20(I8,2X))
3000  FORMAT(1X,A)
3010  FORMAT(1X,A,F8.4)
3020  FORMAT(1X,A,I8)
3030  FORMAT(1X,A,A)
      CLOSE (123)


      WRITE(*,3030)'WROTE BOUNDARY IN VTK FILE : ',FILENAME
      RETURN


      END SUBROUTINE WRITE_CUT_SURFACE_VTK



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GATHER_DATA                                            C
!  Purpose: Gather data from all processes in preparation of           C
!           Writing VTK files                                          C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GATHER_DATA

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE functions

      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L
      INTEGER :: IJK_OFFSET

      INTEGER :: iproc,IERR
      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SHIFTED_CONNECTIVITY

!======================================================================
!  parallel processing
!======================================================================

      CALL allgather_1i (NUMBER_OF_NEW_POINTS,rcount,IERR)

      IF (myPE == 0) THEN
         IJK_OFFSET = 0
      ELSE
         IJK_OFFSET = 0
         DO iproc=0,myPE-1
            IJK_OFFSET = IJK_OFFSET + rcount(iproc)
         ENDDO
      ENDIF

      CALL allgather_1i (IJK_OFFSET,disp,IERR)

      IF(.NOT.GLOBAL_VAR_ALLOCATED) THEN

         IF (myPE == PE_IO) THEN
            allocate (GLOBAL_I_OF(ijkmax3))
            allocate (GLOBAL_J_OF(ijkmax3))
            allocate (GLOBAL_K_OF(ijkmax3))
            allocate (GLOBAL_CONNECTIVITY(ijkmax3,15))
            allocate (GLOBAL_NUMBER_OF_NODES(ijkmax3))
            allocate (GLOBAL_INTERIOR_CELL_AT(ijkmax3))
            allocate (GLOBAL_BLOCKED_CELL_AT(ijkmax3))
            allocate (GLOBAL_STANDARD_CELL_AT(ijkmax3))
            allocate (GLOBAL_CUT_CELL_AT(ijkmax3))
            allocate (GLOBAL_SNAP(ijkmax3))
            allocate (GLOBAL_X_NEW_POINT(ijkmax3))
            allocate (GLOBAL_Y_NEW_POINT(ijkmax3))
            allocate (GLOBAL_Z_NEW_POINT(ijkmax3))
            allocate (GLOBAL_F_AT(ijkmax3))

         ELSE
            allocate (GLOBAL_I_OF(1))
            allocate (GLOBAL_J_OF(1))
            allocate (GLOBAL_K_OF(1))
            allocate (GLOBAL_CONNECTIVITY(1,15))
            allocate (GLOBAL_NUMBER_OF_NODES(1))
            allocate (GLOBAL_INTERIOR_CELL_AT(1))
            allocate (GLOBAL_BLOCKED_CELL_AT(1))
            allocate (GLOBAL_STANDARD_CELL_AT(1))
            allocate (GLOBAL_CUT_CELL_AT(1))
            allocate (GLOBAL_SNAP(1))
            allocate (GLOBAL_X_NEW_POINT(1))
            allocate (GLOBAL_Y_NEW_POINT(1))
            allocate (GLOBAL_Z_NEW_POINT(1))
            allocate (GLOBAL_F_AT(1))
         ENDIF

         GLOBAL_VAR_ALLOCATED = .TRUE.

      ENDIF

      IF(numPEs==1) THEN  ! Serial run
         GLOBAL_X_NEW_POINT(1:NUMBER_OF_NEW_POINTS) =  X_NEW_POINT(1:NUMBER_OF_NEW_POINTS)
         GLOBAL_Y_NEW_POINT(1:NUMBER_OF_NEW_POINTS) =  Y_NEW_POINT(1:NUMBER_OF_NEW_POINTS)
         IF(DO_K) GLOBAL_Z_NEW_POINT(1:NUMBER_OF_NEW_POINTS) =  Z_NEW_POINT(1:NUMBER_OF_NEW_POINTS)
      ELSE !Parallel run
         call gatherv_1d( X_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_X_NEW_POINT, rcount, disp, PE_IO, ierr )
         call gatherv_1d( Y_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_Y_NEW_POINT, rcount, disp, PE_IO, ierr )
         call gatherv_1d( Z_NEW_POINT, NUMBER_OF_NEW_POINTS, GLOBAL_Z_NEW_POINT, rcount, disp, PE_IO, ierr )
      ENDIF

      call global_sum(NUMBER_OF_NEW_POINTS, GLOBAL_NUMBER_OF_NEW_POINTS,  PE_IO, ierr )

      Allocate(  SHIFTED_CONNECTIVITY  (DIMENSION_3,15) )

      SHIFTED_CONNECTIVITY = CONNECTIVITY

      WHERE (SHIFTED_CONNECTIVITY > IJKEND3)
         SHIFTED_CONNECTIVITY = SHIFTED_CONNECTIVITY - IJKEND3 + IJKMAX3 + disp(myPE)
      END WHERE

      DO IJK = IJKSTART3,IJKEND3
         DO L=1,NUMBER_OF_NODES(IJK)
            IF(CONNECTIVITY(IJK,L) <= IJKEND3) THEN
               I = I_OF(CONNECTIVITY(IJK,L))
               J = J_OF(CONNECTIVITY(IJK,L))
               K = K_OF(CONNECTIVITY(IJK,L))
               SHIFTED_CONNECTIVITY(IJK,L) = funijk_gl(I,J,K)
            ENDIF
         ENDDO
      ENDDO


      GLOBAL_INTERIOR_CELL_AT = .FALSE.
      GLOBAL_BLOCKED_CELL_AT  = .FALSE.
      GLOBAL_CUT_CELL_AT      = .FALSE.
      call gather (I_OF,GLOBAL_I_OF,root)
      call gather (J_OF,GLOBAL_J_OF,root)
      call gather (K_OF,GLOBAL_K_OF,root)
      call gather (SHIFTED_CONNECTIVITY,GLOBAL_CONNECTIVITY,root)
      call gather (NUMBER_OF_NODES,GLOBAL_NUMBER_OF_NODES,root)
      call gather (INTERIOR_CELL_AT,GLOBAL_INTERIOR_CELL_AT,root)
      call gather (BLOCKED_CELL_AT,GLOBAL_BLOCKED_CELL_AT,root)
      call gather (STANDARD_CELL_AT,GLOBAL_STANDARD_CELL_AT,root)
      call gather (CUT_CELL_AT,GLOBAL_CUT_CELL_AT,root)
      call gather (SNAP,GLOBAL_SNAP,root)
      call gather (F_AT,GLOBAL_F_AT,root)

      Deallocate(  SHIFTED_CONNECTIVITY )

      IF (myPE == PE_IO) THEN

         POLY_COUNTER = 0

         NUMBER_OF_CELLS = 0

         NUMBER_OF_CUT_CELLS = 0

         NUMBER_OF_BLOCKED_CELLS = 0

         NUMBER_OF_STANDARD_CELLS = 0

         DO IJK = 1, IJKMAX3

            IF (GLOBAL_INTERIOR_CELL_AT(IJK)) THEN

               NUMBER_OF_CELLS = NUMBER_OF_CELLS + 1

               IF (GLOBAL_BLOCKED_CELL_AT(IJK))  NUMBER_OF_BLOCKED_CELLS  = NUMBER_OF_BLOCKED_CELLS + 1
               IF (GLOBAL_STANDARD_CELL_AT(IJK)) NUMBER_OF_STANDARD_CELLS = NUMBER_OF_STANDARD_CELLS + 1
               IF (GLOBAL_CUT_CELL_AT(IJK))      NUMBER_OF_CUT_CELLS = NUMBER_OF_CUT_CELLS + 1

               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   POLY_COUNTER = POLY_COUNTER + GLOBAL_NUMBER_OF_NODES(IJK) + 1

            ENDIF

         END DO


         NUMBER_OF_POINTS = IJKMAX3 + GLOBAL_NUMBER_OF_NEW_POINTS

      ENDIF

      RETURN


      END SUBROUTINE GATHER_DATA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SETUP_VTK_NO_CUTCELL                                   C
!  Purpose: Setup VTK data for the regular grid (no cut cells)         C
!           This i scalled when CARTESIAN_GRID is .FALSE.              C
!           and WRITE>VTK_FILES is .TRUE.              .               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 14-Jan-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SETUP_VTK_NO_CUTCELL

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE functions

      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,L,NODE

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: SHIFTED_CONNECTIVITY

! Only a few arrays need to be allocated here simce we do not need
! all Cartesian-grid arrays
      IF(.NOT.ALLOCATED(XG_E)) Allocate( XG_E(0:DIMENSION_I) )
      IF(.NOT.ALLOCATED(YG_N)) Allocate( YG_N(0:DIMENSION_J) )
      IF(.NOT.ALLOCATED(ZG_T)) Allocate( ZG_T(0:DIMENSION_K) )

      IF(.NOT.ALLOCATED(INTERIOR_CELL_AT)) THEN
         Allocate(  INTERIOR_CELL_AT  (DIMENSION_3) )
         INTERIOR_CELL_AT = .FALSE.
      ENDIF

      IF(.NOT.ALLOCATED(BLOCKED_CELL_AT)) THEN
         Allocate(  BLOCKED_CELL_AT  (DIMENSION_3) )
         BLOCKED_CELL_AT = .FALSE.
      ENDIF

      IF(.NOT.ALLOCATED(STANDARD_CELL_AT)) THEN
         Allocate(  STANDARD_CELL_AT  (DIMENSION_3) )
         STANDARD_CELL_AT = .TRUE.
      ENDIF

      IF(.NOT.ALLOCATED(CUT_CELL_AT)) THEN
         Allocate(  CUT_CELL_AT  (DIMENSION_3) )
         CUT_CELL_AT = .FALSE.
      ENDIF

      IF(.NOT.ALLOCATED(NUMBER_OF_NODES))  Allocate(  NUMBER_OF_NODES  (DIMENSION_3) )
      NUMBER_OF_NODES= 0

      IF(.NOT.ALLOCATED(CONNECTIVITY))     Allocate(  CONNECTIVITY  (DIMENSION_3,15) )


! This is a shoter version of Get_cut_cell_Flags
      DO IJK = IJKSTART3, IJKEND3

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

         CALL GET_CELL_NODE_COORDINATES(IJK,'SCALAR')

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF(NO_K) THEN   ! 2D case

            INTERIOR_CELL_AT(IJK) = (     (I >= ISTART1 ).AND.(I <= IEND1 )  &
                                     .AND.(J >= JSTART1 ).AND.(J <= JEND1 ) )

         ELSE            ! 3D case

            INTERIOR_CELL_AT(IJK) = (     (I >= ISTART1 ).AND.(I <= IEND1 )  &
                                     .AND.(J >= JSTART1 ).AND.(J <= JEND1 )  &
                                     .AND.(K >= KSTART1 ).AND.(K <= KEND1 ) )

         ENDIF


         IF(INTERIOR_CELL_AT(IJK)) THEN

! Set up the connectivity: 4 nodes in 2D, 8 nodes in 3D
            IF(NO_K) THEN
               NUMBER_OF_NODES(IJK) = 4
               CONNECTIVITY(IJK,1) = IJK_OF_NODE(5)
               CONNECTIVITY(IJK,2) = IJK_OF_NODE(6)
               CONNECTIVITY(IJK,3) = IJK_OF_NODE(8)
               CONNECTIVITY(IJK,4) = IJK_OF_NODE(7)
            ELSE
               NUMBER_OF_NODES(IJK) = 8
               DO NODE = 1,8
                  CONNECTIVITY(IJK,NODE) = IJK_OF_NODE(NODE)
               END DO
            ENDIF

! If obstacles are defined, they will be flagged as blocked cells
! and will not be visible in the VTK files
            IF(WALL_AT(IJK)) THEN
               BLOCKED_CELL_AT(IJK)  = .TRUE.
               CUT_CELL_AT(IJK)      = .FALSE.
               STANDARD_CELL_AT(IJK) = .FALSE.
            ENDIF

         ENDIF

      ENDDO


!======================================================================
!  parallel processing
!======================================================================
      call SEND_RECEIVE_1D_LOGICAL(STANDARD_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_CELL_AT,2)

      Allocate(  SHIFTED_CONNECTIVITY  (DIMENSION_3,15) )

      SHIFTED_CONNECTIVITY = CONNECTIVITY

! Replace local node index by global node index before gathering the array
      DO IJK = IJKSTART3,IJKEND3
         DO L=1,NUMBER_OF_NODES(IJK)
            IF(CONNECTIVITY(IJK,L) <= IJKEND3) THEN
               I = I_OF(CONNECTIVITY(IJK,L))
               J = J_OF(CONNECTIVITY(IJK,L))
               K = K_OF(CONNECTIVITY(IJK,L))
               SHIFTED_CONNECTIVITY(IJK,L) = funijk_gl(I,J,K)
            ENDIF
         ENDDO
      ENDDO

! Allocate, initialize and gather arrays
      IF(.NOT.GLOBAL_VAR_ALLOCATED) THEN

         IF (myPE == PE_IO) THEN
            allocate (GLOBAL_I_OF(ijkmax3))
            allocate (GLOBAL_J_OF(ijkmax3))
            allocate (GLOBAL_K_OF(ijkmax3))
            allocate (GLOBAL_CONNECTIVITY(ijkmax3,15))
            allocate (GLOBAL_NUMBER_OF_NODES(ijkmax3))
            allocate (GLOBAL_INTERIOR_CELL_AT(ijkmax3))
            allocate (GLOBAL_BLOCKED_CELL_AT(ijkmax3))
            allocate (GLOBAL_STANDARD_CELL_AT(ijkmax3))
            allocate (GLOBAL_CUT_CELL_AT(ijkmax3))

         ELSE
            allocate (GLOBAL_I_OF(1))
            allocate (GLOBAL_J_OF(1))
            allocate (GLOBAL_K_OF(1))
            allocate (GLOBAL_CONNECTIVITY(1,15))
            allocate (GLOBAL_NUMBER_OF_NODES(1))
            allocate (GLOBAL_INTERIOR_CELL_AT(1))
            allocate (GLOBAL_BLOCKED_CELL_AT(1))
            allocate (GLOBAL_STANDARD_CELL_AT(1))
            allocate (GLOBAL_CUT_CELL_AT(1))
         ENDIF

         GLOBAL_VAR_ALLOCATED = .TRUE.

      ENDIF


      GLOBAL_INTERIOR_CELL_AT = .FALSE.
      GLOBAL_BLOCKED_CELL_AT  = .FALSE.
      GLOBAL_CUT_CELL_AT      = .FALSE.
      GLOBAL_STANDARD_CELL_AT = .TRUE.

      call gather (I_OF,GLOBAL_I_OF,root)
      call gather (J_OF,GLOBAL_J_OF,root)
      call gather (K_OF,GLOBAL_K_OF,root)
      call gather (SHIFTED_CONNECTIVITY,GLOBAL_CONNECTIVITY,root)
      call gather (NUMBER_OF_NODES,GLOBAL_NUMBER_OF_NODES,root)
      call gather (INTERIOR_CELL_AT,GLOBAL_INTERIOR_CELL_AT,root)
      call gather (BLOCKED_CELL_AT,GLOBAL_BLOCKED_CELL_AT,root)
      call gather (STANDARD_CELL_AT,GLOBAL_STANDARD_CELL_AT,root)
      call gather (CUT_CELL_AT,GLOBAL_CUT_CELL_AT,root)

      deAllocate(  SHIFTED_CONNECTIVITY   )



! Count the number of cells
      GLOBAL_NUMBER_OF_NEW_POINTS = 0

      IF (myPE == PE_IO) THEN

         POLY_COUNTER = 0

         NUMBER_OF_CELLS = 0

         NUMBER_OF_CUT_CELLS = 0

         NUMBER_OF_BLOCKED_CELLS = 0

         NUMBER_OF_STANDARD_CELLS = 0

         DO IJK = 1, IJKMAX3

            IF (GLOBAL_INTERIOR_CELL_AT(IJK)) THEN

               NUMBER_OF_CELLS = NUMBER_OF_CELLS + 1

               IF (GLOBAL_BLOCKED_CELL_AT(IJK))  NUMBER_OF_BLOCKED_CELLS  = NUMBER_OF_BLOCKED_CELLS + 1
               IF (GLOBAL_STANDARD_CELL_AT(IJK)) NUMBER_OF_STANDARD_CELLS = NUMBER_OF_STANDARD_CELLS + 1
               IF (GLOBAL_CUT_CELL_AT(IJK))      NUMBER_OF_CUT_CELLS = NUMBER_OF_CUT_CELLS + 1

               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK))   POLY_COUNTER = POLY_COUNTER + GLOBAL_NUMBER_OF_NODES(IJK) + 1

            ENDIF

         END DO

! There are no new points since there a no cut cells
         NUMBER_OF_POINTS = IJKMAX3

      ENDIF

      RETURN


      END SUBROUTINE SETUP_VTK_NO_CUTCELL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PRINT_GRID_STATISTICS                                  C
!  Purpose: PRINT_GRID_STATISTICS ON SCREEN                            C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE PRINT_GRID_STATISTICS

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE functions

      IMPLICIT NONE

      INTEGER :: IJK

      INTEGER :: IERR

      DOUBLE PRECISION :: MIN_VOL, MAX_VOL, GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      DOUBLE PRECISION :: MIN_AYZ, MAX_AYZ, GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
      DOUBLE PRECISION :: MIN_AXZ, MAX_AXZ, GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
      DOUBLE PRECISION :: MIN_AXY, MAX_AXY, GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
      DOUBLE PRECISION :: MIN_CUT, MAX_CUT, GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
      DOUBLE PRECISION :: LOCAL_MIN_Q,LOCAL_MAX_Q, GLOBAL_MIN_Q,GLOBAL_MAX_Q

      IF (myPE == PE_IO) THEN

         IF(.NOT.GRID_INFO_PRINTED_ON_SCREEN) THEN
            WRITE(*,5) 'GRID STATISTICS:'
            WRITE(*,5) 'NUMBER OF CELLS          = ', NUMBER_OF_CELLS
            WRITE(*,10)'NUMBER OF STANDARD CELLS = ', &
                        NUMBER_OF_STANDARD_CELLS,DBLE(NUMBER_OF_STANDARD_CELLS) / DBLE(NUMBER_OF_CELLS) * 100.0D0
            WRITE(*,10)'NUMBER OF CUT CELLS      = ', &
                        NUMBER_OF_CUT_CELLS,DBLE(NUMBER_OF_CUT_CELLS) / DBLE(NUMBER_OF_CELLS) * 100.0D0
            WRITE(*,10)'NUMBER OF BLOCKED CELLS  = ', &
                        NUMBER_OF_BLOCKED_CELLS,DBLE(NUMBER_OF_BLOCKED_CELLS) / DBLE(NUMBER_OF_CELLS) * 100.0D0

5           FORMAT(1X,A,I8)
10          FORMAT(1X,A,I8,' (',F6.2,' % of Total)')

         ENDIF

         GRID_INFO_PRINTED_ON_SCREEN = .TRUE.

      ENDIF


!======================================================================
!  Scalar Cell volumes and areas
!======================================================================

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER
      MAX_AYZ = - LARGE_NUMBER
      MIN_AXZ =   LARGE_NUMBER
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(STANDARD_CELL_AT(IJK)) THEN              ! STANDARD CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
         WRITE(UNIT_CUT_CELL_LOG,1000)  '                       CELLS STATISTICS                         '
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'SCALAR STANDARD CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      ENDIF


      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER
      MAX_AYZ = - LARGE_NUMBER
      MIN_AXZ =   LARGE_NUMBER
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_CELL_AT(IJK)) THEN                   ! CUT CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'SCALAR CUT CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF SMALL SCALAR CELLS   = ', NUMBER_OF_SMALL_CELLS
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
      ENDIF

1000 FORMAT(A,E14.7,2X,E14.7)
1010 FORMAT(A,I8)

!======================================================================
!  U-Momentum Cell volumes and areas
!======================================================================

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER
      MAX_AYZ = - LARGE_NUMBER
      MIN_AXZ =   LARGE_NUMBER
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(STANDARD_U_CELL_AT(IJK)) THEN              ! STANDARD CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_U(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_U(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_U(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_U(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_U(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_U(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_U(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_U(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'U-MOMENTUM STANDARD CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      ENDIF

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER
      MAX_AYZ = - LARGE_NUMBER
      MIN_AXZ =   LARGE_NUMBER
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER
      MIN_CUT =   LARGE_NUMBER
      MAX_CUT = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_U_CELL_AT(IJK).AND.(.NOT.WALL_U_AT(IJK))) THEN      ! CUT CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_U(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_U(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_U(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_U(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_U(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_U(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_U(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_U(IJK))
            MIN_CUT =   DMIN1(MIN_CUT,AREA_U_CUT(IJK))
            MAX_CUT =   DMAX1(MAX_CUT,AREA_U_CUT(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )
      call global_min(MIN_CUT, GLOBAL_MIN_CUT,  PE_IO, ierr )
      call global_max(MAX_CUT, GLOBAL_MAX_CUT,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'U-MOMENTUM CUT CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF CUT AREA              = ', GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF U WALL CELLS         = ', NUMBER_OF_U_WALL_CELLS
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
      ENDIF
!======================================================================
!  V-Momentum Cell volumes and areas
!======================================================================


      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER
      MAX_AYZ = - LARGE_NUMBER
      MIN_AXZ =   LARGE_NUMBER
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(STANDARD_V_CELL_AT(IJK)) THEN              ! STANDARD CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_V(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_V(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_V(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_V(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_V(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_V(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_V(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_V(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'V-MOMENTUM STANDARD CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
      ENDIF

      MIN_VOL =   LARGE_NUMBER
      MAX_VOL = - LARGE_NUMBER
      MIN_AYZ =   LARGE_NUMBER
      MAX_AYZ = - LARGE_NUMBER
      MIN_AXZ =   LARGE_NUMBER
      MAX_AXZ = - LARGE_NUMBER
      MIN_AXY =   LARGE_NUMBER
      MAX_AXY = - LARGE_NUMBER
      MIN_CUT =   LARGE_NUMBER
      MAX_CUT = - LARGE_NUMBER

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_V_CELL_AT(IJK).AND.(.NOT.WALL_V_AT(IJK))) THEN      ! CUT CELLS

            MIN_VOL =   DMIN1(MIN_VOL,VOL_V(IJK))
            MAX_VOL =   DMAX1(MAX_VOL,VOL_V(IJK))
            MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_V(IJK))
            MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_V(IJK))
            MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_V(IJK))
            MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_V(IJK))
            MIN_AXY =   DMIN1(MIN_AXY,AXY_V(IJK))
            MAX_AXY =   DMAX1(MAX_AXY,AXY_V(IJK))
            MIN_CUT =   DMIN1(MIN_CUT,AREA_V_CUT(IJK))
            MAX_CUT =   DMAX1(MAX_CUT,AREA_V_CUT(IJK))

         ENDIF
      END DO

      call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
      call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
      call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
      call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
      call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
      call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
      call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
      call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )
      call global_min(MIN_CUT, GLOBAL_MIN_CUT,  PE_IO, ierr )
      call global_max(MAX_CUT, GLOBAL_MAX_CUT,  PE_IO, ierr )

      IF (myPE == PE_IO) THEN
         WRITE(UNIT_CUT_CELL_LOG,1000)  'V-MOMENTUM CUT CELLS:'
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF CUT AREA              = ', GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
         WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF V WALL CELLS         = ', NUMBER_OF_V_WALL_CELLS
         WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
      ENDIF

!======================================================================
!  W-Momentum Cell volumes and areas
!======================================================================


      IF(DO_K) THEN

         MIN_VOL =   LARGE_NUMBER
         MAX_VOL = - LARGE_NUMBER
         MIN_AYZ =   LARGE_NUMBER
         MAX_AYZ = - LARGE_NUMBER
         MIN_AXZ =   LARGE_NUMBER
         MAX_AXZ = - LARGE_NUMBER
         MIN_AXY =   LARGE_NUMBER
         MAX_AXY = - LARGE_NUMBER

         DO IJK = IJKSTART3, IJKEND3
            IF(STANDARD_W_CELL_AT(IJK)) THEN              ! STANDARD CELLS

               MIN_VOL =   DMIN1(MIN_VOL,VOL_W(IJK))
               MAX_VOL =   DMAX1(MAX_VOL,VOL_W(IJK))
               MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_W(IJK))
               MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_W(IJK))
               MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_W(IJK))
               MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_W(IJK))
               MIN_AXY =   DMIN1(MIN_AXY,AXY_W(IJK))
               MAX_AXY =   DMAX1(MAX_AXY,AXY_W(IJK))

            ENDIF
         END DO

         call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
         call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
         call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
         call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
         call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
         call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
         call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
         call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )

         IF (myPE == PE_IO) THEN
            WRITE(UNIT_CUT_CELL_LOG,1000)  'W-MOMENTUM STANDARD CELLS:'
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
         ENDIF

         MIN_VOL =   LARGE_NUMBER
         MAX_VOL = - LARGE_NUMBER
         MIN_AYZ =   LARGE_NUMBER
         MAX_AYZ = - LARGE_NUMBER
         MIN_AXZ =   LARGE_NUMBER
         MAX_AXZ = - LARGE_NUMBER
         MIN_AXY =   LARGE_NUMBER
         MAX_AXY = - LARGE_NUMBER
         MIN_CUT =   LARGE_NUMBER
         MAX_CUT = - LARGE_NUMBER

         DO IJK = IJKSTART3, IJKEND3
            IF(CUT_W_CELL_AT(IJK).AND.(.NOT.WALL_W_AT(IJK))) THEN      ! CUT CELLS

               MIN_VOL =   DMIN1(MIN_VOL,VOL_W(IJK))
               MAX_VOL =   DMAX1(MAX_VOL,VOL_W(IJK))
               MIN_AYZ =   DMIN1(MIN_AYZ,AYZ_W(IJK))
               MAX_AYZ =   DMAX1(MAX_AYZ,AYZ_W(IJK))
               MIN_AXZ =   DMIN1(MIN_AXZ,AXZ_W(IJK))
               MAX_AXZ =   DMAX1(MAX_AXZ,AXZ_W(IJK))
               MIN_AXY =   DMIN1(MIN_AXY,AXY_W(IJK))
               MAX_AXY =   DMAX1(MAX_AXY,AXY_W(IJK))
               MIN_CUT =   DMIN1(MIN_CUT,AREA_W_CUT(IJK))
               MAX_CUT =   DMAX1(MAX_CUT,AREA_W_CUT(IJK))

            ENDIF
         END DO

         call global_min(MIN_VOL, GLOBAL_MIN_VOL,  PE_IO, ierr )
         call global_max(MAX_VOL, GLOBAL_MAX_VOL,  PE_IO, ierr )
         call global_min(MIN_AYZ, GLOBAL_MIN_AYZ,  PE_IO, ierr )
         call global_max(MAX_AYZ, GLOBAL_MAX_AYZ,  PE_IO, ierr )
         call global_min(MIN_AXZ, GLOBAL_MIN_AXZ,  PE_IO, ierr )
         call global_max(MAX_AXZ, GLOBAL_MAX_AXZ,  PE_IO, ierr )
         call global_min(MIN_AXY, GLOBAL_MIN_AXY,  PE_IO, ierr )
         call global_max(MAX_AXY, GLOBAL_MAX_AXY,  PE_IO, ierr )
         call global_min(MIN_CUT, GLOBAL_MIN_CUT,  PE_IO, ierr )
         call global_max(MAX_CUT, GLOBAL_MAX_CUT,  PE_IO, ierr )

         IF (myPE == PE_IO) THEN
            WRITE(UNIT_CUT_CELL_LOG,1000)  'W-MOMENTUM CUT CELLS:'
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXY                   = ', GLOBAL_MIN_AXY,GLOBAL_MAX_AXY
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AXZ                   = ', GLOBAL_MIN_AXZ,GLOBAL_MAX_AXZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF AYZ                   = ', GLOBAL_MIN_AYZ,GLOBAL_MAX_AYZ
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF CUT AREA              = ', GLOBAL_MIN_CUT,GLOBAL_MAX_CUT
            WRITE(UNIT_CUT_CELL_LOG,1000)  'RANGE OF VOLUME                = ', GLOBAL_MIN_VOL,GLOBAL_MAX_VOL
            WRITE(UNIT_CUT_CELL_LOG,1010)  'NUMBER OF W WALL CELLS         = ', NUMBER_OF_W_WALL_CELLS
            WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'
         ENDIF

      ENDIF



      LOCAL_MIN_Q = MINVAL(Alpha_Ue_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Ue_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO)  WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Alpha_Ue_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Alpha_Un_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Un_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Alpha_Un_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Alpha_Ut_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Ut_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Alpha_Ut_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)

      LOCAL_MIN_Q = MINVAL(Theta_Ue)
      LOCAL_MAX_Q = MAXVAL(Theta_Ue)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_Ue   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Theta_Un)
      LOCAL_MAX_Q = MAXVAL(Theta_Un)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_Un   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Theta_Ut)
      LOCAL_MAX_Q = MAXVAL(Theta_Ut)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_Ut   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)

      LOCAL_MIN_Q = MINVAL(Theta_U_ne)
      LOCAL_MAX_Q = MAXVAL(Theta_U_ne)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_U_ne = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(Theta_U_te)
      LOCAL_MAX_Q = MAXVAL(Theta_U_te)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM Theta_U_te = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)

      LOCAL_MIN_Q = MINVAL(NOC_U_E)
      LOCAL_MAX_Q = MAXVAL(NOC_U_E)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM NOC_U_E    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(NOC_U_N)
      LOCAL_MAX_Q = MAXVAL(NOC_U_N)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM NOC_U_N    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q

      LOCAL_MIN_Q = MINVAL(NOC_U_T)
      LOCAL_MAX_Q = MAXVAL(NOC_U_T)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM NOC_U_T    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)

      LOCAL_MIN_Q = MINVAL(DELH_U)
      LOCAL_MAX_Q = MAXVAL(DELH_U)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF U-MOMENTUM DELH_U     = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'



      LOCAL_MIN_Q = MINVAL(Alpha_Ve_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Ve_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Alpha_Ve_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Alpha_Vn_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Vn_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Alpha_Vn_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Alpha_Vt_c)
      LOCAL_MAX_Q = MAXVAL(Alpha_Vt_c)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Alpha_Vt_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
      LOCAL_MIN_Q = MINVAL(Theta_Ve)
      LOCAL_MAX_Q = MAXVAL(Theta_Ve)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_Ve   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Theta_Vn)
      LOCAL_MAX_Q = MAXVAL(Theta_Vn)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_Vn   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Theta_Vt)
      LOCAL_MAX_Q = MAXVAL(Theta_Vt)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_Vt   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
      LOCAL_MIN_Q = MINVAL(Theta_V_ne)
      LOCAL_MAX_Q = MAXVAL(Theta_V_ne)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_V_ne = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(Theta_V_nt)
      LOCAL_MAX_Q = MAXVAL(Theta_V_nt)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM Theta_V_nt = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
      LOCAL_MIN_Q = MINVAL(NOC_V_E)
      LOCAL_MAX_Q = MAXVAL(NOC_V_E)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM NOC_V_E    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(NOC_V_N)
      LOCAL_MAX_Q = MAXVAL(NOC_V_N)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM NOC_V_N    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      LOCAL_MIN_Q = MINVAL(NOC_V_T)
      LOCAL_MAX_Q = MAXVAL(NOC_V_T)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM NOC_V_T    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
      LOCAL_MIN_Q = MINVAL(DELH_V)
      LOCAL_MAX_Q = MAXVAL(DELH_V)
      call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
      call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF V-MOMENTUM DELH_V     = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
      IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'


      IF(DO_K) THEN

         LOCAL_MIN_Q = MINVAL(Alpha_We_c)
         LOCAL_MAX_Q = MAXVAL(Alpha_We_c)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Alpha_We_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Alpha_Wn_c)
         LOCAL_MAX_Q = MAXVAL(Alpha_Wn_c)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Alpha_Wn_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Alpha_Wt_c)
         LOCAL_MAX_Q = MAXVAL(Alpha_Wt_c)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Alpha_Wt_c = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
         LOCAL_MIN_Q = MINVAL(Theta_We)
         LOCAL_MAX_Q = MAXVAL(Theta_We)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_We   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Theta_Wn)
         LOCAL_MAX_Q = MAXVAL(Theta_Wn)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_Wn   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Theta_Wt)
         LOCAL_MAX_Q = MAXVAL(Theta_Wt)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_Wt   = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
         LOCAL_MIN_Q = MINVAL(Theta_W_te)
         LOCAL_MAX_Q = MAXVAL(Theta_W_te)
         call global_min(LOCAL_MAX_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MIN_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_W_te = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(Theta_W_tn)
         LOCAL_MAX_Q = MAXVAL(Theta_W_tn)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM Theta_W_tn = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
         LOCAL_MIN_Q = MINVAL(NOC_W_E)
         LOCAL_MAX_Q = MAXVAL(NOC_W_E)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM NOC_W_E    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(NOC_W_N)
         LOCAL_MAX_Q = MAXVAL(NOC_W_N)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM NOC_W_N    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         LOCAL_MIN_Q = MINVAL(NOC_W_T)
         LOCAL_MAX_Q = MAXVAL(NOC_W_T)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM NOC_W_T    = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)
         LOCAL_MIN_Q = MINVAL(DELH_W)
         LOCAL_MAX_Q = MAXVAL(DELH_W)
         call global_min(LOCAL_MIN_Q, GLOBAL_MIN_Q,  PE_IO, ierr )
         call global_max(LOCAL_MAX_Q, GLOBAL_MAX_Q,  PE_IO, ierr )
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  ' RANGE OF W-MOMENTUM DELH_W     = ', GLOBAL_MIN_Q, GLOBAL_MAX_Q
         IF (myPE == PE_IO) WRITE(UNIT_CUT_CELL_LOG,1000)  '################################################################'

      ENDIF

      RETURN

      END SUBROUTINE PRINT_GRID_STATISTICS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLEAN_GEOMETRY                                         C
!  Purpose: Clean-up the list of point and only keep points            C
!           that are used in the connectivity list.                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 19-Dec-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CLEAN_GEOMETRY

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk
      USE cdist
      USE functions

      IMPLICIT NONE

      INTEGER :: IJK,L

      INTEGER ::POINT_ID,IJKC
      INTEGER , ALLOCATABLE        ::  POINT_NEW_ID(:)
      INTEGER , ALLOCATABLE        ::  NEW_POINT_NEW_ID(:)
      LOGICAL , ALLOCATABLE        ::  KEEP_POINT(:)
      LOGICAL , ALLOCATABLE        ::  KEEP_NEW_POINT(:)


      IF (myPE == PE_IO.AND.(.NOT.BDIST_IO)) THEN

         IF(ALLOCATED(GLOBAL_CLEANED_CONNECTIVITY)) DEALLOCATE(GLOBAL_CLEANED_CONNECTIVITY)
         IF(ALLOCATED(KEEP_NEW_POINT))              DEALLOCATE (KEEP_NEW_POINT)
         IF(ALLOCATED(POINT_NEW_ID))                DEALLOCATE (POINT_NEW_ID)
         IF(ALLOCATED(NEW_POINT_NEW_ID))            DEALLOCATE (NEW_POINT_NEW_ID)
         IF(ALLOCATED(KEEP_POINT))                  DEALLOCATE (KEEP_POINT)

         ALLOCATE (GLOBAL_CLEANED_CONNECTIVITY(ijkmax3,15))
         ALLOCATE (KEEP_NEW_POINT(GLOBAL_NUMBER_OF_NEW_POINTS))

         ALLOCATE (POINT_NEW_ID(IJKMAX3))
         ALLOCATE (NEW_POINT_NEW_ID(IJKMAX3))
         ALLOCATE (KEEP_POINT(IJKMAX3))

! Step 1: Go through connectivity list and only keep points that are used.
!         For background cell corners, assign KEEP_POINT = .TRUE.
!         For cut cells, the new intersection points were called NEW_POINTS,
!         so assign KEEP_NEW_POINT = .TRUE.
!         A NEW_POINT had an IJK index larger than IJKMAX3

         KEEP_POINT = .FALSE.
         KEEP_NEW_POINT = .FALSE.

         DO IJK = 1,IJKMAX3
            IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
               DO L=1,GLOBAL_NUMBER_OF_NODES(IJK)
                  IJKC = GLOBAL_CONNECTIVITY(IJK,L)
                  IF(IJKC<=IJKMAX3) KEEP_POINT(IJKC) = .TRUE.
                  IF(IJKC>IJKMAX3) KEEP_NEW_POINT(IJKC-IJKMAX3) = .TRUE.
               ENDDO
            ENDIF
         END DO


! Step 2: Clean-up list of used points and store cleaned connectivity
         POINT_NEW_ID = -1
         NEW_POINT_NEW_ID = -1
         POINT_ID = 1
! This is for the background grid cell corners
         DO IJK = 1,IJKMAX3
            IF(KEEP_POINT(IJK)) THEN
               POINT_NEW_ID(IJK) = POINT_ID
               POINT_ID = POINT_ID + 1
            ENDIF
         END DO
! This is for the cut cell new corners
         DO IJK = 1,GLOBAL_NUMBER_OF_NEW_POINTS
            IF(KEEP_NEW_POINT(IJK)) THEN
               NEW_POINT_NEW_ID(IJK) = POINT_ID
               POINT_ID = POINT_ID + 1
            ENDIF
         END DO

! Update the true (clean) number of points
         NUMBER_OF_POINTS = POINT_ID - 1

! Now, store a list of coordinates for all used points
         IF(ALLOCATED(GLOBAL_COORDS_OF_POINTS)) DEALLOCATE(GLOBAL_COORDS_OF_POINTS)

         ALLOCATE(GLOBAL_COORDS_OF_POINTS(3,NUMBER_OF_POINTS))

         POINT_ID = 1
! This is for the background grid cell corners
         DO IJK = 1,IJKMAX3
            IF(KEEP_POINT(IJK)) THEN
               GLOBAL_COORDS_OF_POINTS(1:3,POINT_ID) = &
                    (/REAL(XG_E(GLOBAL_I_OF(IJK))),REAL(YG_N(GLOBAL_J_OF(IJK))),REAL(ZG_T(GLOBAL_K_OF(IJK)))/)
               POINT_ID = POINT_ID + 1
            ENDIF
         END DO
! This is for the cut cell new corners
         DO IJK = 1,GLOBAL_NUMBER_OF_NEW_POINTS
            IF(KEEP_NEW_POINT(IJK)) THEN
               NEW_POINT_NEW_ID(IJK) = POINT_ID
               GLOBAL_COORDS_OF_POINTS(1:3,POINT_ID) = &
                    (/REAL(GLOBAL_X_NEW_POINT(IJK)),REAL(GLOBAL_Y_NEW_POINT(IJK)),REAL(GLOBAL_Z_NEW_POINT(IJK))/)
               POINT_ID = POINT_ID + 1
            ENDIF
         END DO


! Step 3: Shift connectivity with new point indices
         DO IJK = 1,IJKMAX3
            IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
               DO L=1,GLOBAL_NUMBER_OF_NODES(IJK)
                  IF(GLOBAL_CONNECTIVITY(IJK,L)<=IJKMAX3) THEN
                     GLOBAL_CLEANED_CONNECTIVITY(IJK,L) = POINT_NEW_ID(GLOBAL_CONNECTIVITY(IJK,L))
                  ELSE
                     GLOBAL_CLEANED_CONNECTIVITY(IJK,L) = NEW_POINT_NEW_ID(GLOBAL_CONNECTIVITY(IJK,L)-IJKMAX3)
                  ENDIF
               ENDDO
            ENDIF
         END DO



       ELSEIF(BDIST_IO) THEN


          IF(ALLOCATED(CLEANED_CONNECTIVITY))  DEALLOCATE (CLEANED_CONNECTIVITY)
          IF(ALLOCATED(KEEP_NEW_POINT))        DEALLOCATE (KEEP_NEW_POINT)
          IF(ALLOCATED(POINT_NEW_ID))          DEALLOCATE (POINT_NEW_ID)
          IF(ALLOCATED(NEW_POINT_NEW_ID))      DEALLOCATE (NEW_POINT_NEW_ID)
          IF(ALLOCATED(KEEP_POINT))            DEALLOCATE (KEEP_POINT)

          ALLOCATE (CLEANED_CONNECTIVITY(IJKEND3,15))
          ALLOCATE (KEEP_NEW_POINT(NUMBER_OF_NEW_POINTS))

          ALLOCATE (POINT_NEW_ID(IJKEND3))
          ALLOCATE (NEW_POINT_NEW_ID(IJKEND3))
          ALLOCATE (KEEP_POINT(IJKEND3))

! Step 1: Go through connectivity list and only keep points that are used.
!         For background cell corners, assign KEEP_POINT = .TRUE.
!         For cut cells, the new intersection points were called NEW_POINTS,
!         so assign KEEP_NEW_POINT = .TRUE.
!         A NEW_POINT had an IJK index larger than IJKMAX3

          KEEP_POINT = .FALSE.
          KEEP_NEW_POINT = .FALSE.

          DO IJK = 1,IJKEND3
             IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                DO L=1,NUMBER_OF_NODES(IJK)
                   IJKC = CONNECTIVITY(IJK,L)
                   IF(IJKC<=IJKEND3) KEEP_POINT(IJKC) = .TRUE.
                   IF(IJKC>IJKEND3) KEEP_NEW_POINT(IJKC-IJKEND3) = .TRUE.
                ENDDO
             ENDIF
          END DO


! Step 2: Clean-up list of used points and store cleaned connectivity
          POINT_NEW_ID = -1
          NEW_POINT_NEW_ID = -1
          POINT_ID = 1
! This is for the background grid cell corners
          DO IJK = 1,IJKEND3
             IF(KEEP_POINT(IJK)) THEN
                POINT_NEW_ID(IJK) = POINT_ID
                POINT_ID = POINT_ID + 1
             ENDIF
          END DO
! This is for the cut cell new corners
          DO IJK = 1,NUMBER_OF_NEW_POINTS
             IF(KEEP_NEW_POINT(IJK)) THEN
                NEW_POINT_NEW_ID(IJK) = POINT_ID
                POINT_ID = POINT_ID + 1
             ENDIF
          END DO

! Update the true (clean) number of points
          NUMBER_OF_POINTS = POINT_ID - 1

! Now, store a list of coordinates for all used points
          IF(ALLOCATED(COORDS_OF_POINTS)) DEALLOCATE(COORDS_OF_POINTS)

          ALLOCATE(COORDS_OF_POINTS(NUMBER_OF_POINTS,3))

          POINT_ID = 1
! This is for the background grid cell corners
          DO IJK = 1,IJKEND3
             IF(KEEP_POINT(IJK)) THEN
                COORDS_OF_POINTS(POINT_ID,1:3) = &
                     (/REAL(XG_E(I_OF(IJK))),REAL(YG_N(J_OF(IJK))),REAL(ZG_T(K_OF(IJK)))/)
                POINT_ID = POINT_ID + 1
             ENDIF
          END DO
! This is for the cut cell new corners
          DO IJK = 1,NUMBER_OF_NEW_POINTS
             IF(KEEP_NEW_POINT(IJK)) THEN
                NEW_POINT_NEW_ID(IJK) = POINT_ID
                COORDS_OF_POINTS(POINT_ID,1:3) = &
                     (/REAL(X_NEW_POINT(IJK)),REAL(Y_NEW_POINT(IJK)),REAL(Z_NEW_POINT(IJK))/)
                POINT_ID = POINT_ID + 1
             ENDIF
          END DO


! Step 3: Shift connectivity with new point indices
          DO IJK = 1,IJKEND3
             IF (BELONGS_TO_VTK_SUBDOMAIN(IJK)) THEN
                DO L=1,NUMBER_OF_NODES(IJK)
                   IF(CONNECTIVITY(IJK,L)<=IJKEND3) THEN
                      CLEANED_CONNECTIVITY(IJK,L) = POINT_NEW_ID(CONNECTIVITY(IJK,L))
                   ELSE
                      CLEANED_CONNECTIVITY(IJK,L) = NEW_POINT_NEW_ID(CONNECTIVITY(IJK,L)-IJKEND3)
                   ENDIF
                ENDDO
             ENDIF
          END DO

       ENDIF

      RETURN

      END SUBROUTINE CLEAN_GEOMETRY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SETUP_VTK_REGION                                       C
!                                                                      C
!  Purpose: Filter the cells based on the VTK region bounds and        C
!           set the flag BELONGS_TO_VTK_SUBDOMAIN(IJK) to .TRUE.       C
!           to keep the cell.                                          C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 19-Dec-14  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SETUP_VTK_REGION

      USE cdist
      USE compar, only: mype, pe_io, ijkend3
      USE cutcell
      USE geometry
      USE indices, only: i_of, j_of, k_of
      USE vtk

      IMPLICIT NONE

      INTEGER :: IJK,I,J,K,I_E,I_W,J_N,J_S,K_T,K_B
      INTEGER :: NXS,NYS,NZS,NS,I_TMP,J_TMP,K_TMP
      INTEGER :: I_SLICE(DIM_I),J_SLICE(DIM_J),K_SLICE(DIM_K)
      DOUBLE PRECISION :: XE,XW,YS,YN,ZB,ZT
      DOUBLE PRECISION :: XSLICE,YSLICE,ZSLICE
      LOGICAL :: KEEP_XDIR,KEEP_YDIR,KEEP_ZDIR

! Get VTK region bounds
      XE = VTK_X_E(VTK_REGION)
      XW = VTK_X_W(VTK_REGION)
      YS = VTK_Y_S(VTK_REGION)
      YN = VTK_Y_N(VTK_REGION)
      ZB = VTK_Z_B(VTK_REGION)
      ZT = VTK_Z_T(VTK_REGION)

      NXS = VTK_NXS(VTK_REGION)
      NYS = VTK_NYS(VTK_REGION)
      NZS = VTK_NZS(VTK_REGION)

      CALL CALC_CELL (XMIN, VTK_X_W(VTK_REGION), DX, IMAX, I_W)
      I_W = I_W !+ 1
      CALL CALC_CELL (XMIN, VTK_X_E(VTK_REGION), DX, IMAX, I_E)


      CALL CALC_CELL (ZERO, VTK_Y_S(VTK_REGION), DY, JMAX, J_S)
      J_S = J_S !+ 1
      CALL CALC_CELL (ZERO, VTK_Y_N(VTK_REGION), DY, JMAX, J_N)

      IF (NO_K) THEN
         K_B = 1
         K_T = 1
      ELSE
         CALL CALC_CELL (ZERO, VTK_Z_B(VTK_REGION), DZ, KMAX, K_B)
         K_B = K_B !+ 1
         CALL CALC_CELL (ZERO, VTK_Z_T(VTK_REGION), DZ, KMAX, K_T)
      ENDIF

! get slice(s) location
      DO NS = 1,NXS
         XSLICE = XW + (XE-XW)/(NXS-1)*(NS-1)
         CALL CALC_CELL (XMIN, XSLICE, DX, IMAX, I_TMP)
         I_SLICE(NS) = MAX(MIN(I_TMP,IMAX1),IMIN1)
      ENDDO

      DO NS = 1,NYS
         YSLICE = YS + (YN-YS)/(NYS-1)*(NS-1)
         CALL CALC_CELL (ZERO, YSLICE, DY, JMAX, J_TMP)
         J_SLICE(NS) = MAX(MIN(J_TMP,JMAX1),JMIN1)
      ENDDO

      DO NS = 1,NZS
         ZSLICE = ZB + (ZT-ZB)/(NZS-1)*(NS-1)
         CALL CALC_CELL (ZERO, ZSLICE, DZ, KMAX, K_TMP)
         K_SLICE(NS) = MAX(MIN(K_TMP,KMAX1),KMIN1)
      ENDDO

      IF (myPE == PE_IO.AND.(.NOT.BDIST_IO)) THEN

         IF(ALLOCATED(BELONGS_TO_VTK_SUBDOMAIN)) DEALLOCATE(BELONGS_TO_VTK_SUBDOMAIN)

         ALLOCATE (BELONGS_TO_VTK_SUBDOMAIN(ijkmax3))

! Filter the cells based on the VTK region bounds and set the
! flag BELONGS_TO_VTK_SUBDOMAIN(IJK) to .TRUE. to keep the cell.

         BELONGS_TO_VTK_SUBDOMAIN = .FALSE.
         NUMBER_OF_VTK_CELLS      = 0

         DO IJK = 1,IJKMAX3
            IF (GLOBAL_INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.GLOBAL_BLOCKED_CELL_AT(IJK)) THEN
                  I = GLOBAL_I_OF(IJK)
                  J = GLOBAL_J_OF(IJK)
                  K = GLOBAL_K_OF(IJK)

                  IF(VTK_CUTCELL_ONLY(VTK_REGION)) THEN
                     IF(I==IMIN1.OR.I==IMAX1.OR. &
                        J==JMIN1.OR.J==JMAX1.OR. &
                        K==KMIN1.OR.K==KMAX1.OR. &
                        GLOBAL_CUT_CELL_AT(IJK)) THEN

                        BELONGS_TO_VTK_SUBDOMAIN(IJK) = .TRUE.
                        NUMBER_OF_VTK_CELLS = NUMBER_OF_VTK_CELLS + 1
                     ENDIF
                     CYCLE
                  ENDIF


! X-direction
                  KEEP_XDIR=.FALSE.
                  IF(NXS==0) THEN
                     IF(I_W<=I.AND.I<=I_E) KEEP_XDIR=.TRUE.
                  ELSE
                     DO NS = 1,NXS
                        IF(I==I_SLICE(NS)) KEEP_XDIR=.TRUE.
                     ENDDO
                  ENDIF

! Y-direction
                  KEEP_YDIR=.FALSE.
                  IF(NYS==0) THEN
                     IF(J_S<=J.AND.J<=J_N) KEEP_YDIR=.TRUE.
                  ELSE
                     DO NS = 1,NYS
                        IF(J==J_SLICE(NS)) KEEP_YDIR=.TRUE.
                     ENDDO
                  ENDIF

! Z-direction
                  KEEP_ZDIR=.FALSE.
                  IF(NZS==0) THEN
                     IF(K_B<=K.AND.K<=K_T) KEEP_ZDIR=.TRUE.
                  ELSE
                     DO NS = 1,NZS
                        IF(K==K_SLICE(NS)) KEEP_ZDIR=.TRUE.
                     ENDDO
                  ENDIF

! Now combine
                  IF(KEEP_XDIR.AND.KEEP_YDIR.AND.KEEP_ZDIR) THEN
                     BELONGS_TO_VTK_SUBDOMAIN(IJK) = .TRUE.
                     NUMBER_OF_VTK_CELLS = NUMBER_OF_VTK_CELLS + 1
                  ENDIF
               ENDIF
            ENDIF
         END DO

      ELSE  ! BDIST_IO

         IF(ALLOCATED(BELONGS_TO_VTK_SUBDOMAIN)) DEALLOCATE(BELONGS_TO_VTK_SUBDOMAIN)

         ALLOCATE (BELONGS_TO_VTK_SUBDOMAIN(ijkend3))

! Filter the cells based on the VTK region bounds and set the
! flag BELONGS_TO_VTK_SUBDOMAIN(IJK) to .TRUE. to keep the cell.

         BELONGS_TO_VTK_SUBDOMAIN = .FALSE.
         NUMBER_OF_VTK_CELLS      = 0

         DO IJK = 1,IJKEND3
            IF (INTERIOR_CELL_AT(IJK))      THEN
               IF (.NOT.BLOCKED_CELL_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

! X-direction
                  KEEP_XDIR=.FALSE.
                  IF(NXS==0) THEN
                     IF(I_W<=I.AND.I<=I_E) KEEP_XDIR=.TRUE.
                  ELSE
                     DO NS = 1,NXS
                        IF(I==I_SLICE(NS)) KEEP_XDIR=.TRUE.
                     ENDDO
                  ENDIF

! Y-direction
                  KEEP_YDIR=.FALSE.
                  IF(NYS==0) THEN
                     IF(J_S<=J.AND.J<=J_N) KEEP_YDIR=.TRUE.
                  ELSE
                     DO NS = 1,NYS
                        IF(J==J_SLICE(NS)) KEEP_YDIR=.TRUE.
                     ENDDO
                  ENDIF

! Z-direction
                  KEEP_ZDIR=.FALSE.
                  IF(NZS==0) THEN
                     IF(K_B<=K.AND.K<=K_T) KEEP_ZDIR=.TRUE.
                  ELSE
                     DO NS = 1,NZS
                        IF(K==K_SLICE(NS)) KEEP_ZDIR=.TRUE.
                     ENDDO
                  ENDIF

! Now combine
                  IF(KEEP_XDIR.AND.KEEP_YDIR.AND.KEEP_ZDIR) THEN
                     BELONGS_TO_VTK_SUBDOMAIN(IJK) = .TRUE.
                     NUMBER_OF_VTK_CELLS = NUMBER_OF_VTK_CELLS + 1
                  ENDIF
               ENDIF
            ENDIF
         END DO

      ENDIF

      RETURN

      END SUBROUTINE SETUP_VTK_REGION
