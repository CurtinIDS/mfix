!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: stl_dbg_des                                            !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Random functions for debugging STLs with DES.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STL_DBG_DES

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: STL_DBG_DG_REPORT                                       !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Reports the total number of facets in each DES grid cell.  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_DG_REPORT

      use desgrid, only: DG_IJKSTART2, DG_IJKEND2

      use stl, only: FACETS_AT_DG

      use compar, only: myPE, numPEs

      IMPLICIT NONE

      INTEGER :: IJK, TOTAL_FACETS, LC

      CHARACTER(LEN=100) :: FN

      IF(numPEs == 1) THEN
         WRITE(FN,'("FACETS_DG_GRID.DAT")')
      ELSE
         WRITE(FN,'("FACETS_DG_GRID_",I5.5,".DAT")') myPE
      ENDIF

      OPEN(1001, file=TRIM(FN))

      DO IJK=DG_IJKSTART2, DG_IJKEND2
         TOTAL_FACETS = FACETS_AT_DG(IJK)%COUNT  
         IF(TOTAL_FACETS < 1) CYCLE
         WRITE(1001,2000) IJK, TOTAL_FACETS
         DO LC=1, TOTAL_FACETS
            WRITE(1001,'(2x,I10)') FACETS_AT_DG(IJK)%ID(LC)
         ENDDO
      ENDDO

      CLOSE(1001, STATUS = "keep")

 2000 FORMAT(2/2x,'DG CELL: ',I10,3x,'Total STLs: ',I4)

      RETURN
      END SUBROUTINE STL_DBG_DG_REPORT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: STL_DBG_WRITE_FACETS                                    !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: Write back out the STL files read from input files.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_WRITE_FACETS(STL_TYPE)

! Number of facets 
      use stl, only: N_FACETS, N_FACETS_DES
! Facet Vertices and normal
      use stl, only: VERTEX, NORM_FACE
! Processor rank and rank of IO
      use compar, only: myPE, PE_IO
! Start/End position of different STLs
      use stl, only: STL_START, STL_END
! All STLS
      use stl, only: ALL_STL
! STLs read from geometry files
      use stl, only: BASE_STL
! STLs for user specified walls (NSW, PSW, FSW)
      use stl, only: BCWALLS_STL
! STLs for impermeable surfaces
      use stl, only: IMPRMBL_STL
! STLs for default walls
      use stl, only: DEFAULT_STL

      use error_manager

      IMPLICIT NONE

! Type of STL to output
      INTEGER, INTENT(IN) :: STL_TYPE

      INTEGER :: LC, lSTART, lEND
      CHARACTER(len=128) :: FNAME

      IF(myPE /= PE_IO) RETURN

      SELECT CASE(STL_TYPE)
      CASE(BASE_STL)
         lSTART = STL_START(BASE_STL) 
         lEND=STL_END(BASE_STL)
         FNAME='BASE_FACETS.stl'
      CASE(BCWALLS_STL)
         lSTART = STL_START(BCWALLS_STL) 
         lEND=STL_END(BCWALLS_STL)
         FNAME='BCWALLS_FACETS.stl'
      CASE(IMPRMBL_STL)
         lSTART = STL_START(IMPRMBL_STL) 
         lEND=STL_END(IMPRMBL_STL)
         FNAME='IMPRMBL_FACETS.stl'
      CASE(DEFAULT_STL)
         lSTART = STL_START(DEFAULT_STL) 
         lEND=STL_END(DEFAULT_STL)
         FNAME='DEFAULT_FACETS.stl'
      CASE(ALL_STL)
         lSTART = 1
         lEND=N_FACETS_DES
         FNAME='ALL_FACETS.stl'
      END SELECT

      IF(lEND < lSTART) THEN
         WRITE(ERR_MSG,"('No FACETS to report: ',A)") trim(FNAME)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         RETURN
      ENDIF

      OPEN(UNIT=444, FILE=trim(FNAME))

      WRITE(444,*) 'solid vcg'
      DO LC = lSTART, lEND
         WRITE(444,*) '   facet normal ', NORM_FACE(:,LC)
         WRITE(444,*) '      outer loop'
         WRITE(444,*) '         vertex ', VERTEX(1,1:3,LC)
         WRITE(444,*) '         vertex ', VERTEX(2,1:3,LC)
         WRITE(444,*) '         vertex ', VERTEX(3,1:3,LC)
         WRITE(444,*) '      endloop'
         WRITE(444,*) '   endfacet'
      ENDDO
      WRITE(444,*)'endsolid vcg'

      CLOSE(555)

      RETURN
      END SUBROUTINE STL_DBG_WRITE_FACETS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_stl_from_grid_facet                         !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE STL_DBG_WRITE_STL_FROM_DG(WRITE_EACH_CELL,STL_TYPE)

! DES grid loop bounds
      use desgrid, only: DG_IJKSTART2, DG_IJKEND2
! Max numer of STLs
      use stl, only: DIM_STL
! Facets binned to DES grid
      use stl, only: FACETS_AT_DG
! Facet normal and vertex data
      use stl, only: NORM_FACE, VERTEX
! Processor ID and total proc count
      USE compar, only: myPE, numPEs
! Start/End position of different STLs
      use stl, only: STL_START, STL_END
! STLs read from geometry files
      use stl, only: BASE_STL
! STLs for user specified walls (NSW, PSW, FSW)
      use stl, only: BCWALLS_STL
! STLs for impermeable surfaces
      use stl, only: IMPRMBL_STL
! STLs for default walls
      use stl, only: DEFAULT_STL
! All STLs
      use stl, only: ALL_STL
! Total number of STLs for DES
      use stl, only: N_FACETS_DES

      IMPLICIT NONE

      LOGICAL, INTENT(IN), OPTIONAL :: WRITE_EACH_CELL
      INTEGER, INTENT(IN), OPTIONAL :: STL_TYPE

      INTEGER :: IJK, LC1, LC2
      INTEGER :: lSTART, lEND
      CHARACTER(LEN=8) :: FID, IDX
      LOGICAL :: EACH_CELL

      LOGICAL, ALLOCATABLE :: WRITE_FACET(:)

! Initialize flag.
      EACH_CELL = .FALSE.
      IF(present(WRITE_EACH_CELL)) EACH_CELL = WRITE_EACH_CELL

      ALLOCATE (WRITE_FACET(DIM_STL))
      WRITE_FACET = .TRUE.


      IF(present(STL_TYPE)) THEN
         SELECT CASE(STL_TYPE)
         CASE(BASE_STL)
            lSTART = STL_START(BASE_STL) 
            lEND=STL_END(BASE_STL)
            FID='BASE'
         CASE(BCWALLS_STL)
            LSTART = STL_START(BCWALLS_STL) 
            LEND=STL_END(BCWALLS_STL)
            FID='BCWALLS'
         CASE(IMPRMBL_STL)
            LSTART = STL_START(IMPRMBL_STL) 
            LEND=STL_END(IMPRMBL_STL)
            FID='IMPRMBL'
         CASE(DEFAULT_STL)
            LSTART = STL_START(DEFAULT_STL) 
            LEND=STL_END(DEFAULT_STL)
            FID='DEFAULT'
         CASE(ALL_STL)
            LSTART = 1
            LEND=N_FACETS_DES
            FID='ALL'
         END SELECT
      ELSE
         LSTART = 1
         LEND=N_FACETS_DES
         FID='ALL'
      ENDIF

      IF(numPEs == 1) THEN
         OPEN(UNIT=444,FILE='DG_FACETS_'//trim(FID)//&
            '.stl', STATUS='UNKNOWN')
      ELSE
         WRITE(IDX,"(I8.8)") myPE
         OPEN(UNIT=444,FILE='DG_FACETS_'//trim(FID)//&
            '_'//IDX//'.stl', STATUS='UNKNOWN')
      ENDIF

      write(444,*)'solid vcg'
      DO IJK=DG_IJKSTART2,DG_IJKEND2
         IF(FACETS_AT_DG(IJK)%COUNT< 1) CYCLE

         IF(EACH_CELL) CALL WRITE_STLS_THIS_DG(IJK, STL_TYPE)

         DO LC1 = 1, FACETS_AT_DG(IJK)%COUNT
            LC2 = FACETS_AT_DG(IJK)%ID(LC1)

            IF(LC2 < lSTART .OR. LC2 > lEND) &
               WRITE_FACET(LC2) = .FALSE.

            IF(WRITE_FACET(LC2)) THEN
               write(444,*) '   facet normal ', NORM_FACE(:,LC2)
               write(444,*) '      outer loop'
               write(444,*) '         vertex ', VERTEX(1,:,LC2)
               write(444,*) '         vertex ', VERTEX(2,:,LC2)
               write(444,*) '         vertex ', VERTEX(3,:,LC2)
               write(444,*) '      endloop'
               write(444,*) '   endfacet'
               WRITE_FACET(LC2) = .FALSE.
            ENDIF
         ENDDO
      ENDDO
      write(444,*)'endsolid vcg'

      close(444)

      DEALLOCATE (WRITE_FACET)

      RETURN
      END SUBROUTINE STL_DBG_WRITE_STL_FROM_DG




!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_STLS_THIS_DG(DG, STL_TYPE)

! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
! Facets binned to DES grid
      use stl, only: FACETS_AT_DG
! Start/End position of different STLs
      use stl, only: STL_START, STL_END
! STLs read from geometry files
      use stl, only: BASE_STL
! STLs for user specified walls (NSW, PSW, FSW)
      use stl, only: BCWALLS_STL
! STLs for impermeable surfaces
      use stl, only: IMPRMBL_STL
! STLs for default walls
      use stl, only: DEFAULT_STL
! All STLs
      use stl, only: ALL_STL
! Total number of STLs for DES
      use stl, only: N_FACETS_DES


      IMPLICIT NONE
!-----------------------------------------------
      INTEGER, INTENT(IN) :: DG
      INTEGER, INTENT(IN), OPTIONAL :: STL_TYPE

      INTEGER :: ID, FACET, lCOUNT
      INTEGER :: lSTART, lEND

      LOGICAL :: EXISTS
      CHARACTER(LEN=8) :: IDX, FID

      IF(present(STL_TYPE)) THEN
         SELECT CASE(STL_TYPE)
         CASE(BASE_STL)
            lSTART = STL_START(BASE_STL) 
            lEND=STL_END(BASE_STL)
            FID='base'
         CASE(BCWALLS_STL)
            lSTART = STL_START(BCWALLS_STL) 
            lEND=STL_END(BCWALLS_STL)
            FID='bcwalls'
         CASE(IMPRMBL_STL)
            lSTART = STL_START(IMPRMBL_STL) 
            lEND=STL_END(IMPRMBL_STL)
            FID='imprmbl'
         CASE(DEFAULT_STL)
            lSTART = STL_START(DEFAULT_STL) 
            lEND=STL_END(DEFAULT_STL)
            FID='default'
         CASE(ALL_STL)
            lSTART = 1
            lEND=N_FACETS_DES
            FID='all'
         END SELECT
      ELSE
         lSTART = 1
         lEND=N_FACETS_DES
         FID='all'
      ENDIF

      lCOUNT = 0
      DO FACET=1, FACETS_AT_DG(DG)%COUNT
         ID = FACETS_AT_DG(DG)%ID(FACET)
         IF(ID >= lSTART .AND. ID <= lEND) lCOUNT = lCOUNT+1
      ENDDO

      IF(FACETS_AT_DG(DG)%COUNT < 1) RETURN

      write(idx,"(I8.8)") dg
      open(unit=555,file='dg_'//idx//'_'//trim(FID)//&
         '.stl',status='UNKNOWN')

      write(555,*) 'solid vcg'

      DO FACET=1, FACETS_AT_DG(DG)%COUNT

         ID = FACETS_AT_DG(DG)%ID(FACET)
         IF(ID < lSTART .OR. ID > lEND) CYCLE 

         write(555,*) '   facet normal ', NORM_FACE(:,ID)
         write(555,*) '      outer loop'
         write(555,*) '         vertex ', VERTEX(1,1:3,ID)
         write(555,*) '         vertex ', VERTEX(2,1:3,ID)
         write(555,*) '         vertex ', VERTEX(3,1:3,ID)
         write(555,*) '      endloop'
         write(555,*) '   endfacet'
      ENDDO
      CLOSE(555)

      RETURN
      END SUBROUTINE write_stls_this_dg


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_this_stl(this)


! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
      use compar, only: myPE

      IMPLICIT NONE
!-----------------------------------------------
      integer, intent(in) :: this

      logical :: EXISTS
      character(len=4) :: IDX
      character(len=4) :: IPE


      write(idx,"(I4.4)") this
      write(ipe,"(I4.4)") myPE
      open(unit=555, file='idv_'//idx//'_'//IPE//'.stl',&
         status='UNKNOWN')
      write(555,*) 'solid vcg'
      write(555,*) '   facet normal ', NORM_FACE(:,this)
      write(555,*) '      outer loop'
      write(555,*) '         vertex ', VERTEX(1,1:3,this)
      write(555,*) '         vertex ', VERTEX(2,1:3,this)
      write(555,*) '         vertex ', VERTEX(3,1:3,this)
      write(555,*) '      endloop'
      write(555,*) '   endfacet'
      close(555)


      RETURN
      END SUBROUTINE write_this_stl

      END MODULE STL_DBG_DES


