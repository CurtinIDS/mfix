!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: stl_functions_des                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STL_PREPROC_DES

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_STL_PREPROCESSING                                   !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_STL_PREPROCESSING

! Flag to for STL defined geometry
      use cutcell, only: use_stl
! Number of facets from STL files, (plus DES generated)
      use stl, only: N_FACETS, N_FACETS_DES
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

      use stl_dbg_des
      use error_manager

      IMPLICIT NONE

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('Pre-Processing geometry for DES.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

! Process the STL files
      N_FACETS_DES = merge(N_FACETS, 0, USE_STL)
! Store the Start/End of the base STLs from geometry files
      STL_START(BASE_STL)=1;   STL_END(BASE_STL)=N_FACETS_DES

! Process stair-step geometries
      CALL CONVERT_BC_WALLS_TO_STL
! Process stair-step geometries
      CALL CONVERT_IMPERMEABLE_IS_TO_STL
! Process default walls
      CALL CONVERT_DEFAULT_WALLS_TO_STL

! Bin the STL to the DES grid.
      CALL BIN_FACETS_TO_DG

! Some functions for debugging.
!      CALL STL_DBG_WRITE_FACETS(BASE_STL)
!      CALL STL_DBG_WRITE_FACETS(BCWALLS_STL)
!      CALL STL_DBG_WRITE_FACETS(IMPRMBL_STL)
!      CALL STL_DBG_WRITE_FACETS(DEFAULT_STL)
!      CALL STL_DBG_WRITE_FACETS(ALL_STL)
!      CALL STL_DBG_WRITE_STL_FROM_DG(STL_TYPE=BASE_STL)

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('DES geometry pre-processing complete.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      RETURN
      END SUBROUTINE DES_STL_PREPROCESSING


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: BIN_FACETS_TO_DG                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine BIN_FACETS_TO_DG

      use desgrid, only: DG_IJKSIZE2
      use desgrid, only: DG_IEND2, DG_ISTART2
      use desgrid, only: DG_JEND2, DG_JSTART2
      use desgrid, only: DG_KEND2, DG_KSTART2
      use desgrid, only: dg_dxinv, dg_dyinv, dg_dzinv

      use stl, only: FACETS_AT_DG

      use geometry, only: XLENGTH, YLENGTH, ZLENGTH, DO_K
      use stl, only: N_FACETS_DES
      use stl, only: VERTEX

      use stl, only: TOL_STL
      use param1, only: ZERO, ONE

      use desgrid, only: DG_FUNIJK
      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_plus1layers

      IMPLICIT NONE

! DES Grid cell index.
      INTEGER :: IJK, IJK2
! Loop counters:
      INTEGER :: I1, I2, II  ! X-axis
      INTEGER :: J1, J2, JJ  ! Y-axis
      INTEGER :: K1, K2, KK  ! Z-axis
      INTEGER :: NN          ! STLs
! Generic accumulator
      INTEGER :: COUNT_FAC

! Maximum and minimum extents of the indexed STL
      DOUBLE PRECISION:: X1,Y1,Z1
      DOUBLE PRECISION:: X2,Y2,Z2

! Allocate the data storage array.
      IF(.not.allocated(FACETS_AT_DG)) &
         allocate(FACETS_AT_DG(DG_IJKSIZE2))

      FACETS_AT_DG(:)%COUNT = 0

      DO NN = 1,N_FACETS_DES

         X1 = minval(VERTEX(1:3,1,NN))
         X2 = maxval(VERTEX(1:3,1,NN))
         Y1 = minval(VERTEX(1:3,2,NN))
         Y2 = maxval(VERTEX(1:3,2,NN))
         Z1 = minval(VERTEX(1:3,3,NN))
         Z2 = maxval(VERTEX(1:3,3,NN))

         I1 = DG_IEND2
         I2 = DG_ISTART2
         IF(X2>=-TOL_STL .AND. X1<=XLENGTH+TOL_STL) THEN
            I1 = max(iofpos(X1)-1, dg_istart2)
            I2 = min(iofpos(X2)+1, dg_iend2)
         ENDIF

         J1 = DG_JEND2
         J2 = DG_JSTART2
         IF(Y2>=-TOL_STL .AND. Y1<=YLENGTH+TOL_STL) THEN
            J1 = max(jofpos(Y1)-1, dg_jstart2)
            J2 = min(jofpos(Y2)+1, dg_jend2)
         ENDIF

         K1 = DG_KEND2
         K2 = DG_KSTART2
         IF(DO_K) THEN
            IF(Z2>=-TOL_STL .AND. Z1<=ZLENGTH+TOL_STL) THEN
               K1 = max(kofpos(Z1)-1, dg_kstart2)
               K2 = min(kofpos(Z2)+1, dg_kend2)
            ENDIF
         ENDIF

         DO KK=K1,K2
         DO JJ=J1,J2
         DO II=I1,I2
            IF(dg_is_ON_myPE_plus1layers(II,JJ,KK)) THEN
               IJK = DG_FUNIJK(II,JJ,KK)
               CALL ADD_FACET_FOR_DES(II,JJ,KK,IJK,NN)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE BIN_FACETS_TO_DG




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET_FOR_DES                                       !
!  Author: Rahul Garg                                  Date: 24-Oct-13 !
!                                                                      !
!  Purpose: Add facets to DES grid cells.                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET_FOR_DES(I,J,K,IJK,N)

      use geometry, only: DO_K

      use desgrid, only: dg_dxinv, dg_xstart, dg_istart1
      use desgrid, only: dg_dyinv, dg_ystart, dg_jstart1
      use desgrid, only: dg_dzinv, dg_zstart, dg_kstart1

      use discretelement, only: MAX_RADIUS

      use stl, only: VERTEX

      use stl_functions_des, only: TRI_BOX_OVERLAP

      use param1, only: ZERO, HALF, ONE
      use error_manager

      IMPLICIT NONE

! DES grid index and facet index
      INTEGER, INTENT(IN) :: I,J,K,IJK, N

! Center of DES grid cell and half size. Note that a buffer is added to
! the half size to make the cell appear a little larger. This ensures
! that paricles near the edge 'see' STLs that are nearby but do not
! directly intersect the DES grid cell contain the particle center.
      DOUBLE PRECISION :: CENTER(3), HALFSIZE(3)
! Flag: STL intersects the DES grid cell
      LOGICAL :: OVERLAP
! DES grid cell dimensions
      DOUBLE PRECISION :: lDX, lDY, lDZ
! Buffer to ensure all particle-STL collisions are captured.
      DOUBLE PRECISION :: BUFFER
! Legacy variable - should be removed
      INTEGER ::  CURRENT_COUNT

      BUFFER = 1.1d0*MAX_RADIUS

      lDX = ONE/DG_DXINV
      lDY = ONE/DG_DYINV
      lDZ = ONE/DG_DZINV

      CENTER(1) = dg_xstart + (dble(I-dg_istart1)+HALF)*lDX
      HALFSIZE(1) = HALF*lDX + BUFFER

      CENTER(2) = dg_ystart + (dble(J-dg_jstart1)+HALF)*lDY
      HALFSIZE(2) = HALF*lDY + BUFFER

      IF(DO_K)THEN
         CENTER(3) = dg_zstart + (dble(K-dg_kstart1)+HALF)*lDZ
         HALFSIZE(3) = HALF*lDZ + BUFFER
      ELSE
         CENTER(3) = HALF*lDZ
         HALFSIZE(3) = HALF*lDZ
      ENDIF

      CALL TRI_BOX_OVERLAP(CENTER, HALFSIZE, VERTEX(:,:,N), OVERLAP)

      IF(OVERLAP) CALL ADD_FACET(IJK, N)

      RETURN
      END SUBROUTINE ADD_FACET_FOR_DES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET                                               !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET(IJK, FACET_ID)

      use stl, only: VERTEX
      use stl, only: FACETS_AT_DG
      use param1, only: ZERO

      implicit none

      INTEGER, INTENT(IN) :: IJK, facet_id

      INTEGER, ALLOCATABLE :: int_tmp(:)
      DOUBLE PRECISION, ALLOCATABLE :: real_tmp(:)

      INTEGER :: lSIZE, II
      DOUBLE PRECISION :: smallest_extent, min_temp, max_temp


      IF(FACETS_AT_DG(IJK)%COUNT > 0) THEN

         DO II=1, FACETS_AT_DG(IJK)%COUNT
            IF(FACET_ID == FACETS_AT_DG(IJK)%ID(II)) RETURN
         ENDDO

         FACETS_AT_DG(IJK)%COUNT = FACETS_AT_DG(IJK)%COUNT+1

         lSIZE = size(FACETS_AT_DG(IJK)%ID)
         IF(FACETS_AT_DG(IJK)%COUNT +1> lSIZE) THEN
            allocate(int_tmp(2*lSIZE)); int_tmp=0
            int_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%ID(1:lSIZE)
            call move_alloc(int_tmp,FACETS_AT_DG(IJK)%ID)

            allocate(int_tmp(2*lSIZE)); int_tmp=0
            int_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%DIR(1:lSIZE)
            call move_alloc(int_tmp, FACETS_AT_DG(IJK)%DIR)

            allocate(real_tmp(2*lSIZE)); real_tmp=ZERO
            real_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%MIN(1:lSIZE)
            call move_alloc(real_tmp, FACETS_AT_DG(IJK)%MIN)

            allocate(real_tmp(2*lSIZE)); real_tmp=ZERO
            real_tmp(1:lSIZE) = FACETS_AT_DG(IJK)%MAX(1:lSIZE)
            call move_alloc(real_tmp, FACETS_AT_DG(IJK)%MAX)
         ENDIF

      ELSE
         FACETS_AT_DG(IJK)%COUNT = 1
         IF(.not.allocated(FACETS_AT_DG(IJK)%ID)) &
            allocate(FACETS_AT_DG(IJK)%ID(4))
         IF(.not.allocated(FACETS_AT_DG(IJK)%DIR)) &
            allocate(FACETS_AT_DG(IJK)%DIR(4))
         IF(.not.allocated(FACETS_AT_DG(IJK)%MIN)) &
            allocate(FACETS_AT_DG(IJK)%MIN(4))
         IF(.not.allocated(FACETS_AT_DG(IJK)%MAX)) &
            allocate(FACETS_AT_DG(IJK)%MAX(4))
      ENDIF

      FACETS_AT_DG(IJK)%ID(FACETS_AT_DG(IJK)%COUNT) = FACET_ID

      SMALLEST_EXTENT = HUGE(0.0)

      DO II=1,3
         MIN_TEMP = MINVAL(VERTEX(:,II,FACET_ID))
         MAX_TEMP = MAXVAL(VERTEX(:,II,FACET_ID))
         IF(ABS(MAX_TEMP - MIN_TEMP) < SMALLEST_EXTENT ) THEN
            FACETS_AT_DG(IJK)%DIR(FACETS_AT_DG(IJK)%COUNT) = II
            FACETS_AT_DG(IJK)%MIN(FACETS_AT_DG(IJK)%COUNT) = MIN_TEMP
            FACETS_AT_DG(IJK)%MAX(FACETS_AT_DG(IJK)%COUNT) = MAX_TEMP
            SMALLEST_EXTENT = ABS(MAX_TEMP - MIN_TEMP)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE ADD_FACET



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONVERT_BC_WALLS_TO_STL                                 !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Convert user specified walls to STLs for particle-wall     !
!  collision detection.                                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine CONVERT_BC_WALLS_TO_STL

      use geometry, only: ZLENGTH, DO_K

      use bc, only: BC_DEFINED, BC_TYPE_ENUM, FREE_SLIP_WALL, NO_SLIP_WALL, PAR_SLIP_WALL
      use bc, only: BC_I_w, BC_I_e
      use bc, only: BC_J_s, BC_J_n
      use bc, only: BC_K_b, BC_K_t

      use stl, only: N_FACETS_DES
      use stl, only: STL_START, STL_END, BCWALLS_STL

      use discretelement, only: XE, YN, ZT

      use param, only: DIMENSION_BC
      USE param1, only: ZERO

      IMPLICIT NONE

! Loop counter.
      INTEGER :: BCV
! Extents of the BC region with respect to the fluid grid.
      DOUBLE PRECISION :: lXw, lXe, lYs, lYn, lZb, lZt

      STL_START(BCWALLS_STL)=N_FACETS_DES+1

      DO BCV=1, DIMENSION_BC
         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(BC_TYPE_ENUM(BCV) == FREE_SLIP_WALL .OR.   &
            BC_TYPE_ENUM(BCV) == NO_SLIP_WALL   .OR.   &
            BC_TYPE_ENUM(BCV) == PAR_SLIP_WALL) THEN

            lXw = XE(BC_I_w(BCV)-1); lXe = XE(BC_I_e(BCV))
            lYs = YN(BC_J_s(BCV)-1); lYn = YN(BC_J_n(BCV))
            IF(DO_K) THEN
               lZb = ZT(BC_K_b(BCV)-1); lZt = ZT(BC_K_t(BCV))
            ELSE
               lZb = ZERO; lZt = ZLENGTH
            ENDIF
            CALL GENERATE_STL_BOX(lXw, lXe, lYs, lYn, lZb, lZt)
         ENDIF
      ENDDO
      STL_END(BCWALLS_STL)=N_FACETS_DES

      RETURN
      END SUBROUTINE CONVERT_BC_WALLS_TO_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONVERT_IMPERMEABLE_IS_TO_STL                           !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Convert user specified impermeable surfaces to STLs for    !
!  particle-wall collision detection.                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CONVERT_IMPERMEABLE_IS_TO_STL

      use geometry, only: DO_K, ZLENGTH

      use is, only: IS_DEFINED, IS_TYPE
      use is, only: IS_I_w, IS_I_e
      use is, only: IS_J_s, IS_J_n
      use is, only: IS_K_b, IS_K_t

      use stl, only: N_FACETS_DES
      use stl, only: STL_START, STL_END, IMPRMBL_STL
      use discretelement, only: XE, YN, ZT

      use param, only: DIMENSION_IS
      USE param1, only: ZERO

      use error_manager

      IMPLICIT NONE

! Loop counter.
      INTEGER :: ISV
! Extents of the BC region with respect to the fluid grid.
      DOUBLE PRECISION :: lXw, lXe, lYs, lYn, lZb, lZt

      STL_START(IMPRMBL_STL)=N_FACETS_DES+1

      DO ISV=1, DIMENSION_IS
         IF(.NOT.IS_DEFINED(ISV)) CYCLE

         IF(trim(IS_TYPE(ISV)) == 'IMPERMEABLE') THEN

            lXw = XE(IS_I_w(ISV)-1); lXe = XE(IS_I_e(ISV))
            lYs = YN(IS_J_s(ISV)-1); lYn = YN(IS_J_n(ISV))
            IF(DO_K) THEN
               lZb = ZT(IS_K_b(ISV)-1); lZt = ZT(IS_K_t(ISV))
            ELSE
               lZb = ZERO; lZt = ZLENGTH
            ENDIF

            CALL GENERATE_STL_BOX(lXw, lXe, lYs, lYn, lZb, lZt)
         ELSE
            CALL INIT_ERR_MSG('CONVERT_IMPERMEABLE_IS_TO_STL')
            WRITE(ERR_MSG,1000) ISV, IS_TYPE(ISV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      STL_END(IMPRMBL_STL)=N_FACETS_DES

 1000 FORMAT("Error 1000: DES simulations do not support the ",/       &
         'specified IS TYPE:',/3x,'IS: ',I3,/3x,'IS_TYPE=',A)

      RETURN
      END SUBROUTINE CONVERT_IMPERMEABLE_IS_TO_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONVERT_DEFAULT_WALLS_TO_STL                            !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Convert user specified walls to STLs for particle-wall     !
!  collision detection.                                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine CONVERT_DEFAULT_WALLS_TO_STL

      USE geometry, only: DO_K
      USE geometry, only: XLENGTH, YLENGTH, ZLENGTH
      use stl, only: VERTEX, NORM_FACE
      use stl, only: N_FACETS_DES
      use stl, only: STL_START, STL_END, DEFAULT_STL

      use discretelement, only: DES_PERIODIC_WALLS_X
      use discretelement, only: DES_PERIODIC_WALLS_Y
      use discretelement, only: DES_PERIODIC_WALLS_Z

      USE param1, only: ZERO, ONE

      IMPLICIT NONE

      STL_START(DEFAULT_STL)=N_FACETS_DES+1

! West Face
      IF(.NOT.DES_PERIODIC_WALLS_X)THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/ZERO, ZERO, ZERO/)
         VERTEX(2,:,N_FACETS_DES) = (/ZERO, 2*YLENGTH, ZERO/)
         VERTEX(3,:,N_FACETS_DES) = (/ZERO, ZERO, 2*ZLENGTH/)
         NORM_FACE(:,N_FACETS_DES) = (/ONE, ZERO, ZERO/)

! East Face
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/XLENGTH, ZERO, ZERO/)
         VERTEX(2,:,N_FACETS_DES) = (/XLENGTH, 2*YLENGTH, ZERO/)
         VERTEX(3,:,N_FACETS_DES) = (/XLENGTH, ZERO, 2*ZLENGTH/)
         NORM_FACE(:,N_FACETS_DES) = (/-ONE, ZERO, ZERO/)
      ENDIF

! South Face
      IF(.NOT.DES_PERIODIC_WALLS_Y)THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/ZERO, ZERO, ZERO/)
         VERTEX(2,:,N_FACETS_DES) = (/2*XLENGTH, ZERO, ZERO/)
         VERTEX(3,:,N_FACETS_DES) = (/ZERO, ZERO, 2*ZLENGTH/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ONE, ZERO/)

! North Face
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/ZERO, YLENGTH, ZERO/)
         VERTEX(2,:,N_FACETS_DES) = (/2*XLENGTH, YLENGTH, ZERO/)
         VERTEX(3,:,N_FACETS_DES) = (/ZERO, YLENGTH, 2*ZLENGTH/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, -ONE, ZERO/)
      ENDIF

! Bottom Face
      IF(.NOT.DES_PERIODIC_WALLS_Z .AND. DO_K) THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/ZERO, ZERO, ZERO/)
         VERTEX(2,:,N_FACETS_DES) = (/2*XLENGTH, ZERO, ZERO/)
         VERTEX(3,:,N_FACETS_DES) = (/ZERO, 2*YLENGTH, ZERO/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ZERO, ONE/)

! Top Face
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/ZERO, ZERO, ZLENGTH/)
         VERTEX(2,:,N_FACETS_DES) = (/2*XLENGTH, ZERO, ZLENGTH/)
         VERTEX(3,:,N_FACETS_DES) = (/ZERO, 2*YLENGTH, ZLENGTH/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ZERO, -ONE/)
      ENDIF

      STL_END(DEFAULT_STL)=N_FACETS_DES

      RETURN
      END SUBROUTINE CONVERT_DEFAULT_WALLS_TO_STL

      END MODULE STL_PREPROC_DES




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_STL_BOX                                        !
!  Author: J.Musser                                   Date: 03-Nov-15  !
!                                                                      !
!  Purpose: Given the six corners of a box, create the 12 STLs needed  !
!  to define the geometry.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_STL_BOX(pXw, pXe, pYs, pYn, pZb, pZt)

      use stl, only: VERTEX, NORM_FACE
      use stl, only: N_FACETS_DES

      use geometry, only: DO_K

      use param1, only: ZERO, ONE

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: pXw, pXe, pYs, pYn, pZb, pZt

! West Face
      N_FACETS_DES = N_FACETS_DES+1
      VERTEX(1,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
      VERTEX(2,:,N_FACETS_DES) = (/pXw, pYn, pZb/)
      VERTEX(3,:,N_FACETS_DES) = (/pXw, pYn, pZt/)
      NORM_FACE(:,N_FACETS_DES) = (/-ONE, ZERO, ZERO/)

      IF(DO_K)THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/pXw, pYs, pZt/)
         VERTEX(2,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
         VERTEX(3,:,N_FACETS_DES) = (/pXw, pYn, pZt/)
         NORM_FACE(:,N_FACETS_DES) = (/-ONE, ZERO, ZERO/)
      ENDIF

! East Face
      N_FACETS_DES = N_FACETS_DES+1
      VERTEX(3,:,N_FACETS_DES) = (/pXe, pYs, pZb/)
      VERTEX(2,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
      VERTEX(1,:,N_FACETS_DES) = (/pXe, pYn, pZt/)
      NORM_FACE(:,N_FACETS_DES) = (/ONE, ZERO, ZERO/)

      IF(DO_K) THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(3,:,N_FACETS_DES) = (/pXe, pYs, pZt/)
         VERTEX(2,:,N_FACETS_DES) = (/pXe, pYs, pZb/)
         VERTEX(1,:,N_FACETS_DES) = (/pXe, pYn, pZt/)
         NORM_FACE(:,N_FACETS_DES) = (/ONE, ZERO, ZERO/)
      ENDIF

! South Face
      N_FACETS_DES = N_FACETS_DES+1
      VERTEX(1,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
      VERTEX(2,:,N_FACETS_DES) = (/pXe, pYs, pZb/)
      VERTEX(3,:,N_FACETS_DES) = (/pXw, pYs, pZt/)
      NORM_FACE(:,N_FACETS_DES) = (/ZERO, -ONE, ZERO/)

      IF(DO_K) THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/pXe, pYs, pZt/)
         VERTEX(2,:,N_FACETS_DES) = (/pXe, pYs, pZb/)
         VERTEX(3,:,N_FACETS_DES) = (/pXw, pYs, pZt/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, -ONE, ZERO/)
      ENDIF

! North Face
      N_FACETS_DES = N_FACETS_DES+1
      VERTEX(3,:,N_FACETS_DES) = (/pXw, pYn, pZb/)
      VERTEX(2,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
      VERTEX(1,:,N_FACETS_DES) = (/pXw, pYn, pZt/)
      NORM_FACE(:,N_FACETS_DES) = (/ZERO, ONE, ZERO/)

      IF(DO_K) THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(3,:,N_FACETS_DES) = (/pXe, pYn, pZt/)
         VERTEX(2,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
         VERTEX(1,:,N_FACETS_DES) = (/pXw, pYn, pZt/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ONE, ZERO/)
      ENDIF

! Bottom Face
      IF(DO_K)THEN
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
         VERTEX(2,:,N_FACETS_DES) = (/pXe, pYs, pZb/)
         VERTEX(3,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ZERO, -ONE/)

         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(1,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
         VERTEX(2,:,N_FACETS_DES) = (/pXw, pYn, pZb/)
         VERTEX(3,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ZERO, -ONE/)

! Top Face
         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(3,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
         VERTEX(2,:,N_FACETS_DES) = (/pXe, pYs, pZb/)
         VERTEX(1,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ZERO, ONE/)

         N_FACETS_DES = N_FACETS_DES+1
         VERTEX(3,:,N_FACETS_DES) = (/pXe, pYn, pZb/)
         VERTEX(2,:,N_FACETS_DES) = (/pXw, pYn, pZb/)
         VERTEX(1,:,N_FACETS_DES) = (/pXw, pYs, pZb/)
         NORM_FACE(:,N_FACETS_DES) = (/ZERO, ZERO, ONE/)
      ENDIF

      RETURN
      END SUBROUTINE GENERATE_STL_BOX
