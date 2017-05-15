!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_GEOMETRY_DES                                        !
!  Author:   R.Garg                                   Date: 19-Mar-14  !
!                                                                      !
!  Purpose: Allocate des arrays that are based on Eulerian grid.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_GEOMETRY_DES


! Global Variables:
!---------------------------------------------------------------------//
! Arrays for DEM simulations delineating cell edges.
      use discretelement, only: XE, YN, ZT, DIMN
! Domain bounds (max/min).
      use discretelement, only: EX2, TY2, NZ2, WX1, BY1, SZ1
! Fluid grid cell dimensions and mesh size
      USE geometry, only: DX, IMIN2, IMAX2
      USE geometry, only: DY, JMIN2, JMAX2
      USE geometry, only: DZ, KMIN2, KMAX2
! Number of particles in the I/J/K direction
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO

! Module proceedures.
!---------------------------------------------------------------------//
      use mpi_utility
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop indices
      INTEGER :: I, J, K
! Error Flag
      INTEGER :: IER
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_GEOMETRY_DES")

      Allocate( XE (0:DIMENSION_I), STAT=IER )
      Allocate( YN (0:DIMENSION_J), STAT=IER )
      Allocate( ZT (0:DIMENSION_K), STAT=IER )

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
      CALL GLOBAL_ALL_SUM(IER)


! Set boundary edges.
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used. the code should be modified for
! consistency throughout
      EX2 = XLENGTH;    WX1 = ZERO  ! East/West
      TY2 = YLENGTH;    BY1 = ZERO  ! North/South
      NZ2 = ZLENGTH;    SZ1 = ZERO  ! Top/Bottom

! Initialize arrays.
      XE(:) = ZERO
      YN(:) = ZERO
      ZT(:) = ZERO

! Each loop starts at 2 and goes to max+2 (i.e., imin1=2, imax2=imax+2)
! However, the indices range to include ghost cells (0-imax2) to avoid
! multiple if statements in particles_in_cell
      XE(IMIN2-1) = ZERO-DX(IMIN2)
      DO I = IMIN2, IMAX2
         XE(I) = XE(I-1) + DX(I)
      ENDDO

      YN(JMIN2-1) = ZERO-DY(JMIN2)
      DO J  = JMIN2, JMAX2
         YN(J) = YN(J-1) + DY(J)
      ENDDO

      IF(DIMN.EQ.3) THEN
         ZT(KMIN2-1) = ZERO-DZ(KMIN2)
         DO K = KMIN2, KMAX2
            ZT(K) = ZT(K-1) + DZ(K)
         ENDDO
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_GEOMETRY_DES
