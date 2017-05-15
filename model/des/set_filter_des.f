!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_FILTER_DES                                          !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_FILTER_DES

! Runtime Flag: Invoke gas/solids coupled simulation.
      use geometry, only: DX, IMIN1, IMAX1
      use geometry, only: DY, JMIN1, JMAX1
      use geometry, only: DZ, KMIN1, KMAX1, DO_K

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS

      use particle_filter, only: DES_INTERP_WIDTH
      use particle_filter, only: FILTER_WIDTH_INTERP

      use particle_filter, only: OoFILTER_VOL
      use particle_filter, only: FILTER_WIDTH_INTERPx3

      use desgrid, only: DG_DXinv, DG_DYinv, DG_DZinv

      use param1, only: ONE, UNDEFINED

      use sendrecvnode, only: DES_SETNODEINDICES
      use sendrecvnode, only: INIT_DES_COLLECT_gDATA
      use mpi_utility, only: GLOBAL_ALL_MIN

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION :: DXYZ_MIN, DG_DXYZ_MIN


!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_FILTER_DES")

! Get the minimum fluid cell dimension
      DXYZ_MIN = min(minval(DX(IMIN1:IMAX1)),minval(DY(JMIN1:JMAX1)))
      IF(DO_K) DXYZ_MIN = min(DXYZ_MIN,minval(DZ(KMIN1:KMAX1)))

! Get the minimum DES grid cell dimension
      DG_DXYZ_MIN = min(ONE/DG_DXinv, ONE/DG_DYinv)
      IF(DO_K) DG_DXYZ_MIN = min(DG_DXYZ_MIN, ONE/DG_DZinv)
      CALL GLOBAL_ALL_MIN(DG_DXYZ_MIN)

! Verify that the interpolation scheme doesn't exceed the grid.
      IF(DES_INTERP_WIDTH /= UNDEFINED) THEN

         IF(0.5d0*DES_INTERP_WIDTH > DXYZ_MIN .OR.                     &
            0.5d0*DES_INTERP_WIDTH > DG_DXYZ_MIN) THEN

            WRITE(ERR_MSG,2130) DXYZ_MIN, DG_DXYZ_MIN, DES_INTERP_WIDTH
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2130 FORMAT('Error 2130: The specified DES_INTERP_WIDTH is too ',     &
         'large. The',/'interpolation half-width should not exceed ',  &
         'either the minimum fluid',/'cell dimension [DES_INTERP_',    &
         'WIDTH/2.0 <= min(DX,DY,DZ)] nor should it',/'exceed the ',   &
         'minimum DES grid cell.',2/3x,'Minimum fluid cell dimension:',&
         5x,g12.4,/3x,'Minimum DES grid cell dimension:',2x,g12.4/3x,  &
         'DES_INTERP_WIDTH:',17x,g12.4,2/,'By default, the DES grid ', &
         'dimensions are calculated as three times the',/'maximum ',   &
         'particle diameter. This can be altered by specifying',/'DES',&
         'GRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, and DESGRIDSEARH_KMAX ',&
         'in the',/'mfix.dat file.')

         ELSE
            FILTER_WIDTH_INTERP = 0.500d0*DES_INTERP_WIDTH
         ENDIF


      ENDIF


! Calculate reused quanties
      SELECT CASE(DES_INTERP_SCHEME_ENUM)

      CASE(DES_INTERP_GARG)
! Compute the volume of nodes needed in drag_fgs_des0.f
         CALL COMPUTE_VOLUME_OF_NODES
! Setup MPI exchange arrys for nodes
         CALL DES_SETNODEINDICES
      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS)
         OoFILTER_VOL = 0.25d0/(FILTER_WIDTH_INTERP**3)
         FILTER_WIDTH_INTERPx3 = FILTER_WIDTH_INTERP*3
         CALL INIT_DES_COLLECT_gDATA
      END SELECT


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_FILTER_DES
