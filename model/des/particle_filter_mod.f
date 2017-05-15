!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!   Module name: DISCRETELEMENT                                        !
!   Purpose: DES mod file                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE PARTICLE_FILTER

      IMPLICIT NONE

! Switch to decide whether to call drag_gs or to call des_drag_gs via
! drag_fgs to calculate the gas-solids drag coefficient.  if false
! then drag_gs is used, otherwise des_drag_gs via drag_fgs is used
! At a more fundamental level, when des_interp_on is true, then the
! drag on particle is obtained by interpolating the carrier flow
! velocity at the particle location. Similarly the drag force on the
! Eulerian grid is obtained by backward of interpolaiton of the above
! calculated drag force. See the DEM doc for more details.
      LOGICAL DES_INTERP_ON

! Switch to decide if the mean fields (such as solids volume fraction
! and mean solids velocity) are obtained by interpolation or by the more
! crude cell artihmetic averages. For MPPIC, this will switch will always
! be true.
      LOGICAL DES_INTERP_MEAN_FIELDS

! Flag to check if mass is conserved between discrete and continuum
! representations when mean fields are computed by backward interpolation.
! Critical for cut-cells.
      LOGICAL DES_REPORT_MASS_INTERP


      LOGICAL :: DES_DIFFUSE_MEAN_FIELDS
      DOUBLE PRECISION :: DES_DIFFUSE_WIDTH

      DOUBLE PRECISION :: DES_INTERP_WIDTH
      DOUBLE PRECISION :: FILTER_WIDTH_INTERP

      CHARACTER(len=32) :: DES_INTERP_SCHEME

      INTEGER :: DES_INTERP_SCHEME_ENUM
      INTEGER, PARAMETER :: DES_INTERP_NONE  = 0
      INTEGER, PARAMETER :: DES_INTERP_GARG  = 1
      INTEGER, PARAMETER :: DES_INTERP_DPVM  = 2
      INTEGER, PARAMETER :: DES_INTERP_GAUSS = 3
      INTEGER, PARAMETER :: DES_INTERP_LHAT  = 4


      DOUBLE PRECISION :: FILTER_WIDTH_INTERPx3, OoFILTER_VOL

      INTEGER :: FILTER_SIZE = 0
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: FILTER_CELL
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FILTER_WEIGHT


      END MODULE PARTICLE_FILTER
