! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Module name: DERIVED_TYPES                                         C
!   Purpose: contains derived type definitions and enum definitions    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE DERIVED_TYPES

!-----------------------------------------------
! Modules
!-----------------------------------------------
  USE multi_sweep_and_prune, only: multisap_t, boxhandlelist_t
  IMPLICIT NONE
!-----------------------------------------------
  ! the global multisap
  type(multisap_t) multisap

  type(boxhandlelist_t), DIMENSION(:),  ALLOCATABLE :: boxhandle         !(PARTICLES)

! Dynamic information related to computational (eulerian) fluid grid
!----------------------------------------------------------------->>>
! Dynamic variable. for each ijk computational fluid cell store the
! total number of particles and the id's of the particles in that cell
  TYPE iap1
     INTEGER, DIMENSION(:), POINTER:: p
  END TYPE iap1

  ! in order to facilitate the parallel processing the PIC is defined
  ! as single array IJK
  TYPE(iap1), DIMENSION(:), ALLOCATABLE:: pic  ! (DIMENSION_3)

! particle in cell related variable
  type iap2
     integer :: isize
     integer, dimension(:), pointer:: p
  end type iap2

  type(iap2), dimension(:),allocatable:: dg_pic

! Drag model options (see drag_gs for full details)
! default is syam_obrien (may enforce a corrected Umf by defining
! drag_c1 and drag_d1 accordingly)
  CHARACTER(64) :: DRAG_TYPE
  INTEGER :: DRAG_TYPE_ENUM

  ENUM, BIND(C)
     ENUMERATOR :: SYAM_OBRIEN=0
     ENUMERATOR :: GIDASPOW=1
     ENUMERATOR :: GIDASPOW_PCF=2
     ENUMERATOR :: GIDASPOW_BLEND=3
     ENUMERATOR :: GIDASPOW_BLEND_PCF=4
     ENUMERATOR :: WEN_YU=5
     ENUMERATOR :: WEN_YU_PCF=6
     ENUMERATOR :: KOCH_HILL=7
     ENUMERATOR :: KOCH_HILL_PCF=8
     ENUMERATOR :: BVK=9
     ENUMERATOR :: HYS=10
     ENUMERATOR :: USER_DRAG=11
  END ENUM

! filtered/subgrid corrections to the drag coefficient & granular
! stress terms including granular viscosity and solids pressure
! current options are 'igci' and 'milioli'
  CHARACTER(64) :: SUBGRID_TYPE

  INTEGER :: SUBGRID_TYPE_ENUM
  ENUM, BIND(C)
     ENUMERATOR :: UNDEFINED_SUBGRID_TYPE=0
     ENUMERATOR :: IGCI=1
     ENUMERATOR :: MILIOLI=2
  END ENUM

  ! Kinetic theory model options (see calc_mu_s for details)
  ! for m > 1 : IA_nonep, GHD, LUN_1984
  ! for m = 1 : LUN_1984, simonin, ahmadi, or
  !             GD_99 for granular flow or GTSH for gas-solids flow
  CHARACTER(64) :: KT_TYPE
  INTEGER :: KT_TYPE_ENUM
  ENUM, BIND(C)
     ENUMERATOR :: LUN_1984=0
     ENUMERATOR :: SIMONIN_1996=1
     ENUMERATOR :: AHMADI_1995=2
     ENUMERATOR :: GD_1999=3
     ENUMERATOR :: GTSH_2012=4
     ENUMERATOR :: IA_2005=5
     ENUMERATOR :: GHD_2007=6
  END ENUM

  ! Radial distribution function options (see g_0 for details)
  ! for m > 1 options are lebowitz, modified_lebowitz,
  ! mansoori, modified_mansoori.  default = lebowitz
  ! for m = 1 then carnahan and starling rdf used
  CHARACTER(64) :: RDF_TYPE
  INTEGER :: RDF_TYPE_ENUM
  ENUM, BIND(C)
     ENUMERATOR :: LEBOWITZ=0
     ENUMERATOR :: MODIFIED_LEBOWITZ=1
     ENUMERATOR :: MANSOORI=2
     ENUMERATOR :: MODIFIED_MANSOORI=3
     ENUMERATOR :: CARNAHAN_STARLING=4
  END ENUM

 END MODULE DERIVED_TYPES
