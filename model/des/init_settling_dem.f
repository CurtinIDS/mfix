!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#include "version.inc"

      SUBROUTINE INIT_SETTLING_DEM

      USE desgrid, ONLY: desgrid_pic
      USE derived_types, only: multisap, boxhandle
      USE discretelement
      USE error_manager
      USE mpi_funs_des, ONLY: DES_PAR_EXCHANGE
      USE run
      use functions, only: is_nonexistent
      use multi_sweep_and_prune, only: aabb_t, init_multisap, multisap_add, multisap_quicksort, multisap_sweep
      use geometry

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: FACTOR, nn

      type(aabb_t) :: aabb

      real :: mins(3), maxs(3), rad

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------


! Skip this routine if there are no particles.
      IF(PARTICLES == 0) RETURN
! Skip this routine if not a new run.
      IF(RUN_TYPE /= 'NEW') RETURN

! Skip if not coupled.
      IF(.NOT.DES_CONTINUUM_COUPLED) RETURN

! Write the initial configuration before settling
      IF(PRINT_DES_DATA .AND. NFACTOR>0) CALL WRITE_DES_DATA

      WRITE(ERR_MSG, 1100) trim(iVal(NFACTOR))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 1100 FORMAT('Beginning DEM settling period: ',A,' steps.')

! Disable the coupling flag.
      DES_CONTINUUM_COUPLED = .FALSE.

      mins(1) = 0
      mins(2) = 0
      mins(3) = 0
      maxs(1) = XLENGTH
      maxs(2) = YLENGTH
      maxs(3) = ZLENGTH

#ifdef do_sap
      rad = 100*maxval(des_radius)
      print *,"rad = ",rad
      print *,"XLENGTH = ",XLENGTH
      print *,"YLENGTH = ",YLENGTH
      print *,"ZLENGTH = ",ZLENGTH
         call init_multisap(multisap,floor(XLENGTH/rad),floor(YLENGTH/rad),floor(ZLENGTH/rad),mins,maxs)
         ! initialize SAP
         do nn=1, MAX_PIP
            if(is_nonexistent(nn)) cycle
            aabb%minendpoint(:) = DES_POS_NEW(nn,:)-DES_RADIUS(nn)
            aabb%maxendpoint(:) = DES_POS_NEW(nn,:)+DES_RADIUS(nn)

            if ( any(DES_RADIUS(nn)*multisap%one_over_cell_length(1:merge(2,3,NO_K)) > 0.5 ) ) then
               print *,"BAD RADIUS..grid too fine, need to have radius=",des_radius(nn),"  less than half cell length= ",0.5/multisap%one_over_cell_length(:)
               ERROR_STOP __LINE__
            endif

            call multisap_add(multisap,aabb,nn,boxhandle(nn))
         enddo

         call multisap_quicksort(multisap)
         call multisap_sweep(multisap)
#endif

      DO FACTOR = 1, NFACTOR
! calculate forces

         CALL CALC_FORCE_DEM
! update particle position/velocity

         CALL CFNEWVALUES
! set the flag do_nsearch before calling particle in cell (for mpi)
         DO_NSEARCH = (MOD(FACTOR,NEIGHBOR_SEARCH_N)==0)

! Bin the particles to the DES grid.
         CALL DESGRID_PIC(.TRUE.)
! exchange particle crossing boundaries and updates ghost particles
         CALL DES_PAR_EXCHANGE
! find particles on grid
         CALL PARTICLES_IN_CELL
! perform neighbor search
         IF(DO_NSEARCH) CALL NEIGHBOUR
      ENDDO

! Reset the comoupling flag.
      DES_CONTINUUM_COUPLED = .TRUE.

      WRITE(ERR_MSG, 1200)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 1200 FORMAT('DEM settling period complete.')

! this write_des_data is needed to properly show the initial state of
! the simulation (granular or coupled). In the coupled case, the
! particles may have 'settled' according to above.  In the granular
! case, the initial state won't be written until after the particles
! have moved without this call.
!      IF(PRINT_DES_DATA) CALL WRITE_DES_DATA

      RETURN
      END SUBROUTINE INIT_SETTLING_DEM
