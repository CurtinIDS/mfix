!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD                                            !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that P_force is evaluated as -dp/dx                   !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PG_GRAD

! Model B momentum equation
      use run, only: MODEL_B
! Particle volume.
      use discretelement, only: PVOL
! Gas pressure force by fluid cell
      use discretelement, only: P_FORCE
! Particle drag force
      use discretelement, only: DRAG_FC
! Flag for 3D simulatoins.
      use geometry, only: DO_K
! Loop bounds for fluid grid
      USE compar, only: IJKSTART3, IJKEND3
! Flags for cyclic BC with pressure drop 
      use geometry, only: CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD
! Specified pressure drop
      use bc, only: DELP_X, DELP_Y, DELP_Z
! Domain length
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH
! Gas phase pressure
      use fldvar, only: P_G

      use discretelement, only: MAX_PIP, PIJK, DES_EXPLICITLY_COUPLED
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: DES_INTERP_ON

      use functions, only: FLUID_AT

      use functions, only: IS_NORMAL

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO

      implicit none

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lPF(3)
! Loop bound for
      INTEGER :: LP_BND
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: cPG(3)
!......................................................................!

! Calculate the gas phase pressure gradient. (dP/dx)
      CALL CALC_GRAD_DES(P_G, P_FORCE)


! Add in cyclic BC pressure drop.
      cPG(1) = merge(DELP_X/XLENGTH, ZERO, CYCLIC_X_PD)
      cPG(2) = merge(DELP_Y/YLENGTH, ZERO, CYCLIC_Y_PD)
      cPG(3) = merge(DELP_Z/ZLENGTH, ZERO, CYCLIC_Z_PD)

      DO IJK=IJKSTART3, IJKEND3
         P_FORCE(:,IJK) = cPG - P_FORCE(:,IJK)
      ENDDO


      IF(DES_EXPLICITLY_COUPLED .AND. .NOT.MODEL_B) THEN

! Loop bounds for interpolation.
         LP_BND = merge(27,9,DO_K)

! Calculate the gas phase forces acting on each particle.

!$omp  parallel do default(none) &
!$omp  private(NP,lPF,lc,ijk,weight) &
!$omp  shared(MAX_PIP,DES_INTERP_ON,LP_BND,p_force,drag_fc, &
!$omp     filter_cell,filter_weight,pijk,pvol)
         DO NP=1,MAX_PIP

            IF(IS_NORMAL(NP)) THEN
               IF(.NOT.FLUID_AT(PIJK(NP,4))) CYCLE

               IF(DES_INTERP_ON) THEN
                  lPF = ZERO
                  DO LC=1,LP_BND
                     IJK = FILTER_CELL(LC,NP)
                     WEIGHT = FILTER_WEIGHT(LC,NP)
                     lPF = lPF + P_FORCE(:,IJK)*WEIGHT
                  ENDDO
               ELSE
                  lPF = P_FORCE(:,PIJK(NP,4))
               ENDIF

! Include gas pressure and gas-solids drag
               DRAG_FC(NP,:) = DRAG_FC(NP,:) + lPF*PVOL(NP)
            ENDIF

         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_PG_GRAD
