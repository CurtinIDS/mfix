!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLUIDBED_P                                          C
!  Purpose: Set the pressure field inside the bed assuming a fluidized C
!           bed with gravity acting the -ve y-direction                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for including cylindrical geometry           C
!  Author: M. Syamlal                                 Date:  6-MAR-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Set pressure drop for cyclic boundary condition w/         C
!           pressure drop                                              C
!  Author: M. Syamlal                                 Date: 29-APR-94  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, IC_P_g, BC_P_g           C
!                        EP_g, MW_MIX_G, RO_g0, T_g,                   C
!                        SMAX, ROP_s,                                  C
!                        DX, DY, DZ, BFY_G, DELP_X, DELP_Y, DELP_Z,    C
!                        DO_I, DO_J, DO_K, IMIN1, KMIN1, JMIN1, IMAX1, C
!                        IMAX2, JMAX1, JMAX2, KMAX1, KMAX2             C
!  Variables modified: P_g                                             C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLUIDBED_P

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE bodyforce
      USE compar
      USE constant
      USE discretelement
      USE eos, ONLY: EOSG
      USE exit, only: mfix_exit
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE ic
      USE indices
      USE machine, only: start_log, end_log
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE scales
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, IJK, M
! Local loop counter
      INTEGER :: L
! Gas pressure at the axial location j
      DOUBLE PRECISION :: PJ
! Bed weight per unit area
      DOUBLE PRECISION :: BED_WEIGHT
! Total area of a x-z plane
      DOUBLE PRECISION :: AREA
! x-z plane area of one cell
      DOUBLE PRECISION :: dAREA
! Average pressure drop per unit length
      DOUBLE PRECISION :: DPoDX, DPoDY, DPoDZ
!-----------------------------------------------

! If any initial pressures are unspecified skip next section
! calculations.
      DO L = 1, DIMENSION_IC
         IF (IC_DEFINED(L)) THEN
            IF (IC_P_G(L) == UNDEFINED) GOTO 60
            PJ = IC_P_G(L)
         ENDIF
      ENDDO

! Here the pressure in each cell is determined from a specified pressure
! drop across the domain length. This section requires that the pressure
! is already defined in all initial condition regions (otherwise this
! section would be skipped)
! ---------------------------------------------------------------->>>
      IF (DO_I .AND. DELP_X/=UNDEFINED) THEN
         DPODX = DELP_X/XLENGTH
         PJ = PJ - DPODX*HALF*(DX(IMAX1)+DX(IMAX2))
         DO I = IMAX1, IMIN1, -1
            PJ = PJ + DPODX*HALF*(DX(I)+DX(I+1))
            DO K = KMIN1, KMAX1
               DO J = JMIN1, JMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (DO_J .AND. DELP_Y/=UNDEFINED) THEN
         DPODY = DELP_Y/YLENGTH
         PJ = PJ - DPODY*HALF*(DY(JMAX1)+DY(JMAX2))
         DO J = JMAX1, JMIN1, -1
            PJ = PJ + DPODY*HALF*(DY(J)+DY(J+1))
            DO K = KMIN1, KMAX1
               DO I = IMIN1, IMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (DO_K .AND. DELP_Z/=UNDEFINED) THEN
         DPODZ = DELP_Z/ZLENGTH
         PJ = PJ - DPODZ*HALF*(DZ(KMAX1)+DZ(KMAX2))
         DO K = KMAX1, KMIN1, -1
            PJ = PJ + DPODZ*HALF*(DZ(K)+DZ(K+1))
            DO J = JMIN1, JMAX1
               DO I = IMIN1, IMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
! ----------------------------------------------------------------<<<
      GOTO 100   ! pressure in all intial condition region cells was defined

   60 CONTINUE   ! pressure in an initial condition region cell was undefined


! ---------------------------------------------------------------->>>
! Search for an outflow boundary condition where pressure is specified
      PJ = UNDEFINED
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L) .AND. BC_TYPE_ENUM(L)==P_OUTFLOW) PJ = BC_P_G(L)
      ENDDO

      IF (PJ == UNDEFINED) THEN
! either a PO was not specified and/or a PO was specified but not the
! pressure at the outlet
         IF (RO_G0 /= UNDEFINED) THEN
! If incompressible flow set P_g to zero
            DO IJK = IJKSTART3, IJKEND3
               IF (FLUID_AT(IJK)) P_G(IJK) = ZERO
            ENDDO
            GOTO 100

         ELSE   ! compressible case

! Error condition -- no pressure outflow boundary condition is specified
! if a case is compressible and pressure in any of the initial
! conditions regions is unspecified, then a PO is effectively required
! (i.e., is specifies a bc_p_g).
            CALL START_LOG
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified
      DO J = JMAX2, JMIN1, -1

! Find the average weight per unit area over an x-z slice
         BED_WEIGHT = 0.0
         AREA = 0.0
         DO K = KMIN1, KMAX1
            DO I = IMIN1, IMAX1
               IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF (FLUID_AT(IJK)) THEN
                  IF (COORDINATES == 'CARTESIAN') THEN
                     DAREA = DX(I)*DZ(K)
                  ELSE IF (CYLINDRICAL) THEN
                     DAREA = DX(I)*X(I)*DZ(K)
                  ENDIF
                  AREA = AREA + DAREA
                  IF (RO_G0 == UNDEFINED) THEN
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_G(IJK)*EP_G(IJK)*EOSG(&
                        MW_MIX_G(IJK),PJ,T_G(IJK))*DAREA
                  ELSE
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_G(IJK)*EP_G(IJK)*RO_G0&
                        *DAREA
                  ENDIF
! This code is turned off for DEM runs until the value of rop_s can be
! ensured valid values for a DEM run at this point in the code.
                  IF (.NOT.DISCRETE_ELEMENT) THEN
                     DO M = 1, SMAX
                        BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_S(IJK,M)*ROP_S(IJK,M)*&
                           DAREA
                     ENDDO
                  ENDIF  ! end if (.not.discrete_element)
               ENDIF  ! end if (fluid_at(ijk))
            ENDDO    ! end do loop (i=imin1,imax1)
         ENDDO    ! end do loop (k=kmin1,kmax1)

! Global Sum
         call global_all_sum(bed_weight)
         call global_all_sum(area)
         IF (AREA /= 0.0) BED_WEIGHT = BED_WEIGHT/AREA

         PJ = PJ + BED_WEIGHT
         DO K = KMIN1, KMAX1
            DO I = IMIN1, IMAX1
               IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK).AND.P_G(IJK)==UNDEFINED)P_G(IJK)=SCALE_PRESSURE(PJ)
            ENDDO    ! end do (i=imin1,imax1)
         ENDDO   ! end do (k = kmin1,kmax1)
      ENDDO   ! end do (j=jmax2,jimn1, -1)
! end setting an undefined pressure in an initial condition region
! ----------------------------------------------------------------<<<

  100 CONTINUE

      call send_recv(P_G,2)

      RETURN

 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/)

      END SUBROUTINE SET_FLUIDBED_P


