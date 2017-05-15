!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS0

! Modules
      USE bc
      USE compar
      USE constant
      USE cutcell
      USE derived_types, only: pic
      USE desmpi
      USE discretelement
      USE drag
      USE fldvar
      USE functions, only: FLUID_AT
      USE functions, only: FUNIJK
      USE functions, only: IS_ON_myPE_wobnd
      USE geometry
      USE indices
      USE interpolation
      USE mfix_pic
      USE mpi_node_des, only: des_addnodevalues_mean_fields
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE particle_filter, only: DES_REPORT_MASS_INTERP
      USE physprop, only: MMAX
      USE run, only: solids_model
      USE sendrecv

      IMPLICIT NONE
! Local variables
!---------------------------------------------------------------------//
! general i, j, k indices
      INTEGER :: I, J, K, IJK, &
                 II, JJ, KK
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: IDIM, IJK2
      INTEGER :: CUR_IJK
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3):: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to
! set_interpolation_stencil
      INTEGER :: ONEW
! constant whose value depends on dimension of system
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
! index of solid phase that particle NP belongs to
      INTEGER :: M
! particle number index, used for looping
      INTEGER :: NP, NINDX
! Statistical weight of the particle. Equal to one for DEM
      DOUBLE PRECISION :: WTP

      DOUBLE PRECISION :: MASS_SOL1, MASS_SOL2
! sum of mass_sol1 and mass_sol2 across all processors
      DOUBLE PRECISION :: MASS_SOL1_ALL, MASS_SOL2_ALL

      DOUBLE PRECISION :: TEMP1, TEMP2

      DOUBLE PRECISION, DIMENSION(3) :: DES_VEL_DENSITY
      DOUBLE PRECISION :: DES_ROP_DENSITY

      INTEGER :: COUNT_NODES_OUTSIDE, COUNT_NODES_INSIDE, &
                 COUNT_NODES_INSIDE_MAX

      DOUBLE PRECISION :: RESID_ROPS(DIMENSION_M)
      DOUBLE PRECISION :: RESID_VEL(3, DIMENSION_M)
      DOUBLE PRECISION :: NORM_FACTOR

!Handan Liu added on Jan 17 2013
      DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp
      DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft

! total number of 'solids' phases in simulation
      INTEGER :: MMAX_TOT
!......................................................................!

! initializing
      MASS_SOL1 = ZERO
      MASS_SOL2 = ZERO
      MASS_SOL1_ALL = ZERO
      MASS_SOL2_ALL = ZERO
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

! cartesian_grid related quantities
      COUNT_NODES_INSIDE_MAX = merge(4, 8, NO_K)

      MMAX_TOT = DES_MMAX+MMAX
! initialize only information related to the discrete 'phases' of
! these continuous variables
      ROP_S(:,MMAX+1:MMAX_TOT) = zero
      U_S(:,MMAX+1:MMAX_TOT) = ZERO
      V_S(:,MMAX+1:MMAX_TOT) = ZERO
      IF(DO_K) W_S(:,MMAX+1:MMAX_TOT) = ZERO
! these are not 'continuous' variables but used only within this routine...
      DES_VEL_NODE = ZERO
      DES_ROPS_NODE = ZERO

! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      CALL SET_INTERPOLATION_SCHEME(2)

!$omp parallel default(shared) &
!$omp private(IJK, I, J, K, PCELL, IW, IE, JS, JN, KB, KTP, ONEW, &
!$omp         GST_TMP, COUNT_NODES_INSIDE, II, JJ, KK, CUR_IJK, &
!$omp         NINDX, NP, WTP, M, WEIGHT_FT, I1, I2, J1, J2, K1, K2, &
!$omp         IDIM, IJK2, NORM_FACTOR, RESID_ROPS, RESID_VEL, &
!$omp         COUNT_NODES_OUTSIDE, TEMP1)
!$omp do reduction(+:MASS_SOL1) reduction(+:DES_ROPS_NODE,DES_VEL_NODE)
      DO IJK = IJKSTART3,IJKEND3

! Cycle this cell if not in the fluid domain or if it contains no
! particle/parcel
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF( PINC(IJK) == 0) CYCLE

         PCELL(1) = I_OF(IJK)-1
         PCELL(2) = J_OF(IJK)-1
         PCELL(3) = merge(K_OF(IJK)-1, 1, DO_K)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         CALL SET_INTERPOLATION_STENCIL(PCELL, IW, IE, JS, JN, KB, KTP,&
            INTERP_SCHEME, DIMN, ORDERNEW=ONEW)

         COUNT_NODES_OUTSIDE = 0
! Computing/setting the geometric stencil
         DO K=1, merge(1, ONEW, NO_K)
         DO J=1, ONEW
         DO I=1, ONEW
            II = IW + I-1
            JJ = JS + J-1
            KK = KB + K-1
            CUR_IJK = funijk_map_c(II,JJ,KK)
            GST_TMP(I,J,K,1) = XE(II)
            GST_TMP(I,J,K,2) = YN(JJ)
            GST_TMP(I,J,K,3) = merge(DZ(1), ZT(KK), NO_K)
            IF(CARTESIAN_GRID) THEN
               IF(SCALAR_NODE_ATWALL(CUR_IJK))                   &
                  COUNT_NODES_OUTSIDE = COUNT_NODES_OUTSIDE + 1
            ENDIF
         ENDDO
         ENDDO
         ENDDO


! Calculate des_rops_node so rop_s, and in turn, ep_g can be updated
!----------------------------------------------------------------->>>

! looping through particles in the cell
         DO NINDX=1, PINC(IJK)
            NP = PIC(IJK)%P(NINDX)
            call DRAG_WEIGHTFACTOR(gst_tmp,des_pos_new(np,:),weight_ft)
            M = PIJK(NP,5)
            WTP = ONE

            IF(MPPIC) WTP = DES_STAT_WT(NP)

            MASS_SOL1 = MASS_SOL1 + PMASS(NP)*WTP
            TEMP2 = PMASS(NP)*WTP
            DO K = 1, merge(1, ONEW, NO_K)
            DO J = 1, ONEW
            DO I = 1, ONEW
! shift loop index to new variables for manipulation
               II = IW + I-1
               JJ = JS + J-1
               KK = KB + K-1
               CUR_IJK = FUNIJK_MAP_C(II,JJ,KK)
               TEMP1 = WEIGHT_FT(I,J,K)*TEMP2

               DES_ROPS_NODE(CUR_IJK,M) = DES_ROPS_NODE(CUR_IJK,M) +   &
                  TEMP1
               DES_VEL_NODE(CUR_IJK,:,M) = DES_VEL_NODE(CUR_IJK,:,M) + &
                  TEMP1*DES_VEL_NEW(NP,:)
            ENDDO
            ENDDO
            ENDDO
         ENDDO   ! end do (nindx=1,pinc(ijk))
!-----------------------------------------------------------------<<<


! Only for cutcell cases may count_nodes_inside become less than its
! original set value. In such an event, the contribution of scalar nodes
! that do not reside in the domain is added to a residual array. This
! array is then redistribited equally to the nodes that are in the fluid
! domain. These steps are done to conserve mass.
!----------------------------------------------------------------->>>
! initializing
         RESID_ROPS = ZERO
         RESID_VEL = ZERO
         IF (CARTESIAN_GRID) THEN

! only for cartesian_grid will count_nodes_outside be modified from zero
            COUNT_NODES_INSIDE = COUNT_NODES_INSIDE_MAX - &
                                 COUNT_NODES_OUTSIDE

            IF(COUNT_NODES_INSIDE.LT.COUNT_NODES_INSIDE_MAX) THEN

! Convention used to number node numbers
! i=1, j=2           i=2, j=2
!   _____________________
!   |                   |
!   |  I = 2, J = 2     |
!   |___________________|
! i=1, j=1           i=2, j=1
! setting indices based on convention
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               I1 = I-1
               I2 = I
               J1 = J-1
               J2 = J
               K1 = merge(K, K-1, NO_K)
               K2 = K

! first calculate the residual des_rops_node and des_vel_node that was
! computed on nodes that do not belong to the domain
               DO KK = K1, K2
               DO JJ = J1, J2
               DO II = I1, I2
                  IJK2 = funijk(II, JJ, KK)
                  IF(SCALAR_NODE_ATWALL(IJK2)) THEN
                     DO M = MMAX+1,MMAX_TOT
                        RESID_ROPS(M) = RESID_ROPS(M) + &
                           DES_ROPS_NODE(IJK2,M)
                        DES_ROPS_NODE(IJK2,M) = ZERO

                        DO IDIM = 1, merge(2,3,NO_K)
                           RESID_VEL(IDIM,M) = RESID_VEL(IDIM, M) + &
                              DES_VEL_NODE(IJK2,IDIM, M)
                           DES_VEL_NODE(IJK2,IDIM, M) = ZERO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               ENDDO
               ENDDO

! now add this residual equally to the remaining nodes
               NORM_FACTOR = ONE/REAL(COUNT_NODES_INSIDE)
               DO KK = K1, K2
               DO JJ = J1, J2
               DO II = I1, I2
                  IJK2 = funijk(II, JJ, KK)
                  IF(.NOT.SCALAR_NODE_ATWALL(IJK2)) THEN
                     DO M = MMAX+1,MMAX_TOT
                        DES_ROPS_NODE(IJK2,M) = &
                           DES_ROPS_NODE(IJK2,M) + &
                           RESID_ROPS(M)*NORM_FACTOR
                        DO IDIM = 1, merge(2,3,NO_K)
                           DES_VEL_NODE(IJK2,IDIM, M) = &
                              DES_VEL_NODE(IJK2,IDIM, M) + &
                              RESID_VEL(IDIM, M)*NORM_FACTOR
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               ENDDO
               ENDDO

            ENDIF
         ENDIF   ! end if (cartesian_grid)
      ENDDO
!$omp end parallel


! At the interface des_rops_node has to be added since particles
! across the processors will contribute to the same scalar node.
! sendrecv will be called and the node values will be added
! at the junction. des_rops_node is altered by the routine when
! periodic boundaries are invoked
      CALL DES_ADDNODEVALUES_MEAN_FIELDS


! Now go from node to scalar center. Same convention as sketched
! earlier
!----------------------------------------------------------------->>>
! Explanation by RG: 08/17/2012
! the approach used here is to make it general enough for cutcells to be
! included as well. The new changes do not alter earlier calculations
! but make the technique general as to include cartesian grid (cut-cell)
! simulations.
! Previously, the volume of the node (by array des_vol_node) was used to
! first scale the nodal values. Subsequently, these nodal values were
! equally added to compute the cell centered values for the scalar cell.

! Consider an internal node next to an edge node (a node adjacent to a
! boundary). In 2D, the volume of an edge node will be half that of an
! internal node. And an edge node will contribute double compared to
! an internal node to the value of the scalar cell they share. These
! calculations were previously accomplished via the variable volume of
! node.  Now this is accomplished by the ratio vol(ijk2)/vol_sur, where
! vol(ijk2) is the volume of the scalar cell in consideration and
! vol_sur is the sum of all the scalar cell volumes that have this node
! as the common node.

!$omp parallel do default(none) collapse (3) &
!$omp shared(KSTART2, KEND1, JSTART2, JEND1, ISTART2, IEND1, DO_K, &
!$omp        VOL, DEAD_CELL_AT, FUNIJK_MAP_C, VOL_SURR, MMAX_TOT, &
!$omp        MMAX, DES_ROPS_NODE, DES_VEL_NODE) &
!$omp private(I, J, K, IJK, M, II, JJ, KK, IJK2, DES_ROP_DENSITY, &
!$omp         DES_VEL_DENSITY) &
!$omp reduction(+:ROP_S, U_S, V_S, W_S)
      DO K = KSTART2, KEND1
      DO J = JSTART2, JEND1
      DO I = ISTART2, IEND1
         IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
         IJK = funijk(I,J,K)
         if (VOL_SURR(IJK).eq.ZERO) CYCLE ! no FLUID_AT any of the stencil points have

! looping over stencil points (NODE VALUES)
         DO M = MMAX+1, MMAX_TOT
            DES_ROP_DENSITY = DES_ROPS_NODE(IJK, M)/VOL_SURR(IJK)
            DES_VEL_DENSITY(:) = DES_VEL_NODE(IJK, :, M)/VOL_SURR(IJK)

            DO KK = K, merge(K+1, K, DO_K)
            DO JJ = J, J+1
            DO II = I, I+1
               IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells

               IJK2 = funijk_map_c(II, JJ, KK)
               IF(FLUID_AT(IJK2).and.(IS_ON_myPE_wobnd(II, JJ, KK))) THEN
! Since the data in the ghost cells is spurious anyway and overwritten during
! subsequent send receives, do not compute any value here as this will
! mess up the total mass value that is computed below to ensure mass conservation
! between Lagrangian and continuum representations
                  ROP_S(IJK2, M) = ROP_S(IJK2, M) + DES_ROP_DENSITY*VOL(IJK2)
                  U_S(IJK2, M) = U_S(IJK2, M) + DES_VEL_DENSITY(1)*VOL(IJK2)
                  V_S(IJK2, M) = V_S(IJK2, M) + DES_VEL_DENSITY(2)*VOL(IJK2)
                  IF(DO_K) W_S(IJK2, M) = W_S(IJK2, M) + DES_VEL_DENSITY(3)*VOL(IJK2)
               ENDIF
            ENDDO  ! end do (ii=i1,i2)
            ENDDO  ! end do (jj=j1,j2)
            ENDDO  ! end do (kk=k1,k2)
         ENDDO
      ENDDO   ! end do (i=istart2,iend1)
      ENDDO   ! end do (j=jstart2,jend1)
      ENDDO   ! end do (k=kstart2,kend1)
!omp end parallel do
!-----------------------------------------------------------------<<<


!$omp parallel do default(none) private(IJK, M) &
!$omp shared(IJKSTART3, IJKEND3, DO_K, MMAX_TOT, MMAX, ROP_s, U_S, &
!$omp        V_S, W_S, VOL)
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         DO M = MMAX+1, MMAX_TOT
            IF(ROP_S(IJK, M).GT.ZERO) THEN
               U_S(IJK, M) = U_S(IJK,M)/ROP_S(IJK, M)
               V_S(IJK, M) = V_S(IJK,M)/ROP_S(IJK, M)
               IF(DO_K) W_S(IJK, M) = W_S(IJK,M)/ROP_S(IJK, M)
! Divide by scalar cell volume to obtain the bulk density
               ROP_S(IJK, M) = ROP_S(IJK, M)/VOL(IJK)
            ENDIF
         ENDDO   ! end loop over M=1,MMAX_TOT
      ENDDO  ! end loop over IJK=ijkstart3,ijkend3
!omp end parallel do


      IF (MPPIC) CALL SEND_RECV(ROP_S,2)
      CALL CALC_EPG_DES

      IF(MPPIC) THEN
! Now calculate Eulerian mean velocity fields for discrete phases
         CALL SEND_RECV(U_S,2)
         CALL SEND_RECV(V_S,2)
         IF(DO_K) CALL SEND_RECV(W_S,2)

! The Eulerian velocity field is used to set up the stencil to interpolate
! mean solid velocity at the parcel's location. U_S could have also been
! used, but that also would have require the communication at this stage.
! The final interpolated value does not change if the stencil is formed by
! first obtaining face centered Eulerian velocities (U_S, etc.)
! and then computing the node velocities from them or directly computing
! the node velocities from cell centered average velocity field (U_S,
! etc.). We are using the first approach as it is more natural to set
! BC's on solid velocity field in the face centered represenation (U_S,
! etc.)
         IF(.NOT.CARTESIAN_GRID) THEN
            CALL MPPIC_COMP_EULERIAN_VELS_NON_CG
         ELSE
            CALL MPPIC_COMP_EULERIAN_VELS_CG
         ENDIF
      ENDIF   ! end if (.not.mppic)

! turn on the below statements to check if the mass is conserved
! between discrete and continuum representations. Should be turned to
! false for any production runs.
      IF(DES_REPORT_MASS_INTERP) THEN
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.FLUID_AT(IJK)) CYCLE
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
! It is important to check both FLUID_AT and IS_ON_MYPE_WOBND.
            IF(IS_ON_myPE_wobnd(I,J,K)) THEN
               DO M=MMAX+1,MMAX_TOT
                  MASS_SOL2 = MASS_SOL2 + &
                     ROP_S(IJK,M)*VOL(IJK)
               ENDDO
            ENDIF
         ENDDO
         CALL GLOBAL_SUM(MASS_SOL1, MASS_SOL1_ALL)
         CALL GLOBAL_SUM(MASS_SOL2, MASS_SOL2_ALL)
         if(myPE.eq.pe_IO) THEN
            WRITE(*,'(/,5x,A,4(2x,g17.8),/)') &
              'SOLIDS MASS DISCRETE AND CONTINUUM =  ', &
              MASS_SOL1_ALL, MASS_SOL2_ALL
         ENDIF
      ENDIF

      END SUBROUTINE COMP_MEAN_FIELDS0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Subroutine: DRAG_WEIGHTFACTOR                                        C
!  Purpose: DES - Calculate the fluid velocity interpolated at the      C
!           particle's location and weights. Replace 'interpolator'     C
!                       interface for OpenMP implementation.            C
!                                                                       C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_WEIGHTFACTOR(GSTEN,DESPOS,WEIGHTFACTOR)

      use geometry, only: NO_K

        IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: GSTEN
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: DESPOS
        DOUBLE PRECISION, DIMENSION(2,2,2), INTENT(OUT) :: WEIGHTFACTOR
        INTEGER :: II, JJ, KK

        DOUBLE PRECISION, DIMENSION(2) :: XXVAL, YYVAL, ZZVAL
        DOUBLE PRECISION :: DXX, DYY, DZZ
        DOUBLE PRECISION, DIMENSION(3) :: ZETAA

        DXX = GSTEN(2,1,1,1) - GSTEN(1,1,1,1)
        DYY = GSTEN(1,2,1,2) - GSTEN(1,1,1,2)

        ZETAA(1:2) = DESPOS(1:2) - GSTEN(1,1,1,1:2)

        ZETAA(1) = ZETAA(1)/DXX
        ZETAA(2) = ZETAA(2)/DYY

        XXVAL(1)=1-ZETAA(1)
        YYVAL(1)=1-ZETAA(2)
        XXVAL(2)=ZETAA(1)
        YYVAL(2)=ZETAA(2)

        IF(NO_K) THEN
           DO JJ=1,2
              DO II=1,2
                 WEIGHTFACTOR(II,JJ,1) = XXVAL(II)*YYVAL(JJ)
              ENDDO
           ENDDO
        ELSE
           DZZ = GSTEN(1,1,2,3) - GSTEN(1,1,1,3)
           ZETAA(3) = DESPOS(3) - GSTEN(1,1,1,3)
           ZETAA(3) = ZETAA(3)/DZZ
           ZZVAL(1)=1-ZETAA(3)
           ZZVAL(2)=ZETAA(3)
           DO KK=1,2
              DO JJ=1,2
                 DO II=1,2
                    WEIGHTFACTOR(II,JJ,KK) = XXVAL(II)*YYVAL(JJ)*ZZVAL(KK)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

      END SUBROUTINE DRAG_WEIGHTFACTOR
