!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS1

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param, only: dimension_3
      USE param1, only: zero, one
      USE fldvar, only: u_s, v_s, w_s, rop_s, ro_s
      USE geometry
      USE indices
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE mfix_pic
      USE functions
      use particle_filter, only: DES_INTERP_ON
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_SIZE
      use sendrecvnode, only: DES_COLLECT_gDATA
      use physprop, only: mmax
      use param1, only: small_number

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: particles, filter cells, phases
      INTEGER NP, LC, M, mLB, mUB
! Fluid cell index
      INTEGER IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,MMAX+DES_MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight, and times filter weight
      DOUBLE PRECISION :: VOL_WT, VOLxWEIGHT
! Loop bound for filter
      INTEGER :: LP_BND


!-----------------------------------------------


      SOLVOLINC(:,:) = ZERO

      mLB = MMAX+1
      mUB = DES_MMAX+MMAX

      IF(MPPIC) THEN
! initialize only information related to the discrete 'phases' of these 
! continuous variables
         U_S(:,mLB:mUB) = ZERO
         V_S(:,mLB:mUB) = ZERO
         IF(DO_K) W_S(:,mLB:mUB) = ZERO
      ENDIF

! Calculate the gas phase forces acting on each particle.
!$omp parallel default(none) &
!$omp private(NP, VOL_WT, M, LC, IJK, VOLXWEIGHT) &
!$omp shared(MAX_PIP, PVOL, DES_STAT_WT, PIJK, LP_BND, MPPIC, &
!$omp       FILTER_WEIGHT, SOLVOLINC, U_S, V_S, W_S, DO_K, &
!$omp       FILTER_CELL, FILTER_SIZE, DES_VEL_NEW)
!$omp do
      do NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP)) CYCLE

         VOL_WT = PVOL(NP)
         IF(MPPIC) VOL_WT = VOL_WT*DES_STAT_WT(NP)
! Particle phase for data binning.
         M = PIJK(NP,5)

         DO LC=1,FILTER_SIZE
            IJK = FILTER_CELL(LC,NP)
! Particle volume times the weight for this cell.
            VOLxWEIGHT = VOL_WT*FILTER_WEIGHT(LC,NP)
! Accumulate total solids volume (by phase)
!$omp atomic
            SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + VOLxWEIGHT
            IF(MPPIC) THEN
! Accumulate total solids momentum (by phase)
!$omp atomic
               U_S(IJK,1) = U_S(IJK,1) + &
                  DES_VEL_NEW(NP,1)*VOLxWEIGHT
!$omp atomic
               V_S(IJK,1) = V_S(IJK,1) + &
                  DES_VEL_NEW(NP,2)*VOLxWEIGHT
               IF(DO_K) THEN
!$omp atomic
                  W_S(IJK,1) = W_S(IJK,1) + &
                     DES_VEL_NEW(NP,3)*VOLxWEIGHT
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel


! Summ data interpolted into ghost cells into physical cells
!---------------------------------------------------------------------//
      IF(DES_INTERP_ON) THEN
         CALL DES_COLLECT_gDATA(SOLVOLINC(:,mLB:mUB))
         IF(MPPIC) THEN
            CALL DES_COLLECT_gDATA(U_s(:,1))
            CALL DES_COLLECT_gDATA(V_s(:,1))
            IF(DO_K) CALL DES_COLLECT_gDATA(W_s(:,1))
         ENDIF
      ENDIF


! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!---------------------------------------------------------------------//
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)                &
!$omp private(IJK,M,OoSOLVOL)
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

! calculating the bulk density of solids phase
         ROP_S(IJK,mLB:mUB) = RO_S(IJK,mLB:mUB)*&
            SOLVOLINC(IJK,mLB:mUB)/VOL(IJK)

! calculating the cell average solids velocity for each solids phase
         IF(MPPIC) THEN
            OoSOLVOL = sum(SOLVOLINC(IJK,:))
            IF(OoSOLVOL > SMALL_NUMBER) THEN
               OoSOLVOL = ONE/OoSOLVOL
               U_s(IJK,1) = U_s(IJK,1)*OoSOLVOL
               V_s(IJK,1) = V_s(IJK,1)*OoSOLVOL
               IF(DO_K) W_s(IJK,1) = W_s(IJK,1)*OoSOLVOL
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel do


! Halo exchange of solids volume fraction data.
      calL SEND_RECV(ROP_S,2)

      end SUBROUTINE COMP_MEAN_FIELDS1
