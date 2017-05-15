!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MPPIC_MI_BC                                             C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE MASS_INFLOW_PIC

! Modules
!---------------------------------------------------------------------//
      USE bc
      USE constant
      USE cutcell
      USE discretelement
      use error_manager
      USE functions
      USE geometry
      USE mfix_pic
      USE mpi_utility
      USE param, only: dimension_m
      USE param1, only: half, zero
      use physprop, only: mmax, D_p0, ro_s0
      USE pic_bc
      USE randomno
      USE run
      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
      INTEGER :: WDIR, IDIM, IPCOUNT, LAST_EMPTY_SPOT, NEW_SPOT
      INTEGER :: BCV, BCV_I, L, LC, PIP_ADD_COUNT, IPROC
      INTEGER ::  IFLUID, JFLUID, KFLUID, IJK, M
      DOUBLE PRECISION :: CORD_START(3), DOML(3),WALL_NORM(3)
      DOUBLE PRECISION :: AREA_INFLOW, VEL_INFLOW(DIMN), STAT_WT

      LOGICAL :: DELETE_PART
      DOUBLE PRECISION :: EPS_INFLOW(DIMENSION_M)
      DOUBLE PRECISION :: REAL_PARTS(DIMENSION_M)
      DOUBLE PRECISION :: COMP_PARTS(DIMENSION_M)
      DOUBLE PRECISION :: VEL_NORM_MAG, VOL_INFLOW, VOL_IJK

! Temp logical variables for checking constant npc and statwt specification
      LOGICAL :: CONST_NPC, CONST_STATWT

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RANDPOS
      INTEGER :: CNP_CELL_COUNT
      DOUBLE PRECISION :: RNP_CELL_COUNT
      integer :: lglobal_id
      integer, dimension(0:numpes-1) :: add_count_all

      LOGICAL, parameter :: setDBG = .true.
      LOGICAL :: dFlag

      type :: ty_spotlist
         integer spot
         type(ty_spotlist),pointer :: next => NULL()
      end type ty_spotlist

      type(ty_spotlist),pointer :: &
           cur_spotlist => NULL(), &
           prev_spotlist => NULL(), &
           temp_spotlist => NULL()
!......................................................................!


      CALL INIT_ERR_MSG("PIC_MI_BC")

      dFlag = (DMP_LOG .AND. setDBG)
      PIP_ADD_COUNT = 0
      LAST_EMPTY_SPOT = 0

      ALLOCATE(cur_spotlist); cur_spotlist%spot = -1

      DO BCV_I = 1, PIC_BCMI
         BCV = PIC_BCMI_MAP(BCV_I)

         WALL_NORM(1:3) =  PIC_BCMI_NORMDIR(BCV_I,1:3)

         !Find the direction of the normal for this wall
         WDIR = 0
         DO IDIM = 1, DIMN
            WDIR = WDIR + ABS(WALL_NORM(IDIM))*IDIM
         end DO

         DO LC=PIC_BCMI_IJKSTART(BCV_I), PIC_BCMI_IJKEND(BCV_I)
            IJK = PIC_BCMI_IJK(LC)

            IF(.NOT.FLUID_AT(IJK)) CYCLE

            IFLUID = I_OF(IJK)
            JFLUID = J_OF(IJK)
            KFLUID = K_OF(IJK)

            CORD_START(1) = XE(IFLUID) - PIC_BCMI_OFFSET (BCV_I,1)*DX(IFLUID)

            CORD_START(2) = YN(JFLUID) - PIC_BCMI_OFFSET (BCV_I,2)*DY(JFLUID)


            CORD_START(3) = merge(zero, ZT(KFLUID) - PIC_BCMI_OFFSET (BCV_I,3)*DZ(KFLUID), no_k)

            DOML(1) = DX(IFLUID)
            DOML(2) = DY(JFLUID)
            DOML(3) = MERGE(DZ(1), DZ(KFLUID), NO_K)

            AREA_INFLOW = DOML(1)*DOML(2)*DOML(3)/DOML(WDIR)

            VOL_IJK = DOML(1)*DOML(2)*DOML(3)

            DOML(WDIR) = ZERO
            !set this to zero as the particles will
            !be seeded only on the BC plane

            DO M = 1, DES_MMAX+MMAX

               IF(SOLIDS_MODEL(M) /= 'PIC') CYCLE

               EPS_INFLOW(M) = BC_ROP_S(BCV, M)/RO_S0(M)
               VEL_INFLOW(1) = BC_U_S(BCV, M)
               VEL_INFLOW(2) = BC_V_S(BCV, M)
               VEL_INFLOW(3) = BC_W_S(BCV, M)

               VEL_NORM_MAG = ABS(DOT_PRODUCT(VEL_INFLOW(1:DIMN), WALL_NORM(1:DIMN)))
               VOL_INFLOW   = AREA_INFLOW*VEL_NORM_MAG*DTSOLID

               REAL_PARTS(M) = 6.d0*EPS_INFLOW(M)*VOL_INFLOW/&
                               (PI*(D_p0(M)**3.d0))
               COMP_PARTS(M) = zero

               CONST_NPC    = (BC_PIC_MI_CONST_NPC   (BCV, M) .ne. 0)
               CONST_STATWT = (BC_PIC_MI_CONST_STATWT(BCV, M) .ne. ZERO)
               IF(CONST_NPC) THEN
                  IF(EPS_INFLOW(M).GT.ZERO) &
                  COMP_PARTS(M) = REAL(BC_PIC_MI_CONST_NPC(BCV, M))* &
                  VOL_INFLOW/VOL_IJK
               ELSEIF(CONST_STATWT) THEN
                  COMP_PARTS(M) = REAL_PARTS(M)/ &
                  BC_PIC_MI_CONST_STATWT(BCV, M)
               ENDIF

               pic_bcmi_rnp(LC,M) = pic_bcmi_rnp(LC,M) + REAL_PARTS(M)

               pic_bcmi_cnp(LC,M) = pic_bcmi_cnp(LC,M) + COMP_PARTS(M)


               IF(pic_bcmi_cnp(LC,M).GE.1.d0) THEN
                  CNP_CELL_COUNT = INT(pic_bcmi_cnp(LC,M))
                  pic_bcmi_cnp(LC,M) = pic_bcmi_cnp(LC,M) - CNP_CELL_COUNT

                  RNP_CELL_COUNT = pic_bcmi_rnp(LC,M)
                  pic_bcmi_rnp(LC,M)  = zero

                  !set pic_bcmi_rnp to zero to reflect that all real particles have been seeded
                  STAT_WT = RNP_CELL_COUNT/REAL(CNP_CELL_COUNT)
                  ALLOCATE(RANDPOS(CNP_CELL_COUNT*DIMN))
                  CALL UNI_RNO(RANDPOS(:))

                  DO IPCOUNT = 1, CNP_CELL_COUNT


                     CALL PIC_FIND_EMPTY_SPOT(LAST_EMPTY_SPOT, NEW_SPOT)

                     DES_POS_OLD( NEW_SPOT,1:DIMN) =  CORD_START(1:DIMN) &
                          & + RANDPOS((IPCOUNT-1)*DIMN+1: &
                          & (IPCOUNT-1)*DIMN+DIMN)*DOML(1:DIMN)
                     DES_POS_NEW( NEW_SPOT,:) = DES_POS_OLD(NEW_SPOT,:)
                     DES_VEL_OLD(NEW_SPOT,1:DIMN) = VEL_INFLOW(1:DIMN)

                     DES_VEL_NEW( NEW_SPOT,:) = DES_VEL_OLD( NEW_SPOT,:)

                     DES_RADIUS(NEW_SPOT) = D_p0(M)*HALF

                     RO_Sol(NEW_SPOT) =  RO_S0(M)

                     DES_STAT_WT(NEW_SPOT) = STAT_WT

                     PIJK(NEW_SPOT, 1) = IFLUID
                     PIJK(NEW_SPOT, 2) = JFLUID
                     PIJK(NEW_SPOT, 3) = KFLUID
                     PIJK(NEW_SPOT, 4) = IJK
                     PIJK(NEW_SPOT, 5) = M

                     PVOL(NEW_SPOT) = (4.0d0/3.0d0)*Pi*DES_RADIUS(NEW_SPOT)**3
                     PMASS(NEW_SPOT) = PVOL(NEW_SPOT)*RO_SOL(NEW_SPOT)

                     FC( NEW_SPOT,:) = zero
                     DELETE_PART = .false.
                     IF(PIC_BCMI_INCL_CUTCELL(BCV_I)) &
                          CALL CHECK_IF_PARCEL_OVERLAPS_STL &
                          (des_pos_new( NEW_SPOT,1:dimn), &
                          DELETE_PART)

                     IF(.NOT.DELETE_PART) THEN

                        PIP = PIP+1
                        PIP_ADD_COUNT = PIP_ADD_COUNT + 1
                        CALL SET_NORMAL(NEW_SPOT)
                        ! add to the list
                        ALLOCATE(temp_spotlist)
                        temp_spotlist%spot = new_spot
                        temp_spotlist%next => cur_spotlist
                        cur_spotlist => temp_spotlist
                        nullify(temp_spotlist)

                     ELSE
                        CALL SET_NONEXISTENT(NEW_SPOT)
                        LAST_EMPTY_SPOT = NEW_SPOT - 1
                     ENDIF

                     !WRITE(*,'(A,2(2x,i5), 2x, A,2x, 3(2x,i2),2x, A, 3(2x,g17.8))') &
                     !   'NEW PART AT ', NEW_SPOT, MAX_PIP, 'I, J, K = ', IFLUID, JFLUID, KFLUID, 'POS =', DES_POS_NEW(NEW_SPOT,:)
                     !IF(DMP_LOG) WRITE(UNIT_LOG,'(A,2x,i5, 2x, A,2x, 3(2x,i2),2x, A, 3(2x,g17.8))') &
                     !    'NEW PART AT ', NEW_SPOT, 'I, J, K = ', IFLUID, JFLUID, KFLUID, 'POS =', DES_POS_NEW(NEW_SPOT,:)

                     !WRITE(*,*) 'WDIR, DOML = ', WDIR, DOML(:)
                  END DO
                  DEALLOCATE(RANDPOS)

               end IF
            end DO
         end DO
      end DO

!Now assign global id to new particles added
      add_count_all(:) = 0
      add_count_all(mype) = pip_add_count
      call global_all_sum(add_count_all(0:numpes-1))
      lglobal_id = imax_global_id + sum(add_count_all(0:mype-1))

      do l = 1,pip_add_count
         lglobal_id = lglobal_id + 1
         iglobal_id(cur_spotlist%spot)= lglobal_id
         prev_spotlist=> cur_spotlist
         cur_spotlist => cur_spotlist%next
         deallocate(prev_spotlist)
      end do
      deallocate(cur_spotlist)
      imax_global_id = imax_global_id + sum(add_count_all(0:numpes-1))

      IF(PIC_REPORT_SEEDING_STATS) then

         IF(SUM(ADD_COUNT_ALL(:)).GT.0) THEN
            WRITE(err_msg,'(/,2x,A,2x,i10)') &
               'TOTAL NUMBER OF PARCELS ADDED GLOBALLY = ',&
                SUM(ADD_COUNT_ALL(:))

            call flush_err_msg(header = .false., footer = .false.)

            DO IPROC = 0, NUMPES-1
               IF(DMP_LOG) WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
                  'PARCELS ADDED ON PROC:', IPROC,&
                  ' EQUAL TO', ADD_COUNT_ALL(IPROC)
            ENDDO
         ENDIF

      ENDIF

      CALL FINL_ERR_MSG
      RETURN
      END SUBROUTINE MASS_INFLOW_PIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_EMPTY_SPOT                                     C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PIC_FIND_EMPTY_SPOT(LAST_INDEX, EMPTY_SPOT)

! Modules
!---------------------------------------------------------------------//
      USE funits
      USE mpi_utility
      USE error_manager
      USE discretelement, only: max_pip
      USE functions, only: is_nonexistent
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(INOUT) :: LAST_INDEX
      INTEGER, INTENT(OUT) :: EMPTY_SPOT

! Local Variables
!---------------------------------------------------------------------//
      LOGICAL :: SPOT_FOUND
      INTEGER :: LL
!......................................................................!

      CALL INIT_ERR_MSG("PIC_FIND_EMPTY_SPOT")

      IF(LAST_INDEX.EQ.MAX_PIP) THEN
         WRITE(ERR_MSG,2001)

         CALL FLUSH_ERR_MSG(abort = .true.)
         call mfix_exit(mype)
      ENDIF
      SPOT_FOUND = .false.

      DO LL = LAST_INDEX+1, MAX_PIP

         if(IS_NONEXISTENT(LL)) THEN
            EMPTY_SPOT = LL
            LAST_INDEX = LL
            SPOT_FOUND = .true.
            EXIT
         ENDIF
      ENDDO

      IF(.NOT.SPOT_FOUND) THEN
         WRITE(ERR_MSG,2002)
         CALL FLUSH_ERR_MSG(abort = .true.)
      ENDIF

 2001 FORMAT(/,5X,  &
      & 'ERROR IN PIC_FIND_EMPTY_SPOT', /5X, &
      & 'NO MORE EMPTY SPOT IN THE PARTICLE ARRAY TO ADD A NEW PARTICLE',/5X, &
      & 'TERMINAL ERROR: STOPPING')

 2002 FORMAT(/,5X,  &
      & 'ERROR IN PIC_FIND_EMPTY_SPOT', /5X, &
      & 'COULD NOT FIND A SPOT FOR ADDING NEW PARTICLE',/5X, &
      & 'INCREASE THE SIZE OF THE INITIAL ARRAYS', 5X, &
      & 'TERMINAL ERROR: STOPPING')

      CALL FINL_ERR_MSG
      END SUBROUTINE PIC_FIND_EMPTY_SPOT





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PIC_FIND_NEW_CELL                                     C
!  Purpose:                                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PIC_FIND_NEW_CELL(LL)

! Modules
!---------------------------------------------------------------------//
      USE discretelement
      use mpi_utility
      USE cutcell
      USE functions
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: LL

! Local Variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: XPOS, YPOS, ZPOS
      INTEGER :: I, J, K
!......................................................................!

      I = PIJK(LL,1)
      J = PIJK(LL,2)
      K = PIJK(LL,3)
      XPOS = DES_POS_NEW(LL,1)
      YPOS = DES_POS_NEW(LL,2)
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(LL,3)
      ENDIF

      IF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN
         PIJK(LL,1) = I
      ELSEIF(XPOS >= XE(I)) THEN
         PIJK(LL,1) = I+1
      ELSE
         PIJK(LL,1) = I-1
      END IF

      IF(YPOS >= YN(J-1) .AND. YPOS < YN(J)) THEN
         PIJK(LL,2) = J
      ELSEIF(YPOS >= YN(J))THEN
         PIJK(LL,2) = J+1
      ELSE
         PIJK(LL,2) = J-1
      END IF

      IF(DIMN.EQ.2) THEN
         PIJK(LL,3) = 1
      ELSE
         IF(ZPOS >= ZT(K-1) .AND. ZPOS < ZT(K)) THEN
            PIJK(LL,3) = K
         ELSEIF(ZPOS >= ZT(K)) THEN
            PIJK(LL,3) = K+1
         ELSE
            PIJK(LL,3) = K-1
         END IF
      ENDIF

      I = PIJK(LL,1)
      J = PIJK(LL,2)
      K = PIJK(LL,3)

      IF(I.EQ.IEND1+1) then
         IF(XPOS>=XE(IEND1-1) .AND. XPOS<=XE(IEND1)) PIJK(LL,1) = IEND1
      ENDIF

      IF(J.EQ.JEND1+1) then
         IF(YPOS>=YN(JEND1-1) .AND. YPOS<=YN(JEND1)) PIJK(LL,2) = JEND1
      ENDIF

      IF(DIMN.EQ.3.AND.K.EQ.KEND1+1) THEN
         IF(ZPOS>=ZT(KEND1-1) .AND. ZPOS<=ZT(KEND1)) PIJK(LL,3) = KEND1
      ENDIF

      PIJK(LL,4) = FUNIJK(PIJK(LL,1), PIJK(LL,2),PIJK(LL,3))

      END SUBROUTINE PIC_FIND_NEW_CELL

