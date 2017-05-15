!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_CELL_FLAGS                                  C
!  Purpose: Set flags for scalar cut cells, based on intersection      C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_CELL_FLAGS

      USE compar, ONLY: mype, pe_io, ijkstart3, ijkend3, istart1, iend1, jstart1, jend1, kstart1, kend1
      USE cutcell
      USE functions, ONLY: funijk, ip_of, jp_of, kp_of, bottom_of, south_of, west_of, fluid_at, is_on_mype_wobnd
      USE geometry, ONLY: DO_I, DO_J, DO_K, IMIN1, IMAX3, JMIN1, JMAX3, KMIN1, KMAX3, no_k, vol, axy, axz, ayz, dx, dy, dz, flag
      USE indices, ONLY: i_of, j_of, k_of
      USE param, only: dimension_3
      USE param1, only: zero, half, one, undefined
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: tol_f
      USE sendrecv
      USE cut_cell_preproc, ONLY: cad_intersect, clean_intersect, clean_intersect_scalar, eval_f, intersect
      USE vtk, ONLY: GLOBAL_VAR_ALLOCATED, GRID_INFO_PRINTED_ON_SCREEN

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,N_N1,N_N2,Q_ID
      INTEGER :: MIN_INTERSECTIONS,MAX_INTERSECTIONS
      LOGICAL :: CLIP_FLAG
      INTEGER :: BCID

      INTEGER, DIMENSION(DIMENSION_3,15) :: OLD_CONNECTIVITY
      DOUBLE PRECISION, DIMENSION(:), allocatable :: X_OLD_POINT,Y_OLD_POINT,Z_OLD_POINT

      INTEGER, DIMENSION(6) :: NB

      INTEGER :: IJK_NB,NODE_NB,NN,N_NB,NC

      DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,D,TOT_VOL_NODE,TOT_VOL_CELLS

      DOUBLE PRECISION, DIMENSION(:,:), allocatable :: SCALAR_NODE_XYZ_TEMP

      DOUBLE PRECISION :: DIST, NORM1, NORM2, NORM3,Diagonal
      INTEGER :: IJK2, I1, I2, J1, J2, K1, K2, II, JJ, KK
      LOGICAL :: COND_1, COND_2

      allocate(X_OLD_POINT(DIMENSION_MAX_CUT_CELL))
      allocate(Y_OLD_POINT(DIMENSION_MAX_CUT_CELL))
      allocate(Z_OLD_POINT(DIMENSION_MAX_CUT_CELL))
      allocate(SCALAR_NODE_XYZ_TEMP(DIMENSION_3, 3))

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'INTERSECTING GEOMETRY WITH SCALAR CELLS...'
      ENDIF
10    FORMAT(1X,A)

      GLOBAL_VAR_ALLOCATED = .FALSE.
      GRID_INFO_PRINTED_ON_SCREEN = .FALSE.

!  Set arrays for computing indices
!

!      CALL SET_INCREMENTS

      DX(IMAX3+1) = DX(IMAX3)
      DY(JMAX3+1) = DY(JMAX3)
      DZ(KMAX3+1) = DZ(KMAX3)

!     x-Location of U-momentum cells for original (uncut grid)
      IF (DO_I) THEN
         XG_E(1) = ZERO
         DO I = IMIN1, IMAX3
            XG_E(I) = XG_E(I-1) + DX(I)
         END DO
      ENDIF


!     y-Location of V-momentum cells for original (uncut grid)
      IF (DO_J) THEN
         YG_N(1) = ZERO
         DO J = JMIN1, JMAX3
            YG_N(J) = YG_N(J-1) + DY(J)
         END DO
      ENDIF


!     z-Location of W-momentum cells for original (uncut grid)
      IF (DO_K) THEN
         ZG_T(1) = ZERO
         DO K = KMIN1, KMAX3
            ZG_T(K) = ZG_T(K-1) + DZ(K)
         END DO
      ENDIF

      PARTITION = DBLE(myPE)       ! ASSIGN processor ID (for vizualisation)

      CALL SET_GEOMETRY1   ! INITIALIZE ALL VOLUMES AND AREAS


!======================================================================
!  Grid coordinates:
!======================================================================

      SMALL_CELL_AT = .FALSE.
      SMALL_CELL_FLAG = 0
      NUMBER_OF_SMALL_CELLS = 0

!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF(NO_K) THEN   ! 2D case

            INTERIOR_CELL_AT(IJK) = (     (I >= ISTART1 ).AND.(I <= IEND1 )  &
                                     .AND.(J >= JSTART1 ).AND.(J <= JEND1 ) )

         ELSE            ! 3D case

            INTERIOR_CELL_AT(IJK) = (     (I >= ISTART1 ).AND.(I <= IEND1 )  &
                                     .AND.(J >= JSTART1 ).AND.(J <= JEND1 )  &
                                     .AND.(K >= KSTART1 ).AND.(K <= KEND1 ) )

         ENDIF

      END DO

       NUMBER_OF_NODES = 0

      CALL GET_POTENTIAL_CUT_CELLS

!       POTENTIAL_CUT_CELL_AT=.TRUE.

      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
            IF(POTENTIAL_CUT_CELL_AT(IJK))  CALL INTERSECT(IJK,'SCALAR',Xn_int(IJK),Ye_int(IJK),Zt_int(IJK))
!            CALL INTERSECT(IJK,'SCALAR',Xn_int(IJK),Ye_int(IJK),Zt_int(IJK))
         END DO

!======================================================================
!  Clean-up intersection flags by snapping intersection points to close
!  corner (within tolerance TOL_SNAP)
!======================================================================
         CALL CLEAN_INTERSECT_SCALAR

      ELSE
         CALL CAD_INTERSECT('SCALAR',Xn_int,Ye_int,Zt_int)
      ENDIF
!        NUMBER_OF_NODES = 0

      DO IJK = IJKSTART3, IJKEND3

          IF(POTENTIAL_CUT_CELL_AT(IJK))  THEN

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

         CALL GET_CELL_NODE_COORDINATES(IJK,'SCALAR')

         SCALAR_NODE_XYZ_TEMP(IJK,1) = X_NODE(8)
         SCALAR_NODE_XYZ_TEMP(IJK,2) = Y_NODE(8)
         SCALAR_NODE_XYZ_TEMP(IJK,3) = Z_NODE(8)
!======================================================================
!  Initialize location of velocity nodes
!======================================================================

         X_U(IJK) = X_NODE(8)
         Y_U(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
         Z_U(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

         X_V(IJK) = HALF * (X_NODE(7) + X_NODE(8))
         Y_V(IJK) = Y_NODE(8)
         Z_V(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

         X_W(IJK) = HALF * (X_NODE(7) + X_NODE(8))
         Y_W(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
         Z_W(IJK) = Z_NODE(8)

         IF(INTERIOR_CELL_AT(IJK)) THEN
!======================================================================
!  Get Connectivity
!======================================================================

            IF(NO_K) THEN
               MIN_INTERSECTIONS = 2
               MAX_INTERSECTIONS = 2
               N_N1 = 5
               N_N2 = 8

            ELSE
               MIN_INTERSECTIONS = 3
               MAX_INTERSECTIONS = 6
               N_N1 = 1
               N_N2 = 8
            ENDIF

            CALL GET_CONNECTIVITY(IJK,'SCALAR',NUMBER_OF_NEW_POINTS,NUMBER_OF_NODES(IJK),CONNECTIVITY,&
            X_NEW_POINT,Y_NEW_POINT,Z_NEW_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_int,Ye_int,Zt_int)

            IF(TOTAL_NUMBER_OF_INTERSECTIONS < MIN_INTERSECTIONS ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

!               CALL EVAL_STL_FCT_AT('SCALAR',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(USE_MSH.OR.USE_STL) THEN
!                  CALL EVAL_STL_FCT_AT('SCALAR',IJK,1,F_NODE(1),CLIP_FLAG,BCID)
!                  CALL EVAL_STL_FCT_AT('SCALAR',IJK,8,F_NODE(8),CLIP_FLAG,BCID)
!                  F_NODE(0) = HALF*(F_NODE(1)+F_NODE(8))
! Pick the first defined value of F_AT at one of the cell corners
! This should avoid wrongly setting a cell next to the global ghost layer as blocked cell
! when the average between undefined and negative is a positive number
! This occurs for external flows only due to F_AT being initialized as UNDEFINED
! which is a positive number (considered outside the fluid region).
                  DO NODE=1,8
                     CALL EVAL_STL_FCT_AT('SCALAR',IJK,NODE,F_NODE(NODE),CLIP_FLAG,BCID)
                     IF(F_NODE(NODE)/=UNDEFINED.AND.F_NODE(NODE)/=ZERO) THEN
                        F_NODE(0) = F_NODE(NODE)
                        EXIT
                     ENDIF

                  ENDDO
               ENDIF


!  Uncomment the line below to force the cell to be a fluid cell
!  Currently, F_NODE(0) > 0 ==> Blocked cell
!             F_NODE(0) < 0 ==> Fluid cell

!               IF(TOTAL_NUMBER_OF_INTERSECTIONS==-1) F_NODE(0) = -ONE

               BC_ID(IJK) = 0

!               IF(F_NODE(0) < TOL_F) THEN
               IF(F_NODE(0) < ZERO) THEN
                  IF((FLAG(IJK)>=100).AND.(FLAG(IJK)<=102)) THEN
                     BLOCKED_CELL_AT(IJK) = .TRUE.
                  ELSE
                     BLOCKED_CELL_AT(IJK) = .FALSE.
                     STANDARD_CELL_AT(IJK) = .TRUE.           ! Regular fluid cell
                  ENDIF
               ELSE
                  FLAG(IJK) = 100
                  BLOCKED_CELL_AT(IJK) = .TRUE.               ! Blocked fluid cell
                  STANDARD_CELL_AT(IJK) = .FALSE.
                  AXY(IJK) = ZERO
                  AXZ(IJK) = ZERO
                  AYZ(IJK) = ZERO
                  VOL(IJK) = ZERO

                  AXY(BOTTOM_OF(IJK)) = ZERO
                  AXZ(SOUTH_OF(IJK)) = ZERO
                  AYZ(WEST_OF(IJK)) = ZERO
               ENDIF

               IF(NO_K) THEN
                  NUMBER_OF_NODES(IJK) = 4
                  CONNECTIVITY(IJK,1) = IJK_OF_NODE(5)
                  CONNECTIVITY(IJK,2) = IJK_OF_NODE(6)
                  CONNECTIVITY(IJK,3) = IJK_OF_NODE(8)
                  CONNECTIVITY(IJK,4) = IJK_OF_NODE(7)
               ELSE
                  NUMBER_OF_NODES(IJK) = 8
                  DO NODE = N_N1,N_N2
                     CONNECTIVITY(IJK,NODE) = IJK_OF_NODE(NODE)
                  END DO
               ENDIF


            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > MAX_INTERSECTIONS ) THEN

               IF(PRINT_WARNINGS) THEN

                  WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN CELL IJK = ',IJK
                  WRITE(*,*) 'MAXIMUM NUMBER OF INTERSECTIONS = ',MAX_INTERSECTIONS
                  WRITE(*,*) 'TOTAL NUMBER OF INTERSECTIONS = ',TOTAL_NUMBER_OF_INTERSECTIONS
                  WRITE(*,*) 'REMOVING SCALAR CELL FROM COMPUTATION.'
!                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
!                  CALL MFIX_EXIT(MYPE)

               ENDIF

               BLOCKED_CELL_AT(IJK) = .TRUE.               ! Blocked fluid cell

            ELSE                                         ! Cut cell

               CUT_CELL_AT(IJK) = .TRUE.
               BLOCKED_CELL_AT(IJK) = .FALSE.
               STANDARD_CELL_AT(IJK) = .FALSE.

               DO NODE = N_N1,N_N2
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     IF(.NOT.SNAP(IJK_OF_NODE(NODE))) THEN
                        NUMBER_OF_NODES(IJK) = NUMBER_OF_NODES(IJK) + 1
                        CONNECTIVITY(IJK,NUMBER_OF_NODES(IJK)) = IJK_OF_NODE(NODE)
                     ENDIF
                  ENDIF
               END DO

               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'SCALAR',NUMBER_OF_NODES(IJK),CONNECTIVITY,&
               X_NEW_POINT,Y_NEW_POINT,Z_NEW_POINT)

            ENDIF

         ENDIF        ! Interior cell

      ENDIF
      END DO          ! IJK Loop


      call SEND_RECEIVE_1D_LOGICAL(SMALL_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(STANDARD_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_CELL_AT,2)

! Consolidate blocked cells
      DO IJK = IJKSTART3, IJKEND3
         IF(BLOCKED_CELL_AT(IJK)) THEN
            STANDARD_CELL_AT(IJK) = .FALSE.
            CUT_CELL_AT(IJK)      = .FALSE.

            FLAG(IJK)             = 100

            VOL(IJK)              = ZERO
            AXY(IJK)              = ZERO
            AXZ(IJK)              = ZERO
            AYZ(IJK)              = ZERO

            AXY(BOTTOM_OF(IJK))   = ZERO
            AXZ(SOUTH_OF(IJK))    = ZERO
            AYZ(WEST_OF(IJK))     = ZERO
         ENDIF
      ENDDO

      call SEND_RECEIVE_1D_LOGICAL(SMALL_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(STANDARD_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_CELL_AT,2)

      call send_recv(FLAG,2)
      call send_recv(VOL,2)
      call send_recv(AXY,2)
      call send_recv(AXZ,2)
      call send_recv(AYZ,2)

!      print*,'Before removing duplicate nodes:'
      DO IJK = IJKSTART3, IJKEND3,-1
         IF(CUT_CELL_AT(IJK)) THEN
            print*,'==========================================================='
            print*,'IJK,  I,J=',IJK,I_OF(IJK),J_OF(IJK)
            print*,'NUMBER_OF_NODES=',NUMBER_OF_NODES(IJK)
            DO NODE = 1,NUMBER_OF_NODES(IJK)
               IF(CONNECTIVITY(IJK,NODE)>IJKEND3) THEN
                  print*,'CNCT=',NODE,CONNECTIVITY(IJK,NODE), &
                       X_NEW_POINT(CONNECTIVITY(IJK,NODE)-IJKEND3),Y_NEW_POINT(CONNECTIVITY(IJK,NODE)-IJKEND3)
               ELSE
                  print*,'CNCT=',NODE,CONNECTIVITY(IJK,NODE)
               ENDIF
            ENDDO
            print*,''
         ENDIF
      ENDDO



      OLD_CONNECTIVITY = CONNECTIVITY
!      X_OLD_POINT      = X_NEW_POINT
!      Y_OLD_POINT      = Y_NEW_POINT
!      Z_OLD_POINT      = Z_NEW_POINT


!======================================================================
!  Removing duplicate new points
!  Up to here, a cut face node that is shared among neighbor cells had
!  a different index, and appeared as a duplicate node.
!  The list will be updated and duplicate nodes will be removed from
!  the connectivity list
!  This is done to ensure that we can assign a volume surrounding each
!  node, that will be used to compute solids volume fraction for MPPIC
!======================================================================

!      print*,'Removing duplicate nodes...'

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_CELL_AT(IJK)) THEN          ! for each cut cell, identify neibhor cells that are also cut cells
                                             ! Look east and north in 2D, and also Top in 3D

         Diagonal = DSQRT(DX(I_OF(IJK))**2 + DY(J_OF(IJK))**2 + DZ(K_OF(IJK))**2)

            NB(1) = IP_OF(IJK)
            NB(2) = JP_OF(IJK)

            IF(NO_K) THEN
               N_NB = 2
            ELSE
               N_NB = 3
               NB(3) = KP_OF(IJK)
            ENDIF


            DO NN = 1,N_NB
               IF(CUT_CELL_AT(NB(NN))) THEN   ! For two neighbor cut cells, compare each cut-face node to remove duplicates

                  IJK_NB = NB(NN)

!                  print*,'comparing:',IJK,' and',IJK_NB

                  DO NODE = 1,NUMBER_OF_NODES(IJK)
                     IF(OLD_CONNECTIVITY(IJK,NODE)>IJKEND3) THEN  ! node belongs to the cut-face

                        X1 = X_NEW_POINT(OLD_CONNECTIVITY(IJK,NODE)-IJKEND3)
                        Y1 = Y_NEW_POINT(OLD_CONNECTIVITY(IJK,NODE)-IJKEND3)
                        Z1 = Z_NEW_POINT(OLD_CONNECTIVITY(IJK,NODE)-IJKEND3)


                        DO NODE_NB = 1,NUMBER_OF_NODES(IJK_NB)
                           IF(OLD_CONNECTIVITY(IJK_NB,NODE_NB)>IJKEND3) THEN  ! node belongs to the cut-face

                              X2 = X_NEW_POINT(OLD_CONNECTIVITY(IJK_NB,NODE_NB)-IJKEND3)
                              Y2 = Y_NEW_POINT(OLD_CONNECTIVITY(IJK_NB,NODE_NB)-IJKEND3)
                              Z2 = Z_NEW_POINT(OLD_CONNECTIVITY(IJK_NB,NODE_NB)-IJKEND3)

                              ! compare coordinates of cut-face nodes
                              D = (X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2

                              ! Duplicate nodes have identical coordinates (within tolerance TOL_MERGE times diagonal)
                              IF(D<TOL_MERGE*Diagonal) THEN ! keep the smallest node ID
                                 !  print*,'DULICATE NODES:', &
                                 ! NODE,NODE_NB,OLD_CONNECTIVITY(IJK,NODE),OLD_CONNECTIVITY(IJK_NB,NODE_NB)
                                 NC = MIN(OLD_CONNECTIVITY(IJK,NODE),OLD_CONNECTIVITY(IJK_NB,NODE_NB))
                                 CONNECTIVITY(IJK   ,NODE   ) = NC
                                 CONNECTIVITY(IJK_NB,NODE_NB) = NC

                              ENDIF
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDDO

               ENDIF
            ENDDO

         ENDIF
      ENDDO

      ALLOCATE(SCALAR_NODE_XYZ(DIMENSION_3 + NUMBER_OF_NEW_POINTS,3))
      ALLOCATE(Ovol_around_node(DIMENSION_3 + NUMBER_OF_NEW_POINTS))
      ALLOCATE(SCALAR_NODE_ATWALL(DIMENSION_3 + NUMBER_OF_NEW_POINTS))

      !first fill with standard nodes
      DO IJK = IJKSTART3, IJKEND3
         SCALAR_NODE_XYZ(IJK,1:3)  = SCALAR_NODE_XYZ_TEMP(IJK,1:3)
      ENDDO

      !now fill with cut-face nodes
      DO IJK = 1,NUMBER_OF_NEW_POINTS
         SCALAR_NODE_XYZ(IJKEND3+IJK,1) = X_NEW_POINT(IJK)
         SCALAR_NODE_XYZ(IJKEND3+IJK,2) = Y_NEW_POINT(IJK)
         SCALAR_NODE_XYZ(IJKEND3+IJK,3) = Z_NEW_POINT(IJK)
      ENDDO



      SCALAR_NODE_ATWALL(:)  = .true.
      !Rahul:
      !One could either set all scalar_node_atwall to false and
      !then set to true the nodes that are found as on the wall or outside
      !the doamin
      !In some situations, a node can be deemed to be both inside
      !and outside the domain depending upon which cell you look from.
      !Think of a small cell. According to the small cells, all the nodes
      !are outside the domain, but from the perspective or surrounding
      !cut-cells, some of the nodes of this small cell are within the domain.

      !so I'm setting all the points as being outside the domain.
      !if a point is ever found to be in the domain, then it will stay away
      !and further tests will not be able to revert it.
      IJKLOOP: DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IF(.not. IS_ON_myPE_wobnd(I,J,K)) cycle

         I1 = I-1
         I2 = I
         J1 = J-1
         J2 = J

         IF(NO_K) THEN
            K1 = K
            K2 = K
         ELSE
            K1 = K-1
            K2 = K
         ENDIF
         !Convention used to number node numbers is described below

         ! i=1, j=2           i=2, j=2
         !   _____________________
         !   |                   |
         !   |  I = 2, J = 2     |
         !   |___________________|
         ! i=1, j=1           i=2, j=1

         !Let's say the scalar cell with I = 2 and J = 2, i.e., the
         !first scalar cell in either direction.
         !then this scalar cell's node numbering in x- direction
         !will go from 1 to 2 and likewise in other directions.
         DO KK = K1, K2
            DO JJ = J1, J2
               DO II = I1, I2
                  IJK2 = funijk(II, JJ, KK)
                  COND_1 = .false.
                  COND_2 = .false.
                  !if it was already found to be inside the domain, then don't bother
                  IF(.not.SCALAR_NODE_ATWALL(IJK2))  CYCLE

                  IF(.not.FLUID_AT(IJK)) THEN
                     !do nothing
                  ELSE
                     !IJK is a fluid scalar cell. Check if it is
                     !a cut-cell or not
                     IF(CUT_CELL_AT(IJK)) THEN
                        CALL GET_DEL_H_DES(IJK,'SCALAR', &
                        & SCALAR_NODE_XYZ(IJK2,1), SCALAR_NODE_XYZ(IJK2,2), &
                        & SCALAR_NODE_XYZ(IJK2,3), &
                        & DIST, NORM1, NORM2, NORM3, .true.)
                        IF(DIST.GE.ZERO) THEN
                           SCALAR_NODE_ATWALL(IJK2)  = .false.
                           COND_1 =  .true.
                        ENDIF
                     ELSE
                        SCALAR_NODE_ATWALL(IJK2)  = .false.
                        COND_2 = .true.
                     ENDIF
                  ENDIF
                  !if(II == 1) then
                  ! write(*,'(10x, A, 4(2x,i10),3(2x,L2))') &
                  ! 'I1,J1, I, J, ATWALL, COND1, COND2 =  ', II, JJ, I, J,  SCALAR_NODE_ATWALL(IJK2), COND_1, COND_2
                  !endif
               ENDDO
            ENDDO
         ENDDO

      ENDDO IJKLOOP

!      print*,'After removing duplicate nodes:'
      DO IJK = IJKSTART3, IJKEND3,-1
         print*,'==========================================================='
         print*,'IJK,  I,J=',IJK,I_OF(IJK),J_OF(IJK)
         print*,'NUMBER_OF_NODES=',NUMBER_OF_NODES(IJK)
         DO NODE = 1,NUMBER_OF_NODES(IJK)
            print*,'CNCT=',NODE,CONNECTIVITY(IJK,NODE), &
                 SCALAR_NODE_XYZ(CONNECTIVITY(IJK,NODE),1),SCALAR_NODE_XYZ(CONNECTIVITY(IJK,NODE),2)
         ENDDO
         print*,''
      ENDDO



      Ovol_around_node = ZERO !UNDEFINED

      TOT_VOL_CELLS = ZERO

      DO IJK = IJKSTART3, IJKEND3

         DO NODE = 1,NUMBER_OF_NODES(IJK)
            NC = CONNECTIVITY(IJK,NODE)
            Ovol_around_node(NC) = Ovol_around_node(NC) + VOL(IJK)/NUMBER_OF_NODES(IJK)
         ENDDO

!         print*,'IJK,VOL=',IJK,VOL(IJK)


         IF(INTERIOR_CELL_AT(IJK)) TOT_VOL_CELLS = TOT_VOL_CELLS + VOL(IJK)


      ENDDO


      TOT_VOL_NODE= ZERO

      DO IJK = IJKSTART3, IJKEND3 + NUMBER_OF_NEW_POINTS    ! Loop over all nodes

         IF(Ovol_around_node(IJK)>ZERO) THEN
!             print*,'NODE,VOL=',IJK,Ovol_around_node(IJK)
             TOT_VOL_NODE = TOT_VOL_NODE + Ovol_around_node(IJK)
             Ovol_around_node(IJK) = ONE / Ovol_around_node(IJK)    ! Store One/volume
         ENDIF


      ENDDO

      IF(TOT_VOL_CELLS>ZERO) THEN
         IF(DABS(ONE-DABS(TOT_VOL_NODE/TOT_VOL_CELLS))>1.0D-6) THEN
            WRITE(*,*)'==========================================================='
            WRITE(*,*)'Total volume around nodes=',TOT_VOL_NODE
            WRITE(*,*)'Total volume =',TOT_VOL_CELLS
            WRITE(*,*)'The two volumes above should be the same.'
            WRITE(*,*)'==========================================================='
            CALL MFIX_EXIT(MYPE)
         ENDIF
      ENDIF

      deallocate(X_OLD_POINT)
      deallocate(Y_OLD_POINT)
      deallocate(Z_OLD_POINT)
      deallocate(SCALAR_NODE_XYZ_TEMP)

      RETURN
      END SUBROUTINE SET_3D_CUT_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_U_CELL_FLAGS                                C
!  Purpose: Set flags for U-Momentum cut cells, based on intersection  C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_U_CELL_FLAGS

      USE compar, ONLY: mype, pe_io, ijkstart3, ijkend3
      USE cut_cell_preproc, ONLY: cad_intersect, clean_intersect, clean_intersect_scalar, eval_f, intersect
      USE cutcell
      USE geometry, ONLY: no_k, axy_u, ayz_u, vol_u, axz_u
      USE param1, ONLY: half, one, zero
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: tol_f

      IMPLICIT NONE
      INTEGER :: IJK
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,N_N1,N_N2,Q_ID
      INTEGER :: MIN_INTERSECTIONS,MAX_INTERSECTIONS
      LOGICAL :: CLIP_FLAG
      INTEGER :: BCID

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'INTERSECTING GEOMETRY WITH U-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)
!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.
      TOL_SNAP = ZERO

      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
!             IF(POTENTIAL_CUT_CELL_AT(IJK))  CALL INTERSECT(IJK,'U_MOMENTUM',Xn_U_int(IJK),Ye_U_int(IJK),Zt_U_int(IJK))
            CALL INTERSECT(IJK,'U_MOMENTUM',Xn_U_int(IJK),Ye_U_int(IJK),Zt_U_int(IJK))
         END DO

!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================
         DO IJK = IJKSTART3, IJKEND3
            IF(INTERIOR_CELL_AT(IJK)) THEN
!               IF(POTENTIAL_CUT_CELL_AT(IJK))  CALL CLEAN_INTERSECT(IJK,'U_MOMENTUM',Xn_U_int(IJK),Ye_U_int(IJK),Zt_U_int(IJK))
               CALL CLEAN_INTERSECT(IJK,'U_MOMENTUM',Xn_U_int(IJK),Ye_U_int(IJK),Zt_U_int(IJK))
            ENDIF
         END DO

      ELSE
         CALL CAD_INTERSECT('U_MOMENTUM',Xn_U_int,Ye_U_int,Zt_U_int)
      ENDIF


      NUMBER_OF_NEW_U_POINTS = 0

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'U_MOMENTUM')

!======================================================================
!  Initialize location of velocity nodes at center of E,W,T faces
!======================================================================

                  X_U_ec(IJK) = X_NODE(8)
                  Y_U_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_U_nc(IJK) = Y_NODE(8)
                  Z_U_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_U_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_tc(IJK) = Z_NODE(8)

!======================================================================
!  Get Connectivity
!======================================================================

            IF(NO_K) THEN
               MIN_INTERSECTIONS = 2
               MAX_INTERSECTIONS = 2
               N_N1 = 5
               N_N2 = 8

            ELSE
               MIN_INTERSECTIONS = 3
               MAX_INTERSECTIONS = 6
               N_N1 = 1
               N_N2 = 8
            ENDIF

            CALL GET_CONNECTIVITY(IJK,'U_MOMENTUM',NUMBER_OF_NEW_U_POINTS,NUMBER_OF_U_NODES(IJK),CONNECTIVITY_U,&
            X_NEW_U_POINT,Y_NEW_U_POINT,Z_NEW_U_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_U_int,Ye_U_int,Zt_U_int)


            IF(TOTAL_NUMBER_OF_INTERSECTIONS < MIN_INTERSECTIONS ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('U_MOMENTUM',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(TOTAL_NUMBER_OF_INTERSECTIONS==-1) F_NODE(0) = -ONE
               BC_U_ID(IJK) = 0

               IF(F_NODE(0) <= ZERO) THEN
                  BLOCKED_U_CELL_AT(IJK) = .FALSE.
                  STANDARD_U_CELL_AT(IJK) = .TRUE.          ! Regular fluid cell

                  X_U_ec(IJK) = X_NODE(8)
                  Y_U_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_U_nc(IJK) = Y_NODE(8)
                  Z_U_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_U_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_tc(IJK) = Z_NODE(8)

               ELSE
                  BLOCKED_U_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
                  STANDARD_U_CELL_AT(IJK) = .FALSE.
                  AXY_U(IJK) = ZERO
                  AXZ_U(IJK) = ZERO
                  AYZ_U(IJK) = ZERO
                  VOL_U(IJK) = ZERO
               ENDIF

               IF(NO_K) THEN
                  NUMBER_OF_U_NODES(IJK) = 4
                  CONNECTIVITY_U(IJK,1) = IJK_OF_NODE(5)
                  CONNECTIVITY_U(IJK,2) = IJK_OF_NODE(6)
                  CONNECTIVITY_U(IJK,3) = IJK_OF_NODE(8)
                  CONNECTIVITY_U(IJK,4) = IJK_OF_NODE(7)
               ELSE
                  NUMBER_OF_U_NODES(IJK) = 8
                  DO NODE = N_N1,N_N2
                     CONNECTIVITY_U(IJK,NODE) = IJK_OF_NODE(NODE)
                  END DO
               ENDIF

            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > MAX_INTERSECTIONS) THEN

               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN U-CELL IJK = ',IJK
                  WRITE(*,*) 'MAXIMUM NUMBER OF INTERSECTIONS = ',MAX_INTERSECTIONS
                  WRITE(*,*) 'TOTAL NUMBER OF INTERSECTIONS = ',TOTAL_NUMBER_OF_INTERSECTIONS
                  WRITE(*,*) 'REMOVING U-CELL FROM COMPUTATION.'
!                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
!                  CALL MFIX_EXIT(MYPE)

               ENDIF

               BLOCKED_U_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
               STANDARD_U_CELL_AT(IJK) = .FALSE.
               CUT_U_CELL_AT(IJK) = .FALSE.
               AXY_U(IJK) = ZERO
               AXZ_U(IJK) = ZERO
               AYZ_U(IJK) = ZERO
               VOL_U(IJK) = ZERO

            ELSE                                         ! Cut cell

               CUT_U_CELL_AT(IJK) = .TRUE.
               BLOCKED_U_CELL_AT(IJK) = .FALSE.
               STANDARD_U_CELL_AT(IJK) = .FALSE.

               DO NODE = N_N1,N_N2
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     NUMBER_OF_U_NODES(IJK) = NUMBER_OF_U_NODES(IJK) + 1
                     CONNECTIVITY_U(IJK,NUMBER_OF_U_NODES(IJK)) = IJK_OF_NODE(NODE)
                  ENDIF
               END DO

               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'U_MOMENTUM',NUMBER_OF_U_NODES(IJK),CONNECTIVITY_U,&
               X_NEW_U_POINT,Y_NEW_U_POINT,Z_NEW_U_POINT)

            ENDIF

         ENDIF

      END DO

      DO IJK = IJKSTART3, IJKEND3

         CALL GET_CELL_NODE_COORDINATES(IJK,'U_MOMENTUM')

         DELX_Ue(IJK) = X_NODE(8) - X_U(IJK)
         DELX_Uw(IJK) = X_U(IJK) - X_NODE(1)

         DELY_Un(IJK) = Y_NODE(8) - Y_U(IJK)
         DELY_Us(IJK) = Y_U(IJK) - Y_NODE(1)

         DELZ_Ut(IJK) = Z_NODE(8) - Z_U(IJK)
         DELZ_Ub(IJK) = Z_U(IJK) - Z_NODE(1)

      ENDDO

      RETURN

      END SUBROUTINE SET_3D_CUT_U_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_V_CELL_FLAGS                                C
!  Purpose: Set flags for V-Momentum cut cells, based on intersection  C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_V_CELL_FLAGS

      USE compar, ONLY: mype, pe_io, ijkstart3, ijkend3
      USE cut_cell_preproc, ONLY: cad_intersect, clean_intersect, clean_intersect_scalar, eval_f, intersect
      USE cutcell
      USE geometry, ONLY: no_k, axy_v, axz_v, ayz_v, vol_v
      USE param1, ONLY: half, one, zero
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: tol_f

      IMPLICIT NONE
      INTEGER :: IJK
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,N_N1,N_N2,Q_ID
      INTEGER :: MIN_INTERSECTIONS,MAX_INTERSECTIONS
      LOGICAL :: CLIP_FLAG
      INTEGER :: BCID

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)'INTERSECTING GEOMETRY WITH V-MOMENTUM CELLS...'
      ENDIF
!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.

      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
!             IF(POTENTIAL_CUT_CELL_AT(IJK))  CALL INTERSECT(IJK,'V_MOMENTUM',Xn_V_int(IJK),Ye_V_int(IJK),Zt_V_int(IJK))
            CALL INTERSECT(IJK,'V_MOMENTUM',Xn_V_int(IJK),Ye_V_int(IJK),Zt_V_int(IJK))
         END DO

!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================
         DO IJK = IJKSTART3, IJKEND3
            IF(INTERIOR_CELL_AT(IJK)) THEN
!                IF(POTENTIAL_CUT_CELL_AT(IJK)) CALL CLEAN_INTERSECT(IJK,'V_MOMENTUM',Xn_V_int(IJK),Ye_V_int(IJK),Zt_V_int(IJK))
               CALL CLEAN_INTERSECT(IJK,'V_MOMENTUM',Xn_V_int(IJK),Ye_V_int(IJK),Zt_V_int(IJK))
            ENDIF
         END DO

      ELSE
         CALL CAD_INTERSECT('V_MOMENTUM',Xn_V_int,Ye_V_int,Zt_V_int)
      ENDIF


      NUMBER_OF_NEW_V_POINTS = 0

      DO IJK = IJKSTART3, IJKEND3
         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'V_MOMENTUM')

!======================================================================
!  Initialize location of velocity nodes at center of E,W,T faces
!======================================================================

            X_V_ec(IJK) = X_NODE(8)
            Y_V_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_V_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_V_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_V_nc(IJK) = Y_NODE(8)
            Z_V_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_V_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_V_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_V_tc(IJK) = Z_NODE(8)

!======================================================================
!  Get Connectivity
!======================================================================

            IF(NO_K) THEN
               MIN_INTERSECTIONS = 2
               MAX_INTERSECTIONS = 2
               N_N1 = 5
               N_N2 = 8

            ELSE
               MIN_INTERSECTIONS = 3
               MAX_INTERSECTIONS = 6
               N_N1 = 1
               N_N2 = 8
            ENDIF


            CALL GET_CONNECTIVITY(IJK,'V_MOMENTUM',NUMBER_OF_NEW_V_POINTS,NUMBER_OF_V_NODES(IJK),CONNECTIVITY_V,&
            X_NEW_V_POINT,Y_NEW_V_POINT,Z_NEW_V_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_V_int,Ye_V_int,Zt_V_int)


            IF(TOTAL_NUMBER_OF_INTERSECTIONS < MIN_INTERSECTIONS ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('V_MOMENTUM',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(TOTAL_NUMBER_OF_INTERSECTIONS==-1) F_NODE(0) = -ONE
               BC_V_ID(IJK) = 0

               IF(F_NODE(0) <= ZERO) THEN
                  BLOCKED_V_CELL_AT(IJK) = .FALSE.
                  STANDARD_V_CELL_AT(IJK) = .TRUE.          ! Regular fluid cell

                  X_V_ec(IJK) = X_NODE(8)
                  Y_V_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_V_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_V_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_V_nc(IJK) = Y_NODE(8)
                  Z_V_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_V_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_V_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_V_tc(IJK) = Z_NODE(8)

               ELSE
                  BLOCKED_V_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
                  STANDARD_V_CELL_AT(IJK) = .FALSE.
                  AXY_V(IJK) = ZERO
                  AXZ_V(IJK) = ZERO
                  AYZ_V(IJK) = ZERO
                  VOL_V(IJK) = ZERO
               ENDIF

               IF(NO_K) THEN
                  NUMBER_OF_V_NODES(IJK) = 4
                  CONNECTIVITY_V(IJK,1) = IJK_OF_NODE(5)
                  CONNECTIVITY_V(IJK,2) = IJK_OF_NODE(6)
                  CONNECTIVITY_V(IJK,3) = IJK_OF_NODE(8)
                  CONNECTIVITY_V(IJK,4) = IJK_OF_NODE(7)
               ELSE
                  NUMBER_OF_V_NODES(IJK) = 8
                  DO NODE = N_N1,N_N2
                     CONNECTIVITY_V(IJK,NODE) = IJK_OF_NODE(NODE)
                  END DO
               ENDIF

            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > MAX_INTERSECTIONS ) THEN

               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN V-CELL IJK = ',IJK
                  WRITE(*,*) 'MAXIMUM NUMBER OF INTERSECTIONS = ',MAX_INTERSECTIONS
                  WRITE(*,*) 'TOTAL NUMBER OF INTERSECTIONS = ',TOTAL_NUMBER_OF_INTERSECTIONS
                  WRITE(*,*) 'REMOVING V-CELL FROM COMPUTATION.'
!                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
!                  CALL MFIX_EXIT(MYPE)

               ENDIF

               BLOCKED_V_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
               STANDARD_V_CELL_AT(IJK) = .FALSE.
               CUT_V_CELL_AT(IJK) = .FALSE.
               AXY_V(IJK) = ZERO
               AXZ_V(IJK) = ZERO
               AYZ_V(IJK) = ZERO
               VOL_V(IJK) = ZERO

            ELSE                                         ! Cut cell

               CUT_V_CELL_AT(IJK) = .TRUE.
               BLOCKED_V_CELL_AT(IJK) = .FALSE.
               STANDARD_V_CELL_AT(IJK) = .FALSE.

               DO NODE = N_N1,N_N2
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     NUMBER_OF_V_NODES(IJK) = NUMBER_OF_V_NODES(IJK) + 1
                     CONNECTIVITY_V(IJK,NUMBER_OF_V_NODES(IJK)) = IJK_OF_NODE(NODE)
                  ENDIF
               END DO

               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'V_MOMENTUM',NUMBER_OF_V_NODES(IJK),CONNECTIVITY_V,&
               X_NEW_V_POINT,Y_NEW_V_POINT,Z_NEW_V_POINT)

            ENDIF

         ENDIF
      END DO

      DO IJK = IJKSTART3, IJKEND3

         CALL GET_CELL_NODE_COORDINATES(IJK,'V_MOMENTUM')

         DELX_Ve(IJK) = X_NODE(8) - X_V(IJK)
         DELX_Vw(IJK) = X_V(IJK) - X_NODE(1)

         DELY_Vn(IJK) = Y_NODE(8) - Y_V(IJK)
         DELY_Vs(IJK) = Y_V(IJK) - Y_NODE(1)

         DELZ_Vt(IJK) = Z_NODE(8) - Z_V(IJK)
         DELZ_Vb(IJK) = Z_V(IJK) - Z_NODE(1)

      ENDDO

      RETURN

      END SUBROUTINE SET_3D_CUT_V_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_W_CELL_FLAGS                                C
!  Purpose: Set flags for W-Momentum cut cells, based on intersection  C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_W_CELL_FLAGS

      USE compar, ONLY: mype, pe_io, ijkstart3, ijkend3
      USE cut_cell_preproc, ONLY: cad_intersect, clean_intersect, clean_intersect_scalar, eval_f, intersect
      USE cutcell
      USE geometry, ONLY: axy_w, axz_w, ayz_w, vol_w
      USE param1, ONLY: half, one, zero
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: tol_f

      IMPLICIT NONE
      INTEGER :: IJK
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,Q_ID
      LOGICAL :: CLIP_FLAG
      INTEGER :: BCID

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'INTERSECTING GEOMETRY WITH W-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)
!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.

      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
!             IF(POTENTIAL_CUT_CELL_AT(IJK)) CALL INTERSECT(IJK,'W_MOMENTUM',Xn_W_int(IJK),Ye_W_int(IJK),Zt_W_int(IJK))
            CALL INTERSECT(IJK,'W_MOMENTUM',Xn_W_int(IJK),Ye_W_int(IJK),Zt_W_int(IJK))
         END DO

!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================
         DO IJK = IJKSTART3, IJKEND3
            IF(INTERIOR_CELL_AT(IJK)) THEN
!               IF(POTENTIAL_CUT_CELL_AT(IJK)) CALL CLEAN_INTERSECT(IJK,'W_MOMENTUM',Xn_W_int(IJK),Ye_W_int(IJK),Zt_W_int(IJK))
               CALL CLEAN_INTERSECT(IJK,'W_MOMENTUM',Xn_W_int(IJK),Ye_W_int(IJK),Zt_W_int(IJK))
            ENDIF
         END DO

      ELSE
         CALL CAD_INTERSECT('W_MOMENTUM',Xn_W_int,Ye_W_int,Zt_W_int)
      ENDIF


      NUMBER_OF_NEW_W_POINTS = 0

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'W_MOMENTUM')

!======================================================================
!  Initialize location of velocity nodes at center of E,W,T faces
!======================================================================

            X_W_ec(IJK) = X_NODE(8)
            Y_W_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_W_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_W_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_W_nc(IJK) = Y_NODE(8)
            Z_W_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_W_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_W_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_W_tc(IJK) = Z_NODE(8)

!======================================================================
!  Get Connectivity
!======================================================================

            CALL GET_CONNECTIVITY(IJK,'W_MOMENTUM',NUMBER_OF_NEW_W_POINTS,NUMBER_OF_W_NODES(IJK),CONNECTIVITY_W,&
            X_NEW_W_POINT,Y_NEW_W_POINT,Z_NEW_W_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_W_int,Ye_W_int,Zt_W_int)


            IF(TOTAL_NUMBER_OF_INTERSECTIONS < 3 ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('W_MOMENTUM',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(TOTAL_NUMBER_OF_INTERSECTIONS==-1) F_NODE(0) = -ONE
               BC_W_ID(IJK) = 0

               IF(F_NODE(0) <= ZERO) THEN
                  BLOCKED_W_CELL_AT(IJK) = .FALSE.
                  STANDARD_W_CELL_AT(IJK) = .TRUE.          ! Regular fluid cell

                  X_W_ec(IJK) = X_NODE(8)
                  Y_W_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_W_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_W_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_W_nc(IJK) = Y_NODE(8)
                  Z_W_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_W_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_W_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_W_tc(IJK) = Z_NODE(8)

               ELSE
                  BLOCKED_W_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
                  STANDARD_W_CELL_AT(IJK) = .FALSE.
                  AXY_W(IJK) = ZERO
                  AXZ_W(IJK) = ZERO
                  AYZ_W(IJK) = ZERO
                  VOL_W(IJK) = ZERO
               ENDIF

               NUMBER_OF_W_NODES(IJK) = 8
               DO NODE = 1,8
                  CONNECTIVITY_W(IJK,NODE) = IJK_OF_NODE(NODE)
               END DO

            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > 6 ) THEN

               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN W-CELL IJK = ',IJK
                  WRITE(*,*) 'MAXIMUM NUMBER OF INTERSECTIONS = ',6
                  WRITE(*,*) 'TOTAL NUMBER OF INTERSECTIONS = ',TOTAL_NUMBER_OF_INTERSECTIONS
                  WRITE(*,*) 'REMOVING W-CELL FROM COMPUTATION.'
!                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
!                  CALL MFIX_EXIT(MYPE)

               ENDIF

               BLOCKED_W_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
               STANDARD_W_CELL_AT(IJK) = .FALSE.
               CUT_W_CELL_AT(IJK) = .FALSE.
               AXY_W(IJK) = ZERO
               AXZ_W(IJK) = ZERO
               AYZ_W(IJK) = ZERO
               VOL_W(IJK) = ZERO

            ELSE                                         ! Cut cell

               CUT_W_CELL_AT(IJK) = .TRUE.

               DO NODE = 1,8
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     NUMBER_OF_W_NODES(IJK) = NUMBER_OF_W_NODES(IJK) + 1
                     CONNECTIVITY_W(IJK,NUMBER_OF_W_NODES(IJK)) = IJK_OF_NODE(NODE)
                  ENDIF
               END DO


               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'W_MOMENTUM',NUMBER_OF_W_NODES(IJK),CONNECTIVITY_W,&
               X_NEW_W_POINT,Y_NEW_W_POINT,Z_NEW_W_POINT)

            ENDIF

         ENDIF
      END DO

      DO IJK = IJKSTART3, IJKEND3

         CALL GET_CELL_NODE_COORDINATES(IJK,'W_MOMENTUM')

         DELX_We(IJK) = X_NODE(8) - X_W(IJK)
         DELX_Ww(IJK) = X_W(IJK) - X_NODE(1)

         DELY_Wn(IJK) = Y_NODE(8) - Y_W(IJK)
         DELY_Ws(IJK) = Y_W(IJK) - Y_NODE(1)

         DELZ_Wt(IJK) = Z_NODE(8) - Z_W(IJK)
         DELZ_Wb(IJK) = Z_W(IJK) - Z_NODE(1)

      ENDDO

      RETURN

      END SUBROUTINE SET_3D_CUT_W_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_CELL_TREATMENT_FLAGS                        C
!  Purpose: Set flags for scalar cut cells, based on intersection      C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_CELL_TREATMENT_FLAGS

      USE compar, only: mype, pe_io, ijkstart3, ijkend3
      USE cutcell
      USE functions, ONLY: funijk
      USE geometry, ONLY: do_k
      USE indices, ONLY: i_of, j_of, k_of

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,IM,IP,JM,JP,KM,KP
      INTEGER :: IMJK,IPJK,IJMK,IJPK,IJKM,IJKP

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)'SETTING CUT CELL TREATMENT FLAGS...'
      ENDIF
!======================================================================
!  Set flags identifying cells requiring cut cell treatment:
!  These are the cut cells and their neighbours
!======================================================================

      CUT_TREATMENT_AT   = .FALSE.
      CUT_U_TREATMENT_AT = .FALSE.
      CUT_V_TREATMENT_AT = .FALSE.
      CUT_W_TREATMENT_AT = .FALSE.

      call SEND_RECEIVE_1D_LOGICAL(CUT_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_U_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_V_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_W_CELL_AT,2)

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IP = I + 1
            JP = J + 1
            KP = K + 1

            IM = I - 1
            JM = J - 1
            KM = K - 1

            IMJK   = FUNIJK(IM,J,K)
            IJMK   = FUNIJK(I,JM,K)
            IJKM   = FUNIJK(I,J,KM)

            IPJK   = FUNIJK(IP,J,K)
            IJPK   = FUNIJK(I,JP,K)
            IJKP   = FUNIJK(I,J,KP)

            CUT_TREATMENT_AT(IJK) = (CUT_CELL_AT(IJK ).OR.   &
                                     CUT_CELL_AT(IMJK).OR.   &
                                     CUT_CELL_AT(IPJK).OR.   &
                                     CUT_CELL_AT(IJMK).OR.   &
                                     CUT_CELL_AT(IJPK))

            CUT_U_TREATMENT_AT(IJK) = (CUT_U_CELL_AT(IJK ).OR.   &
                                       CUT_U_CELL_AT(IMJK).OR.   &
                                       CUT_U_CELL_AT(IPJK).OR.   &
                                       CUT_U_CELL_AT(IJMK).OR.   &
                                       CUT_U_CELL_AT(IJPK))

            CUT_V_TREATMENT_AT(IJK) = (CUT_V_CELL_AT(IJK ).OR.   &
                                       CUT_V_CELL_AT(IMJK).OR.   &
                                       CUT_V_CELL_AT(IPJK).OR.   &
                                       CUT_V_CELL_AT(IJMK).OR.   &
                                       CUT_V_CELL_AT(IJPK))

            IF(DO_K) THEN

               CUT_TREATMENT_AT(IJK) = (CUT_TREATMENT_AT(IJK).OR.   &
                                            CUT_CELL_AT(IJKM).OR.   &
                                            CUT_CELL_AT(IJKP))

               CUT_U_TREATMENT_AT(IJK) = (CUT_U_TREATMENT_AT(IJK).OR.   &
                                              CUT_U_CELL_AT(IJKM).OR.   &
                                              CUT_U_CELL_AT(IJKP))

               CUT_V_TREATMENT_AT(IJK) = (CUT_V_TREATMENT_AT(IJK).OR.   &
                                              CUT_V_CELL_AT(IJKM).OR.   &
                                              CUT_V_CELL_AT(IJKP))

               CUT_W_TREATMENT_AT(IJK) = (CUT_W_CELL_AT(IJK ).OR.   &
                                          CUT_W_CELL_AT(IMJK).OR.   &
                                          CUT_W_CELL_AT(IPJK).OR.   &
                                          CUT_W_CELL_AT(IJMK).OR.   &
                                          CUT_W_CELL_AT(IJPK).OR.   &
                                          CUT_W_CELL_AT(IJKM).OR.   &
                                          CUT_W_CELL_AT(IJKP))

            ENDIF

         ENDIF

         IF(.NOT.RE_INDEXING) THEN
            IF(BLOCKED_CELL_AT(IJK)) CUT_TREATMENT_AT(IJK) = .FALSE.
            IF(BLOCKED_U_CELL_AT(IJK)) CUT_U_TREATMENT_AT(IJK) = .FALSE.
            IF(BLOCKED_V_CELL_AT(IJK)) CUT_V_TREATMENT_AT(IJK) = .FALSE.
            IF(BLOCKED_W_CELL_AT(IJK)) CUT_W_TREATMENT_AT(IJK) = .FALSE.
         ENDIF

      END DO

      IF(CG_SAFE_MODE(1)==1) CUT_TREATMENT_AT   = .FALSE.
      IF(CG_SAFE_MODE(3)==1) CUT_U_TREATMENT_AT = .FALSE.
      IF(CG_SAFE_MODE(4)==1) CUT_V_TREATMENT_AT = .FALSE.
      IF(CG_SAFE_MODE(5)==1) CUT_W_TREATMENT_AT = .FALSE.

      RETURN

      END SUBROUTINE SET_3D_CUT_CELL_TREATMENT_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_GHOST_CELL_FLAGS                                   C
!  Purpose: Set flags for ghost cell flags                             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_GHOST_CELL_FLAGS

      USE compar, ONLY: ijkstart3, ijkend3
      USE compar, ONLY: iend1, iend2, iend3, jstart1, jstart2, jstart3, kstart1, kstart2, kstart3
      USE compar, only: istart1, istart2, istart3, jend1, jend2, jend3, kend1, kend2, kend3
      USE cutcell
      USE functions, ONLY: funijk
      USE functions, ONLY: WEST_OF, EAST_OF, SOUTH_OF, NORTH_OF, BOTTOM_OF, TOP_OF
      USE functions, ONLY: IM_OF, IP_OF, JM_OF, JP_OF, KM_OF, KP_OF
      USE geometry, ONLY: imax1, imin1, jmax1, jmin1, kmax1, kmin1, vol, vol_u, vol_v, vol_w, axy, axz, ayz, ayz_u, ayz_v, ayz_w, axy_u, axy_v, axy_w, axz_u, axz_v, axz_w, flag, do_k, flag_e, flag_n, flag_t

      USE bc
      USE sendrecv

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,I23,J23,K23
      INTEGER :: BCID
      INTEGER :: IPJK,IMJK,IJPK,IJMK,IJKP,IJKM

!     EAST BOUNDARY
      I = IEND1

      IF(I==IMAX1) THEN

         DO I23 = IEND2,IEND3
            DO K = KSTART3,KEND3
               DO J = JSTART3, JEND3

                  IJK  = FUNIJK(I,J,K)
                  IPJK = FUNIJK(I23,J,K)

                  BLOCKED_CELL_AT(IPJK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IPJK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IPJK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IPJK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IPJK)) FLAG(IPJK) = 100

                  VOL(IPJK)   = VOL(IJK)
                  VOL_U(IPJK) = VOL_U(IJK)
                  VOL_V(IPJK) = VOL_V(IJK)
                  VOL_W(IPJK) = VOL_W(IJK)

                  AYZ(IPJK)   = AYZ(IJK)
                  AYZ_U(IPJK) = AYZ_U(IJK)
                  AYZ_V(IPJK) = AYZ_V(IJK)
                  AYZ_W(IPJK) = AYZ_W(IJK)

                  AXY(IPJK)   = AXY(IJK)
                  AXY_U(IPJK) = AXY_U(IJK)
                  AXY_V(IPJK) = AXY_V(IJK)
                  AXY_W(IPJK) = AXY_W(IJK)

                  AXZ(IPJK)   = AXZ(IJK)
                  AXZ_U(IPJK) = AXZ_U(IJK)
                  AXZ_V(IPJK) = AXZ_V(IJK)
                  AXZ_W(IPJK) = AXZ_W(IJK)

                  BC_ID(IPJK)   = BC_ID(IJK)
                  BC_U_ID(IPJK) = BC_U_ID(IJK)
                  BC_V_ID(IPJK) = BC_V_ID(IJK)
                  BC_W_ID(IPJK) = BC_W_ID(IJK)

               END DO
            END DO
         END DO

      ENDIF

!     WEST BOUNDARY
      I = ISTART1

      IF(I==IMIN1) THEN

         DO I23 = ISTART3,ISTART2
            DO K = KSTART3,KEND3
               DO J = JSTART3, JEND3

                  IJK  = FUNIJK(I,J,K)
                  IMJK = FUNIJK(I23,J,K)

                  BLOCKED_CELL_AT(IMJK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IMJK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IMJK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IMJK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IMJK)) FLAG(IMJK) = 100

                  VOL(IMJK)   = VOL(IJK)
                  VOL_U(IMJK) = VOL_U(IJK)
                  VOL_V(IMJK) = VOL_V(IJK)
                  VOL_W(IMJK) = VOL_W(IJK)

                  AYZ(IMJK)   = AYZ(IJK)
                  AYZ_U(IMJK) = AYZ_U(IJK)
                  AYZ_V(IMJK) = AYZ_V(IJK)
                  AYZ_W(IMJK) = AYZ_W(IJK)

                  AXY(IMJK)   = AXY(IJK)
                  AXY_U(IMJK) = AXY_U(IJK)
                  AXY_V(IMJK) = AXY_V(IJK)
                  AXY_W(IMJK) = AXY_W(IJK)

                  AXZ(IMJK)   = AXZ(IJK)
                  AXZ_U(IMJK) = AXZ_U(IJK)
                  AXZ_V(IMJK) = AXZ_V(IJK)
                  AXZ_W(IMJK) = AXZ_W(IJK)

                  BC_ID(IMJK)   = BC_ID(IJK)
                  BC_U_ID(IMJK) = BC_U_ID(IJK)
                  BC_V_ID(IMJK) = BC_V_ID(IJK)
                  BC_W_ID(IMJK) = BC_W_ID(IJK)

               END DO
            END DO
         END DO

      ENDIF

!     NORTH BOUNDARY
      J = JEND1

      IF(J==JMAX1) THEN

         DO J23 = JEND2,JEND3
            DO K = KSTART3,KEND3
               DO I = ISTART3, IEND3

                  IJK  = FUNIJK(I,J,K)
                  IJPK = FUNIJK(I,J23,K)

                  BLOCKED_CELL_AT(IJPK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IJPK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IJPK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IJPK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IJPK)) FLAG(IJPK) = 100

                  VOL(IJPK)   = VOL(IJK)
                  VOL_U(IJPK) = VOL_U(IJK)
                  VOL_V(IJPK) = VOL_V(IJK)
                  VOL_W(IJPK) = VOL_W(IJK)

                  AYZ(IJPK)   = AYZ(IJK)
                  AYZ_U(IJPK) = AYZ_U(IJK)
                  AYZ_V(IJPK) = AYZ_V(IJK)
                  AYZ_W(IJPK) = AYZ_W(IJK)

                  AXY(IJPK)   = AXY(IJK)
                  AXY_U(IJPK) = AXY_U(IJK)
                  AXY_V(IJPK) = AXY_V(IJK)
                  AXY_W(IJPK) = AXY_W(IJK)

                  AXZ(IJPK)   = AXZ(IJK)
                  AXZ_U(IJPK) = AXZ_U(IJK)
                  AXZ_V(IJPK) = AXZ_V(IJK)
                  AXZ_W(IJPK) = AXZ_W(IJK)

                  BC_ID(IJPK)   = BC_ID(IJK)
                  BC_U_ID(IJPK) = BC_U_ID(IJK)
                  BC_V_ID(IJPK) = BC_V_ID(IJK)
                  BC_W_ID(IJPK) = BC_W_ID(IJK)

               END DO
            END DO
         END DO

      ENDIF
!     SOUTH BOUNDARY
      J = JSTART1

      IF(J==JMIN1) THEN

         DO J23 = JSTART3,JSTART2
            DO K = KSTART3,KEND3
               DO I = ISTART3, IEND3

                  IJK  = FUNIJK(I,J,K)
                  IJMK = FUNIJK(I,J23,K)

                  BLOCKED_CELL_AT(IJMK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IJMK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IJMK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IJMK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IJMK)) FLAG(IJMK) = 100

                  VOL(IJMK)   = VOL(IJK)
                  VOL_U(IJMK) = VOL_U(IJK)
                  VOL_V(IJMK) = VOL_V(IJK)
                  VOL_W(IJMK) = VOL_W(IJK)

                  AYZ(IJMK)   = AYZ(IJK)
                  AYZ_U(IJMK) = AYZ_U(IJK)
                  AYZ_V(IJMK) = AYZ_V(IJK)
                  AYZ_W(IJMK) = AYZ_W(IJK)

                  AXY(IJMK)   = AXY(IJK)
                  AXY_U(IJMK) = AXY_U(IJK)
                  AXY_V(IJMK) = AXY_V(IJK)
                  AXY_W(IJMK) = AXY_W(IJK)

                  AXZ(IJMK)   = AXZ(IJK)
                  AXZ_U(IJMK) = AXZ_U(IJK)
                  AXZ_V(IJMK) = AXZ_V(IJK)
                  AXZ_W(IJMK) = AXZ_W(IJK)

                  BC_ID(IJMK)   = BC_ID(IJK)
                  BC_U_ID(IJMK) = BC_U_ID(IJK)
                  BC_V_ID(IJMK) = BC_V_ID(IJK)
                  BC_W_ID(IJMK) = BC_W_ID(IJK)

               END DO
            END DO
         END DO

      ENDIF

      IF(DO_K) THEN

!        TOP BOUNDARY
         K = KEND1

         IF(K==KMAX1) THEN

            DO K23=KEND2,KEND3
               DO J = JSTART3,JEND3
                  DO I = ISTART3, IEND3

                     IJK  = FUNIJK(I,J,K)
                     IJKP = FUNIJK(I,J,K23)

                     BLOCKED_CELL_AT(IJKP) = BLOCKED_CELL_AT(IJK)
                     BLOCKED_U_CELL_AT(IJKP) = BLOCKED_U_CELL_AT(IJK)
                     BLOCKED_V_CELL_AT(IJKP) = BLOCKED_V_CELL_AT(IJK)
                     BLOCKED_W_CELL_AT(IJKP) = BLOCKED_W_CELL_AT(IJK)

                     IF(BLOCKED_CELL_AT(IJKP)) FLAG(IJKP) = 100

                     VOL(IJKP)   = VOL(IJK)
                     VOL_U(IJKP) = VOL_U(IJK)
                     VOL_V(IJKP) = VOL_V(IJK)
                     VOL_W(IJKP) = VOL_W(IJK)

                     AYZ(IJKP)   = AYZ(IJK)
                     AYZ_U(IJKP) = AYZ_U(IJK)
                     AYZ_V(IJKP) = AYZ_V(IJK)
                     AYZ_W(IJKP) = AYZ_W(IJK)

                     AXY(IJKP)   = AXY(IJK)
                     AXY_U(IJKP) = AXY_U(IJK)
                     AXY_V(IJKP) = AXY_V(IJK)
                     AXY_W(IJKP) = AXY_W(IJK)

                     AXZ(IJKP)   = AXZ(IJK)
                     AXZ_U(IJKP) = AXZ_U(IJK)
                     AXZ_V(IJKP) = AXZ_V(IJK)
                     AXZ_W(IJKP) = AXZ_W(IJK)

                     BC_ID(IJKP)   = BC_ID(IJK)
                     BC_U_ID(IJKP) = BC_U_ID(IJK)
                     BC_V_ID(IJKP) = BC_V_ID(IJK)
                     BC_W_ID(IJKP) = BC_W_ID(IJK)

                  END DO
               END DO
            END DO

         ENDIF

!        BOTTOM BOUNDARY
         K = KSTART1

         IF(K==KMIN1) THEN

            DO K23 = KSTART3,KSTART2
               DO J = JSTART3,JEND3
                  DO I = ISTART3, IEND3

                     IJK  = FUNIJK(I,J,K)
                     IJKM = FUNIJK(I,J,K23)

                     BLOCKED_CELL_AT(IJKM) = BLOCKED_CELL_AT(IJK)
                     BLOCKED_U_CELL_AT(IJKM) = BLOCKED_U_CELL_AT(IJK)
                     BLOCKED_V_CELL_AT(IJKM) = BLOCKED_V_CELL_AT(IJK)
                     BLOCKED_W_CELL_AT(IJKM) = BLOCKED_W_CELL_AT(IJK)

                     IF(BLOCKED_CELL_AT(IJKM)) FLAG(IJKM) = 100

                     VOL(IJKM)   = VOL(IJK)
                     VOL_U(IJKM) = VOL_U(IJK)
                     VOL_V(IJKM) = VOL_V(IJK)
                     VOL_W(IJKM) = VOL_W(IJK)

                     AYZ(IJKM)   = AYZ(IJK)
                     AYZ_U(IJKM) = AYZ_U(IJK)
                     AYZ_V(IJKM) = AYZ_V(IJK)
                     AYZ_W(IJKM) = AYZ_W(IJK)

                     AXY(IJKM)   = AXY(IJK)
                     AXY_U(IJKM) = AXY_U(IJK)
                     AXY_V(IJKM) = AXY_V(IJK)
                     AXY_W(IJKM) = AXY_W(IJK)

                     AXZ(IJKM)   = AXZ(IJK)
                     AXZ_U(IJKM) = AXZ_U(IJK)
                     AXZ_V(IJKM) = AXZ_V(IJK)
                     AXZ_W(IJKM) = AXZ_W(IJK)

                     BC_ID(IJKM)   = BC_ID(IJK)
                     BC_U_ID(IJKM) = BC_U_ID(IJK)
                     BC_V_ID(IJKM) = BC_V_ID(IJK)
                     BC_W_ID(IJKM) = BC_W_ID(IJK)

                  END DO
               END DO
            END DO

         ENDIF

      ENDIF

      DO IJK = IJKSTART3, IJKEND3

         IF(BLOCKED_U_CELL_AT(IJK)) FLAG_E(IJK)=100
         IF(BLOCKED_V_CELL_AT(IJK)) FLAG_N(IJK)=100
         IF(BLOCKED_W_CELL_AT(IJK)) FLAG_T(IJK)=100

      ENDDO

!     BLOCKED CELLS WERE ASSIGNED THE FLAG 100 DURING PRE_PROCESSING
!     THIS IS INCORRECT FOR FREEE-SLIP AND PARTIAL-SLIP
!     THE FLAG IS OVERWRITTEN HERE TO ACCOUNT FOR NSW,FSW AND PSW BCs

      DO IJK = IJKSTART3, IJKEND3

         IF(FLAG(IJK)==100) THEN

            IMJK   = IM_OF(IJK)
            IJMK   = JM_OF(IJK)
            IJKM   = KM_OF(IJK)

            IPJK   = IP_OF(IJK)
            IJPK   = JP_OF(IJK)
            IJKP   = KP_OF(IJK)

            IF(CUT_CELL_AT(IMJK)) THEN
               BCID=BC_ID(IMJK)
               IF(BCID>0) THEN
                  IF(IS_NSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 100
                  IF(IS_FSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 101
                  IF(IS_PSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 102
               ENDIF

            ELSEIF(CUT_CELL_AT(IPJK)) THEN
               BCID=BC_ID(IPJK)
               IF(BCID>0) THEN
                  IF(IS_NSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 100
                  IF(IS_FSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 101
                  IF(IS_PSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 102
               ENDIF

            ELSEIF(CUT_CELL_AT(IJMK)) THEN
               BCID=BC_ID(IJMK)
               IF(BCID>0) THEN
                  IF(IS_NSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 100
                  IF(IS_FSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 101
                  IF(IS_PSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 102
               ENDIF

            ELSEIF(CUT_CELL_AT(IJPK)) THEN
               BCID=BC_ID(IJPK)
               IF(BCID>0) THEN
                  IF(IS_NSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 100
                  IF(IS_FSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 101
                  IF(IS_PSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 102
               ENDIF

            ELSEIF(CUT_CELL_AT(IJKM)) THEN
               BCID=BC_ID(IJKM)
               IF(BCID>0) THEN
                  IF(IS_NSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 100
                  IF(IS_FSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 101
                  IF(IS_PSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 102
               ENDIF

            ELSEIF(CUT_CELL_AT(IJKP)) THEN
               BCID=BC_ID(IJKP)
               IF(BCID>0) THEN
                  IF(IS_NSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 100
                  IF(IS_FSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 101
                  IF(IS_PSW(BC_TYPE_ENUM(BCID))) FLAG(IJK) = 102
               ENDIF

            ENDIF
         ENDIF
      ENDDO

      call send_recv(FLAG,2)
      RETURN

      END SUBROUTINE SET_GHOST_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_POTENTIAL_CUT_CELLS                                  C
!  Purpose: Set flags for scalar cut cells, based on intersection      C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_POTENTIAL_CUT_CELLS

      USE compar, only: ijkstart3, ijkend3, mype, pe_io
      USE cutcell
      USE functions, only: funijk, bottom_of, south_of, west_of
      USE geometry, ONLY: dx, dy, dz, do_k, imin3, imax3, jmin3, jmax3, kmin3, kmax3, flag, axy, axz, ayz, vol, no_k
      USE indices, only: i_of, j_of, k_of
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: tol_f
      USE param, ONLY: DIMENSION_3
      USE param1, ONLY: HALF, ZERO

      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,II,JJ,KK
      INTEGER :: I1,I2,J1,J2,K1,K2
      INTEGER :: NODE,N_N1,N_N2
      LOGICAL :: ALL_NEGATIVE,ALL_POSITIVE,CLIP_FLAG
      INTEGER :: Q_ID,BCID

      INTEGER :: IJK_NB,NUMBER_OF_POTENTIAL_CUT_CELLS

      DOUBLE PRECISION :: xc,yc,zc,fc

      LOGICAL, DIMENSION(DIMENSION_3) ::POSITIVE_F_AT

      POTENTIAL_CUT_CELL_AT=.TRUE.

      RETURN  ! This subroutine is currently disabled

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'ESTIMATING POTENTIAL SCALAR CUT CELLS...'
      ENDIF
10    FORMAT(1X,A)

!======================================================================
!  Evaluate f at cell center and store where f>0
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         xc = XG_E(I) - HALF*DX(I)
         yc = YG_N(J) - HALF*DY(J)

         IF(DO_K) THEN
            zc = ZG_T(K) - HALF*DZ(K)
         ELSE
            zc = zero
         ENDIF

          Q_ID = 1
          CALL EVAL_F('QUADRIC',xc,yc,zc,Q_ID,fc,CLIP_FLAG)

          CALL EVAL_F('POLYGON',xc,yc,zc,N_POLYGON,fc,CLIP_FLAG)

          CALL EVAL_F('USR_DEF',xc,yc,zc,N_USR_DEF,Fc,CLIP_FLAG)

          X_NODE(15) = xc
          Y_NODE(15) = yc
          Z_NODE(15) = zc
          CALL EVAL_STL_FCT_AT('SCALAR',IJK,15,fc,CLIP_FLAG,BCID)

          IF(fc>TOL_F) THEN
             POSITIVE_F_AT(IJK)=.TRUE.
          ELSE
             POSITIVE_F_AT(IJK)=.FALSE.
          ENDIF

       ENDDO

      DO IJK = IJKSTART3, IJKEND3

         IF(INTERIOR_CELL_AT(IJK)) THEN

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            I1 = MAX(I - 2,IMIN3)
            I2 = MIN(I + 2,IMAX3)
            J1 = MAX(J - 2,JMIN3)
            J2 = MIN(J + 2,JMAX3)

            IF(DO_K) THEN
               K1 = MAX(K - 2,KMIN3)
               K2 = MIN(K + 2,KMAX3)
            ELSE
               K1=K
               K2=K
            ENDIF

            IF(POSITIVE_F_AT(IJK)) THEN
               ALL_POSITIVE=.TRUE.
               DO KK=K1,K2
                  DO JJ=J1,J2
                     DO II=I1,I2
                        IJK_NB = FUNIJK(II,JJ,KK)
                        IF(.NOT.POSITIVE_F_AT(IJK_NB)) THEN
                           ALL_POSITIVE=.FALSE.
                        ENDIF
                     ENDDO
                  ENDDO
                ENDDO

                IF(ALL_POSITIVE) THEN
                   POTENTIAL_CUT_CELL_AT(IJK)=.FALSE.
                   FLAG(IJK) = 100
                   BLOCKED_CELL_AT(IJK) = .TRUE.               ! Blocked fluid cell
                   STANDARD_CELL_AT(IJK) = .FALSE.
                   AXY(IJK) = ZERO
                   AXZ(IJK) = ZERO
                   AYZ(IJK) = ZERO
                   VOL(IJK) = ZERO

                   AXY(BOTTOM_OF(IJK)) = ZERO
                   AXZ(SOUTH_OF(IJK)) = ZERO
                   AYZ(WEST_OF(IJK)) = ZERO
                ENDIF

            ELSE
               ALL_NEGATIVE=.TRUE.
               DO KK=K1,K2
                  DO JJ=J1,J2
                     DO II=I1,I2
                        IJK_NB = FUNIJK(II,JJ,KK)
                        IF(POSITIVE_F_AT(IJK_NB)) THEN
                           ALL_NEGATIVE=.FALSE.
                        ENDIF
                     ENDDO
                  ENDDO
                ENDDO

                IF(ALL_NEGATIVE) THEN
                   POTENTIAL_CUT_CELL_AT(IJK)=.FALSE.
                   BLOCKED_CELL_AT(IJK) = .FALSE.
                   STANDARD_CELL_AT(IJK) = .TRUE.           ! Regular fluid cell
                ENDIF

            ENDIF

            IF((FLAG(IJK)>=100).AND.(FLAG(IJK)<=102)) THEN
               POTENTIAL_CUT_CELL_AT(IJK)=.FALSE.
               BLOCKED_CELL_AT(IJK) = .TRUE.
               STANDARD_CELL_AT(IJK) = .FALSE.          ! Blocked cell = wall cell
            ENDIF

         ENDIF

      END DO

      NUMBER_OF_POTENTIAL_CUT_CELLS = 0

      IF(NO_K) THEN
         N_N1 = 5
         N_N2 = 8
      ELSE
         N_N1 = 1
         N_N2 = 8
      ENDIF

      DO IJK=IJKSTART3,IJKEND3
         IF(POTENTIAL_CUT_CELL_AT(IJK)) THEN
            NUMBER_OF_POTENTIAL_CUT_CELLS = NUMBER_OF_POTENTIAL_CUT_CELLS + 1
         ELSE
            CALL GET_CELL_NODE_COORDINATES(IJK,'SCALAR')

            IF(NO_K) THEN
               NUMBER_OF_NODES(IJK) = 4
               CONNECTIVITY(IJK,1) = IJK_OF_NODE(5)
               CONNECTIVITY(IJK,2) = IJK_OF_NODE(6)
               CONNECTIVITY(IJK,3) = IJK_OF_NODE(8)
               CONNECTIVITY(IJK,4) = IJK_OF_NODE(7)
            ELSE
               NUMBER_OF_NODES(IJK) = 8
               DO NODE = N_N1,N_N2
                  CONNECTIVITY(IJK,NODE) = IJK_OF_NODE(NODE)
               END DO
            ENDIF

            X_U(IJK) = X_NODE(8)
            Y_U(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_U(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_V(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_V(IJK) = Y_NODE(8)
            Z_V(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_W(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_W(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_W(IJK) = Z_NODE(8)

         ENDIF
      ENDDO

!      call SEND_RECEIVE_1D_LOGICAL(SNAP,2)
      IF(MyPE == PE_IO) THEN
         WRITE(*,*)'DONE ESTIMATING POTENTIAL SCALAR CUT CELLS.',NUMBER_OF_POTENTIAL_CUT_CELLS,IJKEND3
      ENDIF

      RETURN
      END SUBROUTINE GET_POTENTIAL_CUT_CELLS
