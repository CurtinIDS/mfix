
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_CUT_CELL_VOLUME_AND_AREAS                          C
!  Purpose: Compute volume and face areas of a cut cell                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_CUT_CELL_VOLUME_AND_AREAS(IJK,TYPE_OF_CELL,N_NODES,CONNECT,X_NP,Y_NP,Z_NP)

      USE bc
      USE compar, ONLY: ijkend3, istart1, jstart1, kstart1, mype
      USE cut_cell_preproc, only: eval_f
      USE cutcell
      USE exit, only: mfix_exit
      USE functions, only: funijk
      USE geometry, ONLY: do_k, no_k, ayz, axz, axy, ayz_u, axy_u, axz_u, ayz_v, axz_v, axy_v, ayz_w, axz_w, axy_w, vol, vol_u, vol_v, vol_w, dx, dy, dz, zlength, flag, flag_n, flag_t
      USE indices, only: i_of, j_of, k_of
      USE param, only: dimension_3
      USE param1, only: undefined
      USE param1, only: zero
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: n_quadric, tol_f, bc_id_q, dim_quadric
      USE stl, ONLY: n_facet_at, list_facet_at, bc_id_stl_face

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: I,J,K,IJK,L,NODE
      INTEGER :: IM,JM,KM,IMJK,IJMK,IJKM
      INTEGER :: NN,N_N1,N_N2,Q_ID
      INTEGER :: N_CUT_FACE_NODES
      INTEGER :: N_EAST_FACE_NODES,N_NORTH_FACE_NODES,N_TOP_FACE_NODES
      INTEGER :: N_WEST_FACE_NODES,N_SOUTH_FACE_NODES,N_BOTTOM_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_CUT_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_EAST_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_NORTH_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_TOP_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_WEST_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_SOUTH_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_BOTTOM_FACE_NODES
      DOUBLE PRECISION :: X_COPY,Y_COPY,Z_COPY,F_COPY,F_Q,F_MIN
      DOUBLE PRECISION :: X_MEAN,Y_MEAN,Z_MEAN
      DOUBLE PRECISION :: AREA_EAST,AREA_NORTH,AREA_TOP
      DOUBLE PRECISION :: AREA_WEST,AREA_SOUTH,AREA_BOTTOM
      DOUBLE PRECISION :: CUT_AREA
      DOUBLE PRECISION, DIMENSION(3) :: CENTROID_EAST,CENTROID_NORTH,CENTROID_TOP
      DOUBLE PRECISION, DIMENSION(3) :: CENTROID_WEST,CENTROID_SOUTH,CENTROID_BOTTOM
      DOUBLE PRECISION, DIMENSION(3) :: CENTROID_CUT
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz,VOLUME,NORM
      DOUBLE PRECISION :: DIFFERENCE
      INTEGER :: N_NODES
      INTEGER :: NDIFF,SPACE
      INTEGER, DIMENSION(DIMENSION_3,15) :: CONNECT
      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) :: X_NP,Y_NP,Z_NP

      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: F_CUT_FACE

      INTEGER, DIMENSION(15) :: TEMP_CONNECTIVITY

      LOGICAL :: CLIP_FLAG,INSIDE_FACET
      INTEGER :: BCID,N_BC,QID_FMIN,BCID2,NF

      LOGICAL :: CORNER_POINT
      INTEGER :: NODE_OF_CORNER,IERROR

!======================================================================
!  Filter the connectivity to identify nodes belonging to
!  East, North, Top, and cut faces
!======================================================================

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      IM = I - 1
      JM = J - 1
      KM = K - 1

      IMJK   = FUNIJK(IM,J,K)
      IJMK   = FUNIJK(I,JM,K)
      IJKM   = FUNIJK(I,J,KM)

      N_EAST_FACE_NODES = 0
      N_NORTH_FACE_NODES = 0
      N_TOP_FACE_NODES = 0

      N_CUT_FACE_NODES = 0

      N_WEST_FACE_NODES = 0
      N_SOUTH_FACE_NODES = 0
      N_BOTTOM_FACE_NODES = 0

      X_MEAN = ZERO
      Y_MEAN = ZERO
      Z_MEAN = ZERO

      IF(NO_K) THEN
         N_N1 = 5
         N_N2 = 8
      ELSE
         N_N1 = 1
         N_N2 = 8
      ENDIF

      DO L = 1, N_NODES
         IF(CONNECT(IJK,L)>IJKEND3) THEN   ! One of the new point
            X_COPY = X_NP(CONNECT(IJK,L)-IJKEND3)
            Y_COPY = Y_NP(CONNECT(IJK,L)-IJKEND3)
            Z_COPY = Z_NP(CONNECT(IJK,L)-IJKEND3)
            CORNER_POINT = .FALSE.

         ELSE                                   ! An existing point
            CORNER_POINT = .TRUE.
            DO NODE = N_N1,N_N2
               IF(CONNECT(IJK,L) == IJK_OF_NODE(NODE)) THEN
                  NODE_OF_CORNER = NODE
                  X_COPY = X_NODE(NODE)
                  Y_COPY = Y_NODE(NODE)
                  Z_COPY = Z_NODE(NODE)

                  IF(SNAP(IJK_OF_NODE(NODE))) THEN
                     N_CUT_FACE_NODES = N_CUT_FACE_NODES + 1
                     COORD_CUT_FACE_NODES(1,N_CUT_FACE_NODES) = X_COPY
                     COORD_CUT_FACE_NODES(2,N_CUT_FACE_NODES) = Y_COPY
                     COORD_CUT_FACE_NODES(3,N_CUT_FACE_NODES) = Z_COPY

                     IF(N_QUADRIC>0) THEN
                        F_MIN = UNDEFINED
                        QID_FMIN = 0
                        DO Q_ID = 1, N_QUADRIC
                           CALL GET_F_QUADRIC(X_COPY,Y_COPY,Z_COPY,Q_ID,F_Q,CLIP_FLAG)
                           IF(DABS(F_Q) < F_MIN) THEN
                              F_MIN = DABS(F_Q)
                              QID_FMIN = Q_ID
                           ENDIF
                        ENDDO
                        BC_ID(IJK) = BC_ID_Q(QID_FMIN)
                     ENDIF

                 ENDIF

               ENDIF
            END DO
         ENDIF

         X_MEAN = X_MEAN + X_COPY
         Y_MEAN = Y_MEAN + Y_COPY
         Z_MEAN = Z_MEAN + Z_COPY

         IF(CORNER_POINT) THEN
            Q_ID = 1
            CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

            CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

            CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,NODE_OF_CORNER,F_COPY,CLIP_FLAG,BCID2)
         ELSE
            F_COPY = ZERO
         ENDIF

         IF (DABS(F_COPY) < TOL_F ) THEN ! belongs to cut face
            N_CUT_FACE_NODES = N_CUT_FACE_NODES + 1
            COORD_CUT_FACE_NODES(1,N_CUT_FACE_NODES) = X_COPY
            COORD_CUT_FACE_NODES(2,N_CUT_FACE_NODES) = Y_COPY
            COORD_CUT_FACE_NODES(3,N_CUT_FACE_NODES) = Z_COPY
         ENDIF

         IF (X_COPY > (X_NODE(8) - TOL_F) ) THEN ! belongs to East face
            N_EAST_FACE_NODES = N_EAST_FACE_NODES + 1
            COORD_EAST_FACE_NODES(1,N_EAST_FACE_NODES) = X_COPY
            COORD_EAST_FACE_NODES(2,N_EAST_FACE_NODES) = Y_COPY
            COORD_EAST_FACE_NODES(3,N_EAST_FACE_NODES) = Z_COPY
         ENDIF

         IF (Y_COPY > (Y_NODE(8) - TOL_F) ) THEN ! belongs to North face
            N_NORTH_FACE_NODES = N_NORTH_FACE_NODES + 1
            COORD_NORTH_FACE_NODES(1,N_NORTH_FACE_NODES) = X_COPY
            COORD_NORTH_FACE_NODES(2,N_NORTH_FACE_NODES) = Y_COPY
            COORD_NORTH_FACE_NODES(3,N_NORTH_FACE_NODES) = Z_COPY
         ENDIF

         IF (Z_COPY > (Z_NODE(8) - TOL_F) ) THEN ! belongs to Top face
            N_TOP_FACE_NODES = N_TOP_FACE_NODES + 1
            COORD_TOP_FACE_NODES(1,N_TOP_FACE_NODES) = X_COPY
            COORD_TOP_FACE_NODES(2,N_TOP_FACE_NODES) = Y_COPY
            COORD_TOP_FACE_NODES(3,N_TOP_FACE_NODES) = Z_COPY
         ENDIF

         IF(I>=ISTART1) THEN
            IF (X_COPY < (X_NODE(1) + TOL_F) ) THEN ! belongs to West face
               N_WEST_FACE_NODES = N_WEST_FACE_NODES + 1
               COORD_WEST_FACE_NODES(1,N_WEST_FACE_NODES) = X_COPY
               COORD_WEST_FACE_NODES(2,N_WEST_FACE_NODES) = Y_COPY
               COORD_WEST_FACE_NODES(3,N_WEST_FACE_NODES) = Z_COPY
            ENDIF
         ENDIF

         IF(J>=JSTART1) THEN
            IF (Y_COPY < (Y_NODE(1) + TOL_F) ) THEN ! belongs to South face
               N_SOUTH_FACE_NODES = N_SOUTH_FACE_NODES + 1
               COORD_SOUTH_FACE_NODES(1,N_SOUTH_FACE_NODES) = X_COPY
               COORD_SOUTH_FACE_NODES(2,N_SOUTH_FACE_NODES) = Y_COPY
               COORD_SOUTH_FACE_NODES(3,N_SOUTH_FACE_NODES) = Z_COPY
            ENDIF
         ENDIF

         IF(DO_K) THEN
            IF(K>=KSTART1) THEN
               IF (Z_COPY < (Z_NODE(1) + TOL_F) ) THEN ! belongs to Bottom face
                  N_BOTTOM_FACE_NODES = N_BOTTOM_FACE_NODES + 1
                  COORD_BOTTOM_FACE_NODES(1,N_BOTTOM_FACE_NODES) = X_COPY
                  COORD_BOTTOM_FACE_NODES(2,N_BOTTOM_FACE_NODES) = Y_COPY
                  COORD_BOTTOM_FACE_NODES(3,N_BOTTOM_FACE_NODES) = Z_COPY
               ENDIF
            ENDIF
         ENDIF

      END DO

      X_MEAN = X_MEAN / N_NODES
      Y_MEAN = Y_MEAN / N_NODES
      Z_MEAN = Z_MEAN / N_NODES


      IF( N_CUT_FACE_NODES == N_EAST_FACE_NODES) THEN
         DIFFERENCE = ZERO
         DO NDIFF = 1,N_CUT_FACE_NODES
            DO SPACE = 1,3
               DIFFERENCE = DIFFERENCE +DABS(COORD_CUT_FACE_NODES(SPACE,NDIFF) &
                                           - COORD_EAST_FACE_NODES(SPACE,NDIFF))
            END DO
         END DO
         IF( DIFFERENCE < TOL_F ) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'WARNING: CUT FACE IS IDENTICAL TO EAST FACE AT IJK = ',IJK
               WRITE(*,*)'TYPE OF CELL = ',TYPE_OF_CELL
               WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
               WRITE(*,*)'REMOVING EAST FACE.'
            ENDIF
            N_EAST_FACE_NODES = 0
         ENDIF
      ENDIF

      IF( N_CUT_FACE_NODES == N_WEST_FACE_NODES) THEN
         DIFFERENCE = ZERO
         DO NDIFF = 1,N_CUT_FACE_NODES
            DO SPACE = 1,3
               DIFFERENCE = DIFFERENCE +DABS(COORD_CUT_FACE_NODES(SPACE,NDIFF) &
                                           - COORD_WEST_FACE_NODES(SPACE,NDIFF))
            END DO
         END DO
         IF( DIFFERENCE < TOL_F ) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'WARNING: CUT FACE IS IDENTICAL TO WEST FACE AT IJK = ',IJK
               WRITE(*,*)'TYPE OF CELL = ',TYPE_OF_CELL
               WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
               WRITE(*,*)'REMOVING WEST FACE.'
            ENDIF
            N_WEST_FACE_NODES = 0
         ENDIF
      ENDIF

      IF( N_CUT_FACE_NODES == N_NORTH_FACE_NODES) THEN
         DIFFERENCE = ZERO
         DO NDIFF = 1,N_CUT_FACE_NODES
            DO SPACE = 1,3
               DIFFERENCE = DIFFERENCE +DABS(COORD_CUT_FACE_NODES(SPACE,NDIFF) &
                                           - COORD_NORTH_FACE_NODES(SPACE,NDIFF))
            END DO
         END DO
         IF( DIFFERENCE < TOL_F ) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'WARNING: CUT FACE IS IDENTICAL TO NORTH FACE AT IJK = ',IJK
               WRITE(*,*)'TYPE OF CELL = ',TYPE_OF_CELL
               WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
               WRITE(*,*)'REMOVING NORTH FACE.'
            ENDIF
            N_NORTH_FACE_NODES = 0
         ENDIF
      ENDIF

      IF( N_CUT_FACE_NODES == N_SOUTH_FACE_NODES) THEN
         DIFFERENCE = ZERO
         DO NDIFF = 1,N_CUT_FACE_NODES
            DO SPACE = 1,3
               DIFFERENCE = DIFFERENCE +DABS(COORD_CUT_FACE_NODES(SPACE,NDIFF) &
                                           - COORD_SOUTH_FACE_NODES(SPACE,NDIFF))
            END DO
         END DO
         IF( DIFFERENCE < TOL_F ) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'WARNING: CUT FACE IS IDENTICAL TO SOUTH FACE AT IJK = ',IJK
               WRITE(*,*)'TYPE OF CELL = ',TYPE_OF_CELL
               WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
               WRITE(*,*)'REMOVING SOUTH FACE.'
            ENDIF
            N_SOUTH_FACE_NODES = 0
         ENDIF
      ENDIF

      IF((NO_K).AND.(N_TOP_FACE_NODES <=2)) THEN
         IF(PRINT_WARNINGS) THEN
            WRITE(*,*)'WARNING: TOP FACE HAS ONLY TWO NODES AT IJK = ',IJK
            WRITE(*,*)'TYPE OF CELL = ',TYPE_OF_CELL
            WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
            WRITE(*,*)'REMOVING CUT CELL.'
         ENDIF

         SELECT CASE (TYPE_OF_CELL)
            CASE('SCALAR')
               BLOCKED_CELL_AT(IJK) = .TRUE.
            CASE('U_MOMENTUM')
               BLOCKED_U_CELL_AT(IJK) = .TRUE.
            CASE('V_MOMENTUM')
               BLOCKED_V_CELL_AT(IJK) = .TRUE.
            CASE('W_MOMENTUM')
               BLOCKED_W_CELL_AT(IJK) = .TRUE.
            CASE DEFAULT
               WRITE(*,*)'SUBROUTINE: GET_CUT_CELL_VOLUME_AND_AREAS'
               WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
               WRITE(*,*)'ACCEPTABLE TYPES ARE:'
               WRITE(*,*)'SCALAR'
               WRITE(*,*)'U_MOMENTUM'
               WRITE(*,*)'V_MOMENTUM'
               WRITE(*,*)'W_MOMENTUM'
               CALL MFIX_EXIT(myPE)
         END SELECT
         RETURN

      ELSE
         IF( N_CUT_FACE_NODES == N_TOP_FACE_NODES) THEN
            DIFFERENCE = ZERO
            DO NDIFF = 1,N_CUT_FACE_NODES
               DO SPACE = 1,3
                  DIFFERENCE = DIFFERENCE +DABS(COORD_CUT_FACE_NODES(SPACE,NDIFF) &
                                              - COORD_TOP_FACE_NODES(SPACE,NDIFF))
               END DO
            END DO
            IF( DIFFERENCE < TOL_F ) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'WARNING: CUT FACE IS IDENTICAL TO TOP FACE AT IJK = ',IJK
                  WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
                  WRITE(*,*)'REMOVING TOP FACE.'
               ENDIF
               N_TOP_FACE_NODES = 0
            ENDIF
         ENDIF

         IF( N_CUT_FACE_NODES == N_BOTTOM_FACE_NODES) THEN
            DIFFERENCE = ZERO
            DO NDIFF = 1,N_CUT_FACE_NODES
               DO SPACE = 1,3
                  DIFFERENCE = DIFFERENCE +DABS(COORD_CUT_FACE_NODES(SPACE,NDIFF) &
                                              - COORD_BOTTOM_FACE_NODES(SPACE,NDIFF))
               END DO
            END DO
            IF( DIFFERENCE < TOL_F ) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'WARNING: CUT FACE IS IDENTICAL TO BOTTOM FACE AT IJK = ',IJK
                  WRITE(*,*)'THIS USUALLY OCCURS WHEN QUADRIC SURFACE IS ALIGNED WITH GRID.'
                  WRITE(*,*)'REMOVING BOTTOM FACE.'
               ENDIF
               N_BOTTOM_FACE_NODES = 0
            ENDIF
         ENDIF
      ENDIF

      IERROR = 0

      CALL GET_POLYGON_AREA_AND_CENTROID(N_EAST_FACE_NODES,COORD_EAST_FACE_NODES,AREA_EAST,CENTROID_EAST,IERROR)

      CALL GET_POLYGON_AREA_AND_CENTROID(N_NORTH_FACE_NODES,COORD_NORTH_FACE_NODES,AREA_NORTH,CENTROID_NORTH,IERROR)

      CALL GET_POLYGON_AREA_AND_CENTROID(N_TOP_FACE_NODES,COORD_TOP_FACE_NODES,AREA_TOP,CENTROID_TOP,IERROR)

      CALL GET_POLYGON_AREA_AND_CENTROID(N_CUT_FACE_NODES,COORD_CUT_FACE_NODES,CUT_AREA,CENTROID_CUT,IERROR)

      CALL STORE_CUT_FACE_INFO(IJK,TYPE_OF_CELL,N_CUT_FACE_NODES,COORD_CUT_FACE_NODES,X_MEAN,Y_MEAN,Z_MEAN)

      CALL GET_DEL_H(IJK,TYPE_OF_CELL,X_MEAN,Y_MEAN,Z_MEAN,Del_H,Nx,Ny,Nz)

      IF(DEL_H == UNDEFINED) DEL_H = ZERO

      IF(I==ISTART1) CALL GET_POLYGON_AREA_AND_CENTROID(N_WEST_FACE_NODES,COORD_WEST_FACE_NODES,AREA_WEST,CENTROID_WEST,IERROR)

      IF(J==JSTART1) CALL GET_POLYGON_AREA_AND_CENTROID(N_SOUTH_FACE_NODES,COORD_SOUTH_FACE_NODES,AREA_SOUTH,CENTROID_SOUTH,IERROR)

      IF(DO_K.AND.K==KSTART1) CALL GET_POLYGON_AREA_AND_CENTROID(N_BOTTOM_FACE_NODES,COORD_BOTTOM_FACE_NODES,AREA_BOTTOM,&
      CENTROID_BOTTOM,IERROR)

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')
            IF(I>ISTART1) AREA_WEST   = AYZ(IMJK)
            IF(J>JSTART1) AREA_SOUTH  = AXZ(IJMK)
            IF(K>KSTART1) AREA_BOTTOM = AXY(IJKM)

            AYZ(IJK) = AREA_EAST
            AXZ(IJK) = AREA_NORTH
            AXY(IJK) = AREA_TOP
            Area_CUT(IJK) = CUT_AREA

! Compute normal and store cut face info (normal vector and reference point)
            NORMAL_S(IJK,1) = AREA_EAST  - AREA_WEST
            NORMAL_S(IJK,2) = AREA_NORTH - AREA_SOUTH

            IF(DO_K) THEN
               NORMAL_S(IJK,3) = AREA_TOP   - AREA_BOTTOM
            ELSE
               NORMAL_S(IJK,3) = ZERO
            ENDIF

            NORM = sqrt(NORMAL_S(IJK,1)**2 + NORMAL_S(IJK,2)**2 + NORMAL_S(IJK,3)**2)

            IF(NORM==ZERO) NORM = UNDEFINED

            NORMAL_S(IJK,:) = NORMAL_S(IJK,:) / NORM

            REFP_S(IJK,:)   = CENTROID_CUT(:)

            CALL GET_DEL_H(IJK,TYPE_OF_CELL,X_MEAN,Y_MEAN,Z_MEAN,Del_H,Nx,Ny,Nz)
            IF(DEL_H == UNDEFINED) DEL_H = ZERO
            DELH_Scalar(IJK) = DEL_H

            IF(NO_K) THEN
               VOL(IJK) = AXY(IJK) * ZLENGTH
            ELSE

!             The cell is divided into pyramids having the same apex (Mean value of nodes)
!             Cell Volume = sum of pyramid volums
               VOLUME =   AREA_EAST   * (X_NODE(8) - X_MEAN)                   &
                        + AREA_NORTH  * (Y_NODE(8) - Y_MEAN)                   &
                        + AREA_TOP    * (Z_NODE(8) - Z_MEAN)                   &
                        + AREA_WEST   * (X_MEAN - X_NODE(1))                   &
                        + AREA_SOUTH  * (Y_MEAN - Y_NODE(1))                   &
                        + AREA_BOTTOM * (Z_MEAN - Z_NODE(1))                   &
                        + CUT_AREA    * DEL_H

               VOL(IJK) = VOLUME / 3.0D0

            ENDIF

            IF(IERROR==1) print*,'scalar IERROR',IJK
            IF(IERROR==1) VOL(IJK) = ZERO

            IF(VOL(IJK)==ZERO) THEN

               IF(PRINT_WARNINGS) THEN
                  PRINT*,'WARNING: ZERO VOLUME CUT CELL DETECTED AT IJK =',IJK
                  PRINT*,'REMOVING CUT CELL FROM COMPUTATION (FLAGGED AS BLOCKED CELL)'
               ENDIF

               FLAG(IJK) = 100
               CUT_CELL_AT(IJK) = .FALSE.
               BLOCKED_CELL_AT(IJK) = .TRUE.
               STANDARD_CELL_AT(IJK) = .FALSE.
               AYZ(IJK) = ZERO
               AXZ(IJK) = ZERO
               AXY(IJK) = ZERO
               VOL(IJK) = ZERO

            ELSEIF(VOL(IJK) < TOL_SMALL_CELL * DX(I)*DY(J)*DZ(K) ) THEN
               SMALL_CELL_AT(IJK) = .TRUE.
               NUMBER_OF_SMALL_CELLS = NUMBER_OF_SMALL_CELLS + 1


               IF(PRINT_WARNINGS) THEN
                  PRINT*,'WARNING: SMALL SCALAR CELL CUT CELL DETECTED AT IJK =',IJK
                  PRINT*,'REMOVING CUT CELL FROM COMPUTATION (FLAGGED AS BLOCKED CELL)'
               ENDIF

               BLOCKED_CELL_AT(IJK) = .TRUE.

            ENDIF

            IF(AREA_EAST > TOL_SMALL_AREA*DY(J)*DZ(K)) THEN
               WALL_U_AT(IJK) = .FALSE.
               X_U(IJK) = CENTROID_EAST(1)
               Y_U(IJK) = CENTROID_EAST(2)
               Z_U(IJK) = CENTROID_EAST(3)
            ENDIF

            IF(AREA_NORTH > TOL_SMALL_AREA*DX(I)*DZ(K)) THEN
               WALL_V_AT(IJK) = .FALSE.
               X_V(IJK) = CENTROID_NORTH(1)
               Y_V(IJK) = CENTROID_NORTH(2)
               Z_V(IJK) = CENTROID_NORTH(3)
            ENDIF

            IF(AREA_TOP > TOL_SMALL_AREA*DX(I)*DY(J)) THEN
               WALL_W_AT(IJK) = .FALSE.
               X_W(IJK) = CENTROID_TOP(1)
               Y_W(IJK) = CENTROID_TOP(2)
               Z_W(IJK) = CENTROID_TOP(3)
            ENDIF

            IF(N_QUADRIC>0) THEN
               F_CUT_FACE = UNDEFINED
               N_BC = 0
!               BC_ID(IJK) = 0
               DO Q_ID = 1, N_QUADRIC
                  DO NODE = 1,N_CUT_FACE_NODES

                     CALL GET_F_QUADRIC(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                     Q_ID,F_CUT_FACE(Q_ID),CLIP_FLAG)

                     IF(DABS(F_CUT_FACE(Q_ID)) < TOL_F) THEN
                        IF(BC_ID(IJK)/=BC_ID_Q(Q_ID)) N_BC = N_BC + 1
                        BC_ID(IJK) = BC_ID_Q(Q_ID)
                     ENDIF
                  ENDDO
                  IF(BC_ID(IJK)>0) THEN
!                     IF(BC_TYPE_ENUM(BC_ID(IJK))  == CG_MI) EXIT
                     IF(BC_TYPE_ENUM(BC_ID(IJK))  == CG_NSW) EXIT
                  ENDIF
               ENDDO
            ENDIF

            IF(N_POLYGON>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_POLY_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  Q_ID,F_CUT_FACE(1),CLIP_FLAG,BCID)
                   IF(F_CUT_FACE(1)==ZERO) THEN
                     BC_ID(IJK) = BCID
                     IF(BC_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            IF(N_USR_DEF>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_USR_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  BCID,F_CUT_FACE(1),CLIP_FLAG)
                  IF(DABS(F_CUT_FACE(1)) < TOL_F) THEN
                     BC_ID(IJK) = BCID
                     IF(BC_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            IF(USE_STL.OR.USE_MSH) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  DO NN=1,N_FACET_AT(IJK)
                     NF=LIST_FACET_AT(IJK,NN)
                     CALL IS_POINT_INSIDE_FACET(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),&
                                                COORD_CUT_FACE_NODES(3,NODE),NF,INSIDE_FACET)
                     IF(INSIDE_FACET) THEN
                        BC_ID(IJK) = BC_ID_STL_FACE(NF)
                        IF(BC_ID(IJK)>0) THEN
                           IF(BC_TYPE_ENUM(BC_ID(IJK))  == CG_MI) EXIT
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

            IF (BC_ID(IJK)>0) THEN

               IF (.not. IS_CG(BC_TYPE_ENUM(BC_ID(IJK))))  THEN
                  WRITE(*,'(/1X,A)')'ERROR: INVALID BOUNDARY CONDITION ASSIGNED TO A CUT CELL.'
                  WRITE(*,'(/1X,A,I9,I6)')'IJK, MyPE = ', IJK,myPE
                  WRITE(*,'(1X,A,I6)')'BC_ID  = ',BC_ID(IJK)
                  WRITE(*,'(1X,A,A)')'BC_TYPE = ',BC_TYPE_ENUM(BC_ID(IJK))
                  WRITE(*,'(1X,A)')'VALID BC_TYPE FOR CUT CELLS ARE: CG_NSW, CG_FSW, CG_PSW, CG_MI, and CG_PO'
                  WRITE(*,'(/1X,A)')'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  call mfix_exit(myPE)
               ENDIF

            ENDIF

!            Reordering connectivity such that polygon is defined appropriately for 2D vtk file

            IF(NO_K) THEN
               DO NN = 1,N_TOP_FACE_NODES
                  TEMP_CONNECTIVITY(NN) = CONNECT(IJK,ORDER(NN))
               END DO
               CONNECT(IJK,:) = TEMP_CONNECTIVITY
             ENDIF

         CASE('U_MOMENTUM')

            IF(I>ISTART1) AREA_WEST   = AYZ_U(IMJK)
            IF(J>JSTART1) AREA_SOUTH  = AXZ_U(IJMK)
            IF(K>KSTART1) AREA_BOTTOM = AXY_U(IJKM)

            AYZ_U(IJK) = AREA_EAST
            AXZ_U(IJK) = AREA_NORTH
            AXY_U(IJK) = AREA_TOP
            Area_U_CUT(IJK) = CUT_AREA

            IF(AREA_EAST > ZERO) THEN
               X_U_ec(IJK) = CENTROID_EAST(1)
               Y_U_ec(IJK) = CENTROID_EAST(2)
               Z_U_ec(IJK) = CENTROID_EAST(3)
            ENDIF

            IF(AREA_NORTH > ZERO) THEN
               X_U_nc(IJK) = CENTROID_NORTH(1)
               Y_U_nc(IJK) = CENTROID_NORTH(2)
               Z_U_nc(IJK) = CENTROID_NORTH(3)
            ENDIF

            IF(AREA_TOP > ZERO) THEN
               X_U_tc(IJK) = CENTROID_TOP(1)
               Y_U_tc(IJK) = CENTROID_TOP(2)
               Z_U_tc(IJK) = CENTROID_TOP(3)
            ENDIF


! Compute normal and store cut face info (normal vector and reference point)
            NORMAL_U(IJK,1) = AREA_EAST  - AREA_WEST
            NORMAL_U(IJK,2) = AREA_NORTH - AREA_SOUTH

            IF(DO_K) THEN
               NORMAL_U(IJK,3) = AREA_TOP   - AREA_BOTTOM
            ELSE
               NORMAL_U(IJK,3) = ZERO
            ENDIF

            NORM = sqrt(NORMAL_U(IJK,1)**2 + NORMAL_U(IJK,2)**2 + NORMAL_U(IJK,3)**2)

            IF(NORM==ZERO) NORM = UNDEFINED

            NORMAL_U(IJK,:) = NORMAL_U(IJK,:) / NORM

            REFP_U(IJK,:)   = CENTROID_CUT(:)

            CALL GET_DEL_H(IJK,TYPE_OF_CELL,X_MEAN,Y_MEAN,Z_MEAN,Del_H,Nx,Ny,Nz)
            IF(DEL_H == UNDEFINED) DEL_H = ZERO

            IF(NO_K) THEN
               VOL_U(IJK) = AXY_U(IJK) * ZLENGTH
            ELSE

!             The cell is divided into pyramids having the same apex (Mean value of nodes)
!             Cell Volume = sum of pyramid volums
               VOLUME =   AREA_EAST   * (X_NODE(8) - X_MEAN)                   &
                        + AREA_NORTH  * (Y_NODE(8) - Y_MEAN)                   &
                        + AREA_TOP    * (Z_NODE(8) - Z_MEAN)                   &
                        + AREA_WEST   * (X_MEAN - X_NODE(1))                   &
                        + AREA_SOUTH  * (Y_MEAN - Y_NODE(1))                   &
                        + AREA_BOTTOM * (Z_MEAN - Z_NODE(1))                   &
                        + CUT_AREA    * DEL_H

               VOL_U(IJK) = VOLUME / 3.0D0
            ENDIF

            IF(IERROR==1) print*,'U_VEL IERROR',IJK
            IF(IERROR==1) AYZ(IJK) = ZERO

            IF(AYZ(IJK) <= TOL_SMALL_AREA*DY(J)*DZ(K)) THEN
               WALL_U_AT(IJK) = .TRUE.

               NUMBER_OF_U_WALL_CELLS = NUMBER_OF_U_WALL_CELLS + 1

               X_U(IJK) = CENTROID_CUT(1)
               Y_U(IJK) = CENTROID_CUT(2)
               Z_U(IJK) = CENTROID_CUT(3)
            ENDIF


            IF(N_QUADRIC>0) THEN
               F_CUT_FACE = UNDEFINED
               N_BC = 0
               BC_U_ID(IJK) = 0
               DO Q_ID = 1, N_QUADRIC
                  DO NODE = 1,N_CUT_FACE_NODES

                     CALL GET_F_QUADRIC(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                     Q_ID,F_CUT_FACE(Q_ID),CLIP_FLAG)

                     IF(DABS(F_CUT_FACE(Q_ID)) < TOL_F) THEN
                        IF(BC_U_ID(IJK)/=BC_ID_Q(Q_ID)) N_BC = N_BC + 1
                        BC_U_ID(IJK) = BC_ID_Q(Q_ID)
                     ENDIF
                  ENDDO
                  IF(BC_U_ID(IJK)>0) THEN
                     IF(BC_TYPE_ENUM(BC_U_ID(IJK))  == CG_MI) EXIT
!                     IF(BC_TYPE_ENUM(BC_U_ID(IJK))  == CG_NSW) EXIT
                  ENDIF

               ENDDO

            ENDIF

            IF(N_POLYGON>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_POLY_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  Q_ID,F_CUT_FACE(1),CLIP_FLAG,BCID)
                   IF(F_CUT_FACE(1)==ZERO) THEN
                     BC_U_ID(IJK) = BCID
                     IF(BC_U_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_U_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF


            IF(N_USR_DEF>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_USR_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  BCID,F_CUT_FACE(1),CLIP_FLAG)
                  IF(DABS(F_CUT_FACE(1)) < TOL_F) THEN
                     BC_U_ID(IJK) = BCID
                     IF(BC_U_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_U_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF


            IF(USE_STL.OR.USE_MSH) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  DO NN=1,N_FACET_AT(IJK)
                     NF=LIST_FACET_AT(IJK,NN)
                     CALL IS_POINT_INSIDE_FACET(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),&
                                                COORD_CUT_FACE_NODES(3,NODE),NF,INSIDE_FACET)
                     IF(INSIDE_FACET) THEN
                        BC_U_ID(IJK) = BC_ID_STL_FACE(NF)
                        IF(BC_U_ID(IJK)>0) THEN
                           IF(BC_TYPE_ENUM(BC_U_ID(IJK))  == CG_MI) EXIT
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

         CASE('V_MOMENTUM')

            IF(I>ISTART1) AREA_WEST   = AYZ_V(IMJK)
            IF(J>JSTART1) AREA_SOUTH  = AXZ_V(IJMK)
            IF(K>KSTART1) AREA_BOTTOM = AXY_V(IJKM)

            AYZ_V(IJK) = AREA_EAST
            AXZ_V(IJK) = AREA_NORTH
            AXY_V(IJK) = AREA_TOP
            Area_V_CUT(IJK) = CUT_AREA

            IF(AREA_EAST > ZERO) THEN
               X_V_ec(IJK) = CENTROID_EAST(1)
               Y_V_ec(IJK) = CENTROID_EAST(2)
               Z_V_ec(IJK) = CENTROID_EAST(3)
            ENDIF

            IF(AREA_NORTH > ZERO) THEN
               X_V_nc(IJK) = CENTROID_NORTH(1)
               Y_V_nc(IJK) = CENTROID_NORTH(2)
               Z_V_nc(IJK) = CENTROID_NORTH(3)
            ENDIF

            IF(AREA_TOP > ZERO) THEN
               X_V_tc(IJK) = CENTROID_TOP(1)
               Y_V_tc(IJK) = CENTROID_TOP(2)
               Z_V_tc(IJK) = CENTROID_TOP(3)
            ENDIF

! Compute normal and store cut face info (normal vector and reference point)
            NORMAL_V(IJK,1) = AREA_EAST  - AREA_WEST
            NORMAL_V(IJK,2) = AREA_NORTH - AREA_SOUTH

            IF(DO_K) THEN
               NORMAL_V(IJK,3) = AREA_TOP   - AREA_BOTTOM
            ELSE
               NORMAL_V(IJK,3) = ZERO
            ENDIF


            NORM = sqrt(NORMAL_V(IJK,1)**2 + NORMAL_V(IJK,2)**2 + NORMAL_V(IJK,3)**2)

            IF(NORM==ZERO) NORM = UNDEFINED

            NORMAL_V(IJK,:) = NORMAL_V(IJK,:) / NORM

            REFP_V(IJK,:)   = CENTROID_CUT(:)

            CALL GET_DEL_H(IJK,TYPE_OF_CELL,X_MEAN,Y_MEAN,Z_MEAN,Del_H,Nx,Ny,Nz)
            IF(DEL_H == UNDEFINED) DEL_H = ZERO

            IF(NO_K) THEN
               VOL_V(IJK) = AXY_V(IJK) * ZLENGTH
            ELSE

!             The cell is divided into pyramids having the same apex (Mean value of nodes)
!             Cell Volume = sum of pyramid volums
               VOLUME =   AREA_EAST   * (X_NODE(8) - X_MEAN)                   &
                        + AREA_NORTH  * (Y_NODE(8) - Y_MEAN)                   &
                        + AREA_TOP    * (Z_NODE(8) - Z_MEAN)                   &
                        + AREA_WEST   * (X_MEAN - X_NODE(1))                   &
                        + AREA_SOUTH  * (Y_MEAN - Y_NODE(1))                   &
                        + AREA_BOTTOM * (Z_MEAN - Z_NODE(1))                   &
                        + CUT_AREA    * DEL_H

               VOL_V(IJK) = VOLUME / 3.0D0
            ENDIF

            IF(IERROR==1) print*,'V_VEL IERROR',IJK
            IF(IERROR==1) AXZ(IJK) = ZERO

            IF(AXZ(IJK) <= TOL_SMALL_AREA*DX(I)*DZ(K)) THEN
               WALL_V_AT(IJK) = .TRUE.
               FLAG_N(IJK) = 0
               NUMBER_OF_V_WALL_CELLS = NUMBER_OF_V_WALL_CELLS + 1
               X_V(IJK) = CENTROID_CUT(1)
               Y_V(IJK) = CENTROID_CUT(2)
               Z_V(IJK) = CENTROID_CUT(3)
            ENDIF

            IF(N_QUADRIC>0) THEN
               F_CUT_FACE = UNDEFINED
               N_BC = 0
               BC_V_ID(IJK) = 0
               DO Q_ID = 1, N_QUADRIC
                  DO NODE = 1,N_CUT_FACE_NODES

                     CALL GET_F_QUADRIC(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                     Q_ID,F_CUT_FACE(Q_ID),CLIP_FLAG)

                     IF(DABS(F_CUT_FACE(Q_ID)) < TOL_F) THEN
                        IF(BC_V_ID(IJK)/=BC_ID_Q(Q_ID)) N_BC = N_BC + 1
                        BC_V_ID(IJK) = BC_ID_Q(Q_ID)
                     ENDIF
                  ENDDO
                  IF(BC_V_ID(IJK)>0) THEN
                      IF(BC_TYPE_ENUM(BC_V_ID(IJK))  == CG_MI) EXIT
!                     IF(BC_TYPE_ENUM(BC_V_ID(IJK))  == CG_NSW) EXIT
                  ENDIF

               ENDDO

            ENDIF

            IF(N_POLYGON>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_POLY_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  Q_ID,F_CUT_FACE(1),CLIP_FLAG,BCID)
                   IF(F_CUT_FACE(1)==ZERO) THEN
                     BC_V_ID(IJK) = BCID
                     IF(BC_V_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_V_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF


            IF(N_USR_DEF>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_USR_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  BCID,F_CUT_FACE(1),CLIP_FLAG)
                  IF(DABS(F_CUT_FACE(1)) < TOL_F) THEN
                     BC_V_ID(IJK) = BCID
                     IF(BC_V_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_V_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF


            IF(USE_STL.OR.USE_MSH) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  DO NN=1,N_FACET_AT(IJK)
                     NF=LIST_FACET_AT(IJK,NN)
                     CALL IS_POINT_INSIDE_FACET(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),&
                                                COORD_CUT_FACE_NODES(3,NODE),NF,INSIDE_FACET)
                     IF(INSIDE_FACET) THEN
                        BC_V_ID(IJK) = BC_ID_STL_FACE(NF)
                        IF(BC_V_ID(IJK)>0) THEN
                           IF(BC_TYPE_ENUM(BC_V_ID(IJK))  == CG_MI) EXIT
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

         CASE('W_MOMENTUM')
            IF(I>ISTART1) AREA_WEST   = AYZ_W(IMJK)
            IF(J>JSTART1) AREA_SOUTH  = AXZ_W(IJMK)
            IF(K>KSTART1) AREA_BOTTOM = AXY_W(IJKM)

            AYZ_W(IJK) = AREA_EAST
            AXZ_W(IJK) = AREA_NORTH
            AXY_W(IJK) = AREA_TOP
            Area_W_CUT(IJK) = CUT_AREA

            IF(AREA_EAST > ZERO) THEN
               X_W_ec(IJK) = CENTROID_EAST(1)
               Y_W_ec(IJK) = CENTROID_EAST(2)
               Z_W_ec(IJK) = CENTROID_EAST(3)
            ENDIF

            IF(AREA_NORTH > ZERO) THEN
               X_W_nc(IJK) = CENTROID_NORTH(1)
               Y_W_nc(IJK) = CENTROID_NORTH(2)
               Z_W_nc(IJK) = CENTROID_NORTH(3)
            ENDIF

            IF(AREA_TOP > ZERO) THEN
               X_W_tc(IJK) = CENTROID_TOP(1)
               Y_W_tc(IJK) = CENTROID_TOP(2)
               Z_W_tc(IJK) = CENTROID_TOP(3)
            ENDIF

! Compute normal and store cut face info (normal vector and reference point)
            NORMAL_W(IJK,1) = AREA_EAST  - AREA_WEST
            NORMAL_W(IJK,2) = AREA_NORTH - AREA_SOUTH
            NORMAL_W(IJK,3) = AREA_TOP   - AREA_BOTTOM


            NORM = sqrt(NORMAL_W(IJK,1)**2 + NORMAL_W(IJK,2)**2 + NORMAL_W(IJK,3)**2)

            IF(NORM==ZERO) NORM = UNDEFINED

            NORMAL_W(IJK,:) = NORMAL_W(IJK,:) / NORM

            REFP_W(IJK,:)   = CENTROID_CUT(:)


            CALL GET_DEL_H(IJK,TYPE_OF_CELL,X_MEAN,Y_MEAN,Z_MEAN,Del_H,Nx,Ny,Nz)
            IF(DEL_H == UNDEFINED) DEL_H = ZERO


!          The cell is divided into pyramids having the same apex (Mean value of nodes)
!          Cell Volume = sum of pyramid volums
            VOLUME =   AREA_EAST   * (X_NODE(8) - X_MEAN)                   &
                     + AREA_NORTH  * (Y_NODE(8) - Y_MEAN)                   &
                     + AREA_TOP    * (Z_NODE(8) - Z_MEAN)                   &
                     + AREA_WEST   * (X_MEAN - X_NODE(1))                   &
                     + AREA_SOUTH  * (Y_MEAN - Y_NODE(1))                   &
                     + AREA_BOTTOM * (Z_MEAN - Z_NODE(1))                   &
                     + CUT_AREA    * DEL_H

            VOL_W(IJK) = VOLUME / 3.0D0

            IF(IERROR==1) print*,'W_VEL IERROR',IJK
            IF(IERROR==1) AXY(IJK) = ZERO


            IF(AXY(IJK) <= TOL_SMALL_AREA*DX(I)*DY(J)) THEN
               WALL_W_AT(IJK) = .TRUE.
               FLAG_T(IJK) = 0
               NUMBER_OF_W_WALL_CELLS = NUMBER_OF_W_WALL_CELLS + 1
               X_W(IJK) = CENTROID_CUT(1)
               Y_W(IJK) = CENTROID_CUT(2)
               Z_W(IJK) = CENTROID_CUT(3)
            ENDIF

            IF(N_QUADRIC>0) THEN
               F_CUT_FACE = UNDEFINED
               N_BC = 0
               BC_W_ID(IJK) = 0
               DO Q_ID = 1, N_QUADRIC
                  DO NODE = 1,N_CUT_FACE_NODES
                     CALL GET_F_QUADRIC(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                     Q_ID,F_CUT_FACE(Q_ID),CLIP_FLAG)
                     IF(DABS(F_CUT_FACE(Q_ID)) < TOL_F) THEN
                        IF(BC_W_ID(IJK)/=BC_ID_Q(Q_ID)) N_BC = N_BC + 1
                        BC_W_ID(IJK) = BC_ID_Q(Q_ID)
                     ENDIF
                  ENDDO
                  IF(BC_W_ID(IJK)>0) THEN
                     IF(BC_TYPE_ENUM(BC_W_ID(IJK))  == CG_MI) EXIT
!                     IF(BC_TYPE_ENUM(BC_W_ID(IJK))  == CG_NSW) EXIT
                  ENDIF
               ENDDO
            ENDIF


            IF(N_POLYGON>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_POLY_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  Q_ID,F_CUT_FACE(1),CLIP_FLAG,BCID)
                   IF(F_CUT_FACE(1)==ZERO) THEN
                     BC_W_ID(IJK) = BCID
                     IF(BC_W_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_W_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF



            IF(N_USR_DEF>0) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  CALL EVAL_USR_FCT(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),COORD_CUT_FACE_NODES(3,NODE),&
                  BCID,F_CUT_FACE(1),CLIP_FLAG)
                  IF(DABS(F_CUT_FACE(1)) < TOL_F) THEN
                     BC_W_ID(IJK) = BCID
                     IF(BC_W_ID(IJK)>0) THEN
                        IF(BC_TYPE_ENUM(BC_W_ID(IJK))  == CG_MI) EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF


            IF(USE_STL.OR.USE_MSH) THEN
               DO NODE = 1,N_CUT_FACE_NODES
                  DO NN=1,N_FACET_AT(IJK)
                     NF=LIST_FACET_AT(IJK,NN)
                     CALL IS_POINT_INSIDE_FACET(COORD_CUT_FACE_NODES(1,NODE),COORD_CUT_FACE_NODES(2,NODE),&
                                                COORD_CUT_FACE_NODES(3,NODE),NF,INSIDE_FACET)
                     IF(INSIDE_FACET) THEN
                        BC_W_ID(IJK) = BC_ID_STL_FACE(NF)
                        IF(BC_W_ID(IJK)>0) THEN
                           IF(BC_TYPE_ENUM(BC_W_ID(IJK))  == CG_MI) EXIT
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF


         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_CUT_CELL_VOLUME_AND_AREAS'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN

      END SUBROUTINE GET_CUT_CELL_VOLUME_AND_AREAS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_POLYGON_AREA_AND_CENTROID                          C
!  Purpose: Compute location of centroid of a 3D polygon               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_POLYGON_AREA_AND_CENTROID(NP,COORD,AREA,CENTROID,IERROR)

      USE compar, only: mype
      USE cutcell
      USE exit, only: mfix_exit
      USE geometry, ONLY: no_k, zlength
      USE quadric, ONLY: cross_product
      USE param1, only: half, undefined, zero

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: NP,IERROR
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(3) :: CENTROID
      DOUBLE PRECISION, INTENT(INOUT) :: AREA
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(3,15) :: COORD

      INTEGER :: NN,NC,NEW_NP
      INTEGER :: LAST,LASTM
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_BCK
      DOUBLE PRECISION, DIMENSION(3) :: NORMAL,CP,SUMCP
      DOUBLE PRECISION, DIMENSION(3) :: LASTM_VEC,LAST_VEC
      DOUBLE PRECISION, DIMENSION(3) :: R,SUM_AR,AV
      DOUBLE PRECISION :: A,SUM_A
      LOGICAL, DIMENSION(15) :: KEEP
      DOUBLE PRECISION, DIMENSION(15) :: D

!======================================================================
!   Remove duplicate points in the list
!======================================================================
      IF(.FALSE.) THEN
      DO NN=1,NP
         D(NN) = sqrt(dot_product(COORD(:,NN),COORD(:,NN)))
         KEEP(NN) = .TRUE.
         COORD_BCK(:,NN) = COORD(:,NN)
      ENDDO

      DO NN=1,NP-1
         DO NC=NN+1,NP
            IF(DABS(D(NN)-D(NC))<TOL_MERGE) KEEP(NC)=.FALSE.
         ENDDO
      ENDDO

      NEW_NP = 0

      DO NN=1,NP
         IF(KEEP(NN)) THEN
            NEW_NP = NEW_NP + 1
            COORD(:,NEW_NP) = COORD_BCK(:,NN)
         ENDIF
      ENDDO

      NP = NEW_NP
      ENDIF

      IF( NP < 2 ) THEN
         AREA = ZERO
         CENTROID = UNDEFINED
         RETURN
      ELSEIF( NP == 2 ) THEN
         IF(NO_K) THEN  ! 2D case
            AREA = sqrt((COORD(1,2)-COORD(1,1))**2 + (COORD(2,2)-COORD(2,1))**2) * ZLENGTH
            CENTROID(1) = HALF * (COORD(1,1)+COORD(1,2))
            CENTROID(2) = HALF * (COORD(2,1)+COORD(2,2))
            CENTROID(3) = ZERO
            RETURN
         ELSE
           AREA = ZERO
           CENTROID = UNDEFINED
           RETURN
!           WRITE(*,*)'CRITICAL ERROR: POLYGON WITH ONLY 2 POINTS in 3D.'
!           WRITE(*,*)'LIST OF COORDINATES:'
!           WRITE(*,*)'NODE 1: (X,Y,Z)=', COORD(1,1),COORD(2,1),COORD(3,1)
!           WRITE(*,*)'NODE 2: (X,Y,Z)=', COORD(1,2),COORD(2,2),COORD(3,2)
!           WRITE(*,*)'MFiX will exit now.'
!           CALL MFIX_EXIT(myPE)
         ENDIF
      ELSEIF( NP > 6 ) THEN
         WRITE(*,*)'CRITICAL ERROR: POLYGON WITH MORE THAN 6 POINTS.',NP
         DO NN=1,NP
            print*,COORD(:,NN)
         ENDDO
         WRITE(*,*)'MFiX will exit now.'
         CALL MFIX_EXIT(myPE)
      ENDIF

! The following instructions are only executed if 3 <= NP < 6

!======================================================================
!   Reorder nodes to create a convex polygon
!======================================================================

      CALL REORDER_POLYGON(NP,COORD,NORMAL,IERROR)

!======================================================================
!   Duplicate last node to close the polygon
!======================================================================

      COORD(:,NP+1) = COORD(:,1)

!======================================================================
!   Accumulate sum of cross products
!======================================================================

      SUMCP = ZERO

      DO NN = 1,NP
         CP = CROSS_PRODUCT(COORD(:,NN),COORD(:,NN+1))
         SUMCP = SUMCP + CP
      ENDDO

!======================================================================
!   Take dot product by Normal unit vector and divide by two
!======================================================================

      AREA = DABS(HALF * DOT_PRODUCT(NORMAL, SUMCP))

!======================================================================
!   Alternate computation: Split the polygon into triangles
!======================================================================

      SUM_AR = ZERO
      SUM_A  = ZERO
      DO NN = 1,NP-2
         LAST  = NN + 2
         LASTM = LAST - 1
         R = (COORD(:,1)+COORD(:,LASTM)+COORD(:,LAST))/3.0D0
         LASTM_VEC = COORD(:,LASTM)-COORD(:,1)
         LAST_VEC  = COORD(:,LAST)-COORD(:,1)
         AV = CROSS_PRODUCT(LASTM_VEC,LAST_VEC)
         A = sqrt(dot_product(av(:),av(:)))
         SUM_AR = SUM_AR + A * R
         SUM_A  = SUM_A  + A
      ENDDO
      CENTROID = SUM_AR / SUM_A

      RETURN

      END SUBROUTINE GET_POLYGON_AREA_AND_CENTROID

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REORDER_POLYGON                                        C
!  Purpose: Re-order polygon nodes to form a convex polygon            C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE REORDER_POLYGON(NP,COORD,NORMAL,IERROR)

      USE compar, ONLY: mype
      USE cutcell
      USE quadric, only: cross_product
      USE param1, only: zero

      IMPLICIT NONE
      INTEGER :: I,K,NN,NP,NC,NEW_NP
      INTEGER :: LONGEST
      INTEGER :: I_SWAP,IERROR

      DOUBLE PRECISION :: Xc,Yc,Zc,NORM,ANGLE,ANGLE_REF
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD,COORD_BCK
      DOUBLE PRECISION, DIMENSION(3) :: NORMAL,SWAP
      DOUBLE PRECISION, DIMENSION(3) :: VECTMP,VECTMP2

      LOGICAL, DIMENSION(15) :: KEEP
      DOUBLE PRECISION, DIMENSION(15) :: D

!======================================================================
!   Remove duplicate points in the list
!======================================================================
      IF(.FALSE.) THEN
      DO NN=1,NP
         D(NN) = sqrt(dot_product(COORD(:,NN),COORD(:,NN)))
         KEEP(NN) = .TRUE.
         COORD_BCK(:,NN) = COORD(:,NN)
      ENDDO

      DO NN=1,NP-1
         DO NC=NN+1,NP
            IF(DABS(D(NN)-D(NC))<TOL_MERGE) KEEP(NC)=.FALSE.
         ENDDO
      ENDDO

      NEW_NP = 0

      DO NN=1,NP
         IF(KEEP(NN)) THEN
            NEW_NP = NEW_NP + 1
            COORD(:,NEW_NP) = COORD_BCK(:,NN)
         ENDIF
      ENDDO

      NP = NEW_NP
      ENDIF

!======================================================================
!   Exit if polygon has less than 3 vertices (no need to reorder)
!======================================================================

      IF(NP<=2) RETURN

!======================================================================
!   Initialize order list
!======================================================================

      DO NN=1,NP
         ORDER(NN) = NN
      ENDDO

!======================================================================
!   Find Mean of polygon points
!======================================================================

      Xc = ZERO
      Yc = ZERO
      Zc = ZERO

      DO NN=1,NP
         Xc = Xc + COORD(1,NN)
         Yc = Yc + COORD(2,NN)
         Zc = Zc + COORD(3,NN)
      ENDDO

      Xc = Xc / DBLE(NP)
      Yc = Yc / DBLE(NP)
      Zc = Zc / DBLE(NP)

!======================================================================
!   Find unit normal vector
!======================================================================

      VECTMP  = COORD(:,2)-COORD(:,1)
      VECTMP2 = COORD(:,3)-COORD(:,1)
      NORMAL = CROSS_PRODUCT(VECTMP,VECTMP2)

      NORM = sqrt(dot_product(NORMAL(:),NORMAL(:)))
      IF(NORM==ZERO) THEN
!         DO NN=1,NP
!            print*,N,COORD(:,N)
!         ENDDO

         IERROR = 1

         RETURN

         WRITE(*,*)'ERROR IN REORDER_POLYGON: NORMAL UNIT VECTOR HAS ZERO NORM'
         WRITE(*,*)'THIS USUALLY HAPPENS WHEN NODES WERE IMPROPERLY MERGED.'
         WRITE(*,*)'TO PREVENT THIS, INCREASE TOL_SNAP, OR DECREASE TOL_MERGE.'
         WRITE(*,*)'CURRENT VALUE OF TOL_SNAP  = ', TOL_SNAP
         WRITE(*,*)'CURRENT VALUE OF TOL_MERGE = ', TOL_MERGE
         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         CALL MFIX_EXIT(MYPE)
      ELSE
         NORMAL = NORMAL / NORM
      ENDIF
!======================================================================
!   Find longest component of normal vector
!======================================================================
      LONGEST = 3
      IF(ABS(NORMAL(1)) > ABS(NORMAL(2))) THEN
         IF(ABS(NORMAL(1)) > ABS(NORMAL(3))) LONGEST = 1
      ELSE
         IF(ABS(NORMAL(2)) > ABS(NORMAL(3))) LONGEST = 2
      ENDIF

!======================================================================
!   Reorder nodes to create a convex polygon
!======================================================================
      SELECT CASE(LONGEST)

         CASE(1)
            DO I = 1, NP-1
               ANGLE_REF = ATAN2 (COORD(3,I) - Zc,COORD(2,I) - Yc)
               DO K = I+1,NP
                  ANGLE = ATAN2 (COORD(3,K) - Zc,COORD(2,K) - Yc)
                  IF(ANGLE < ANGLE_REF) THEN
                     ANGLE_REF = ANGLE
                     SWAP  = COORD(:,K)
                     COORD(:,K) = COORD(:,I)
                     COORD(:,I) = SWAP
                     I_SWAP = ORDER(K)
                     ORDER(K) = ORDER(I)
                     ORDER(I) = I_SWAP
                  ENDIF
               END DO
            END DO

         CASE(2)
            DO I = 1, NP-1
               ANGLE_REF = ATAN2 (COORD(1,I) - Xc,COORD(3,I) - Zc)
               DO K = I+1,NP
                  ANGLE = ATAN2 (COORD(1,K) - Xc,COORD(3,K) - Zc)
                  IF(ANGLE < ANGLE_REF) THEN
                     ANGLE_REF = ANGLE
                     SWAP  = COORD(:,K)
                     COORD(:,K) = COORD(:,I)
                     COORD(:,I) = SWAP
                     I_SWAP = ORDER(K)
                     ORDER(K) = ORDER(I)
                     ORDER(I) = I_SWAP
                  ENDIF
               END DO
            END DO

         CASE(3)
            DO I = 1, NP-1
               ANGLE_REF = ATAN2 (COORD(2,I) - Yc,COORD(1,I) - Xc)
               DO K = I+1,NP
                  ANGLE = ATAN2 (COORD(2,K) - Yc,COORD(1,K) - Xc)
                  IF(ANGLE < ANGLE_REF) THEN
                     ANGLE_REF = ANGLE
                     SWAP  = COORD(:,K)
                     COORD(:,K) = COORD(:,I)
                     COORD(:,I) = SWAP
                     I_SWAP = ORDER(K)
                     ORDER(K) = ORDER(I)
                     ORDER(I) = I_SWAP
                  ENDIF
               END DO
            END DO

      END SELECT

      RETURN

      END SUBROUTINE REORDER_POLYGON

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SORT                                                   C
!  Purpose: Rearrange column COL2 of an array                          C
!           by sorting column COL1 in increasing order                 C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SORT(Array, NUM_ROWS,COL1,COL2)


      IMPLICIT NONE

      INTEGER :: M,N,NUM_ROWS,COL,COL1,COL2,Saved_index
      DOUBLE PRECISION :: Swap,Local_min
      INTEGER, PARAMETER :: ARRAY_SIZE = 10000
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE,10) :: Array


      DO N = 1, NUM_ROWS-1
         Local_min = Array (N,COL1)
         Saved_index = N
         DO M = N,NUM_ROWS
            if (Array(M,COL1) < Local_min ) Then
               Local_min = Array (M,COL1)
               Saved_index = M
            endif
         END DO

         DO COL = COL1,COL2
            Swap = Array(N,COL)
            Array(N,COL) = Array(Saved_index,COL)
            Array(Saved_index,COL) = Swap
         END DO

      END DO

      RETURN

      END SUBROUTINE SORT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_INTERPOLATION_TERMS_G                              C
!  Purpose: gets terms involved in velocity interpolation  (Gas phase) C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_INTERPOLATION_TERMS_G(IJK,TYPE_OF_CELL,ALPHA_CUT,AW,HW,VELW)

      USE bc
      USE compar, ONLY: mype
      USE cutcell
      USE exit, only: mfix_exit
      USE param1, only: one, zero, undefined

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK
      INTEGER :: BCV
      INTEGER :: BCT
      DOUBLE PRECISION :: ALPHA_CUT,AW,HW,VELW

!======================================================================
! The alpha correction term is only used for No-slip walls
! In case of PSW, only the two extreme conditions ( BC_HW_G(BCV)==UNDEFINED, corresponding to NSW
! and BC_HW_G(BCV)==ZERO, corresponding to FSW) are active.
! True partial slip will be implemented in the future.
!======================================================================

     AW   = ALPHA_CUT
     HW   = ZERO
     VELW = ZERO

      SELECT CASE (TYPE_OF_CELL)

         CASE('U_MOMENTUM')

            BCV = BC_U_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
               SELECT CASE(BCT)

                  CASE(CG_NSW)
                     AW   = ALPHA_CUT
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_FSW)
                     AW   = ONE
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_PSW)

                     IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                        AW   = ALPHA_CUT
                        HW   = ZERO
                        VELW = ZERO
                     ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                        AW   = ONE
                        HW   = ZERO
                        VELW = ZERO
                     ELSE                              ! partial slip
                        AW   = ONE !ALPHA_CUT
                        HW   = BC_HW_G(BCV)
                        VELW = BC_UW_G(BCV)
                     ENDIF

               END SELECT

            ENDIF


         CASE('V_MOMENTUM')

            BCV = BC_V_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
               SELECT CASE(BCT)

                  CASE(CG_NSW)
                     AW   = ALPHA_CUT
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_FSW)
                     AW   = ONE
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_PSW)

                     IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                        AW   = ALPHA_CUT
                        HW   = ZERO
                        VELW = ZERO
                     ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                        AW   = ONE
                        HW   = ZERO
                        VELW = ZERO
                     ELSE                              ! partial slip
                        AW   = ONE !ALPHA_CUT
                        HW   = BC_HW_G(BCV)
                        VELW = BC_VW_G(BCV)
                     ENDIF

               END SELECT

            ENDIF


         CASE('W_MOMENTUM')

            BCV = BC_W_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
               SELECT CASE(BCT)

                  CASE(CG_NSW)
                     AW   = ALPHA_CUT
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_FSW)
                     AW   = ONE
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_PSW)
                     IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                        AW   = ALPHA_CUT
                        HW   = ZERO
                        VELW = ZERO
                     ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                        AW   = ONE
                        HW   = ZERO
                        VELW = ZERO
                     ELSE                              ! partial slip
                        AW   = ONE !ALPHA_CUT
                        HW   = BC_HW_G(BCV)
                        VELW = BC_WW_G(BCV)
                     ENDIF

               END SELECT

            ENDIF

         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_INTERPOLATION_TERMS_G'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN

      END SUBROUTINE GET_INTERPOLATION_TERMS_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_INTERPOLATION_TERMS_S                              C
!  Purpose: gets terms involved in velocity interpolation              C
!  (Solids phase)                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_INTERPOLATION_TERMS_S(IJK,M,TYPE_OF_CELL,ALPHA_CUT,AW,HW,VELW)

      USE bc
      USE compar, ONLY: mype
      USE cutcell
      USE exit, only: mfix_exit
      USE param1, ONLY: one, zero, undefined

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,M
      INTEGER :: BCV
      INTEGER :: BCT
      DOUBLE PRECISION :: ALPHA_CUT,AW,HW,VELW

!======================================================================
! The alpha correction term is only used for No-slip walls
! In case of PSW, only the two extreme conditions ( BC_HW_G(BCV)==UNDEFINED, corresponding to NSW
! and BC_HW_G(BCV)==ZERO, corresponding to FSW) are active.
! True partial slip will be implemented in the future.
!======================================================================

     AW   = ALPHA_CUT
     HW   = ZERO
     VELW = ZERO

      SELECT CASE (TYPE_OF_CELL)

         CASE('U_MOMENTUM')

            BCV = BC_U_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
               SELECT CASE(BCT)

                  CASE(CG_NSW)
                     AW   = ALPHA_CUT
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_FSW)
                     AW   = ONE
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_PSW)

                     IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                        AW   = ALPHA_CUT
                        HW   = ZERO
                        VELW = ZERO
                     ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                        AW   = ONE
                        HW   = ZERO
                        VELW = ZERO
                     ELSE                              ! partial slip
                        AW   = ALPHA_CUT
                        HW   = BC_HW_S(BCV,M)
                        VELW = BC_UW_S(BCV,M)
                     ENDIF

               END SELECT

            ENDIF


         CASE('V_MOMENTUM')

            BCV = BC_V_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
               SELECT CASE(BCT)

                  CASE(CG_NSW)
                     AW   = ALPHA_CUT
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_FSW)
                     AW   = ONE
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_PSW)

                     IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                        AW   = ALPHA_CUT
                        HW   = ZERO
                        VELW = ZERO
                     ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                        AW   = ONE
                        HW   = ZERO
                        VELW = ZERO
                     ELSE                              ! partial slip
                        AW   = ALPHA_CUT
                        HW   = BC_HW_S(BCV,M)
                        VELW = BC_VW_S(BCV,M)
                     ENDIF

               END SELECT

            ENDIF


         CASE('W_MOMENTUM')

            BCV = BC_W_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
               SELECT CASE(BCT)

                  CASE(CG_NSW)
                     AW   = ALPHA_CUT
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_FSW)
                     AW   = ONE
                     HW   = ZERO
                     VELW = ZERO

                  CASE(CG_PSW)
                     IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                        AW   = ALPHA_CUT
                        HW   = ZERO
                        VELW = ZERO
                     ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                        AW   = ONE
                        HW   = ZERO
                        VELW = ZERO
                     ELSE                              ! partial slip
                        AW   = ALPHA_CUT
                        HW   = BC_HW_S(BCV,M)
                        VELW = BC_WW_S(BCV,M)
                     ENDIF

               END SELECT

            ENDIF

         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_INTERPOLATION_TERMS_S'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN

      END SUBROUTINE GET_INTERPOLATION_TERMS_S
