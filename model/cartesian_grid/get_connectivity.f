!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_CONNECTIVITY                                       C
!  Purpose: Set flags for saclar cut cells, based on intersection      C
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
  SUBROUTINE GET_CONNECTIVITY(IJK,TYPE_OF_CELL,N_NEW_POINTS,N_NODES,CONNECT,X_NP,Y_NP,Z_NP,TOTAL_NUMBER_OF_INTERSECTIONS,&
             X_intersect,Y_intersect,Z_intersect)

      USE compar, ONLY: ijkend3
      USE cutcell
      USE cut_cell_preproc, ONLY: eval_f
      USE functions, ONLY: FUNIJK
      USE geometry, ONLY: DO_K, NO_K
      USE indices, ONLY: I_OF, J_OF, K_OF
      USE polygon, ONLY: n_polygon
      USE quadric, ONLY: tol_f
      USE param, ONLY: DIMENSION_3

      IMPLICIT NONE
      INTEGER :: I,J,K,IM,JM,KM
      INTEGER :: IJK,IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM
      LOGICAL :: CLIP_FLAG
      LOGICAL, DIMENSION(8) :: CORNER_INTERSECTION
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS,NUMBER_OF_EDGE_INTERSECTIONS
      INTEGER :: NUMBER_OF_CORNER_INTERSECTIONS,MAX_CORNER_INTERSECTIONS
      INTEGER :: N_NODES,N_NEW_POINTS
      INTEGER :: NODE,N_N1,N_N2,Q_ID,BCID
      INTEGER, DIMENSION(DIMENSION_3,15) :: CONNECT
      DOUBLE PRECISION, DIMENSION(DIMENSION_MAX_CUT_CELL) :: X_NP,Y_NP,Z_NP
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: X_intersect,Y_intersect,Z_intersect
      CHARACTER (LEN=*) :: TYPE_OF_CELL

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      IM = I - 1
      JM = J - 1
      KM = K - 1

      IMJK   = FUNIJK(IM,J,K)
      IJMK   = FUNIJK(I,JM,K)
      IJKM   = FUNIJK(I,J,KM)

      IMJMK  = FUNIJK(IM,JM,K)
      IMJKM  = FUNIJK(IM,J,KM)
      IJMKM  = FUNIJK(I,JM,KM)

      IMJMKM = FUNIJK(IM,JM,KM)


!======================================================================
!  Evaluate Quadric at all corners
!======================================================================

         N_NODES = 0

         CORNER_INTERSECTION = .FALSE.
         NUMBER_OF_CORNER_INTERSECTIONS = 0
         NUMBER_OF_EDGE_INTERSECTIONS = 0

         IF(NO_K) THEN
            N_N1 = 5
            N_N2 = 8
         ELSE
            N_N1 = 1
            N_N2 = 8
         ENDIF

         DO NODE = N_N1,N_N2

            Q_ID = 1
            CALL EVAL_F('QUADRIC',X_NODE(NODE),Y_NODE(NODE),Z_NODE(NODE),Q_ID,F_NODE(NODE),CLIP_FLAG)

            CALL EVAL_F('POLYGON',X_NODE(NODE),Y_NODE(NODE),Z_NODE(NODE),N_POLYGON,F_NODE(NODE),CLIP_FLAG)

            CALL EVAL_F('USR_DEF',X_NODE(NODE),Y_NODE(NODE),Z_NODE(NODE),N_USR_DEF,F_NODE(NODE),CLIP_FLAG)

            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,NODE,F_NODE(NODE),CLIP_FLAG,BCID)

            IF (ABS(F_NODE(NODE)) < TOL_F ) THEN
               CORNER_INTERSECTION(NODE) = .TRUE.
               NUMBER_OF_CORNER_INTERSECTIONS = NUMBER_OF_CORNER_INTERSECTIONS + 1
               N_NODES = N_NODES + 1
               CONNECT(IJK,N_NODES) = IJK_OF_NODE(NODE)
            ENDIF


            IF(SNAP(IJK_OF_NODE(NODE))) THEN
               CORNER_INTERSECTION(NODE) = .TRUE.
               NUMBER_OF_CORNER_INTERSECTIONS = NUMBER_OF_CORNER_INTERSECTIONS + 1
               N_NODES = N_NODES + 1
               CONNECT(IJK,N_NODES) = IJK_OF_NODE(NODE)
            ENDIF

         END DO

!======================================================================
!  Count the number of edge intersections (excluding corner intersections)
!  For each new edge intersection found:
!     - Increment the total number of points
!     - Store the location of the additional point
!======================================================================

            IF(DO_K) THEN

               IF(INTERSECT_X(IJMKM)) THEN  ! Edge 1  = Nodes 1-2
                  IF((.NOT.CORNER_INTERSECTION(1)).AND.(.NOT.CORNER_INTERSECTION(2))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_intersect(IJMKM)
                     Y_NP(N_NEW_POINTS) = Y_NODE(1)
                     Z_NP(N_NEW_POINTS) = Z_NODE(1)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

               IF(INTERSECT_Y(IJKM)) THEN  ! Edge 2  = Nodes 2-4
                  IF((.NOT.CORNER_INTERSECTION(2)).AND.(.NOT.CORNER_INTERSECTION(4))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_NODE(2)
                     Y_NP(N_NEW_POINTS) = Y_intersect(IJKM)
                     Z_NP(N_NEW_POINTS) = Z_NODE(2)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

               IF(INTERSECT_X(IJKM)) THEN  ! Edge 3  = Nodes 3-4
                  IF((.NOT.CORNER_INTERSECTION(3)).AND.(.NOT.CORNER_INTERSECTION(4))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_intersect(IJKM)
                     Y_NP(N_NEW_POINTS) = Y_NODE(3)
                     Z_NP(N_NEW_POINTS) = Z_NODE(3)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

               IF(INTERSECT_Y(IMJKM)) THEN  ! Edge 4  = Nodes 1-3
                  IF((.NOT.CORNER_INTERSECTION(1)).AND.(.NOT.CORNER_INTERSECTION(3))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_NODE(1)
                     Y_NP(N_NEW_POINTS) = Y_intersect(IMJKM)
                     Z_NP(N_NEW_POINTS) = Z_NODE(1)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

            ENDIF


            IF(INTERSECT_X(IJMK)) THEN  ! Edge 5  = Nodes 5-6
               IF((.NOT.CORNER_INTERSECTION(5)).AND.(.NOT.CORNER_INTERSECTION(6))) THEN
                  NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                  N_NEW_POINTS = N_NEW_POINTS + 1
                  X_NP(N_NEW_POINTS) = X_intersect(IJMK)
                  Y_NP(N_NEW_POINTS) = Y_NODE(5)
                  Z_NP(N_NEW_POINTS) = Z_NODE(5)
                  N_NODES = N_NODES + 1
                  CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
               ENDIF
            ENDIF

            IF(INTERSECT_Y(IJK)) THEN  ! Edge 6  = Nodes 6-8
               IF((.NOT.CORNER_INTERSECTION(6)).AND.(.NOT.CORNER_INTERSECTION(8))) THEN
                  NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                  N_NEW_POINTS = N_NEW_POINTS + 1
                  X_NP(N_NEW_POINTS) = X_NODE(6)
                  Y_NP(N_NEW_POINTS) = Y_intersect(IJK)
                  Z_NP(N_NEW_POINTS) = Z_NODE(6)
                  N_NODES = N_NODES + 1
                  CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
               ENDIF
            ENDIF

            IF(INTERSECT_X(IJK)) THEN  ! Edge 7  = Nodes 7-8
               IF((.NOT.CORNER_INTERSECTION(7)).AND.(.NOT.CORNER_INTERSECTION(8))) THEN
                  NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                  N_NEW_POINTS = N_NEW_POINTS + 1
                  X_NP(N_NEW_POINTS) = X_intersect(IJK)
                  Y_NP(N_NEW_POINTS) = Y_NODE(7)
                  Z_NP(N_NEW_POINTS) = Z_NODE(7)
                  N_NODES = N_NODES + 1
                  CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
               ENDIF
            ENDIF

            IF(INTERSECT_Y(IMJK)) THEN  ! Edge 8  = Nodes 5-7
               IF((.NOT.CORNER_INTERSECTION(5)).AND.(.NOT.CORNER_INTERSECTION(7))) THEN
                  NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                  N_NEW_POINTS = N_NEW_POINTS + 1
                  X_NP(N_NEW_POINTS) = X_NODE(5)
                  Y_NP(N_NEW_POINTS) = Y_intersect(IMJK)
                  Z_NP(N_NEW_POINTS) = Z_NODE(5)
                  N_NODES = N_NODES + 1
                  CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
               ENDIF
            ENDIF


            IF(DO_K) THEN

               IF(INTERSECT_Z(IMJMK)) THEN  ! Edge 9  = Nodes 1-5
                  IF((.NOT.CORNER_INTERSECTION(1)).AND.(.NOT.CORNER_INTERSECTION(5))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_NODE(1)
                     Y_NP(N_NEW_POINTS) = Y_NODE(1)
                     Z_NP(N_NEW_POINTS) = Z_intersect(IMJMK)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

               IF(INTERSECT_Z(IJMK))  THEN  ! Edge 10 = Nodes 2-6
                  IF((.NOT.CORNER_INTERSECTION(2)).AND.(.NOT.CORNER_INTERSECTION(6))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_NODE(2)
                     Y_NP(N_NEW_POINTS) = Y_NODE(2)
                     Z_NP(N_NEW_POINTS) = Z_intersect(IJMK)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

               IF(INTERSECT_Z(IJK))  THEN  ! Edge 11 = Nodes 4-8
                  IF((.NOT.CORNER_INTERSECTION(4)).AND.(.NOT.CORNER_INTERSECTION(8))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_NODE(4)
                     Y_NP(N_NEW_POINTS) = Y_NODE(4)
                     Z_NP(N_NEW_POINTS) = Z_intersect(IJK)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

               IF(INTERSECT_Z(IMJK)) THEN  ! Edge 12 = Nodes 3-7
                  IF((.NOT.CORNER_INTERSECTION(3)).AND.(.NOT.CORNER_INTERSECTION(7))) THEN
                     NUMBER_OF_EDGE_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + 1
                     N_NEW_POINTS = N_NEW_POINTS + 1
                     X_NP(N_NEW_POINTS) = X_NODE(3)
                     Y_NP(N_NEW_POINTS) = Y_NODE(3)
                     Z_NP(N_NEW_POINTS) = Z_intersect(IMJK)
                     N_NODES = N_NODES + 1
                     CONNECT(IJK,N_NODES) = N_NEW_POINTS + IJKEND3
                  ENDIF
               ENDIF

            ENDIF

!======================================================================
!  Count the total number of intersections (corner and edge intersections)
!======================================================================

            IF(NO_K) THEN
               MAX_CORNER_INTERSECTIONS = 2
            ELSE
               MAX_CORNER_INTERSECTIONS = 4
            ENDIF


           IF(NUMBER_OF_CORNER_INTERSECTIONS > MAX_CORNER_INTERSECTIONS) THEN
              IF(PRINT_WARNINGS) THEN
                 WRITE(*,*)'WARNING:',NUMBER_OF_CORNER_INTERSECTIONS,&
                        ' CORNER INTERSECTIONS DETECTED IN CELL IJK=',IJK
                 WRITE(*,*)'THIS USUALLY INDICATE A FALSE CUT-CELL (CORNER CELL)'
!                 WRITE(*,*)'RESETTING NUMBER_OF_CORNER_INTERSECTIONS TO', MAX_CORNER_INTERSECTIONS
                 WRITE(*,*)'REMOVING CUT CELL'
              ENDIF
!              NUMBER_OF_CORNER_INTERSECTIONS = MAX_CORNER_INTERSECTIONS

              ! Force the total number of intersections to be zero, and therefore, the cell will be considered as a non-cut cell
!              NUMBER_OF_CORNER_INTERSECTIONS = -NUMBER_OF_EDGE_INTERSECTIONS

              ! Force the total number of intersections to be -one, and therefore, the cell will be considered as a non-cut cell
              NUMBER_OF_CORNER_INTERSECTIONS = -NUMBER_OF_EDGE_INTERSECTIONS -1


           ENDIF

            TOTAL_NUMBER_OF_INTERSECTIONS = NUMBER_OF_EDGE_INTERSECTIONS + NUMBER_OF_CORNER_INTERSECTIONS

      RETURN

      END SUBROUTINE GET_CONNECTIVITY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_CELL_NODE_COORDINATES                              C
!  Purpose: Get the cell corners (x,y,z) coordinates                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      USE compar, ONLY: mype
      USE cutcell
      USE exit, only: mfix_exit
      USE functions, ONLY: FUNIJK
      USE geometry, ONLY: NO_K, dx, dy, dz
      USE indices, ONLY: I_OF, J_OF, K_OF
      USE param1, ONLY: HALF, ZERO

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys,Zb,Zt
      INTEGER :: I,J,K,IP,JP,KP,IM,JM,KM
      INTEGER :: IJK,IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM

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

      IMJMK  = FUNIJK(IM,JM,K)
      IMJKM  = FUNIJK(IM,J,KM)
      IJMKM  = FUNIJK(I,JM,KM)

      IMJMKM = FUNIJK(IM,JM,KM)

      IJK_OF_NODE(0) = IJK
      IJK_OF_NODE(1) = IMJMKM
      IJK_OF_NODE(2) = IJMKM
      IJK_OF_NODE(3) = IMJKM
      IJK_OF_NODE(4) = IJKM
      IJK_OF_NODE(5) = IMJMK
      IJK_OF_NODE(6) = IJMK
      IJK_OF_NODE(7) = IMJK
      IJK_OF_NODE(8) = IJK
      IJK_OF_NODE(15) = IJK


      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')
            Xw = XG_E(I) - DX(I)          ! west face location
            Xe = XG_E(I)                  ! east face location

            Ys = YG_N(J) - DY(J)          ! south face location
            Yn = YG_N(J)                  ! north face location

            IF(NO_K) THEN
               Zb = ZERO                  ! bottom face location
               Zt = ZERO                  ! top face location
            ELSE
               Zb = ZG_T(K) - DZ(K)       ! bottom face location
               Zt = ZG_T(K)               ! top face location
            ENDIF

         CASE('U_MOMENTUM')
            Xw = XG_E(I) - HALF * DX(I)   ! west face location
            Xe = XG_E(I) + HALF * DX(IP)  ! east face location

            Ys = YG_N(J) - DY(J)          ! south face location
            Yn = YG_N(J)                  ! north face location

            IF(NO_K) THEN
               Zb = ZERO                  ! bottom face location
               Zt = ZERO                  ! top face location
            ELSE
               Zb = ZG_T(K) - DZ(K)       ! bottom face location
               Zt = ZG_T(K)               ! top face location
            ENDIF

         CASE('V_MOMENTUM')
            Xw = XG_E(I) - DX(I)          ! west face location
            Xe = XG_E(I)                  ! east face location

            Ys = YG_N(J) - HALF * DY(J)   ! south face location
            Yn = YG_N(J) + HALF * DY(JP)  ! north face location

            IF(NO_K) THEN
               Zb = ZERO                  ! bottom face location
               Zt = ZERO                  ! top face location
            ELSE
               Zb = ZG_T(K) - DZ(K)       ! bottom face location
               Zt = ZG_T(K)               ! top face location
            ENDIF

         CASE('W_MOMENTUM')
            Xw = XG_E(I) - DX(I)          ! west face location
            Xe = XG_E(I)                  ! east face location

            Ys = YG_N(J) - DY(J)          ! south face location
            Yn = YG_N(J)                  ! north face location

            Zb = ZG_T(K) - HALF * DZ(K)   ! bottom face location
            Zt = ZG_T(K) + HALF * DZ(KP)  ! top face location


         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_CELL_NODE_COORDINATES'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

!     CELL CENTER :
      X_NODE(0) = HALF * ( xw + xe )
      Y_NODE(0) = HALF * ( ys + yn )
      Z_NODE(0) = HALF * ( zb + zt )

!     NODE 1 : IM,JM,KM
      X_NODE(1) = xw
      Y_NODE(1) = ys
      Z_NODE(1) = zb

!     NODE 2 : I,JM,KM
      X_NODE(2) = xe
      Y_NODE(2) = ys
      Z_NODE(2) = zb

!     NODE 3 : IM,J,KM
      X_NODE(3) = xw
      Y_NODE(3) = yn
      Z_NODE(3) = zb

!     NODE 4 : I,J,KM
      X_NODE(4) = xe
      Y_NODE(4) = yn
      Z_NODE(4) = zb

!     NODE 5 : IM,JM,K
      X_NODE(5) = xw
      Y_NODE(5) = ys
      Z_NODE(5) = zt

!     NODE 6 : I,JM,K
      X_NODE(6) = xe
      Y_NODE(6) = ys
      Z_NODE(6) = zt

!     NODE 7 : IM,J,K
      X_NODE(7) = xw
      Y_NODE(7) = yn
      Z_NODE(7) = zt

!     NODE 8 : I,J,K
      X_NODE(8) = xe
      Y_NODE(8) = yn
      Z_NODE(8) = zt

      RETURN

      END SUBROUTINE GET_CELL_NODE_COORDINATES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_GLOBAL_CELL_NODE_COORDINATES                       C
!  Purpose: Get the cell corners (x,y,z) coordinates                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_GLOBAL_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      USE compar, ONLY: MYPE
      USE cutcell
      USE exit, only: mfix_exit
      USE functions, ONLY: FUNIJK_GL
      USE geometry, ONLY: NO_K, dx, dy, dz
      USE param1, only: half, zero
      USE vtk, ONLY: GLOBAL_I_OF, GLOBAL_J_OF, GLOBAL_K_OF

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: Xw,Xe,Yn,Ys,Zb,Zt
      INTEGER :: I,J,K,IP,JP,KP,IM,JM,KM
      INTEGER :: IJK,IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM

      I = GLOBAL_I_OF(IJK)
      J = GLOBAL_J_OF(IJK)
      K = GLOBAL_K_OF(IJK)

      IP = I + 1
      JP = J + 1
      KP = K + 1

      IM = I - 1
      JM = J - 1
      KM = K - 1

      IMJK   = FUNIJK_GL(IM,J,K)
      IJMK   = FUNIJK_GL(I,JM,K)
      IJKM   = FUNIJK_GL(I,J,KM)

      IMJMK  = FUNIJK_GL(IM,JM,K)
      IMJKM  = FUNIJK_GL(IM,J,KM)
      IJMKM  = FUNIJK_GL(I,JM,KM)

      IMJMKM = FUNIJK_GL(IM,JM,KM)

      IJK_OF_NODE(1) = IMJMKM
      IJK_OF_NODE(2) = IJMKM
      IJK_OF_NODE(3) = IMJKM
      IJK_OF_NODE(4) = IJKM
      IJK_OF_NODE(5) = IMJMK
      IJK_OF_NODE(6) = IJMK
      IJK_OF_NODE(7) = IMJK
      IJK_OF_NODE(8) = IJK

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')
            Xw = XG_E(I) - DX(I)          ! west face location
            Xe = XG_E(I)                  ! east face location

            Ys = YG_N(J) - DY(J)          ! south face location
            Yn = YG_N(J)                  ! north face location

            IF(NO_K) THEN
               Zb = ZERO                  ! bottom face location
               Zt = ZERO                  ! top face location
            ELSE
               Zb = ZG_T(K) - DZ(K)       ! bottom face location
               Zt = ZG_T(K)               ! top face location
            ENDIF

         CASE('U_MOMENTUM')
            Xw = XG_E(I) - HALF * DX(I)   ! west face location
            Xe = XG_E(I) + HALF * DX(IP)  ! east face location

            Ys = YG_N(J) - DY(J)          ! south face location
            Yn = YG_N(J)                  ! north face location

            IF(NO_K) THEN
               Zb = ZERO                  ! bottom face location
               Zt = ZERO                  ! top face location
            ELSE
               Zb = ZG_T(K) - DZ(K)       ! bottom face location
               Zt = ZG_T(K)               ! top face location
            ENDIF

         CASE('V_MOMENTUM')
            Xw = XG_E(I) - DX(I)          ! west face location
            Xe = XG_E(I)                  ! east face location

            Ys = YG_N(J) - HALF * DY(J)   ! south face location
            Yn = YG_N(J) + HALF * DY(JP)  ! north face location

            IF(NO_K) THEN
               Zb = ZERO                  ! bottom face location
               Zt = ZERO                  ! top face location
            ELSE
               Zb = ZG_T(K) - DZ(K)       ! bottom face location
               Zt = ZG_T(K)               ! top face location
            ENDIF

         CASE('W_MOMENTUM')
            Xw = XG_E(I) - DX(I)          ! west face location
            Xe = XG_E(I)                  ! east face location

            Ys = YG_N(J) - DY(J)          ! south face location
            Yn = YG_N(J)                  ! north face location

            Zb = ZG_T(K) - HALF * DZ(K)   ! bottom face location
            Zt = ZG_T(K) + HALF * DZ(KP)  ! top face location


         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_CELL_NODE_COORDINATES'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

!     CELL CENTER :
      X_NODE(0) = HALF * ( xw + xe )
      Y_NODE(0) = HALF * ( ys + yn )
      Z_NODE(0) = HALF * ( zb + zt )

!     NODE 1 : IM,JM,KM
      X_NODE(1) = xw
      Y_NODE(1) = ys
      Z_NODE(1) = zb

!     NODE 2 : I,JM,KM
      X_NODE(2) = xe
      Y_NODE(2) = ys
      Z_NODE(2) = zb

!     NODE 3 : IM,J,KM
      X_NODE(3) = xw
      Y_NODE(3) = yn
      Z_NODE(3) = zb

!     NODE 4 : I,J,KM
      X_NODE(4) = xe
      Y_NODE(4) = yn
      Z_NODE(4) = zb

!     NODE 5 : IM,JM,K
      X_NODE(5) = xw
      Y_NODE(5) = ys
      Z_NODE(5) = zt

!     NODE 6 : I,JM,K
      X_NODE(6) = xe
      Y_NODE(6) = ys
      Z_NODE(6) = zt

!     NODE 7 : IM,J,K
      X_NODE(7) = xw
      Y_NODE(7) = yn
      Z_NODE(7) = zt

!     NODE 8 : I,J,K
      X_NODE(8) = xe
      Y_NODE(8) = yn
      Z_NODE(8) = zt

      RETURN

      END SUBROUTINE GET_GLOBAL_CELL_NODE_COORDINATES
