!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  GET_DEL_H                                             C
!  Purpose: Finds the normal distance and unit vector from a cut face  C
!  to any point (x0,y0,z0)                                             C
!  The unit normal vector points away from the boundary,               C
!  towards the fluid.                                                  C
!  This subroutine must be called from a cut-cell                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_DEL_H(IJK,TYPE_OF_CELL,X0,Y0,Z0,Del_H,Nx,Ny,Nz)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: X0,Y0,Z0,XREF,YREF,ZREF
      INTEGER :: IJK,I,J,K
      DOUBLE PRECISION :: Del_H,Diagonal
      DOUBLE PRECISION :: Nx,Ny,Nz

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            IF(.NOT.CUT_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' SCALAR CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_S(IJK,1)
            Ny = NORMAL_S(IJK,2)
            Nz = NORMAL_S(IJK,3)

            Xref = REFP_S(IJK,1)
            Yref = REFP_S(IJK,2)
            Zref = REFP_S(IJK,3)

         CASE('U_MOMENTUM')

            IF(.NOT.CUT_U_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' U-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_U(IJK,1)
            Ny = NORMAL_U(IJK,2)
            Nz = NORMAL_U(IJK,3)

            Xref = REFP_U(IJK,1)
            Yref = REFP_U(IJK,2)
            Zref = REFP_U(IJK,3)

            IF(WALL_U_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF

         CASE('V_MOMENTUM')

            IF(.NOT.CUT_V_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' V-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_V(IJK,1)
            Ny = NORMAL_V(IJK,2)
            Nz = NORMAL_V(IJK,3)

            Xref = REFP_V(IJK,1)
            Yref = REFP_V(IJK,2)
            Zref = REFP_V(IJK,3)

            IF(WALL_V_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF

         CASE('W_MOMENTUM')

            IF(.NOT.CUT_W_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' W-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_W(IJK,1)
            Ny = NORMAL_W(IJK,2)
            Nz = NORMAL_W(IJK,3)

            Xref = REFP_W(IJK,1)
            Yref = REFP_W(IJK,2)
            Zref = REFP_W(IJK,3)

            IF(WALL_W_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF

         CASE DEFAULT
            WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT


      DEL_H = Nx * (X0 - Xref) + Ny * (Y0 - Yref) + Nz * (Z0 - Zref)

!======================================================================
! Negative values of DEL_H are not accepted
!======================================================================

       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)

       IF(NO_K) THEN
          Diagonal = sqrt(DX(I)**2 + DY(J)**2 )
       ELSE
          Diagonal = sqrt(DX(I)**2 + DY(J)**2 + DZ(K)**2)
       ENDIF

      IF (DEL_H <= TOL_DELH * Diagonal) THEN

         DEL_H = UNDEFINED
         Nx = ZERO
         Ny = ZERO
         Nz = ZERO

      ENDIF

      RETURN


      END SUBROUTINE GET_DEL_H

  SUBROUTINE GET_DEL_H_DES(IJK,TYPE_OF_CELL,X0,Y0,Z0,Del_H,Nx,Ny,Nz, allow_neg_dist)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: X0,Y0,Z0,XREF,YREF,ZREF
      INTEGER :: IJK,I,J,K
      DOUBLE PRECISION :: Del_H,Diagonal
      DOUBLE PRECISION :: Nx,Ny,Nz

      LOGICAL :: ALLOW_NEG_DIST

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            IF(.NOT.CUT_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
               WRITE(*,*)' SCALAR CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' I, J, K =',I_OF(IJK), J_OF(IJK), K_OF(IJK)
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_S(IJK,1)
            Ny = NORMAL_S(IJK,2)
            Nz = NORMAL_S(IJK,3)

            Xref = REFP_S(IJK,1)
            Yref = REFP_S(IJK,2)
            Zref = REFP_S(IJK,3)

         CASE('U_MOMENTUM')

            IF(.NOT.CUT_U_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
               WRITE(*,*)' U-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_U(IJK,1)
            Ny = NORMAL_U(IJK,2)
            Nz = NORMAL_U(IJK,3)

            Xref = REFP_U(IJK,1)
            Yref = REFP_U(IJK,2)
            Zref = REFP_U(IJK,3)

            IF(WALL_U_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF

         CASE('V_MOMENTUM')

            IF(.NOT.CUT_V_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
               WRITE(*,*)' V-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_V(IJK,1)
            Ny = NORMAL_V(IJK,2)
            Nz = NORMAL_V(IJK,3)

            Xref = REFP_V(IJK,1)
            Yref = REFP_V(IJK,2)
            Zref = REFP_V(IJK,3)

            IF(WALL_V_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF

         CASE('W_MOMENTUM')

            IF(.NOT.CUT_W_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
               WRITE(*,*)' W-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'
               CALL MFIX_EXIT(myPE)
            ENDIF

            Nx = NORMAL_W(IJK,1)
            Ny = NORMAL_W(IJK,2)
            Nz = NORMAL_W(IJK,3)

            Xref = REFP_W(IJK,1)
            Yref = REFP_W(IJK,2)
            Zref = REFP_W(IJK,3)

            IF(WALL_W_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF

         CASE DEFAULT
            WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

      DEL_H = Nx * (X0 - Xref) + Ny * (Y0 - Yref) + Nz * (Z0 - Zref)

!======================================================================
! Negative values of DEL_H are not accepted
!======================================================================

       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)

       IF(.NOT.ALLOW_NEG_DIST) THEN
          Diagonal = sqrt(DX(I)**2 + DY(J)**2 + DZ(K)**2)

          IF (DEL_H <= TOL_DELH * Diagonal) THEN
             DEL_H = UNDEFINED
             Nx = ZERO
             Ny = ZERO
             Nz = ZERO

          ENDIF
       ENDIF

      RETURN

      END SUBROUTINE GET_DEL_H_DES
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_CUT_FACE_INFO                                    C
!  Purpose: Compute and store unit normal vector and reference point   C
!           Defining a cut face                                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE STORE_CUT_FACE_INFO(IJK,TYPE_OF_CELL,N_CUT_FACE_NODES,COORD_CUT_FACE_NODES,X_MEAN,Y_MEAN,Z_MEAN)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE functions

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK
      INTEGER :: N_CUT_FACE_NODES
      DOUBLE PRECISION, DIMENSION(3,15) :: COORD_CUT_FACE_NODES
      DOUBLE PRECISION :: X_MEAN,Y_MEAN,Z_MEAN
      DOUBLE PRECISION, DIMENSION(3) :: NN,V,N1,N2
      DOUBLE PRECISION :: NORM,N1_dot_N2
      DOUBLE PRECISION, DIMENSION(3) :: VECTEMP,VECTEMP2

     IF(NO_K) THEN  ! 2D case

        NN(1) = COORD_CUT_FACE_NODES(2,1) - COORD_CUT_FACE_NODES(2,2) ! y1-y2
        NN(2) = COORD_CUT_FACE_NODES(1,2) - COORD_CUT_FACE_NODES(1,1) ! x2-x1
        NN(3) = ZERO

     ELSE

!======================================================================
! Make sure there are a least three points along the plane
!======================================================================

        IF(N_CUT_FACE_NODES < 3) THEN
           WRITE(*,*)' ERROR IN SUBROUTINE STORE_CUT_FACE_INFO:'
           WRITE(*,*)' CUT FACE HAS LESS THAN 3 NODES.'
           WRITE(*,*)' MFIX WILL EXIT NOW.'
           CALL MFIX_EXIT(myPE)
        END IF

!======================================================================
!  Find tentative unit normal vector
!  and reverse sign if necessary
!  (unit vector must be pointing towards the fluid)
!======================================================================

        VECTEMP  = COORD_CUT_FACE_NODES(:,2)-COORD_CUT_FACE_NODES(:,1)
        VECTEMP2 = COORD_CUT_FACE_NODES(:,3)-COORD_CUT_FACE_NODES(:,1)
        NN = CROSS_PRODUCT(VECTEMP,VECTEMP2)

     ENDIF

      NORM = sqrt(dot_product(nn(:),nn(:)))
      NN = NN / NORM

      V(1) = X_MEAN - COORD_CUT_FACE_NODES(1,1)
      V(2) = Y_MEAN - COORD_CUT_FACE_NODES(2,1)
      V(3) = Z_MEAN - COORD_CUT_FACE_NODES(3,1)

      IF (DOT_PRODUCT(NN,V) < ZERO) NN = - NN

      IF(N_CUT_FACE_NODES > 3) THEN     ! FOR 3D geometry, check normal of plane defined by nodes 1,2, and 4

         N1 = NN  ! Keep copy of previous N (nodes 1,2,3)

        VECTEMP  = COORD_CUT_FACE_NODES(:,2)-COORD_CUT_FACE_NODES(:,1)
        VECTEMP2 = COORD_CUT_FACE_NODES(:,4)-COORD_CUT_FACE_NODES(:,1)

         N2 = CROSS_PRODUCT(VECTEMP,VECTEMP2)

         NORM = sqrt(dot_product(n2(:),n2(:)))
         N2 = N2 / NORM

         V(1) = X_MEAN - COORD_CUT_FACE_NODES(1,1)
         V(2) = Y_MEAN - COORD_CUT_FACE_NODES(2,1)
         V(3) = Z_MEAN - COORD_CUT_FACE_NODES(3,1)

         IF (DOT_PRODUCT(N2,V) < ZERO) N2 = - N2

      ENDIF

      N1_dot_N2 = DOT_PRODUCT(N1,N2)
      DEBUG_CG(IJK,1)=N1_dot_N2

      IF(N1_dot_N2<0.99) THEN

!         What should be done when the unit vectors are different ?

      ENDIF

!======================================================================
! Store unit normal vector and reference point
!======================================================================

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            NORMAL_S(IJK,:) = NN
            REFP_S(IJK,:)   = COORD_CUT_FACE_NODES(:,1)

            IF(DO_K) CALL TEST_DEL_H(IJK,'SCALAR') ! test for negative del_H

         CASE('U_MOMENTUM')

            NORMAL_U(IJK,:) = NN
            REFP_U(IJK,:)   = COORD_CUT_FACE_NODES(:,1)

         CASE('V_MOMENTUM')

            NORMAL_V(IJK,:) = NN
            REFP_V(IJK,:)   = COORD_CUT_FACE_NODES(:,1)

         CASE('W_MOMENTUM')

            NORMAL_W(IJK,:) = NN
            REFP_W(IJK,:)   = COORD_CUT_FACE_NODES(:,1)

         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: STORE_CUT_FACE_INFO'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:'
            WRITE(*,*)'SCALAR'
            WRITE(*,*)'U_MOMENTUM'
            WRITE(*,*)'V_MOMENTUM'
            WRITE(*,*)'W_MOMENTUM'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN

      END SUBROUTINE STORE_CUT_FACE_INFO

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TEST_DEL_H                                             C
!  Purpose: tests the computation of wall distance                     C
!           If a negative distance is detected, the normal vector      C
!           is inverted                                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 22-Feb-12  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE TEST_DEL_H(IJK,TYPE_OF_CELL)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell

      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: X_COPY,Y_COPY,Z_COPY
      INTEGER :: IJK
      INTEGER :: NODE,N_N1,N_N2,NN
      DOUBLE PRECISION :: Del_H
      DOUBLE PRECISION :: Nx,Ny,Nz

      LOGICAL :: ALLOW_NEG_DIST = .TRUE.  ! forces GET_DEL_H_DES to output negative delh
                                           ! i.e. do not let the subroutine overwrite negative values

! This subroutine tests values of del_H for nodes defining a cut cell.
! Only nodes that are in the fluid region are tested.
! Nodes belonging to the cut face should return zero (or near zero) values and are not tested.
! If a negative del_H is detected, the unit normal vector is reversed.

      IF(NO_K) THEN
         N_N1 = 5
         N_N2 = 8
      ELSE
         N_N1 = 1
         N_N2 = 8
      ENDIF

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      DO NODE = 1,NUMBER_OF_NODES(IJK)
         IF(CONNECTIVITY(IJK,NODE)<=IJKEND3) THEN  ! node does not belong to the cut-face
                                                    ! i.e. is in the fluid region

            DO NN = N_N1,N_N2                     ! get node coordinate
               IF(CONNECTIVITY(IJK,NODE) == IJK_OF_NODE(NN)) THEN
                  X_COPY = X_NODE(NN)
                  Y_COPY = Y_NODE(NN)
                  Z_COPY = Z_NODE(NN)
                  EXIT
               ENDIF
            ENDDO

!           Compute del_H
            CALL GET_DEL_H_DES(IJK,TYPE_OF_CELL,X_COPY,Y_COPY,Z_COPY,Del_H,Nx,Ny,Nz, ALLOW_NEG_DIST)


            IF(DEL_H<ZERO) THEN


               IF(PRINT_WARNINGS.AND.MyPE==PE_IO) THEN
                  WRITE(*,*) ' Warning: Negative delh detected in scalar cell :',IJK
                  WRITE(*,*) ' Location (X,Y,Z) = ',X_COPY,Y_COPY,Z_COPY
                  WRITE(*,*) ' Reverting unit normal vector.'
               ENDIF

               NORMAL_S(IJK,1) = -NORMAL_S(IJK,1)
               NORMAL_S(IJK,2) = -NORMAL_S(IJK,2)
               NORMAL_S(IJK,3) = -NORMAL_S(IJK,3)

            ENDIF

         ENDIF
      ENDDO

      RETURN

      END SUBROUTINE TEST_DEL_H

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  GET_DISTANCE_TO_WALL                                  C
!  Purpose: Finds the distance fraom any scalar cell to the closest    C
!  wall                                                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 16-May-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_DISTANCE_TO_WALL

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE functions

      USE mpi_utility

      IMPLICIT NONE
!      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: X0,Y0,Z0,XREF,YREF,ZREF
      INTEGER :: IJK,NN

      DOUBLE PRECISION :: D_TO_CUT, D_TO_PE_REF

      INTEGER :: N_CUT_CELLS
      INTEGER :: LIST_OF_CUT_CELLS(DIMENSION_3)

      INTEGER :: iproc,IERR,IJK_OFFSET,nb,n1,n2
      INTEGER :: GLOBAL_N_CUT_CELLS
      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount
      LOGICAL, DIMENSION(0:numPEs-1) :: ALREADY_VISITED
      DOUBLE PRECISION, DIMENSION(0:numPEs-1,3) :: PE_REFP,ALL_PE_REFP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: LOCAL_REFP_S,GLOBAL_REFP_S

      IF(MyPE==PE_IO) WRITE(*,*)'COMPUTING WALL DISTANCE...'

! Count the number of cut cells
      N_CUT_CELLS = 0
      DO IJK = IJKSTART3, IJKEND3
         IF(INTERIOR_CELL_AT(IJK).AND.CUT_CELL_AT(IJK)) THEN
            N_CUT_CELLS = N_CUT_CELLS + 1
            LIST_OF_CUT_CELLS(N_CUT_CELLS) = IJK
         ENDIF
      ENDDO

! Store cut cell reference points (centroid of cut cell face)
! and save a reference point for the entire processor
      ALLOCATE (LOCAL_REFP_S(N_CUT_CELLS,3))

      N_CUT_CELLS = 0
      PE_REFP(MyPE,1) = ZERO
      PE_REFP(MyPE,2) = ZERO
      PE_REFP(MyPE,3) = ZERO

      DO IJK = IJKSTART3, IJKEND3
         IF(INTERIOR_CELL_AT(IJK).AND.CUT_CELL_AT(IJK)) THEN
            N_CUT_CELLS = N_CUT_CELLS + 1
            LOCAL_REFP_S(N_CUT_CELLS,1) = REFP_S(IJK,1)
            LOCAL_REFP_S(N_CUT_CELLS,2) = REFP_S(IJK,2)
            LOCAL_REFP_S(N_CUT_CELLS,3) = REFP_S(IJK,3)

            PE_REFP(MyPE,1) = PE_REFP(MyPE,1) + REFP_S(IJK,1)
            PE_REFP(MyPE,2) = PE_REFP(MyPE,2) + REFP_S(IJK,2)
            PE_REFP(MyPE,3) = PE_REFP(MyPE,3) + REFP_S(IJK,3)
         ENDIF
      ENDDO

      IF(N_CUT_CELLS>0) THEN
         PE_REFP(MyPE,1) = PE_REFP(MyPE,1) / N_CUT_CELLS
         PE_REFP(MyPE,2) = PE_REFP(MyPE,2) / N_CUT_CELLS
         PE_REFP(MyPE,3) = PE_REFP(MyPE,3) / N_CUT_CELLS
      ELSE
         PE_REFP(MyPE,1) = UNDEFINED
         PE_REFP(MyPE,2) = UNDEFINED
         PE_REFP(MyPE,3) = UNDEFINED
      ENDIF

!======================================================================
! Now, gather the local reference points to head node
! to get a global list, and broadcast it to each processor.
!======================================================================

!======================================================================
! First get the offset and build the rcount and disp arrays.
! rcount is the number of elements to be gathered.
! disp is the displacement of the variable size gather,
! i.e. the cumulative sum at a given procesor.
!======================================================================

      CALL allgather_1i (N_CUT_CELLS,rcount,IERR)

      IF (myPE == 0) THEN
         IJK_OFFSET = 0
      ELSE
         IJK_OFFSET = 0
         DO iproc=0,myPE-1
            IJK_OFFSET = IJK_OFFSET + rcount(iproc)
         ENDDO
      ENDIF

      CALL allgather_1i (IJK_OFFSET,disp,IERR)

      CALL allgather_1d (PE_REFP(MyPE,1),ALL_PE_REFP(:,1),IERR)
      CALL allgather_1d (PE_REFP(MyPE,2),ALL_PE_REFP(:,2),IERR)
      CALL allgather_1d (PE_REFP(MyPE,3),ALL_PE_REFP(:,3),IERR)

!======================================================================
! Get the global number of cut cells and allocate the
! global reference point array. Each processor gets its own
! copy of all cut cell reference points !!
!======================================================================
      CALL GLOBAL_ALL_SUM(N_CUT_CELLS, GLOBAL_N_CUT_CELLS,  PE_IO, ierr )
      ALLOCATE (GLOBAL_REFP_S(GLOBAL_N_CUT_CELLS,3))

!======================================================================
! For a serial run, the global and local arrays are the same.
! for a parallel run, first gather on head node, then broadcast to all.
!======================================================================

      IF(numPEs==1) THEN  ! Serial run
         GLOBAL_REFP_S =  LOCAL_REFP_S
      ELSE !Parallel run
         call gatherv_1d( LOCAL_REFP_S(:,1), N_CUT_CELLS, GLOBAL_REFP_S(:,1), rcount, disp, PE_IO, ierr )
         call gatherv_1d( LOCAL_REFP_S(:,2), N_CUT_CELLS, GLOBAL_REFP_S(:,2), rcount, disp, PE_IO, ierr )
         call gatherv_1d( LOCAL_REFP_S(:,3), N_CUT_CELLS, GLOBAL_REFP_S(:,3), rcount, disp, PE_IO, ierr )

         call bcast(GLOBAL_REFP_S(:,1))
         call bcast(GLOBAL_REFP_S(:,2))
         call bcast(GLOBAL_REFP_S(:,3))
      ENDIF


      ALREADY_VISITED(:) = .FALSE.
!======================================================================
!  Loop through all scalar cells, grouped by processor.
!  Compute the distance to each cut cell and take the minimum
!  This is doe in three passes, in order of likelyhood to find a wall:
!  1) Within local processor.
!  2) Within neighboring processors (from send2 schedule)
!  3) Within all other processors
!
!  The idea is that it is very likely that cut cells located
!  on current or neighbor processors will contribute to the
!  wall distance calculation, whereas cut cells located
!  on the other processors will not contribute
!  During the third pass, the entire processor is skipped
!  if a reference point (average of all cut faces reference
!  points) is farther than the current wall distance.
!  This should speed up the brute calculation of going
!  explicitely through all cut cells.
!  However, it is not guaranteed that dwall will be
!  independent of the grid partitionning. This is considered
!  very unlikely at the moment
!  Change DWALL_BRUTE_FORCE to .TRUE. to force brute-force
!  calculation.
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN
            DWALL(IJK) = UNDEFINED

            IF(.NOT.BLOCKED_CELL_AT(IJK)) THEN

! Get coordinates of cell center
               CALL GET_CELL_NODE_COORDINATES(IJK,'SCALAR')

               X0 = X_NODE(0)
               Y0 = Y_NODE(0)
               Z0 = Z_NODE(0)

!======================================================================
! First pass: Loop through local cut cells
!======================================================================

               ALREADY_VISITED(MyPE) = .TRUE.

               DO NN = 1,N_CUT_CELLS

                  Xref = LOCAL_REFP_S(NN,1)
                  Yref = LOCAL_REFP_S(NN,2)
                  Zref = LOCAL_REFP_S(NN,3)

                  D_TO_CUT = sqrt((X0 - Xref)**2 + (Y0 - Yref)**2 + (Z0 - Zref)**2)

                  DWALL(IJK) = MIN(DWALL(IJK),D_TO_CUT)

               ENDDO


!======================================================================
! Second pass: Loop through neighbor processors (use send2 schedule)
!======================================================================

               DO nb=1,nsend2
                  iproc = sendproc2(nb)
                  ALREADY_VISITED(iproc) = .TRUE.
                  n1 = disp(iproc)+1
                  n2 = n1 + rcount(iproc) - 1

                  DO NN = n1,n2

                     Xref = GLOBAL_REFP_S(NN,1)
                     Yref = GLOBAL_REFP_S(NN,2)
                     Zref = GLOBAL_REFP_S(NN,3)

                     D_TO_CUT = sqrt((X0 - Xref)**2 + (Y0 - Yref)**2 + (Z0 - Zref)**2)

                     DWALL(IJK) = MIN(DWALL(IJK),D_TO_CUT)

                  ENDDO

               ENDDO

!======================================================================
! Third pass: Loop through all other processors
! skip already visited processors (myPE and its neigbhors)
! First test if the average reference point is farther than dwall
! if this is true, then skip the entire processor cut cells.
!======================================================================


               DO iproc=0,numPEs-1
                  IF(ALREADY_VISITED(iproc)) CYCLE

                  Xref = ALL_PE_REFP(iproc,1)
                  Yref = ALL_PE_REFP(iproc,2)
                  Zref = ALL_PE_REFP(iproc,3)

                  D_TO_PE_REF = sqrt((X0 - Xref)**2 + (Y0 - Yref)**2 + (Z0 - Zref)**2)

                  IF((DWALL(IJK) < D_TO_PE_REF).AND. &
                     (.NOT.DWALL_BRUTE_FORCE  )) THEN
                     CYCLE
                  ELSE
                     n1 = disp(iproc)+1
                     n2 = n1 + rcount(iproc) - 1

                     DO NN = n1,n2

                        Xref = GLOBAL_REFP_S(NN,1)
                        Yref = GLOBAL_REFP_S(NN,2)
                        Zref = GLOBAL_REFP_S(NN,3)

                        D_TO_CUT = sqrt((X0 - Xref)**2 + (Y0 - Yref)**2 + (Z0 - Zref)**2)

                        DWALL(IJK) = MIN(DWALL(IJK),D_TO_CUT)

                     ENDDO
                  ENDIF
               ENDDO

            ENDIF

         ENDIF

      ENDDO

      call send_recv(DWALL,2)
!======================================================================
!  Deallocate arrays before leaving
!======================================================================
      DEALLOCATE (LOCAL_REFP_S)
      DEALLOCATE (GLOBAL_REFP_S)
      RETURN

      END SUBROUTINE GET_DISTANCE_TO_WALL
