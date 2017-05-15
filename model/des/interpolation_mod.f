! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INTERPOLATION                                          C
!  Purpose: Contains number of generic interfaces that are used if
!           DES_INTERP_ON is TRUE
!                                                                      C
!                                                                      C
!  Author: Chidambaram Narayanan and Rahul Garg       Date: 01-Aug-07  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE interpolation

  USE constant
  USE discretelement
  USE param1, only: half, one, zero
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: set_interpolation_stencil,  set_interpolation_scheme, set_interpolation_pstencil, array_dot_product
  PUBLIC:: INTERPOLATE_CC


  PUBLIC:: interpolator
  INTERFACE interpolator
     MODULE PROCEDURE interp_oned_scalar
     MODULE PROCEDURE interp_oned_vector
     MODULE PROCEDURE interp_threed_scalar
     MODULE PROCEDURE interp_threed_vector
     MODULE PROCEDURE interp_twod_scalar
     MODULE PROCEDURE interp_twod_vector
     MODULE PROCEDURE calc_weightderiv_threed
     MODULE PROCEDURE calc_weightderiv_twod
  END INTERFACE

  INTERFACE array_dot_product
      Module procedure array_dot_product_2d
      Module procedure array_dot_product_3d
  END INTERFACE

  SAVE

! the interpolator subroutine essentially returns the value of interp_scl/interp_vec
! which is the interpolated value of the quantity passed through 'scalar/vector'.
! the value is interpolated based on grid (x,y,z position of grid)  passed
! through coor at the x,y,z position of the particle passed through ppos
! the interpolation order is specified by order and the interpolation scheme is
! specified by isch/fun

! interp_oned_scalar  (coor,scalar,ppos,interp_scl,order, isch,weight_pointer)
!    REAL coor(:), scalar(:), PPOS, INTERP_SCL
!    INTEGER ORDER; CHARACTER ISCH; REAL WEIGHT_POINTER(:)
! interp_oned_vector  (coor,vector,ppos,interp_vec,order, fun, weight_pointer)
!    REAL coor(:), vector(:,:), PPOS, INTERP_VEC(:)
!    INTEGER ORDER; CHARACTER FUN; REAL WEIGHT_POINTER(:)
! interp_twod_scalar  (coor,scalar,ppos,interp_scl,order, isch,weight_pointer)
!    REAL coor(:,:,:), scalar(:,:), PPOS(2), INTERP_SCL
!    INTEGER ORDER; CHARACTER ISCH; REAL WEIGHT_POINTER(:,:,:)
! interp_twod_vector  (coor,vector,ppos,interp_vec,order, fun, weight_pointer)
!    REAL coor(:,:,:), vector(:,:,:), PPOS(2), INTERP_VEC(:)
!    INTEGER ORDER; CHARACTER FUN; REAL WEIGHT_POINTER(:,:,:)
! interp_threed_scalar(coor,scalar,ppos,interp_scl,order, isch,weight_pointer)
!    REAL coor(:,:,:,:), scalar(:,:,:), PPOS(3), INTERP_SCL
!    INTEGER ORDER; CHARACTER ISCH; REAL WEIGHT_POINTER(:,:,:)
! interp_threed_vector(coor,vector,ppos,interp_vec,order, fun, weight_pointer)
!    REAL coor(:,:,:,:), vector(:,:,:,:), PPOS(3), INTERP_VEC(:)
!    INTEGER ORDER; CHARACTER FUN; REAL WEIGHT_POINTER(:,:,:)

  INTEGER, PARAMETER :: prcn=8
  INTEGER, PARAMETER :: iprec=8
  !double precision, Parameter:: zero = 0.0_iprec
  !double precision, Parameter:: one  = 1.0_iprec
  DOUBLE PRECISION, PARAMETER  :: two = 2.0_iprec
  DOUBLE PRECISION, PARAMETER  :: three = 3.0_iprec
  DOUBLE PRECISION, PARAMETER  :: four = 4.0_iprec
  DOUBLE PRECISION, PARAMETER :: six = 6.0_iprec

  !double precision, Parameter:: half = 0.5_iprec
  DOUBLE PRECISION, PARAMETER  :: fourth = 0.25_iprec

  INTEGER, PARAMETER :: maxorder=6
  DOUBLE PRECISION, DIMENSION(maxorder), TARGET :: xval, yval, zval
  DOUBLE PRECISION, DIMENSION(maxorder-1) :: dx, dy, dz
  DOUBLE PRECISION, DIMENSION(maxorder,maxorder,maxorder), TARGET :: weights



 CONTAINS

        !----------
  ! A . B where A and B are arrays of the same dimensions
  !----------
  !----------
  FUNCTION array_dot_product_3d(A,B)  !3D Arrays
    Real(prcn):: array_dot_product_3d

    Integer:: j, k, s1, s2, s3
    Real(prcn), Dimension(:,:,:), Intent(in):: A, B

    s1 = SIZE(A,1)
    s2 = SIZE(A,2)
    s3 = SIZE(A,3)

    array_dot_product_3d = zero
    DO k = 1,s3
       DO j = 1,s2
          array_dot_product_3d = array_dot_product_3d + dot_product(A(:,j,k), B(1:s1,j,k))
       ENDDO
    ENDDO
  END FUNCTION array_dot_product_3d

  FUNCTION array_dot_product_2d(A,B)  !2D Arrays
    Real(prcn):: array_dot_product_2d

    Integer:: j, s2
    Real(prcn), Dimension(:,:), Intent(in):: A, B

    s2 = SIZE(A,2)

    array_dot_product_2d = zero
    DO j = 1,s2
       array_dot_product_2d = array_dot_product_2d + dot_product(A(:,j), B(:,j))
    ENDDO
  END FUNCTION array_dot_product_2d

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE set_interpolation_stencil(PC, IW, IE, JS, JN, KB, KTP,&
      isch,dimprob,ordernew)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       !USE discretelement, ONLY : order,ob2l,ob2r, des_periodic_walls_x,y,z

    USE geometry

    IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    INTEGER, DIMENSION(3), INTENT(in):: pc   ! i,j,k indices of particle - 1
    INTEGER, INTENT(out):: IW, IE, JS, JN, KB, KTP
    CHARACTER(LEN=5), INTENT(in) :: isch   ! interpolation scheme
    INTEGER, INTENT(in) :: dimprob   ! dimension of system = DIMN
    INTEGER, OPTIONAL :: ordernew   ! interpolation order

    INTEGER :: im, jm, km   ! local variables assigned maximum i,j,k fluid cell indices
    INTEGER :: ob2rtmp, ob2ltmp, ordertemp, ordernewtmp
!-----------------------------------------------

! reassigning ob2l and ob2r ?
    ob2l = (order+1)/2
    ob2r = order/2
! local variables for maximum fluid cell indices
    im = imax1
    jm = jmax1
    km = kmax1

    SELECT CASE(isch)
    CASE('csi')

       ob2rtmp = ob2r
       ob2ltmp = ob2l
       ordertemp = order

! lowest IW will be assigned is 1 (ghost cell)
       IW = MAX(1 ,pc(1) - (ob2l - 1))
! highest IE will be assigned is maximum fluid cell index
       IE = MIN(im,pc(1) + ob2r)
! if IW is west ghost cell and/or IE is maximum fluid cell
! reassign IW and/or IE accordingly
       IF(.NOT.DES_PERIODIC_WALLS_X) THEN
          IF (IW.EQ.1 ) IE = IW + order - 1
          IF (IE.EQ.im) IW = IE - order + 1
       ELSE
          IF (IW.EQ.1 ) IW = IE - order + 1
          IF (IE.EQ.im) IE = IW + order - 1
       ENDIF

       JS = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
       JN = MIN(jm,pc(2) + ob2r)
       IF(.NOT.DES_PERIODIC_WALLS_Y) THEN
          IF (JS.EQ.1 ) JN = JS + order - 1
          IF (JN.EQ.jm) JS = JN - order + 1
       ELSE
          IF (JS.EQ.1 ) JS = JN - order + 1
          IF (JN.EQ.jm) JN = JS + order - 1
       ENDIF

       KB = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
       KTP = MIN(km,pc(3) + ob2r)
       IF(.NOT.DES_PERIODIC_WALLS_Z) THEN
          IF (KB.EQ.1 ) KTP = KB + order - 1
          IF (KTP.EQ.km) KB = KTP - order + 1
       ELSE
          IF (KB.EQ.1 ) KB = KTP - order + 1
          IF (KTP.EQ.km) KTP = KB + order - 1
       ENDIF

       ob2r =  ob2rtmp
       ob2l = ob2ltmp
       ordernewtmp = order
       order = ordertemp !reset the order

    CASE('lpi')

       IW = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
       IE = MIN(im,pc(1) + ob2r)
       IF(.NOT.DES_PERIODIC_WALLS_X) THEN
          IF (IW.EQ.1 ) IE = IW + order - 1
          IF (IE.EQ.im) IW = IE - order + 1
       ELSE
          IF (IW.EQ.1 ) IW = IE - order + 1
          IF (IE.EQ.im) IE = IW + order - 1
       ENDIF

       JS = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
       JN = MIN(jm,pc(2) + ob2r)
       IF(.NOT.DES_PERIODIC_WALLS_Y) THEN
          IF (JS.EQ.1 ) JN = JS + order - 1
          IF (JN.EQ.jm) JS = JN - order + 1
       ELSE
          IF (JS.EQ.1 ) JS = JN - order + 1
          IF (JN.EQ.jm) JN = JS + order - 1
       ENDIF

       KB = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
       KTP = MIN(km,pc(3) + ob2r)
       IF(.NOT.DES_PERIODIC_WALLS_Z) THEN
          IF (KB.EQ.1 ) KTP = KB + order - 1
          IF (KTP.EQ.km) KB = KTP - order + 1
       ELSE
          IF (KB.EQ.1 ) KB = KTP - order + 1
          IF (KTP.EQ.km) KTP = KB + order - 1
       ENDIF
       ordernewtmp = order

    END SELECT

    IF(dimprob == 2) THEN
       KB = pc(3)
       KTP = pc(3)
    ENDIF

! for debugging
    !print*, 'order in set stencil = ', order,pc(3)
    !Print*,'IW, IE = ', pc(1), IW, IE

    IF(PRESENT(ordernew)) ordernew = ordernewtmp

  END SUBROUTINE set_interpolation_stencil


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: SET_INTERPOLATION_STENCIL_CC                            !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_INTERPOLATION_STENCIL_CC(NP, PC, IW, JS, KB, FOCUS)

      USE discretelement
      USE geometry
      USE param1

      IMPLICIT NONE

! I, J, K indices of the cell containing the particle.
      INTEGER, INTENT(IN) :: PC(3)
! particle number
      INTEGER, INTENT(IN) :: NP

! These values indicate the I,J,K starting indices used to select the
! nodes for interpolation.
      INTEGER, INTENT(OUT) :: IW, JS, KB
! flag to tell whether local debugging was called in drag_fgs
      LOGICAL, INTENT(IN) :: FOCUS

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER I_SHIFT, J_SHIFT, K_SHIFT, IE, JN, KTP

! Determine the shift index necessary to account for the quadrant of
! the fluid cell where the particle is located.
      IF( DES_POS_NEW(NP,1) <= XE(PC(1))-HALF*DX(PC(1)))THEN
         I_SHIFT = 0
      ELSE
         I_SHIFT = 1
      ENDIF
      IF( DES_POS_NEW(NP,2) <= YN(PC(2))-HALF*DY(PC(2)))THEN
         J_SHIFT = 0
      ELSE
         J_SHIFT = 1
      ENDIF
      IF(DIMN == 2) THEN
         K_SHIFT = 0
      ELSE
         IF( DES_POS_NEW(NP,3) <= ZT(PC(3))-HALF*DZ(PC(3)))THEN
            K_SHIFT = 0
         ELSE
            K_SHIFT = 1
         ENDIF
      ENDIF

! Calculate the west and east I indices to be used.
!-----------------------------------------------------------------------
! Restrict the west I index to a minimum of 1
      IW = MAX( 1, (PC(1) + I_SHIFT-1))
! Restrict the east I index to a maximum of IMAX1
      IE = MIN( IMAX2, PC(1) + I_SHIFT)
      IF(.NOT.DES_PERIODIC_WALLS_X) THEN
! If the YZ walls are not periodic, ensure that the indices used for
! interpolation stay within the domain. This may require a shift in the
! indices for cells near the walls and interpolation schemes of high
! order.
         IF (IE .EQ. IMAX2) IW = IMAX2 - 1
      ELSE
! If the YZ walls are periodic, the starting index can been less than 1.
! This is accounted for when filling GSTENCIL. corrected later.
         IF (IW .EQ. 1) IW = IE - 1
      ENDIF

! Calculate the south and north J indices to be used.
!-----------------------------------------------------------------------
! Restrict the south J index to a minimum of 1
      JS = MAX( 1, PC(2)+ J_SHIFT - 1)
! Restrict the north J index to a maximum of JMAX1
      JN = MIN( JMAX2, PC(2) + J_SHIFT)
      IF(.NOT.DES_PERIODIC_WALLS_Y) THEN
! Shift indices for non periodic boundaries if needed.
         IF (JN .EQ. JMAX2) JS = JMAX2 - 1
      ELSE
! Shift indices for periodic boundaries if needed.
         IF (JS .EQ. 1 ) JS = JN - 1
      ENDIF

! Calculate the bottom and top K indices to be used.
!-----------------------------------------------------------------------
      IF(DIMN == 2) THEN
! If 2 dimension simulation, set the bottom index to 1.
         KB = 1
      ELSE
! Restric the bottom K index to a minimum  of 1
         KB = MAX( 1, PC(3)+K_SHIFT-1)
! Restric the top K index to a maximum  of KMAX1
         KTP = MIN( KMAX2, PC(3) + K_SHIFT)
         IF(.NOT.DES_PERIODIC_WALLS_Z) THEN
! Shift indices for non periodic boundaries if needed.
            IF (KTP .EQ. KMAX1) KB = KTP - 1
         ELSE
! Shift indices for periodic boundaries if needed.
            IF (KB .EQ. 1 ) KB = KTP - 1
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE SET_INTERPOLATION_STENCIL_CC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: INTERPOLATE_CC                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jan-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INTERPOLATE_CC(NP, INTP_IJK, INTP_WEIGHTS, FOCUS)

      USE compar
      USE discretelement
      USE geometry
      USE indices
      USE param
      USE param1
      USE functions

      IMPLICIT NONE

! Particle index number
      INTEGER, INTENT(IN)  :: NP
! Local debugging flag
      LOGICAL, INTENT(IN)  :: FOCUS

! IJK values of cells. These have already been adjusted for periodic
! and non-periodic boundary conditions.
      INTEGER, INTENT(INOUT)  :: INTP_IJK(2**DIMN)
! Weights associated with the IJK cells in INTP_IJK
      DOUBLE PRECISION, INTENT(INOUT)  :: INTP_WEIGHTS(2**DIMN)

!-----------------------------------------------
! Local variables
!-----------------------------------------------

! Starting grid indices used in interpolation
      INTEGER IW, JS, KB

! Weights associated with the nodes used for interpolation. These values
! only need to be returned with interpolating routine is invoked from
! the continuum phase.
      DOUBLE PRECISION, POINTER :: WEIGHTS_CC (:,:,:)
! indices
      INTEGER I, J, K, II, JJ, KK
! Loop counter
      INTEGER LC
! geometery and field stencils for cell-center interpolations
      DOUBLE PRECISION STNCL_GEMTY (2, 2, MAX(1, 2*(DIMN-2)), DIMN)
      DOUBLE PRECISION STNCL_FEILD (2, 2, MAX(1, 2*(DIMN-2)))
! Generic value. In the original interpolation routine, this would be
! the value that was interpolated within the code. However, this routine
! only needs to return the IJK indicies of the cells and the interpolation
! weights. In order to satisfy the function interface of "interpolator"
! is to send it is dummy value.
      DOUBLE PRECISION INTERP_scalar
! Interpolatin scheme. It is 'hard coded' so that only a second order
! Lagrange polynomial is used. Higher order isn't necessary.
      CHARACTER(LEN=5) :: INTP_SCHM = 'LPI'
!-----------------------------------------------

! Obtain the starting cell index values for interpolation around
! particle NP
      CALL SET_INTERPOLATION_STENCIL_CC(NP, PIJK(NP,1:3), IW, JS, KB, &
         FOCUS)

      LC = 0
      DO K = 1, (3-DIMN)*1+(DIMN-2)*2
         DO J = 1, 2
            DO I = 1, 2
! increment the loop counter
               LC = LC + 1
! Shift loop indicies to the correct values
               II = IW + I-1
               JJ = JS + J-1
               KK = KB + K-1
! Adjust for boundary conditions.
               IF (DES_PERIODIC_WALLS_X) THEN
                  IF(II .LT. IMIN1)THEN
! Shift from the east side of the domain to the west side.
                     II = IMAX2 + (II-IMIN1)
                     STNCL_GEMTY(I,J,K,1) = XE(II) - HALF*DX(II) - XLENGTH
                  ELSEIF(II .GT. IMAX1) THEN
! Shift for the west side of the domain to the east side.
                     II = II - IMAX
                     STNCL_GEMTY(I,J,K,1) = XLENGTH + XE(II) - HALF*DX(II)
                  ELSE
! No shift needed. Use actual node location.
                     STNCL_GEMTY(I,J,K,1) = XE(II) - HALF*DX(II)
                  ENDIF
               ELSE
! If the YZ plane is not periodic, calculate the position of the node.
                  STNCL_GEMTY(I,J,K,1) = XE(II) - HALF*DX(II)
! Due to the restrictions on non periodic 1 <= II and II <= IMAX2.
! Shift the II east/west to reflect the cell value at the boundary.
                  IF(II .EQ. 1)THEN
                     II = IMIN1
                  ELSEIF(II .EQ. IMAX2)THEN
                     II = IMAX1
                  ENDIF
               ENDIF
               IF (DES_PERIODIC_WALLS_Y) THEN
                  IF(JJ .LT. JMIN1) THEN
! Shift from the north side of the domain to the south side.
                     JJ = JMAX2 + (JJ-JMIN1)
                     STNCL_GEMTY(I,J,K,2) = YN(JJ) - HALF*DY(JJ) - YLENGTH
                  ELSEIF(JJ .GT. JMAX1) THEN
! Shift from the south side of the domain to the north side.
                     JJ = JJ - JMAX
                     STNCL_GEMTY(I,J,K,2) = YLENGTH + YN(JJ) - HALF*DY(JJ)
                  ELSE
! No shift needed. Use the nodes actual position.
                     STNCL_GEMTY(I,J,K,2) = YN(JJ) - HALF*DY(JJ)
                  ENDIF
               ELSE
! If the XZ plane is not periodic, calculate the position of the node.
                  STNCL_GEMTY(I,J,K,2) = YN(JJ) - HALF*DY(JJ)
! Due to the restrictions on non periodic 1 <= II and II <= IMAX2.
! Shift the II east/west to reflect the cell value at the boundary.
                  IF(JJ .EQ. 1)THEN
                     JJ = JMIN1
                  ELSEIF(JJ .EQ. JMAX2)THEN
                     JJ = JMAX1
                  ENDIF
               ENDIF
               IF (DIMN .EQ. 3) THEN
                  IF (DES_PERIODIC_WALLS_Z) THEN
                     IF(KK .LT. KMIN1) THEN
! Shift from the top side of the domain to the bottom side.
                        KK = KMAX1 + (KK-KMIN1)
                        STNCL_GEMTY(I,J,K,3) = ZT(KK) - HALF*DZ(KK) - ZLENGTH
                     ELSEIF(KK .GT. KMAX1) THEN
! Shift from the bottom side of the domain to the top side.
                        KK = KK - KMAX
                        STNCL_GEMTY(I,J,K,3) = ZLENGTH + ZT(KK) - HALF*DZ(KK)
                     ELSE
! No shift needed. Use actual node location.
                        STNCL_GEMTY(I,J,K,3) = ZT(KK) - HALF*DZ(KK)
                     ENDIF
                  ELSE
! If the XY plane is not periodic, calculate the position of the node.
                     STNCL_GEMTY(I,J,K,3) = ZT(KK) - HALF*DZ(KK)
! Due to the restrictions on non periodic 1 <= II and II <= IMAX2.
! Shift the II east/west to reflect the cell value at the boundary.
                     IF(KK .EQ. 1)THEN
                        KK = KMIN1
                     ELSEIF(KK .EQ. KMAX2)THEN
                        KK = KMAX1
                     ENDIF
                  ENDIF
               ELSE
                  KK = 1
               ENDIF
! Store the IJK in the interp array
               INTP_IJK(LC) = FUNIJK(II,JJ,KK)
! Assign a dummy value for the feild variable
               STNCL_FEILD(I,J,K) = 1.0d0
            ENDDO
         ENDDO
      ENDDO

! Call the interpolation routine. If the call to this routine originates
! from the continuum phase, the interpolation weights are of interest.
! Otherwise, the interpolated scalar value (INTERP_scalar) is desired.
      IF(DIMN.EQ.2) THEN
         CALL interpolator( STNCL_GEMTY(1:2, 1:2, 1, 1:DIMN), &
            STNCL_FEILD(1:2, 1:2, 1), DES_POS_NEW(NP,1:2), INTERP_scalar,&
            2, INTP_SCHM, WEIGHTS_CC )
      ELSE
         CALL interpolator( STNCL_GEMTY(1:2, 1:2, 1:2, 1:DIMN), &
            STNCL_FEILD(1:2, 1:2, 1:2), DES_POS_NEW(NP,:), INTERP_scalar,&
            2, INTP_SCHM, WEIGHTS_CC )
      ENDIF
! Store interpolation weights in an array for calling routine.
      LC = 0
      DO K = 1, (3-DIMN)*1+(DIMN-2)*2
         DO J = 1, 2
            DO I = 1, 2
! increment the loop counter
               LC = LC + 1
               INTP_WEIGHTS(LC) = WEIGHTS_CC(I,J,K)
            ENDDO
         ENDDO
      ENDDO

! Local debugging
      IF(FOCUS .AND. DEBUG_DES)THEN
         WRITE(*,"(3X,A,1X,A3,4X,A,I2)") &
            'Interpolation Scheme:','LPI','Interpolation Order:',2
         WRITE(*,"(3X,A,I5,3X,A,I6)")'Particle: ',NP,'IJK: ',PIJK(NP,4)
         WRITE(*,"(/5X,'|',29('-'),'|')")
         WRITE(*,"(5X,A)")'| Fluid Cell | Interp. Weight |'
         WRITE(*,"(5X,'|',12('-'),'|',16('-'),'|')")
         DO LC = 1, 2**DIMN
            WRITE(*,"(5X,'|',2X,I8,2X,'|',2X,F13.10,2X,'|')")&
               INTP_IJK(LC), INTP_WEIGHTS(LC)
         ENDDO
         WRITE(*,"(5X,'|',29('-'),'|'//)")
      ENDIF ! Local Debugging

      END SUBROUTINE INTERPOLATE_CC




  SUBROUTINE set_interpolation_pstencil(pc, ib,ie,jb,je,kb,ke, isch&
       &,dimprob, ordernew)
       !USE discretelement, ONLY : order,ob2l,ob2r,  intx_per, inty_per, intz_per
    USE geometry
    USE param1, only: half

    IMPLICIT NONE
    INTEGER, DIMENSION(3), INTENT(in):: pc
    INTEGER :: im,jm,km
    INTEGER, INTENT(out):: ib, ie, jb, je, kb, ke
    INTEGER, INTENT(in) :: dimprob
    INTEGER, OPTIONAL :: ordernew
    CHARACTER(LEN=5), INTENT(in) :: isch
    INTEGER :: ob2rtmp, ob2ltmp, ordertemp, ordernewtmp
    ob2l = (order+1)/2
    ob2r = order/2
    im = imax1
    jm = jmax1
    km = kmax1

    SELECT CASE(isch)
    CASE('csi')

       ob2rtmp = ob2r
       ob2ltmp = ob2l
       ordertemp = order
!!$       IF(order.NE.3.AND.(pc(1).EQ.1.OR.pc(1).EQ.im&
!!$            &-1.OR.pc(2).EQ.1.OR.pc(2).EQ.jm&
!!$            &-1.OR.pc(3).EQ.1.OR.pc(3).EQ.km-1)) order = 3
       !To make the order at boundary cells for csi equal to 3.
       !print*,'ob2l = ',ob2l,ob2r
       ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
       ie = MIN(im,pc(1) + ob2r)
          IF (ib.EQ.1 ) ie = ib + order - 1
          IF (ie.EQ.im) ib = ie - order + 1


       jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
       je = MIN(jm,pc(2) + ob2r)

          IF (jb.EQ.1 ) je = jb + order - 1
          IF (je.EQ.jm) jb = je - order + 1


       kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
       ke = MIN(km,pc(3) + ob2r)

          IF (kb.EQ.1 ) ke = kb + order - 1
          IF (ke.EQ.km) kb = ke - order + 1


       ob2r =  ob2rtmp
       ob2l = ob2ltmp
       ordernewtmp = order

       !Print*,'ib,ie .... in processing =', ib,ie,jb,je,kb,ke,&
       !     & ordernewtmp
       !PRINT*, 'pc = ',pc(1),pc(2),pc(3)!
       order = ordertemp !reset the order

    CASE('lpi')
       !print*, 'order in set stencil = ', order,pc(3)
       !Print*,'ib, ie = ', pc(1), ib, ie


       ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
       ie = MIN(im,pc(1) + ob2r)

          IF (ib.EQ.1 ) ie = ib + order - 1
          IF (ie.EQ.im) ib = ie - order + 1


       jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
       je = MIN(jm,pc(2) + ob2r)

          IF (jb.EQ.1 ) je = jb + order - 1
          IF (je.EQ.jm) jb = je - order + 1


       kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
       ke = MIN(km,pc(3) + ob2r)

          IF (kb.EQ.1 ) ke = kb + order - 1
          IF (ke.EQ.km) kb = ke - order + 1

       ordernewtmp = order
    END SELECT
    IF(dimprob == 2) THEN
       kb = pc(3)
       ke = pc(3)
       !Print*,'kb,ke  =', kb,ke
    ENDIF
    IF(PRESENT(ordernew)) ordernew = ordernewtmp
  END SUBROUTINE set_interpolation_pstencil



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE interp_oned_scalar(coor,scalar,ppos,interp_scl,order, &
      isch, weight_pointer)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:), INTENT(in):: coor, scalar
    REAL(prcn), INTENT(in):: ppos
    REAL(prcn), INTENT(out):: interp_scl
    INTEGER, INTENT(in):: order
    CHARACTER(LEN=5), INTENT(in) :: isch
    REAL(prcn), DIMENSION(:), POINTER, OPTIONAL:: weight_pointer

    REAL(prcn), DIMENSION(:), ALLOCATABLE:: zetacsi
    INTEGER ::  iorig, i
    REAL(prcn):: zeta
!-----------------------------------------------

    DO i = 1,order-1
       dx(i) = coor(i+1) - coor(i)
    ENDDO


    SELECT CASE(isch)
    CASE('lpi')

       !order = SIZE(coor)
       iorig = order/2

       zeta = (ppos - coor(iorig))/dx(iorig)

       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (2)
          DO i = 1,order
             xval(i) = shape2(zeta,i)
          ENDDO
       CASE (3)
          DO i = 1,order
             xval(i) = shape3(zeta,i,dx)
          ENDDO
       CASE (4)
          DO i = 1,order
             xval(i) = shape4(zeta,i,dx)
          ENDDO
       CASE (5)
          DO i = 1,order
             xval(i) = shape5(zeta,i,dx)
          ENDDO
       CASE (6)
          DO i = 1,order
             xval(i) = shape6(zeta,i,dx)
          ENDDO
       END SELECT
    CASE('csi')
       ! order = SIZE(coor,1)
       iorig = (order+1)/2
       !-------
       ! Find out center cell widths
       !-------
       ALLOCATE(zetacsi(order))
       !Zetacsi as defined in Yueng and Pope hence the name
       !The defintions for zetacsi are true only for a structured grid
       !if (order.eq.4) then

       !end if
       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (4)
          DO i = 1, order
             zetacsi(i) = (-ppos + coor(i))/dx(1)
          END DO
          DO i = 1,order
             xval(i) = shape4csi(zetacsi(i),i,dx,1)
          ENDDO

       end SELECT

       DEALLOCATE(zetacsi)
    end SELECT !Scheme

    !-------
    ! Calculate interpolated value
    !-------
    interp_scl = 0.0
    DO i = 1,order
       interp_scl = interp_scl + scalar(i)*xval(i)
    ENDDO

    !-------
    ! Return the weights (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => xval
    ENDIF

  END SUBROUTINE interp_oned_scalar



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE interp_oned_vector(coor,vector,ppos,interp_vec,order,&
      fun, weight_pointer)

! Interpolate an arbitrary sized array in one space dimension.
! Uses the scalar interpolation to obtain the weights.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:), INTENT(in):: coor
    REAL(prcn), DIMENSION(:,:), INTENT(in):: vector
    REAL(prcn), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:), INTENT(out):: interp_vec
    INTEGER, INTENT(in) :: order
    CHARACTER(LEN=5) :: fun
    REAL(prcn), DIMENSION(:), POINTER, OPTIONAL:: weight_pointer

    INTEGER:: vec_size, nv, i
    REAL(prcn), DIMENSION(:), POINTER:: weights_scalar
!-----------------------------------------------

    !order    = SIZE(coor)
    !print*,'In Interp_onedvector:order = ',order
    vec_size = SIZE(vector,2)

    !-------
    ! Interpolate first component and get weights
    !-------
    CALL interp_oned_scalar(coor,vector(:,1),ppos             &
         ,interp_vec(1),order, fun, weights_scalar)

    !-------
    ! Interpolate remaining components
    !-------
    DO nv = 2,vec_size
       interp_vec(nv) = 0.0
       DO i = 1,order
          interp_vec(nv) =  interp_vec(nv)  &
               + vector(i,nv)*weights_scalar(i)
       ENDDO
    ENDDO

    !-------
    ! Return the weights (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights_scalar
    ENDIF

  END SUBROUTINE interp_oned_vector



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE interp_twod_scalar(coor,scalar,ppos,interp_scl,order,&
      isch,weight_pointer)

! Interpolate a scalar quantity in two dimensions.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:), INTENT(in):: coor
    REAL(prcn), DIMENSION(:,:), INTENT(in):: scalar
    REAL(prcn), DIMENSION(2), INTENT(in):: ppos
    REAL(prcn), INTENT(out):: interp_scl
    INTEGER, INTENT(in):: order
    CHARACTER(LEN=5), INTENT(in) :: isch
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer

    REAL(prcn), DIMENSION(:,:), ALLOCATABLE:: zetacsi
    INTEGER:: i, j
    INTEGER:: iorig
    REAL(prcn), DIMENSION(2):: zeta
!-----------------------------------------------

    !-------
    ! Get order of interpolation
    !-------
    !
    weights  = zero
    DO i = 1,order-1
       dx(i) = coor(i+1,1,1)-coor(i,1,1)
       dy(i) = coor(1,i+1,2)-coor(1,i,2)
       !dz(i) = coor(1,1,i+1,3)-coor(1,1,i,3)
    ENDDO

    SELECT CASE(isch)

    CASE('lpi')
       !order = SIZE(coor,1)
       iorig = order/2

       !-------
       ! Find out center cell widths
       !-------

       !-------
       ! Zeta as defined in Bala/Maxey
       !-------
       zeta(1:2) = ppos(1:2) - coor(iorig,iorig,1:2)
       zeta(1) = zeta(1)/dx(iorig)
       zeta(2) = zeta(2)/dy(iorig)

       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (2)
          DO i = 1,order
             xval(i) = shape2(zeta(1),i)
             yval(i) = shape2(zeta(2),i)
             !zval(i) = shape2(zeta(3),i)
          ENDDO
       CASE (3)
          DO i = 1,order
             xval(i) = shape3(zeta(1),i,dx)
             yval(i) = shape3(zeta(2),i,dy)
             !zval(i) = shape3(zeta(3),i,dz)
          ENDDO
       CASE (4)
          DO i = 1,order
             ! print*, 'in interp....zetayp(3,i) = ', zetayp(3,i),zval(i),i
             xval(i) = shape4(zeta(1),i,dx)
             yval(i) = shape4(zeta(2),i,dy)
             !zval(i) = shape4(zeta(3),i,dz)
             !Print*, ppos(1:3)
             !Print*,'int',i,xval(i),yval(i),zval(i)
!!$          xval(i) = shape4new(ppos(1),coor(1:order,1,1,1),i)
!!$          yval(i) = shape4new(ppos(2),coor(1,1:order,1,2),i)
!!$          zval(i) = shape4new(ppos(3),coor(1,1,1:order,3),i)
          ENDDO
       CASE (5)
          DO i = 1,order
             xval(i) = shape5(zeta(1),i,dx)
             yval(i) = shape5(zeta(2),i,dy)
             !zval(i) = shape5(zeta(3),i,dz)
          ENDDO
       CASE (6)
          DO i = 1,order
             xval(i) = shape6(zeta(1),i,dx)
             yval(i) = shape6(zeta(2),i,dy)
             !zval(i) = shape6(zeta(3),i,dz)
          ENDDO
       END SELECT

    CASE('csi')
       ! order = SIZE(coor,1)
       iorig = (order+1)/2

       !-------
       ! Find out center cell widths
       !-------
       ALLOCATE(zetacsi(2,order))
       !Zetacsi as defined in Yueng and Pope hence the name
       !The defintions for zetacsi are true only for a structured grid

       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (4)
          DO i = 1, order
             zetacsi(1,i) = (-ppos(1) + coor(i,1,1))/dx(1)
             zetacsi(2,i) = (-ppos(2) + coor(1,i,2))/dy(1)
             !zetacsi(3,i) = (-ppos(3) + coor(1,1,i,3))/dz(1)
          END DO
          DO i = 1,order
             xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
             yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
             !zval(i) = shape4csi(zetacsi(3,i),i,dz,3)

             !Print*,'zetacsi  = ', zetacsi(3,i), coor(1,1,i,3), zval(i)
          ENDDO
       CASE(3)
          DO i = 1, order

             zetacsi(1,i) = ((-ppos(1) + coor(i,1,1))/dx(1))
             zetacsi(2,i) =((-ppos(2) + coor(1,i,2))/dy(1))
             !zetacsi(3,i) = ((-ppos(3) + coor(1,1,i,3))/dz(1))
!!$             zetacsi(1,i) = (ppos(1) - coor(1,1,1,1))/(coor(order,1,1&
!!$                  &,1)-coor(1,1,1,1))
!!$             zetacsi(2,i) = (ppos(2) - coor(1,1,1,2))/(coor(1,order,1&
!!$                  &,2)-coor(1,1,1,2))
!!$             zetacsi(3,i) = (ppos(3) - coor(1,1,1,3))/(coor(1,1,order&
!!$                  &,3)-coor(1,1,1,3))
          END DO
          DO i = 1,order
             if((xval(1)-coor(1,1,1)).lt.dx(1)) then
                xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
             else
                xval(i) = shape3csiright(zetacsi(1,i),i,dx,1)
             endif
             if((yval(1)-coor(1,1,2)).lt.dy(1)) then

                yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
             else

                yval(i) = shape3csiright(zetacsi(2,i),i,dy,2)
             end if

             print*,'zeta = ',zetacsi(1,i), xval(i),i
          ENDDO
       END SELECT

       DEALLOCATE(zetacsi)
    END SELECT !SCHEME
    !-------
    ! Calculate weights for the different nodes
    !-------
    ! DO 10 k = 1,order
!!$    DO 10 j = 1,order
!!$    DO 10 i = 1,order
!!$
!!$    10 CONTINUE
    !If(order.eq.3) Print*,'in interpo...sum wt=,',sum(weights),order

    !-------
    ! Calculate the interpolated value
    !-------
    interp_scl = 0.0
    !DO 20 k = 1,order
    DO  j = 1,order
       DO  i = 1,order
          weights(i,j,1) = xval(i)*yval(j)
          interp_scl = interp_scl + scalar(i,j)*weights(i,j,1)

       end DO
    end DO

    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights
    ENDIF

  END SUBROUTINE interp_twod_scalar



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE interp_twod_vector(coor,vector,ppos,interp_vec,order,&
      fun,weight_pointer)

! Interpolate an arbitrary sized array in two dimensions.
! Uses the scalar interpolation to obtain the weights.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:), INTENT(in):: coor, vector
    REAL(prcn), DIMENSION(2), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:), INTENT(out):: interp_vec
    INTEGER, INTENT(in):: order
    CHARACTER(LEN=5) :: fun
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer

    INTEGER :: vec_size
    INTEGER:: i, j, nv
    REAL(prcn), DIMENSION(:,:,:), POINTER:: weights_scalar
!-----------------------------------------------

    !-------
    ! Get order of interpolation
    !-------
    !order    = SIZE(coor,1)
    vec_size = SIZE(vector,3)
    !print*,'In Interp_threedvector:order = ',order,vec_size

    CALL interp_twod_scalar(coor,vector(:,:,1),ppos  &
         ,interp_vec(1),order,fun,weights_scalar)

    !-------
    ! Calculate the interpolated velocities
    !-------
    DO nv = 2,vec_size
       interp_vec(nv) = 0.0
       !DO  k = 1,order
       DO j = 1,order
          DO  i = 1,order
             interp_vec(nv) =  interp_vec(nv)  &
                  + vector(i,j,nv)*weights_scalar(i,j,1)
          ENDDO
       ENDDO
    ENDDO


    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights_scalar
    ENDIF

  END SUBROUTINE interp_twod_vector



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE interp_threed_scalar(coor,scalar,ppos,interp_scl,order,&
      isch,weight_pointer)

! Interpolate a scalar quantity in three dimensions.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor
    REAL(prcn), DIMENSION(:,:,:), INTENT(in):: scalar
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), INTENT(out):: interp_scl
    INTEGER, INTENT(in):: order
    CHARACTER(LEN=5), INTENT(in) :: isch
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer

    REAL(prcn), DIMENSION(:,:), ALLOCATABLE:: zetacsi
    INTEGER:: i, j, k
    INTEGER:: iorig
    REAL(prcn) :: zeta(3), zetasph, sigma, bandwidth
    LOGICAL :: calcwts
!-----------------------------------------------

    !-------
    ! Get order of interpolation
    !-------
    !
    weights  = zero
    DO i = 1,order-1
       dx(i) = coor(i+1,1,1,1)-coor(i,1,1,1)
       dy(i) = coor(1,i+1,1,2)-coor(1,i,1,2)
       dz(i) = coor(1,1,i+1,3)-coor(1,1,i,3)
    ENDDO

    !Print*,'In interpolator', isch
    SELECT CASE(isch)

    CASE('lpi')
       calcwts = .true.
       !order = SIZE(coor,1)
       iorig = order/2

       !-------
       ! Find out center cell widths
       !-------

       !-------
       ! Zeta as defined in Bala/Maxey
       !-------
       zeta(1:3) = ppos(1:3) - coor(iorig,iorig,iorig,1:3)
       zeta(1) = zeta(1)/dx(iorig)
       zeta(2) = zeta(2)/dy(iorig)
       zeta(3) = zeta(3)/dz(iorig)

       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (2)
          DO i = 1,order
             xval(i) = shape2(zeta(1),i)
             yval(i) = shape2(zeta(2),i)
             zval(i) = shape2(zeta(3),i)
          ENDDO
       CASE (3)
          DO i = 1,order
             xval(i) = shape3(zeta(1),i,dx)
             yval(i) = shape3(zeta(2),i,dy)
             zval(i) = shape3(zeta(3),i,dz)
          ENDDO
       CASE (4)
          DO i = 1,order
             ! print*, 'in interp....zetayp(3,i) = ', zetayp(3,i),zval(i),i
             xval(i) = shape4(zeta(1),i,dx)
             yval(i) = shape4(zeta(2),i,dy)
             zval(i) = shape4(zeta(3),i,dz)
             !Print*, ppos(1:3)
             !Print*,'int',i,xval(i),yval(i),zval(i)
!!$          xval(i) = shape4new(ppos(1),coor(1:order,1,1,1),i)
!!$          yval(i) = shape4new(ppos(2),coor(1,1:order,1,2),i)
!!$          zval(i) = shape4new(ppos(3),coor(1,1,1:order,3),i)
          ENDDO
       CASE (5)
          DO i = 1,order
             xval(i) = shape5(zeta(1),i,dx)
             yval(i) = shape5(zeta(2),i,dy)
             zval(i) = shape5(zeta(3),i,dz)
          ENDDO
       CASE (6)
          DO i = 1,order
             xval(i) = shape6(zeta(1),i,dx)
             yval(i) = shape6(zeta(2),i,dy)
             zval(i) = shape6(zeta(3),i,dz)
          ENDDO
       END SELECT

    CASE('csi')

       calcwts = .true.
       ! order = SIZE(coor,1)
       iorig = (order+1)/2

       !-------
       ! Find out center cell widths
       !-------
       ALLOCATE(zetacsi(3,order))
       !Zetacsi as defined in Yueng and Pope hence the name
       !The defintions for zetacsi are true only for a structured grid

       !-------
       ! Get shape function values
       !-------
       SELECT CASE (order)
       CASE (4)
          DO i = 1, order
             zetacsi(1,i) = (-ppos(1) + coor(i,1,1,1))/dx(1)
             zetacsi(2,i) = (-ppos(2) + coor(1,i,1,2))/dy(1)
             zetacsi(3,i) = (-ppos(3) + coor(1,1,i,3))/dz(1)
          END DO
          DO i = 1,order
             xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
             yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
             zval(i) = shape4csi(zetacsi(3,i),i,dz,3)

             !Print*,'zetacsi  = ', zetacsi(3,i), coor(1,1,i,3), zval(i)
          ENDDO
       CASE(3)
          DO i = 1, order

             zetacsi(1,i) = ((-ppos(1) + coor(i,1,1,1))/dx(1))
             zetacsi(2,i) =((-ppos(2) + coor(1,i,1,2))/dy(1))
             zetacsi(3,i) = ((-ppos(3) + coor(1,1,i,3))/dz(1))
!!$             zetacsi(1,i) = (ppos(1) - coor(1,1,1,1))/(coor(order,1,1&
!!$                  &,1)-coor(1,1,1,1))
!!$             zetacsi(2,i) = (ppos(2) - coor(1,1,1,2))/(coor(1,order,1&
!!$                  &,2)-coor(1,1,1,2))
!!$             zetacsi(3,i) = (ppos(3) - coor(1,1,1,3))/(coor(1,1,order&
!!$                  &,3)-coor(1,1,1,3))
          END DO
          DO i = 1,order
             if((xval(1)-coor(1,1,1,1)).lt.dx(1)) then
                xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
             else
                xval(i) = shape3csiright(zetacsi(1,i),i,dx,1)
             endif
             if((yval(1)-coor(1,1,1,2)).lt.dy(1)) then

                yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
             else

                yval(i) = shape3csiright(zetacsi(2,i),i,dy,2)
             end if
             if((zval(1)-coor(1,1,1,3)).lt.dz(1)) then
                zval(i) = shape3csileft(zetacsi(3,i),i,dz,3)
             else

                zval(i) = shape3csiright(zetacsi(3,i),i,dz,3)
             endif

             print*,'zeta = ',zetacsi(1,i), xval(i),i
          ENDDO
       END SELECT

       DEALLOCATE(zetacsi)

    CASE('sph')

       iorig = (order+1)/2
       calcwts = .false.
       !-------
       SELECT CASE (order)
       CASE (4)
          bandwidth  = one*(dx(1)*dy(1))**(half)
          sigma = one/pi
          do k = 1, order
             do j = 1, order
                DO i = 1, order
                   zetasph = (-ppos(1) + coor(i,j,k,1))**2.0
                   zetasph = zetasph + (-ppos(2) + coor(i,j,k,2))**2.0
                   !zetasph = zetasph + (-ppos(3) + coor(i,j,k,3))**2.0
                   zetasph = sqrt(zetasph)

                   zetasph = zetasph/bandwidth

                   if(zetasph.ge.zero.and.zetasph.lt.one) then
                      weights(i,j,k) = one - (three/two)*zetasph**two&
                           & + (three/four)*zetasph**three
                   elseif(zetasph.ge.one.and.zetasph.lt.two) then
                      weights(i,j,k) = fourth*(two-zetasph)**three
                   else
                      weights(i,j,k) = zero
                   end if
                   weights(i,j,k) = (sigma*weights(i,j,k))!/bandwidth&
                        !&**three

                   !Print*,weights(i,j,k), i,j,k
                END DO
             end do
          end DO
       END SELECT

    END SELECT !SCHEME
    !-------
    ! Calculate weights for the different nodes
    !-------
    if(calcwts) then
       DO  k = 1,order
          DO  j = 1,order
             DO  i = 1,order
                weights(i,j,k) = xval(i)*yval(j)*zval(k)
             end DO
          end DO
       end DO
    end if

    !If(order.eq.3) Print*,'in interpo...sum wt=,',sum(weights),order

    !-------
    ! Calculate the interpolated value
    !-------
    interp_scl = 0.0
    DO  k = 1,order
       DO  j = 1,order
          DO  i = 1,order
             interp_scl = interp_scl + scalar(i,j,k)*weights(i,j,k)
          end DO
       end DO
    end DO


    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights
    ENDIF
  END SUBROUTINE interp_threed_scalar




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE interp_threed_vector(coor,vector,ppos,interp_vec,order,&
      fun,weight_pointer)

! Interpolate an arbitrary sized array in three dimensions.
! Uses the scalar interpolation to obtain the weights.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor, vector
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:), INTENT(out):: interp_vec
    INTEGER, INTENT(in):: order
    CHARACTER(LEN=5) :: fun
    REAL(prcn), DIMENSION(:,:,:), POINTER, OPTIONAL:: weight_pointer

    INTEGER :: vec_size
    INTEGER:: i, j, k, nv
    REAL(prcn), DIMENSION(:,:,:), POINTER:: weights_scalar
!-----------------------------------------------

    !-------
    ! Get order of interpolation
    !-------
    !order    = SIZE(coor,1)
    vec_size = SIZE(vector,4)
    !print*,'In Interp_threedvector:order = ',order,vec_size

    CALL interp_threed_scalar(coor,vector(:,:,:,1),ppos  &
         ,interp_vec(1),order,fun,weights_scalar)

    !-------
    ! Calculate the interpolated velocities
    !-------
    DO nv = 2,vec_size
       interp_vec(nv) = 0.0
       DO  k = 1,order
          DO j = 1,order
             DO  i = 1,order
                interp_vec(nv) =  interp_vec(nv)  &
                     + vector(i,j,k,nv)*weights_scalar(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !-------
    ! Return the weights for the force distribution (optional)
    !-------
    IF (PRESENT(weight_pointer)) THEN
       weight_pointer => weights_scalar
    ENDIF

  END SUBROUTINE interp_threed_vector



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE calc_weightderiv_threed(coor,ppos,weight_pointer,order, isch)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in):: coor
    REAL(prcn), DIMENSION(3), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:,:,:,:) :: weight_pointer
    INTEGER, INTENT(in) :: order
    CHARACTER(len=5), INTENT(in) :: isch

    INTEGER:: i, j, k, kk
    Real(prcn) :: dx1,dy1,dz1
    INTEGER :: iorig
    REAL(prcn):: zeta(3), zetasph, bandwidth, sigma, tmp
    REAL(prcn), DIMENSION(3,order) :: zetacsi !right now only
    ! for cubic spline interp rg 03/14/05
!-----------------------------------------------

    weight_pointer = zero
    !-------
    ! Get order of interpolation
    !-------
    !order = SIZE(coor,1)
    iorig = order/2
    ! print*,'in interp...Iorig = ',iorig !Debug

    !-------
    ! Find out center cell widths
    !-------
    DO i = 1,order-1
       dx(i) = coor(i+1,1,1,1)-coor(i,1,1,1)
       dy(i) = coor(1,i+1,1,2)-coor(1,i,1,2)
       dz(i) = coor(1,1,i+1,3)-coor(1,1,i,3)
    ENDDO

    !Zetacsi as defined in Yueng and Pope hence the name
    !The defintions for zetacsi are true only for a structured grid
    !if (order.eq.4) then

    !end if
    !-------
    ! Zeta as defined in Bala/Maxey
    !-------


    !-------
    ! Get shape function values
    !-------
    SELECT CASE (isch)
    CASE ('csi')
       !Allocate(zetacsi(3,order))

       DO k = 1,3
          SELECT CASE(order)
          CASE(4)
             DO i = 1, order
                zetacsi(1,i) = (-ppos(1) + coor(i,1,1,1))/dx(1)
                zetacsi(2,i) = (-ppos(2) + coor(1,i,1,2))/dy(1)
                zetacsi(3,i) = (-ppos(3) + coor(1,1,i,3))/dz(1)
             END DO
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape4deriv_csi(zetacsi(1,i),i,dx))/dx(1)
                   yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
                   zval(i) = shape4csi(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
                   yval(i) = (shape4deriv_csi(zetacsi(2,i),i,dy))/(dy(1))
                   zval(i) = shape4csi(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.3) THEN
                   xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
                   yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
                   zval(i) = (shape4deriv_csi(zetacsi(3,i),i,dz))/(dz(1))
                ENDIF
             ENDDO
          CASE(3)
             dx1 = coor(order,1,1,1)-coor(1,1,1,1)
             dy1 = coor(1,order,1,2)-coor(1,1,1,2)
             dz1 = coor(1,1,order,3)-coor(1,1,1,3)
             zetacsi(1,i) = (ppos(1) - coor(1,1,1,1))/dx1
             zetacsi(2,i) = (ppos(2) - coor(1,1,1,2))/dy1
             zetacsi(3,i) = (ppos(3) - coor(1,1,1,3))/dz1
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape3deriv_csi(zetacsi(1,i),i,dx))/dx1
                   yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
                   zval(i) = shape3csileft(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
                   yval(i) = (shape3deriv_csi(zetacsi(2,i),i,dy))/(dy1)
                   zval(i) = shape3csileft(zetacsi(3,i),i,dz,3)
                ELSEIF(k.EQ.3) THEN
                   xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
                   yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
                   zval(i) = (shape3deriv_csi(zetacsi(3,i),i,dz))/(dz1)
                ENDIF
             ENDDO
          END SELECT!order
          DO kk = 1,order
             DO j = 1,order
                DO i = 1,order
                   weight_pointer(i,j,kk,k) = xval(i)*yval(j)*zval(kk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO!end loop over the coordinate directions
       !deallocate(zetacsi)
    CASE('lpi')
       zeta(1:3) = ppos(1:3) - coor(iorig,iorig,iorig,1:3)
       zeta(1) = zeta(1)/dx(iorig)
       zeta(2) = zeta(2)/dy(iorig)
       zeta(3) = zeta(3)/dz(iorig)
       DO k = 1,3
          SELECT CASE(order)
          CASE(4)
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape4deriv(zeta(1),i,dx))/dx(iorig)
                   yval(i) = shape4(zeta(2),i,dy)
                   zval(i) = shape4(zeta(3),i,dz)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape4(zeta(1),i,dx)
                   yval(i) = (shape4deriv(zeta(2),i,dy))/(dy(iorig))
                   zval(i) = shape4(zeta(3),i,dz)
                ELSEIF(k.EQ.3) THEN
                   xval(i) = shape4(zeta(1),i,dx)
                   yval(i) = shape4(zeta(2),i,dy)
                   zval(i) = (shape4deriv(zeta(3),i,dz))/(dz(iorig))
                ENDIF
             ENDDO
          CASE(2)
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape2deriv(zeta(1),i))/dx(iorig)
                   yval(i) = shape2(zeta(2),i)
                   zval(i) = shape2(zeta(3),i)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape2(zeta(1),i)
                   yval(i) = (shape2deriv(zeta(2),i))/(dy(iorig))
                   zval(i) = shape2(zeta(3),i)
                ELSEIF(k.EQ.3) THEN
                   xval(i) = shape2(zeta(1),i)
                   yval(i) = shape2(zeta(2),i)
                   zval(i) = (shape2deriv(zeta(3),i))/(dz(iorig))
                ENDIF
             ENDDO
          end SELECT

          DO kk = 1,order
             DO j = 1,order
                DO i = 1,order
                   weight_pointer(i,j,kk,k) = -xval(i)*yval(j)*zval(kk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO!end loop over the coordinate directions

    CASE('sph')
       SELECT CASE (order)

       CASE (4)
                    bandwidth  = one*(dx(1)*dy(1))**(half)
          sigma = one/pi

          do k = 1, order
             do j = 1, order
                DO i = 1, order
                   zetasph = (-ppos(1) + coor(i,j,k,1))**2.0
                   zetasph = zetasph + (-ppos(2) + coor(i,j,k,2))**2.0
                   !zetasph = zetasph + (-ppos(3) + coor(i,j,k,3))**2.0
                   zetasph = sqrt(zetasph)

                   zetasph = zetasph/bandwidth

                   if(zetasph.ge.zero.and.zetasph.lt.one) then
                      tmp = -two*(three/two)*zetasph &
                           & + three*(three/four)*zetasph**two
                   elseif(zetasph.ge.one.and.zetasph.lt.two) then
                      tmp = -three*fourth*(two-zetasph)**two
                   else
                      tmp = zero
                   end if

                   weight_pointer(i,j,k,1) = (tmp/zetasph)*(-ppos(1) &
                        &+ coor(i,j,k,1))
                   weight_pointer(i,j,k,2) = (tmp/zetasph)*(-ppos(2) &
                        &+ coor(i,j,k,2))
                   weight_pointer(i,j,k,3) = (tmp/zetasph)*(-ppos(3) &
                        &+ coor(i,j,k,3))
                   weight_pointer(i,j,k,:) = (sigma*weight_pointer(i,j&
                        &,k,:))/bandwidth&
                        &**two

                   !Print*,weights(i,j,k), i,j,k
                END DO
             end do
          end DO
       END SELECT

    END SELECT
  END SUBROUTINE calc_weightderiv_threed



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE calc_weightderiv_twod(coor,ppos,weight_pointer,order, isch)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:), INTENT(in):: coor
    REAL(prcn), DIMENSION(2), INTENT(in):: ppos
    REAL(prcn), DIMENSION(:,:,:,:) :: weight_pointer
    INTEGER, INTENT(in) :: order
    CHARACTER(len=5), INTENT(in) :: isch

    INTEGER:: i, j, k
    REAL(prcn) :: dx1,dy1
    INTEGER :: iorig
    REAL(prcn):: zeta(3), zetasph, bandwidth, sigma, tmp
    REAL(prcn), DIMENSION(2,order) :: zetacsi !right now only
    ! for cubic spline interp rg 03/14/05
!-----------------------------------------------

    weight_pointer = zero
    !-------
    ! Get order of interpolation
    !-------
    !order = SIZE(coor,1)
    iorig = order/2
    ! print*,'in interp...Iorig = ',iorig !Debug

    !-------
    ! Find out center cell widths
    !-------
    DO i = 1,order-1
       dx(i) = coor(i+1,1,1)-coor(i,1,1)
       dy(i) = coor(1,i+1,2)-coor(1,i,2)
    ENDDO

    !Zetacsi as defined in Yueng and Pope hence the name
    !The defintions for zetacsi are true only for a structured grid
    !if (order.eq.4) then

    !end if
    !-------
    ! Zeta as defined in Bala/Maxey
    !-------


    !-------
    ! Get shape function values
    !-------
    SELECT CASE (isch)
    CASE ('csi')
       !Allocate(zetacsi(3,order))

       DO k = 1,2
          SELECT CASE(order)
          CASE(4)
             DO i = 1, order
                zetacsi(1,i) = (-ppos(1) + coor(i,1,1))/dx(1)
                zetacsi(2,i) = (-ppos(2) + coor(1,i,2))/dy(1)

             END DO
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape4deriv_csi(zetacsi(1,i),i,dx))/dx(1)
                   yval(i) = shape4csi(zetacsi(2,i),i,dy,2)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape4csi(zetacsi(1,i),i,dx,1)
                   yval(i) = (shape4deriv_csi(zetacsi(2,i),i,dy))/(dy(1))
                ENDIF
             ENDDO
          CASE(3)
             dx1 = coor(order,1,1)-coor(1,1,1)
             dy1 = coor(1,order,2)-coor(1,1,2)
             zetacsi(1,i) = (ppos(1) - coor(1,1,1))/dx1
             zetacsi(2,i) = (ppos(2) - coor(1,1,2))/dy1
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape3deriv_csi(zetacsi(1,i),i,dx))/dx1
                   yval(i) = shape3csileft(zetacsi(2,i),i,dy,2)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape3csileft(zetacsi(1,i),i,dx,1)
                   yval(i) = (shape3deriv_csi(zetacsi(2,i),i,dy))/(dy1)
                ENDIF
             ENDDO
          END SELECT!order
             DO j = 1,order
                DO i = 1,order
                   weight_pointer(i,j,1,k) = xval(i)*yval(j)
                ENDDO
             ENDDO
          ENDDO!end loop over the coordinate directions
       !deallocate(zetacsi)
    CASE('lpi')
       zeta(1:2) = ppos(1:2) - coor(iorig,iorig,1:2)
       zeta(1) = zeta(1)/dx(iorig)
       zeta(2) = zeta(2)/dy(iorig)
       DO k = 1,2
          SELECT CASE(order)
          CASE(4)
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape4deriv(zeta(1),i,dx))/dx(iorig)
                   yval(i) = shape4(zeta(2),i,dy)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape4(zeta(1),i,dx)
                   yval(i) = (shape4deriv(zeta(2),i,dy))/(dy(iorig))
                ENDIF
             ENDDO
          CASE(2)
             DO i = 1,order
                IF(k.EQ.1) THEN
                   xval(i) = (shape2deriv(zeta(1),i))/dx(iorig)
                   yval(i) = shape2(zeta(2),i)
                ELSEIF(k.EQ.2) THEN
                   xval(i) = shape2(zeta(1),i)
                   yval(i) = (shape2deriv(zeta(2),i))/(dy(iorig))
                ENDIF
             ENDDO
          end SELECT

          DO j = 1,order
                DO i = 1,order
                   weight_pointer(i,j,1,k) = -xval(i)*yval(j)
                ENDDO
             ENDDO
          ENDDO!end loop over the coordinate directions

    CASE('sph')
       SELECT CASE (order)

       CASE (4)
          bandwidth  = one*(dx(1)*dy(1))**(half)
          sigma = one/pi

          do j = 1, order
             DO i = 1, order
                zetasph = (-ppos(1) + coor(i,j,1))**2.0
                zetasph = zetasph + (-ppos(2) + coor(i,j,2))**2.0
                !zetasph = zetasph + (-ppos(3) + coor(i,j,k,3))**2.0
                zetasph = sqrt(zetasph)

                zetasph = zetasph/bandwidth

                if(zetasph.ge.zero.and.zetasph.lt.one) then
                   tmp = -two*(three/two)*zetasph &
                        & + three*(three/four)*zetasph**two
                elseif(zetasph.ge.one.and.zetasph.lt.two) then
                   tmp = -three*fourth*(two-zetasph)**two
                else
                   tmp = zero
                end if

                weight_pointer(i,j,1,1) = (tmp/zetasph)*(-ppos(1) &
                        &+ coor(i,j,1))
                   weight_pointer(i,j,1,2) = (tmp/zetasph)*(-ppos(2) &
                        &+ coor(i,j,2))
                   !weight_pointer(i,j,1,3) = (tmp/zetasph)*(-ppos(3) &
                   !     &+ coor(i,j,3))
                   weight_pointer(i,j,1,:) = (sigma*weight_pointer(i,j&
                        &,1,:))/bandwidth&
                        &**two

                   !Print*,weights(i,j,k), i,j,k
                END DO
             end do
          END SELECT

       END SELECT

     END SUBROUTINE calc_weightderiv_twod



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION justweights(coor,ppos)
! To calculate just the weights using trilinear interpolation
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn), DIMENSION(:,:,:), POINTER:: justweights
    REAL(prcn), DIMENSION(:,:,:,:), INTENT(IN):: coor
    REAL(prcn), DIMENSION(3), INTENT(IN):: ppos

    INTEGER, PARAMETER:: order=2
    INTEGER:: i, j, k, iorig
    REAL(prcn):: dxl, dyl, dzl
    REAL(prcn), DIMENSION(3):: zeta
!-----------------------------------------------

    iorig = order/2

    dxl = coor(order,1,1,1)-coor(order-1,1,1,1)
    dyl = coor(1,order,1,2)-coor(1,order-1,1,2)
    dzl = coor(1,1,order,3)-coor(1,1,order-1,3)

    zeta(1:3) = ppos(1:3) - coor(iorig,iorig,iorig,1:3)
    zeta(1) = zeta(1)/dxl
    zeta(2) = zeta(2)/dyl
    zeta(3) = zeta(3)/dzl

    DO i = 1,order
       xval(i) = shape2(zeta(1),i)
       yval(i) = shape2(zeta(2),i)
       zval(i) = shape2(zeta(3),i)
    ENDDO

    !-------
    ! Calculate weights for the different nodes
    !-------
    DO  k = 1,order
       DO  j = 1,order
          DO  i = 1,order
             weights(i,j,k) = xval(i)*yval(j)*zval(k)
          end DO
       end DO
    end DO


    justweights => weights
  END FUNCTION justweights



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape2(zeta,i)
! Second-order (linear) shape functions
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape2
    REAL(prcn):: zeta
    INTEGER:: i
!-----------------------------------------------
    SELECT CASE (i)
    CASE (1)
       shape2 = 1 - zeta
    CASE (2)
       shape2 = zeta
    END SELECT
  END FUNCTION shape2



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape3(zeta,i,dx)
! Third-order (quadratic) shape functions
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape3
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn):: zh, num, denom
!-----------------------------------------------

    SELECT CASE (i)
    CASE (1)
       zh    = dx(1)*zeta
       denom = dx(1)*(dx(1)+dx(2))
       num   = (zh - dx(1))*(zh - dx(1) - dx(2))
       shape3 = num/denom
    CASE (2)
       zh    = dx(1)*zeta
       denom = -dx(1)*dx(2)
       num   = zh*(zh - dx(1) - dx(2))
       shape3 = num/denom
    CASE (3)
       zh    = dx(1)*zeta
       denom = dx(2)*(dx(1)+dx(2))
       num   = zh*(zh - dx(1))
       shape3 = num/denom
    END SELECT
  END FUNCTION shape3



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape4(zeta,i,dx)
! Fourth-order (cubic) shape functions
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape4
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn):: r1, r2, c1, c2, c3
!-----------------------------------------------

    r1 = dx(2)/dx(1)
    r2 = dx(3)/dx(1)

    SELECT CASE (i)
    CASE (1)
       c1 = 1.0/(1.0 + r1)
       c3 = 1.0/(1.0 + r1 + r1*r2)
       shape4 = -c1*c3*(r1**3.0)*(zeta)*(zeta-1.0)*(zeta-(1.0+r2))
       !shape4 = -(one/6.)*(zeta)*(zeta-one)*(zeta-two)
    CASE (2)
       c2 = 1.0/(1.0 + r2)
       shape4 = c2*(zeta-1.0)*(zeta*r1+1.0)*(zeta-(1.0+r2))
       ! shape4 = half*(zeta-one)*(zeta+one)*(zeta-two)
    CASE (3)
       c1 = 1.0/(1.0 + r1)
       shape4 = -(c1/r2)*(zeta)*(zeta*r1+1.0)*(zeta-(1.0+r2))
       ! shape4 = -half*(zeta)*(zeta+1)*(zeta-two)
    CASE (4)
       c2 = 1.0/(1.0 + r2)
       c3 = 1.0/(1.0 + r1 + r1*r2)
       shape4 = (c3*c2/r2)*(zeta)*(zeta*r1+1.0)*(zeta-1.0)
       ! shape4 = (one/6.)*(zeta)*(zeta+one)*(zeta-one)
    END SELECT
  END FUNCTION shape4



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape5(zeta,i,dx)
! Fifth-order (power of 4) shape functions
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape5
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn):: d1, d2, d3, d4
    REAL(prcn):: num, denom, zh
!-----------------------------------------------

    SELECT CASE (i)
    CASE (1)
       d1 = -dx(1)
       d2 = d1 - dx(2)
       d3 = d2 - dx(3)
       d4 = d3 - dx(4)
       denom = d1*d2*d3*d4
       zh = zeta*dx(2)
       num = zh*(zh -dx(2))*(zh -dx(2) -dx(3)) &
            *(zh -dx(2) -dx(3) -dx(4))
       shape5 = num/denom
    CASE (2)
       d1 =  dx(1)
       d2 = -dx(2)
       d3 =  d2 - dx(3)
       d4 =  d3 - dx(4)
       denom = d1*d2*d3*d4
       zh = zeta*dx(2)
       num = (zh +dx(1))*(zh -dx(2))*(zh -dx(2) -dx(3)) &
            *(zh -dx(2) -dx(3) -dx(4))
       shape5 = num/denom
    CASE (3)
       d1 =  dx(1) + dx(2)
       d2 =  dx(2)
       d3 = -dx(3)
       d4 =  d3 - dx(4)
       denom = d1*d2*d3*d4
       zh = zeta*dx(2)
       num = (zh +dx(1))*(zh)*(zh -dx(2) -dx(3)) &
            *(zh -dx(2) -dx(3) -dx(4))
       shape5 = num/denom
    CASE (4)
       d1 =  dx(1) + dx(2) + dx(3)
       d2 =  d1 - dx(1)
       d3 =  d2 - dx(2)
       d4 = -dx(4)
       denom = d1*d2*d3*d4
       zh = zeta*dx(2)
       num = (zh +dx(1))*(zh)*(zh -dx(2)) &
            *(zh -dx(2) -dx(3) -dx(4))
       shape5 = num/denom
    CASE (5)
       d1 =  dx(1) + dx(2) + dx(3) + dx(4)
       d2 =  d1 - dx(1)
       d3 =  d2 - dx(2)
       d4 =  d3 - dx(3)
       denom = d1*d2*d3*d4
       zh = zeta*dx(2)
       num = (zh +dx(1))*(zh)*(zh -dx(2)) &
            *(zh -dx(2) -dx(3))
       shape5 = num/denom
    END SELECT
  END FUNCTION shape5



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape6(zeta,i,dx)
! Sixth-order (power of 5) shape functions
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape6
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn):: d1, d2, d3, d4, d5
    REAL(prcn):: num, denom, zh
!-----------------------------------------------

    SELECT CASE (i)
    CASE (1)
       d1 = -dx(1)
       d2 = d1 - dx(2)
       d3 = d2 - dx(3)
       d4 = d3 - dx(4)
       d5 = d4 - dx(5)
       denom = d1*d2*d3*d4*d5
       zh = zeta*dx(3)
       num = (zh + dx(2))*(zh)*(zh - dx(3)) &
            *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
       shape6 = num/denom
    CASE (2)
       d1 =  dx(1)
       d2 = -dx(2)
       d3 =  d2 - dx(3)
       d4 =  d3 - dx(4)
       d5 =  d4 - dx(5)
       denom = d1*d2*d3*d4*d5
       zh = zeta*dx(3)
       num = (zh +dx(1) +dx(2))*(zh)*(zh -dx(3)) &
            *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
       shape6 = num/denom
    CASE (3)
       d1 =  dx(1) + dx(2)
       d2 =  dx(2)
       d3 = -dx(3)
       d4 =  d3 - dx(4)
       d5 =  d4 - dx(5)
       denom = d1*d2*d3*d4*d5
       zh = zeta*dx(3)
       num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh -dx(3)) &
            *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
       shape6 = num/denom
    CASE (4)
       d1 =  dx(1) + dx(2) + dx(3)
       d2 =  d1 - dx(1)
       d3 =  d2 - dx(2)
       d4 = -dx(4)
       d5 =  d4 - dx(5)
       denom = d1*d2*d3*d4*d5
       zh = zeta*dx(3)
       num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh) &
            *(zh -dx(3) -dx(4))*(zh - dx(3) -dx(4) -dx(5))
       shape6 = num/denom
    CASE (5)
       d1 =  dx(1) + dx(2) + dx(3) + dx(4)
       d2 =  d1 - dx(1)
       d3 =  d2 - dx(2)
       d4 =  d3 - dx(3)
       d5 = -dx(5)
       denom = d1*d2*d3*d4*d5
       zh = zeta*dx(3)
       num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh) &
            *(zh -dx(3))*(zh - dx(3) -dx(4) -dx(5))
       shape6 = num/denom
    CASE (6)
       d1 =  dx(1) + dx(2) + dx(3) + dx(4) + dx(5)
       d2 =  d1 - dx(1)
       d3 =  d2 - dx(2)
       d4 =  d3 - dx(3)
       d5 =  d4 - dx(4)
       denom = d1*d2*d3*d4*d5
       zh = zeta*dx(3)
       num = (zh +dx(1) +dx(2))*(zh +dx(2))*(zh) &
            *(zh -dx(3))*(zh - dx(3) -dx(4))
       shape6 = num/denom
    END SELECT
  END FUNCTION shape6


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape3deriv_csi(zeta,i,dx)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape3deriv_csi
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
!-----------------------------------------------

!!$    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN
!!$       shape3deriv_csi = (half)*(two+zeta)**2.0
!!$    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN
!!$       shape3deriv_csi = (one/six)*(-9.0*zeta**2.0-12.0*zeta)
!!$    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN
!!$       shape3deriv_csi = (one/six)*(9.0*zeta**2.0-12.0*zeta)
!!$    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN
!!$       shape3deriv_csi = (-half)*(two-zeta)**2.0
!!$    ELSE
!!$       shape3deriv_csi = zero
!!$    ENDIF
    SELECT CASE (i)
    CASE (1)
       shape3deriv_csi = -two*(1-zeta)
    CASE (2)
       shape3deriv_csi = -two*zeta+two*(1-zeta)
    CASE (3)
       shape3deriv_csi = two*zeta
    END SELECT
  END FUNCTION shape3deriv_csi



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape4deriv_csi(zeta,i,dx)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape4deriv_csi
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
!-----------------------------------------------

    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN
       shape4deriv_csi = (half)*(two+zeta)**2.0
    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN
       shape4deriv_csi = (one/six)*(-9.0*zeta**2.0-12.0*zeta)
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN
       shape4deriv_csi = (one/six)*(9.0*zeta**2.0-12.0*zeta)
    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN
       shape4deriv_csi = (-half)*(two-zeta)**2.0
    ELSE
       shape4deriv_csi = zero
    ENDIF
  END FUNCTION shape4deriv_csi


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape4deriv(zeta,i,dx)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape4deriv
    REAL(prcn):: zeta
    REAL(prcn), DIMENSION(:):: dx
    INTEGER:: i
    REAL(prcn) :: tmp
    REAL(prcn):: r1, r2, c1, c2, c3
!-----------------------------------------------

    r1 = dx(2)/dx(1)
    r2 = dx(3)/dx(1)

    SELECT CASE (i)
    CASE (1)
       c1 = 1.0/(1.0 + r1)
       c3 = 1.0/(1.0 + r1 + r1*r2)
       tmp = -c1*c3*(r1**3.0)
       !shape4deriv = -c1*c3*(r1**3.0)*(zeta)*(zeta-1.0)*(zeta-(1.0+r2))
       shape4deriv = (1.0)*(zeta-1.0)*(zeta-(1.0+r2))&
            + (zeta)*(1.0)*(zeta-(1.0+r2))&
            +(zeta)*(zeta-1.0)*(1.0)
       shape4deriv = tmp*shape4deriv
    CASE (2)
       c2 = 1.0/(1.0 + r2)

       shape4deriv = (1.0)*(zeta*r1+1.0)*(zeta-(1.0+r2)) &
            +(zeta-1.0)*(r1)*(zeta-(1.0+r2)) &
            +(zeta-1.0)*(zeta*r1+1.0)*(1.0)
       shape4deriv = c2*shape4deriv
    CASE (3)
       c1 = 1.0/(1.0 + r1)
       tmp = -c1/r2
       shape4deriv = (1.0)*(zeta*r1+1.0)*(zeta-(1.0+r2))&
            +(zeta)*(r1)*(zeta-(1.0+r2))&
            +(zeta)*(zeta*r1+1.0)*(1.0)
       shape4deriv = tmp*shape4deriv
    CASE (4)
       c2 = 1.0/(1.0 + r2)
       c3 = 1.0/(1.0 + r1 + r1*r2)
       tmp = (c3*c2/r2)
       shape4deriv = (1.0)*(zeta*r1+1.0)*(zeta-1.0)&
            +(zeta)*(r1)*(zeta-1.0)&
            +(zeta)*(zeta*r1+1.0)*(1.0)
       shape4deriv = tmp*shape4deriv
    END SELECT
  END FUNCTION shape4deriv



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape2deriv(zeta,i)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape2deriv
    REAL(prcn):: zeta
    INTEGER:: i
!-----------------------------------------------

    SELECT CASE (i)
    CASE (1)
       shape2deriv = -one!zeta
    CASE (2)
       shape2deriv = one
    END SELECT
  END FUNCTION shape2deriv



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape4csi(zeta,i,dx,dim)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape4csi
    REAL(prcn),INTENT(in):: zeta
    REAL(prcn), DIMENSION(:),INTENT(in):: dx
    INTEGER,INTENT(in):: i,dim
!-----------------------------------------------

    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN
       shape4csi = (one/six)*(two+zeta)**3.0
    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN
       shape4csi = (one/six)*(-3.0*zeta**3.0-six*zeta**2.0+4.0)
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN
       shape4csi = (one/six)*(3.0*zeta**3.0-six*zeta**2.0+4.0)
    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN
       shape4csi = (one/six)*(two-zeta)**3.0
    ELSE
       shape4csi = zero
       !print*,'shape4rg .... zeta=',zeta,dim
    ENDIF

  END FUNCTION shape4csi



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape3csileft(zeta,i,dx,dim)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape3csileft
    REAL(prcn),INTENT(in):: zeta
    REAL(prcn), DIMENSION(:),INTENT(in):: dx
    INTEGER,INTENT(in):: i,dim
!-----------------------------------------------

    IF (zeta.GE.-one.AND.zeta.LE.zero) THEN
       shape3csileft = (half)*(zeta+one)**2.0
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN
       shape3csileft = (half)*(-two*zeta**2.+2.0*zeta+1.0)
    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN
       shape3csileft = (half)*(zeta-two)**2.0
    ELSE
       shape3csileft = zero
       !print*,'shape4rg .... zeta=',zeta,dim
    ENDIF

!!$    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN
!!$       shape3csileft = (half)*(two+zeta)**2.0
!!$    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN
!!$       shape3csileft = (half)*(-two*zeta**2.-two*zeta+one)
!!$    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN
!!$       shape3csileft = (half)*(one-zeta)**2.0
!!$    ELSEIF(zeta.GT.one.AND.zeta.LE.two) THEN
!!$       shape3csileft = (half)*(two-zeta)**2.0
!!$    ELSE
!!$       shape3csileft = zero
!!$       !print*,'shape4rg .... zeta=',zeta,dim
!!$    ENDIF
!!$    SELECT CASE (i)
!!$    CASE (1)
!!$       shape3csileft = (1-zeta)**2.
!!$       !shape4 = -(one/6.)*(zeta)*(zeta-one)*(zeta-two)
!!$    CASE (2)
!!$        shape3csileft = two*zeta*(1-zeta)
!!$     CASE (3)
!!$        shape3csileft = zeta**2.
!!$     END SELECT
  END FUNCTION shape3csileft



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape3csiright(zeta,i,dx,dim)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape3csiright
    REAL(prcn),INTENT(in):: zeta
    REAL(prcn), DIMENSION(:),INTENT(in):: dx
    INTEGER,INTENT(in):: i,dim
!-----------------------------------------------

    IF (zeta.GE.-two.AND.zeta.LE.-one) THEN
       shape3csiright = (half)*(two+zeta)**2.0
    ELSEIF(zeta.GT.-one.AND.zeta.LE.zero) THEN
       shape3csiright = (half)*(-two*zeta**2.-two*zeta+one)
    ELSEIF(zeta.GT.zero.AND.zeta.LE.one) THEN
       shape3csiright = (half)*(one-zeta)**2.0
    ELSE
       shape3csiright = zero
       !print*,'shape4rg .... zeta=',zeta,dim
    ENDIF
  END FUNCTION shape3csiright



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  FUNCTION shape4new(pos,x,i)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    REAL(prcn):: shape4new
    REAL(prcn):: pos
    REAL(prcn), DIMENSION(:):: x
    INTEGER:: i
    REAL(prcn):: num,den
!-----------------------------------------------

    SELECT CASE (i)
    CASE (1)
       num = (pos-x(2))*(pos - x(3))*(pos - x(4))
       den = (x(1) - x(2))*(x(1) - x(3))*(x(1) - x(4))
       shape4new = num/den
    CASE (2)
       num = (pos-x(1))*(pos - x(3))*(pos - x(4))
       den = (x(2) - x(1))*(x(2) - x(3))*(x(2) - x(4))
       shape4new = num/den
    CASE (3)
       num = (pos-x(1))*(pos - x(2))*(pos - x(4))
       den = (x(3) - x(1))*(x(3) - x(2))*(x(3) - x(4))
       shape4new = num/den
    CASE (4)
       num = (pos-x(1))*(pos - x(2))*(pos - x(3))
       den = (x(4) - x(1))*(x(4) - x(2))*(x(4) - x(3))
       shape4new = num/den
    END SELECT
  END FUNCTION shape4new


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  SUBROUTINE set_interpolation_scheme(choice)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !USE discretelement, ONLY : scheme, interp_scheme, order,ob2l,ob2r,&
    !     & gstencil, vstencil, sstencil, wtderivp
      USE param
    IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    INTEGER, INTENT(in) :: choice
    INTEGER :: order_orig
!-----------------------------------------------

    order_orig = order
    IF(choice.EQ.1) THEN
       interp_scheme = 'lpi'
       scheme = '4-order'
    ELSE IF(choice.EQ.2) THEN
       interp_scheme = 'lpi'
       scheme = '2-order'
    ELSE IF(choice.EQ.3) THEN
       interp_scheme = 'csi'
       scheme = '4-order'
    ENDIF
    SELECT CASE(scheme)
    CASE("2-order")
       order = 2
    CASE("3-order")
       order = 3
    CASE("4-order")
       order = 4
    CASE("5-order")
       order = 5
    CASE("6-order")
       order = 6
    END SELECT

! if ob2l is even then ob2l will equal ob2r since results after decimal are
! discarded (truncated)
    ob2l = (order+1)/2
    ob2r = order/2

    IF(.not.allocated(gstencil)) THEN
! max(1*(3-dimn), order*(dimn-2)) =order (in 3D) or =1 (in 2D)
       ALLOCATE(gstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ELSEIF(ALLOCATED(gstencil).AND.order_orig.NE.order) THEN
       DEALLOCATE(gstencil)
       ALLOCATE(gstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ENDIF

    IF(.not.allocated(vstencil)) THEN
       ALLOCATE(vstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ELSEIF(ALLOCATED(vstencil).AND.order_orig.NE.order) THEN
       DEALLOCATE(vstencil)
       ALLOCATE(vstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ENDIF

    IF(.not.allocated(pgradstencil)) THEN
       ALLOCATE(pgradstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ELSEIF(ALLOCATED(pgradstencil).AND.order_orig.NE.order) THEN
       DEALLOCATE(pgradstencil)
       ALLOCATE(pgradstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ENDIF

    IF(.not.allocated(psgradstencil)) THEN
       ALLOCATE(psgradstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ELSEIF(ALLOCATED(psgradstencil).AND.order_orig.NE.order) THEN
       DEALLOCATE(psgradstencil)
       ALLOCATE(psgradstencil  (order,order,max(1*(3-dimn), order*(dimn-2)),3))
    ENDIF

    IF(.not.allocated(VEL_SOL_STENCIL)) THEN
       ALLOCATE(VEL_SOL_STENCIL  (order,order,max(1*(3-dimn), order*(dimn-2)),3, DIM_M))
    ELSEIF(ALLOCATED(VEL_SOL_STENCIL).AND.order_orig.NE.order) THEN
       DEALLOCATE(VEL_SOL_STENCIL)
       ALLOCATE(VEL_SOL_STENCIL  (order,order,max(1*(3-dimn), order*(dimn-2)),3, DIM_M))
    ENDIF

    IF(.not.allocated(sstencil)) THEN
       ALLOCATE(sstencil  (order,order,max(1*(3-dimn), order*(dimn-2))))
    ELSEIF(ALLOCATED(sstencil).AND.order_orig.NE.order) THEN
       DEALLOCATE(sstencil)
       ALLOCATE(sstencil  (order,order,max(1*(3-dimn), order*(dimn-2))))
    ENDIF

  END SUBROUTINE set_interpolation_scheme



END MODULE interpolation
