!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DEFINE_QUADRICS                                        C
!  Purpose: Defines all matrices used to evaluate quadrics             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE DEFINE_QUADRICS

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
      USE vtk

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3,3) :: Rx,Ry,Rz,C_QUADRIC,R_QUADRIC

!======================================================================
! Quadric Normal form : lambda_x * x^2 + lambda_y * y^2 + lambda_z * z^2 + d = 0
!======================================================================

      DO QUADRIC_ID = 1 , N_QUADRIC

!       Build translation matrices
        CALL BUILD_1x3_MATRIX(t_x(QUADRIC_ID),t_y(QUADRIC_ID),t_z(QUADRIC_ID),T_QUADRIC(:,:,QUADRIC_ID))
!       Build characteristic matrices
        CALL BUILD_C_QUADRIC_MATRIX(lambda_x(QUADRIC_ID),lambda_y(QUADRIC_ID),lambda_z(QUADRIC_ID),C_QUADRIC)
!       Build Rotation matrices
        CALL BUILD_X_ROTATION_MATRIX(Theta_x(QUADRIC_ID), Rx)
        CALL BUILD_Y_ROTATION_MATRIX(Theta_y(QUADRIC_ID), Ry)
        CALL BUILD_Z_ROTATION_MATRIX(Theta_z(QUADRIC_ID), Rz)
        R_QUADRIC = MATMUL(Rz,MATMUL(Ry,Rx))
!       Build A-matrices
        A_QUADRIC(:,:,QUADRIC_ID) = MATMUL(TRANSPOSE(R_QUADRIC),MATMUL(C_QUADRIC,R_QUADRIC))

      END DO

      ! Activate Quadric
      QUADRIC_ID = N_QUADRIC

      RETURN

      END SUBROUTINE DEFINE_QUADRICS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_F_QUADRIC                                          C
!  Purpose: Evaluate the function f(x,y,z) defining the quadric        C
!                                                                      C
!  Author: Jeff Dietiker                               Date: 21-FEB-08 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!  Modified: ##                                        Date: ##-###-## C
!  Purpose: ##                                                         C
!                                                                      C
!  Literature/Document References: ##                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
        SUBROUTINE GET_F_QUADRIC(x1,x2,x3,Q_ID,f,CLIP_FLAG)

      USE constant, only: pi
      USE exit, only: mfix_exit
      USE param1, only: half, one, undefined, zero
      USE quadric
      USE sendrecv, only: mype

      IMPLICIT NONE

      DOUBLE PRECISION :: x1,x2,x3,xt,yt,zt,R1,R2,xtr,ytr,ztr
      DOUBLE PRECISION :: THETA,THETA1,THETA2,THETA3m,THETA3
      DOUBLE PRECISION :: THETA1CYL1,THETA2TOR,THETA3CYL2,THETA_MIN
      DOUBLE PRECISION :: Y1,Y2,R
      DOUBLE PRECISION :: c,ss,ytest1,ytest2
      DOUBLE PRECISION :: f
      DOUBLE PRECISION :: fxmin,fxmax,fymin,fymax,fzmin,fzmax
      DOUBLE PRECISION, DIMENSION(1,3) :: X_VECTOR,XMT
      DOUBLE PRECISION, DIMENSION(3,1) :: TXMT
      DOUBLE PRECISION, DIMENSION(1,1) :: TEMP_1x1
      INTEGER :: Q_ID
      LOGICAL :: CLIP_X,CLIP_Y,CLIP_Z,CLIP_FLAG
      LOGICAL :: PIECE_X,PIECE_Y,PIECE_Z,PIECE_FLAG

      DOUBLE PRECISION :: YR1,YR2,RR1,RR2
      DOUBLE PRECISION :: YRR1,YRR2,RC1,RC2,YC1,YC2


!======================================================================

       PIECE_X = (piece_xmin(Q_ID) <= x1).AND.( x1 <= piece_xmax(Q_ID))
       PIECE_Y = (piece_ymin(Q_ID) <= x2).AND.( x2 <= piece_ymax(Q_ID))
       PIECE_Z = (piece_zmin(Q_ID) <= x3).AND.( x3 <= piece_zmax(Q_ID))

       PIECE_FLAG = (PIECE_X.AND.PIECE_Y.AND.PIECE_Z)

       IF(.NOT.PIECE_FLAG) THEN
          f = UNDEFINED
          RETURN
       ENDIF

       CLIP_X = (clip_xmin(Q_ID) <= x1).AND.( x1 <= clip_xmax(Q_ID))
       CLIP_Y = (clip_ymin(Q_ID) <= x2).AND.( x2 <= clip_ymax(Q_ID))
       CLIP_Z = (clip_zmin(Q_ID) <= x3).AND.( x3 <= clip_zmax(Q_ID))

       CLIP_FLAG = (CLIP_X.AND.CLIP_Y.AND.CLIP_Z)


         IF(TRIM(quadric_form(Q_ID))=='PLANE') THEN

            f = lambda_x(Q_ID)*x1 + lambda_y(Q_ID)*x2 +lambda_z(Q_ID)*x3 + dquadric(Q_ID)

         ELSEIF(TRIM(quadric_form(Q_ID))=='TORUS_INT') THEN

            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

            R1 = Torus_R1(Q_ID)
            R2 = Torus_R2(Q_ID)

            f = -(4*(xt**2+zt**2)*R1**2-(xt**2+yt**2+zt**2+R1**2-R2**2)**2)


         ELSEIF(TRIM(quadric_form(Q_ID))=='TORUS_EXT') THEN

            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

            xtr = xt
            ytr = yt
            ztr = zt

            R1 = Torus_R1(Q_ID)
            R2 = Torus_R2(Q_ID)

            f = 4*(xtr**2+ztr**2)*R1**2-(xtr**2+ytr**2+ztr**2+R1**2-R2**2)**2


         ELSEIF(TRIM(quadric_form(Q_ID))=='Y_UCOIL_EXT') THEN
! This shape represents a pair of parallel cylinders (y-direction)
! capped at both ends by half a torus
! to create a U-shaped coil
! UCOIL_Y1 and UCOIL_Y2 are the min and max y-values of the coil
! The coil is translated in the x and z direction (t_x, t_z)
! and can be rotated about the y-direction, centered about (t_x,t_z)
! by the angle THETA_Y

            c = DCOS(THETA_Y(Q_ID))
            ss = DSIN(THETA_Y(Q_ID))

!           Translation
            xt = x1-t_x(Q_ID)
            zt = x3-t_z(Q_ID)

!           Rotation
            xtr =  xt*c + zt*ss
            ztr = -xt*ss + zt*c

            R1 = UCOIL_R1(Q_ID)
            R2 = UCOIL_R2(Q_ID)

! Limits of the cylinders.
! There is half a torus above ytest1 and half a torus below ytest2
            ytest1 = UCOIL_Y1(Q_ID) + R1 + R2
            ytest2 = UCOIL_Y2(Q_ID) - (R1 + R2)

            IF(x2>=ytest2) THEN
               ytr = x2 - ytest2
               ELSEIF(x2<=ytest1) THEN
               ytr = ytest1 - x2
            ELSE
               ytr = ZERO ! setting ytr = zero degenerates a torus into a pair of cylinders
            ENDIF

            f = 4*(xtr**2+ytr**2)*R1**2-(xtr**2+ytr**2+ztr**2+R1**2-R2**2)**2

         ELSEIF(TRIM(quadric_form(Q_ID))=='Y_UCOIL2_EXT') THEN
! This shape represents a pair of parallel cylinders (y-direction)
! capped at both ends by a cylinder at 90 degree angle
! to create a U-shaped coil
! UCOIL_R1 is half the spacing between vertical cylinders
! UCOIL_R2 is the cylinders radius
! UCOIL_Y1 and UCOIL_Y2 are the min and max y-values of the coil
! The coil is translated in the x and z direction (t_x, t_z)
! and can be rotated about the y-direction, centered about (t_x,t_z)
! by the angle THETA_Y

            c = DCOS(THETA_Y(Q_ID))
            ss = DSIN(THETA_Y(Q_ID))

!           Translation
            xt = x1-t_x(Q_ID)
            yt = x2
            zt = x3-t_z(Q_ID)

!           Rotation
            xtr =  xt*c + zt*ss
            ytr =  yt
            ztr = -xt*ss + zt*c

            R1 = UCOIL_R1(Q_ID)
            R2 = UCOIL_R2(Q_ID)




! Limits of the cylinders.
! There is half a torus above ytest1 and half a torus below ytest2
            ytest1 = UCOIL_Y1(Q_ID) + R1 + R2
            ytest2 = UCOIL_Y2(Q_ID) - (R1 + R2)

            IF(ytest1<=ytr.AND.ytr<=ytest2) THEN
               f = 4*(xtr**2)*R1**2-(xtr**2+ztr**2+R1**2-R2**2)**2 ! a pair of cylinders
            ELSEIF(ytr<=ytest1) THEN
! Convert (x,y) into angle theta, and adjust its range from zero to 2*pi
               THETA  = ATAN2(ytr-ytest1,xtr) ! Result is from -pi to pi
               IF(-0.75*PI<=THETA.AND.THETA<=-0.25*PI) THEN
                  f = - ( ((ytr-(UCOIL_Y1(Q_ID) + R2))/R2)**2 + (ztr/R2)**2 -1.0 ) ! horizontal cylinder
               ELSE
                  f = 4*(xtr**2)*R1**2-(xtr**2+ztr**2+R1**2-R2**2)**2 ! a pair of cylinders
               ENDIF
            ELSEIF(ytr>=ytest2) THEN
! Convert (x,y) into angle theta, and adjust its range from zero to 2*pi
               THETA  = ATAN2(ytr-ytest2,xtr) ! Result is from -pi to pi
               IF(0.25*PI<=THETA.AND.THETA<0.75*PI) THEN
                  f = - ( ((ytr-(UCOIL_Y2(Q_ID) - R2))/R2)**2 + (ztr/R2)**2 -1.0 ) ! horizontal cylinder
               ELSE
                  f = 4*(xtr**2)*R1**2-(xtr**2+ztr**2+R1**2-R2**2)**2 ! a pair of cylinders
               ENDIF
            ENDIF



         ELSEIF(TRIM(quadric_form(Q_ID))=='XY_BEND_INT') THEN
! This shape represent a bend between two cylinders in the XY plane
! BEND_R1 is the radius of the bend
! BEND_R2 is the cylinders radius
! BEND_THETA1 is the orientation of the first cylinder (Deg.)
! BEND_THETA2 is the orientation of the second cylinder (Deg.).
! The orientation is measured from the y-axis.
! For example BEND_THETA1=0.0 and BEND_THETA2=90.0
! represents a 90 deg bend between a vertical (first cylinder)
! and a horizontal (second) cylinder.
! The bend itself is represented by a section of a torus
! The translation (t_x,t_y,t_z) defines the center of the bend (torus)
! The shape is defines as three pieces: 2 cylinders and a torus
! The switch between these pieces occur bases on angular position
! There is a complicated definition of intermediate angles to allow
! for a wide range of angles. This may be simplified in the future.
! Both BEND_THETA1 and BEND_THETA2 must be between zero and 360.0 degrees.
! Typically, BEND_THETA2 > BEND_THETA1, unless for a left bend,
! For example, BEND_THETA1=330.0 and BEND_THETA2=30.0
! will define a 60.0 deg. left-bend.
! Specifying a bend angle larger that 180 degrees will likey fail
! unless the bend is clipped somewhere else.

!           Translation
            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

            R1 = BEND_R1(Q_ID)
            R2 = BEND_R2(Q_ID)
! Convert angles from degrees to radians
            THETA1 = BEND_THETA1(Q_ID)*(pi/180.0D0)
            THETA2 = BEND_THETA2(Q_ID)*(pi/180.0D0)

! Convert (x,y) into angle theta, and adjust its range from zero to 2*pi
            THETA  = ATAN2(yt,xt) ! Result is from -pi to pi
            IF(THETA<ZERO) THETA = THETA + 2.0D0*PI

! THETA3 correspond to the point on a unit circle between THETA1 and THETA2
            IF(THETA2>THETA1) THEN
               THETA3m = HALF*(THETA1+THETA2)
               IF(THETA3m<PI) THEN
                  THETA3 = THETA3m + PI
               ELSE
                  THETA3 = THETA3m - PI
               ENDIF
            ELSE
               THETA3 = THETA2 + HALF * (THETA1-THETA2)
            ENDIF

! This angles are adjusted to wrap nicely along the discontinuity where
! Theta=zero or 2*pi. There may be a simpler way to do this...

            IF(THETA1>THETA3) THEN
               THETA1CYL1 = THETA1
            ELSE
               THETA1CYL1 = THETA1 + 2.0*PI
            ENDIF

            IF(THETA2>THETA1) THEN
               THETA2TOR = THETA2
            ELSE
               THETA2TOR = THETA2 + 2.0*PI
            ENDIF

            IF(THETA3>THETA2) THEN
               THETA3CYL2 = THETA3
            ELSE
               THETA3CYL2 = THETA3 + 2.0*PI
            ENDIF

            THETA_MIN = DMIN1(THETA1,THETA2,THETA3,THETA1CYL1,THETA2TOR,THETA3CYL2)

            IF(THETA<THETA_MIN) THETA = THETA + 2.0*PI

! Now join the pieces together:
            IF(THETA3<=THETA.AND.THETA<=THETA1CYL1) THEN  ! cylinder 1
               c = DCOS(THETA1)
               ss = DSIN(THETA1)
!              translation
               xt = x1-(t_x(Q_ID)+R1*c)
               yt = x2-(t_y(Q_ID)+R1*ss)
               zt = x3-t_z(Q_ID)
!              Rotation
               xtr =  xt*c  + yt*ss
               ytr = -xt*ss  + yt*c
               ztr = zt

               f = (xtr/R2)**2 + (ztr/R2)**2 -1.0

            ELSEIF(THETA1<=THETA.AND.THETA<=THETA2TOR) THEN  ! Torus=elbow
!              translation
               xt = x1-t_x(Q_ID)
               yt = x2-t_y(Q_ID)
               zt = x3-t_z(Q_ID)
!              Rotation
               xtr = xt
               ytr = yt
               ztr = zt

               f = -4*(xtr**2+ytr**2)*R1**2+(xtr**2+ytr**2+ztr**2+R1**2-R2**2)**2

            ELSEIF(THETA2<=THETA.AND.THETA<=THETA3CYL2) THEN  ! cylinder 2
               c = DCOS(THETA2)
               ss = DSIN(THETA2)
!              translation
               xt = x1-(t_x(Q_ID)+R1*c)
               yt = x2-(t_y(Q_ID)+R1*ss)
               zt = x3-t_z(Q_ID)
!              Rotation
               xtr =  xt*c  + yt*ss
               ytr = -xt*ss  + yt*c
               ztr = zt

               f = (xtr/R2)**2 + (ztr/R2)**2 -1.0

            ELSE
               WRITE(*,*)' Error in processing elbow.','THETA = ',THETA/PI*180.0
               CALL MFIX_EXIT(myPE)
            ENDIF

         ELSEIF(TRIM(quadric_form(Q_ID))=='Y_C2C_INT') THEN
! This shape connects two vertical cylinders by a conical section
! This is a more convenient way of doing it the traditional way
! with three quadrics (cylinder-cone-cylinder)

!           Translation
            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

!           Rotation
            xtr = xt
            ytr = yt
            ztr = zt

! Radii
            R1 = C2C_R1(Q_ID)
            R2 = C2C_R2(Q_ID)
! Extent of the conical section
            Y1 = C2C_Y1(Q_ID)
            Y2 = C2C_Y2(Q_ID)

            IF(ytr>=Y2) THEN
               R = R2
            ELSEIF(ytr<=Y1) THEN
               R = R1
            ELSE
               IF(Y2/=Y1) THEN ! when Y2=Y1, then R2=R1 (see check_data_cartesian)
                  R = R1 + (R2-R1)/(Y2-Y1)*(ytr-Y1)
               ELSE
                  R = R1
               ENDIF
            ENDIF

            f = (xtr/R)**2 + (ztr/R)**2 - 1.0


         ELSEIF(TRIM(quadric_form(Q_ID))=='REACTOR1') THEN
! This shape defines a reactor (interior flow), made of two vertical cylinders,
! connected by a conical section.
! Each cylinder is rounded and closed by a conical cap.

!           Translation
            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

!           Rotation
            xtr = xt
            ytr = yt
            ztr = zt

! Cylinders Radii
            R1  = REACTOR1_R1(Q_ID)  ! Lower section
            R2  = REACTOR1_R2(Q_ID)  ! Upper section

! Conical transition
            Y1 = REACTOR1_Y1(Q_ID)  ! Conical transition between lower
            Y2 = REACTOR1_Y2(Q_ID)  ! and upper sections

! Rounding
            YR1    = REACTOR1_YR1(Q_ID)  ! Lower rounding
            RR1    = REACTOR1_RR1(Q_ID)  ! Lower section rounding
            THETA1 = REACTOR1_THETA1(Q_ID) ! angle (radians)

            YR2    = REACTOR1_YR2(Q_ID)  ! Upper rounding
            RR2    = REACTOR1_RR2(Q_ID)  ! Lower section rounding
            THETA2 = REACTOR1_THETA2(Q_ID) ! angle (radians)


            YRR1 = YR1 - RR1 * DSIN(THETA1)
            RC1  = R1-RR1*(ONE-DCOS(THETA1))
            YC1  = YRR1 - RC1/DTAN(THETA1)
            YRR2 = YR2 + RR2 * DSIN(THETA2)
            RC2  = R2-RR2*(ONE-DCOS(THETA2))
            YC2  = YRR2 + RC2/DTAN(THETA2)

            IF(ytr>=YC2) THEN  ! Above upper conical cap

               R = ZERO

            ELSEIF(YRR2<=ytr.AND.ytr<=YC2) THEN   ! upper conical cap

               R = RC2/(YRR2-YC2)*(ytr-YC2)

            ELSEIF(YR2<=ytr.AND.ytr<=YRR2) THEN   ! upper rounding

               R = R2 - RR2 + DSQRT(RR2**2 - (YR2-ytr)**2)

            ELSEIF(Y2<=ytr.AND.ytr<=YR2) THEN     ! upper cylinder

               R = R2

            ELSEIF(Y1<=ytr.AND.ytr<=Y2) THEN      ! conical transition

               R = R1 + (R2-R1)/(Y2-Y1)*(ytr-Y1)

            ELSEIF(YR1<=ytr.AND.ytr<=Y1) THEN     ! lower cylinder

               R = R1

            ELSEIF(YRR1<=ytr.AND.ytr<=YR1) THEN   ! lower rounding

               R = R1 - RR1 + DSQRT(RR1**2 - (YR1-ytr)**2)

            ELSEIF(YC1<=ytr.AND.ytr<=YRR1) THEN   ! lower conical cap

               R = RC1/(YRR1-YC1)*(ytr-YC1)

            ELSE

               R = ZERO

            ENDIF


            IF(R>ZERO) THEN
               f = (xtr/R)**2 + (ztr/R)**2 - 1.0
            ELSE
               IF(ytr<=YC1) THEN
                  f = YC1 - ytr     ! below lower conical cap
               ELSE
                  f = ytr - YC2     ! above upper conical cap
               ENDIF
            ENDIF


         ELSE

            CALL BUILD_1x3_MATRIX(x1,x2,x3,X_VECTOR)

            XMT = X_VECTOR - T_QUADRIC(:,:,Q_ID)
            TXMT = TRANSPOSE(XMT)

            TEMP_1x1 = MATMUL(XMT,MATMUL(A_QUADRIC(:,:,Q_ID),TXMT))

            f = TEMP_1x1(1,1) + dquadric(Q_ID)

         ENDIF

! Each clipping limit is treated as a plane. For example, fxmin is
! the equation of the plane describing x=xmin, and a value of fxmin
! is compared with the current value of f to determine if the location
! is part of the computational domain. The comparison (min of max)
! follows the same logis as the 'AND' (max) , or 'OR' (min)
! logic when combining two quadrics.
! The clipping procedure is ignored when CLIP_FLAG is .FALSE.
! This will happen when we are in a 'PIECEWISE' group



         IF(FLUID_IN_CLIPPED_REGION(Q_ID)) THEN

            IF(clip_xmin(Q_ID)/=UNDEFINED) THEN
               fxmin = -(clip_xmin(Q_ID)-x1)
               f = dmin1(f,fxmin)
            ENDIF

            IF(clip_xmax(Q_ID)/=UNDEFINED) THEN
               fxmax = -(x1-clip_xmax(Q_ID))
               f = dmin1(f,fxmax)
            ENDIF

            IF(clip_ymin(Q_ID)/=UNDEFINED) THEN
               fymin = -(clip_ymin(Q_ID)-x2)
               f = dmin1(f,fymin)
            ENDIF

            IF(clip_ymax(Q_ID)/=UNDEFINED) THEN
               fymax = -(x2-clip_ymax(Q_ID))
               f = dmin1(f,fymax)
            ENDIF

            IF(clip_zmin(Q_ID)/=UNDEFINED) THEN
               fzmin = -(clip_zmin(Q_ID)-x3)
               f = dmin1(f,fzmin)
            ENDIF

            IF(clip_zmax(Q_ID)/=UNDEFINED) THEN
               fzmax = -(x3-clip_zmax(Q_ID))
               f = dmin1(f,fzmax)
            ENDIF

         ELSE

            IF(clip_xmin(Q_ID)/=UNDEFINED) THEN
               fxmin = clip_xmin(Q_ID)-x1
               f = dmax1(f,fxmin)
            ENDIF

            IF(clip_xmax(Q_ID)/=UNDEFINED) THEN
               fxmax = x1-clip_xmax(Q_ID)
               f = dmax1(f,fxmax)
            ENDIF

            IF(clip_ymin(Q_ID)/=UNDEFINED) THEN
               fymin = clip_ymin(Q_ID)-x2
               f = dmax1(f,fymin)
            ENDIF

            IF(clip_ymax(Q_ID)/=UNDEFINED) THEN
               fymax = x2-clip_ymax(Q_ID)
               f = dmax1(f,fymax)
            ENDIF

            IF(clip_zmin(Q_ID)/=UNDEFINED) THEN
               fzmin = clip_zmin(Q_ID)-x3
               f = dmax1(f,fzmin)
            ENDIF

            IF(clip_zmax(Q_ID)/=UNDEFINED) THEN
               fzmax = x3-clip_zmax(Q_ID)
               f = dmax1(f,fzmax)
            ENDIF

         ENDIF

      RETURN
      END SUBROUTINE GET_F_QUADRIC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REASSIGN_QUADRIC                                       C
!  Purpose: Reassign the quadric based on location                     C
!                                                                      C
!  Author: Jeff Dietiker                               Date: 21-FEB-08 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!  Modified: ##                                        Date: ##-###-## C
!  Purpose: ##                                                         C
!                                                                      C
!  Literature/Document References: ##                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE REASSSIGN_QUADRIC(x1,x2,x3,GROUP,Q_ID)

      USE compar
      USE exit, only: mfix_exit
      USE quadric

      IMPLICIT NONE

      DOUBLE PRECISION x1,x2,x3
      INTEGER :: I,Q_ID,GROUP,GS,P
      LOGICAL :: PIECE_X,PIECE_Y,PIECE_Z,PIECE_FLAG
      CHARACTER(LEN=9) :: GR

      Q_ID = 0

      GS = GROUP_SIZE(GROUP)
      GR = TRIM(GROUP_RELATION(GROUP))

      IF( GR /= 'PIECEWISE') RETURN

      DO P = 1 , GS

         I = GROUP_Q(GROUP,P)

         PIECE_X = (piece_xmin(I) <= x1).AND.( x1 <= piece_xmax(I))
         PIECE_Y = (piece_ymin(I) <= x2).AND.( x2 <= piece_ymax(I))
         PIECE_Z = (piece_zmin(I) <= x3).AND.( x3 <= piece_zmax(I))

         PIECE_FLAG = (PIECE_X.AND.PIECE_Y.AND.PIECE_Z)

         IF (PIECE_FLAG) Q_ID = I

      ENDDO

      IF(Q_ID == 0 ) THEN
         WRITE(*,*)' No Quadric defined at current location x,y,z=', x1,x2,x3
         WRITE(*,*)' Please Check piecewise limits of quadric(s)'
         WRITE(*,*)' Mfix will exit now.'
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN


      END SUBROUTINE REASSSIGN_QUADRIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_1x3_MATRIX                                       C
!  Purpose: Catesian Grid - Build a (1x3) matrix                       C
!           from 3 scalars                                             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_1x3_MATRIX(scalar1,scalar2,scalar3,M1x3)

      IMPLICIT NONE

      DOUBLE PRECISION:: scalar1,scalar2,scalar3
      DOUBLE PRECISION, DIMENSION(1,3) :: M1x3
      M1x3(1,1) = scalar1
      M1x3(1,2) = scalar2
      M1x3(1,3) = scalar3

      RETURN

  END SUBROUTINE BUILD_1x3_MATRIX


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_C_QUADRIC_MATRIX                                 C
!  Purpose: Catesian Grid - Build a (3x3) diagonal matrix              C
!           whose diagonal elements are the characteristic values of   C
!           the quadric                                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_C_QUADRIC_MATRIX(lambda1,lambda2,lambda3,C_QUADRIC)

      USE param1, only: zero

      IMPLICIT NONE

      DOUBLE PRECISION:: lambda1,lambda2,lambda3
      DOUBLE PRECISION, DIMENSION(3,3) :: C_QUADRIC

!      Transpose is used because matrices are stored column-wise

      C_QUADRIC = TRANSPOSE(RESHAPE((/                                           &
                                          lambda1  ,   ZERO      ,  ZERO        ,&
                                          ZERO     ,   lambda2   ,  ZERO        ,&
                                          ZERO     ,   ZERO      ,  lambda3   /),&
                                                                                  (/3,3/)))
      RETURN

      END SUBROUTINE BUILD_C_QUADRIC_MATRIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_X_ROTATION_MATRIX                                C
!  Purpose: Catesian Grid - Build a (3x3) rotation matrix about x-axis,C
!           given the rotation angle in degrees                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_X_ROTATION_MATRIX(Theta, R)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: Theta
      DOUBLE PRECISION, DIMENSION(3,3) :: R

      theta = theta * (pi/180.0D0) ! Rotation angle about x-axis (radians)

!      Transpose is used because matrices are stored column-wise

      R = TRANSPOSE(RESHAPE((/                                                  &
                                  ONE   ,   ZERO         ,  ZERO               ,&
                                  ZERO  ,   dcos(theta)  ,  dsin(theta)        ,&
                                  ZERO  ,  -dsin(theta)  ,  dcos(theta)      /),&
                                                                                  (/3,3/)))

      RETURN

      END SUBROUTINE BUILD_X_ROTATION_MATRIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_Y_ROTATION_MATRIX                                C
!  Purpose: Catesian Grid - Build a (3x3) rotation matrix about y-axis,C
!           given the rotation angle in degrees                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_Y_ROTATION_MATRIX(Theta, R)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: Theta
      DOUBLE PRECISION, DIMENSION(3,3) :: R

      theta = theta * (pi/180.0D0) ! Rotation angle about x-axis (radians)

!      Transpose is used because matrices are stored column-wise

      R = TRANSPOSE(RESHAPE((/                                                &
                                  dcos(theta)  ,  ZERO  ,  dsin(theta)       ,&
                                  ZERO         ,  ONE   ,  ZERO              ,&
                                 -dsin(theta)  ,  ZERO  ,  dcos(theta)     /),&
                                                                                (/3,3/)))

      RETURN

      END SUBROUTINE BUILD_Y_ROTATION_MATRIX


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_Z_ROTATION_MATRIX                                C
!  Purpose: Catesian Grid - Build a (3x3) rotation matrix about z-axis,C
!           given the rotation angle in degrees                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_Z_ROTATION_MATRIX(Theta, R)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: Theta
      DOUBLE PRECISION, DIMENSION(3,3) :: R

      theta = theta * (pi/180.0D0) ! Rotation angle about x-axis (radians)

!      Transpose is used because matrices are stored column-wise

      R = TRANSPOSE(RESHAPE((/                                                 &
                                   dcos(theta)  ,  dsin(theta)  ,  ZERO       ,&
                                  -dsin(theta)  ,  dcos(theta)  ,  ZERO       ,&
                                   ZERO         ,  ZERO         ,  ONE      /),&
                                                                                 (/3,3/)))

      RETURN

      END SUBROUTINE BUILD_Z_ROTATION_MATRIX
