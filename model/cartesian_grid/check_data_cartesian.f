MODULE CHECK_DATA_CG
   USE exit, only: mfix_exit

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_CARTESIAN                                   C
!  Purpose: check the data related to cartesian grid implementation    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_CARTESIAN

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE cutcell
      USE dashboard
      USE discretelement
      USE funits
      USE indices
      USE iterate, ONLY: max_nit
      USE leqsol
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE polygon
      USE quadric
      USE run
      USE scalars
      USE stl
      USE vtk
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I,J,Q,BCV
      DOUBLE PRECISION :: norm, tan_half_angle
      CHARACTER(LEN=9) :: GR
!-----------------------------------------------

      IF(.NOT.CARTESIAN_GRID) RETURN

      IF(DISCRETE_ELEMENT) THEN
         IF(MyPE == PE_IO) THEN

            WRITE(*,10)'######################################################################'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##  ===>   WARNING: RUNNING CARTESIAN GRID WITH DISCRETE ELEMENT.   ##'
            WRITE(*,10)'##                  THIS HAS NOT BEEN FULLY TESTED.                 ##'
            WRITE(*,10)'##                  PLEASE USE WITH CAUTION.                        ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'######################################################################'

         ENDIF
      ENDIF

10    FORMAT(1X,A)

      IF(COORDINATES=='CYLINDRICAL') THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT AVAILABLE'
            WRITE(*,*)'WITH CYLINDRICAL COORDINATE SYSTEM.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(USE_STL.AND.(.NOT.USE_MSH)) THEN
         IF(DO_K) THEN
            CALL GET_STL_DATA
         ELSE
            IF(MyPE == PE_IO) WRITE(*,*) &
               'ERROR: STL METHOD VALID ONLY IN 3D.'
            CALL MFIX_EXIT(MYPE)
         ENDIF
         IF(N_QUADRIC > 0) THEN
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND ',&
               'STL INPUT ARE SPECIFIED.'
            IF(MyPE == PE_IO) WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE ',&
               'OF SURFACE INPUT.'
            CALL MFIX_EXIT(MYPE)
         ENDIF
      ENDIF

      IF(USE_MSH.AND.(.NOT.USE_STL)) THEN
         IF(DO_K) THEN
            CALL GET_MSH_DATA
         ELSE
            IF(MyPE == PE_IO) WRITE(*,*) &
               'ERROR: MSH METHOD VALID ONLY IN 3D.'
            CALL MFIX_EXIT(MYPE)
         ENDIF
         IF(N_QUADRIC > 0) THEN
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND ',&
               'MSH INPUT ARE SPECIFIED.'
            IF(MyPE == PE_IO) WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE ',&
               'OF SURFACE INPUT.'
            CALL MFIX_EXIT(MYPE)
         ENDIF
      ENDIF

      IF(USE_POLYGON) THEN
         IF(DO_K) THEN
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: POLYGON METHOD ',&
               'VALID ONLY IN 2D.'
            CALL MFIX_EXIT(MYPE)
         ELSE
            CALL GET_POLY_DATA
         ENDIF
      ENDIF

      IF(N_QUADRIC > 0) THEN
         IF(N_POLYGON > 0) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND POLYGON(S) ',&
                  'DEFINED.'
               WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE INPUT.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF
         IF(N_USR_DEF > 0) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND USER-DEFINED ',&
                  'FUNCTION DEFINED.'
               WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF
         IF(QUADRIC_SCALE <= ZERO) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: QUADRIC_SCALE MUST BE POSITIVE.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ELSEIF(QUADRIC_SCALE /= ONE) THEN
            DO Q = 1, N_QUADRIC
               lambda_x(Q)  = lambda_x(Q)  * quadric_scale**2
               lambda_y(Q)  = lambda_y(Q)  * quadric_scale**2
               lambda_z(Q)  = lambda_z(Q)  * quadric_scale**2
               Radius(Q)    = Radius(Q)    * quadric_scale
               t_x(Q)       = t_x(Q)       * quadric_scale
               t_y(Q)       = t_y(Q)       * quadric_scale
               t_z(Q)       = t_z(Q)       * quadric_scale
               clip_xmin(Q) = clip_xmin(Q) * quadric_scale
               clip_xmax(Q) = clip_xmax(Q) * quadric_scale
               clip_ymin(Q) = clip_ymin(Q) * quadric_scale
               clip_ymax(Q) = clip_ymax(Q) * quadric_scale
               clip_zmin(Q) = clip_zmin(Q) * quadric_scale
               clip_zmax(Q) = clip_zmax(Q) * quadric_scale
               piece_xmin(Q) = piece_xmin(Q) * quadric_scale
               piece_xmax(Q) = piece_xmax(Q) * quadric_scale
               piece_ymin(Q) = piece_ymin(Q) * quadric_scale
               piece_ymax(Q) = piece_ymax(Q) * quadric_scale
               piece_zmin(Q) = piece_zmin(Q) * quadric_scale
               piece_zmax(Q) = piece_zmax(Q) * quadric_scale
            ENDDO
         ENDIF
      ELSE
         IF((N_POLYGON > 0).AND.(N_USR_DEF > 0)) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: POLYGON(S) AND USER-DEFINED ',&
                  'FUNCTION DEFINED.'
               WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF
      ENDIF


      IF(N_QUADRIC > DIM_QUADRIC) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF N_QUADRIC =', &
               N_QUADRIC
            WRITE(*,*)'MAXIMUM ACCEPTABLE VALUE IS DIM_QUADRIC =', &
               DIM_QUADRIC
            WRITE(*,*)'CHANGE MAXIMUM VALUE IN QUADRIC_MOD.F ',&
               'IF NECESSARY.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


      DO Q = 1, N_QUADRIC



         SELECT CASE (TRIM(QUADRIC_FORM(Q)))

            CASE ('NORMAL')

               lambda_x(Q) = lambda_x(Q)
               lambda_y(Q) = lambda_y(Q)
               lambda_z(Q) = lambda_z(Q)

               norm = dsqrt(lambda_x(Q)**2 + lambda_y(Q)**2 + &
                            lambda_z(Q)**2)

               IF(norm < TOL_F) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, &
                        ' HAS ZERO COEFFICIENTS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            CASE ('PLANE')   ! The quadric is predefined as a plane

               lambda_x(Q) = n_x(Q)
               lambda_y(Q) = n_y(Q)
               lambda_z(Q) = n_z(Q)

               norm = dsqrt(lambda_x(Q)**2 + lambda_y(Q)**2 + &
                            lambda_z(Q)**2)

               IF( norm > TOL_F) THEN
                  lambda_x(Q) = lambda_x(Q) / norm
                  lambda_y(Q) = lambda_y(Q) / norm
                  lambda_z(Q) = lambda_z(Q) / norm
               ELSE
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: PLANE:', Q, &
                        ' HAS ZERO NORMAL VECTOR.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

              dquadric(Q) = - (lambda_x(Q)*t_x(Q) + lambda_y(Q)*t_y(Q) + &
                 lambda_z(Q)*t_z(Q))

            CASE ('X_CYL_INT')   ! The quadric is predefined as a cylinder, along x-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, &
                        ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = ZERO
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

            CASE ('Y_CYL_INT')   ! The quadric is predefined as a cylinder, along y-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, &
                        ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ZERO
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

            CASE ('Z_CYL_INT')   ! The quadric is predefined as a cylinder, along z-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, &
                        ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ZERO
                  dquadric(Q) = -Radius(Q)**2
               ENDIF


            CASE ('X_CYL_EXT')   ! The quadric is predefined as a cylinder, along x-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, &
                        ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = ZERO
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('Y_CYL_EXT')   ! The quadric is predefined as a cylinder, along y-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, &
                        ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = ZERO
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('Z_CYL_EXT')   ! The quadric is predefined as a cylinder, along z-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, &
                        ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = ZERO
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('SPHERE_INT')   ! The quadric is predefined as a sphere
                                  ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: SPHERE:', Q, &
                     ' HAS INVALID RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

           CASE ('SPHERE_EXT')   ! The quadric is predefined as a sphere
                                  ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: SPHERE:', Q, &
                     ' HAS INVALID RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF


            CASE ('X_CONE')    ! The quadric is predefined as a cone, along x-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CONE:', Q, &
                        ' HAS INCORRECT HALF-ANGLE.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = ONE/(tan_half_angle)**2
                  lambda_z(Q) = ONE/(tan_half_angle)**2
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('Y_CONE')    ! The quadric is predefined as a cone, along y-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CONE:', Q, &
                        ' HAS INCORRECT HALF-ANGLE.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = ONE/(tan_half_angle)**2
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = ONE/(tan_half_angle)**2
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('Z_CONE')    ! The quadric is predefined as a cone, along z-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CONE:', Q, &
                        ' HAS INCORRECT HALF-ANGLE.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = ONE/(tan_half_angle)**2
                  lambda_y(Q) = ONE/(tan_half_angle)**2
                  lambda_z(Q) = -ONE
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('C2C')        ! Cylinder to cylinder junction using cone
                                ! Internal flow

               CALL BUILD_CONE_FOR_C2C(Q)


            CASE ('TORUS_INT','TORUS_EXT')      ! Torus - Hard coded in define_quadrics.f
               IF((Torus_R1(Q) <= ZERO).OR.(Torus_R1(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: TORUS:', Q, &
                        ' HAS INVALID RADIUS R1:',Torus_R1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF
               IF((Torus_R2(Q) <= ZERO).OR.(Torus_R2(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: TORUS:', Q, &
                        ' HAS INVALID RADIUS R2:',Torus_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            CASE ('Y_UCOIL_EXT')      ! UCOIL - Hard coded in define_quadrics.f
               IF((UCOIL_R1(Q) <= ZERO).OR.(UCOIL_R1(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_UCOIL_EXT:', Q, &
                        ' HAS INVALID RADIUS R1:',UCOIL_R1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF
               IF((UCOIL_R2(Q) <= ZERO).OR.(UCOIL_R2(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_UCOIL_EXT:', Q, &
                        ' HAS INVALID RADIUS R2:',UCOIL_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF(UCOIL_Y2(Q)<UCOIL_Y1(Q)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_UCOIL_EXT:', Q, &
                        ' COIL_Y2 < COIL_Y1: Y2,Y1=',UCOIL_Y2(Q),UCOIL_Y1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            CASE ('Y_UCOIL2_EXT')      ! UCOIL - Hard coded in define_quadrics.f
               IF((UCOIL_R1(Q) <= ZERO).OR.(UCOIL_R1(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_UCOIL2_EXT:', Q, &
                        ' HAS INVALID RADIUS R1:',UCOIL_R1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF
               IF((UCOIL_R2(Q) <= ZERO).OR.(UCOIL_R2(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_UCOIL2_EXT:', Q, &
                        ' HAS INVALID RADIUS R2:',UCOIL_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF(UCOIL_Y2(Q)<UCOIL_Y1(Q)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_UCOIL2_EXT:', Q, &
                        ' COIL_Y2 < COIL_Y1: Y2,Y1=',UCOIL_Y2(Q),UCOIL_Y1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            CASE ('XY_BEND_INT')       ! Bend  - Hard coded in define_quadrics.f
               IF((BEND_R1(Q) <= ZERO).OR.(BEND_R1(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: XY_BEND_INT:', Q, &
                        ' HAS INVALID RADIUS R1:',BEND_R1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF((BEND_R2(Q) <= ZERO).OR.(BEND_R2(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: XY_BEND_INT:', Q, &
                        ' HAS INVALID RADIUS R2:',BEND_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF((BEND_THETA1(Q) < ZERO).OR.(BEND_THETA1(Q)>360.0)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: XY_BEND_INT:', Q, &
                        ' HAS INVALID ANGLE THETA1:',BEND_THETA1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF((BEND_THETA2(Q) < ZERO).OR.(BEND_THETA2(Q)>360.0)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: XY_BEND_INT:', Q, &
                        ' HAS INVALID ANGLE THETA2:',BEND_THETA2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            CASE ('Y_C2C_INT')       ! Cylinder-cone-cylinder  - Hard coded in define_quadrics.f
               IF((C2C_R1(Q) <= ZERO).OR.(C2C_R1(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_C2C_INT:', Q, &
                        ' HAS INVALID RADIUS R1:',C2C_R1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF((C2C_R2(Q) <= ZERO).OR.(C2C_R2(Q)==UNDEFINED)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: C2C_XY_INT:', Q, &
                        ' HAS INVALID RADIUS R2:',C2C_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF(C2C_Y2(Q) < C2C_Y1(Q)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_C2C_INT:', Q
                     WRITE(*,*)'MUST HAVE C2C_Y2 >= C2C_Y1.'
                     WRITE(*,*)'C2C_Y1,C2C_Y2 =', C2C_Y1(Q),C2C_Y2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF((C2C_Y1(Q) == C2C_Y2(Q)).AND.(C2C_R1(Q)/=C2C_R2(Q))) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: Y_C2C_INT:', Q, &
                        ' C2C_Y1=C2C_Y2 BUT C2C_R1/=C2C_R2:'
                     WRITE(*,*)'C2C_Y1,C2C_Y2 =', C2C_Y1(Q),C2C_Y2(Q)
                     WRITE(*,*)'C2C_R1,C2C_R2 =', C2C_R1(Q),C2C_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            CASE ('REACTOR1')       ! Cylinder-cone-cylinder  - Hard coded in define_quadrics.f

               IF(REACTOR1_Y2(Q) < REACTOR1_Y1(Q)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: REACTOR1:', Q
                     WRITE(*,*)'MUST HAVE REACTOR1_Y2 >= REACTOR1_Y1.'
                     WRITE(*,*)'REACTOR1_Y1,REACTOR1_Y2 =', REACTOR1_Y1(Q),REACTOR1_Y2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF((REACTOR1_Y1(Q) == REACTOR1_Y2(Q)).AND.(REACTOR1_R1(Q)/=REACTOR1_R2(Q))) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: REACTOR1:', Q, &
                        ' REACTOR1_Y1=REACTOR1_Y2 BUT REACTOR1_R1/=REACTOR1_R2:'
                     WRITE(*,*)'REACTOR1_Y1,REACTOR1_Y2 =', REACTOR1_Y1(Q),REACTOR1_Y2(Q)
                     WRITE(*,*)'REACTOR1_R1,REACTOR1_R2 =', REACTOR1_R1(Q),REACTOR1_R2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF


               IF(REACTOR1_YR2(Q) <= REACTOR1_Y2(Q)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: REACTOR1:', Q
                     WRITE(*,*)'MUST HAVE REACTOR1_YR2 > REACTOR1_Y2.'
                     WRITE(*,*)'REACTOR1_YR2,REACTOR1_Y2 =', REACTOR1_YR2(Q),REACTOR1_Y2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF(REACTOR1_YR1(Q) >= REACTOR1_Y1(Q)) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: REACTOR1:', Q
                     WRITE(*,*)'MUST HAVE REACTOR1_YR1 < REACTOR1_Y1.'
                     WRITE(*,*)'REACTOR1_YR1,REACTOR1_Y1 =', REACTOR1_YR1(Q),REACTOR1_Y1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF(REACTOR1_THETA1(Q) <= ZERO.OR.REACTOR1_THETA1(Q) > 90.0D0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: REACTOR1:', Q
                     WRITE(*,*)'MUST HAVE 0.0 < REACTOR1_THETA1 <= 90 DEGREES.'
                     WRITE(*,*)'REACTOR1_THETA1 =', REACTOR1_THETA1(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF

               IF(REACTOR1_THETA2(Q) <= ZERO.OR.REACTOR1_THETA2(Q) > 90.0D0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: REACTOR1:', Q
                     WRITE(*,*)'MUST HAVE 0.0 < REACTOR1_THETA2 <= 90 DEGREES.'
                     WRITE(*,*)'REACTOR1_THETA2 =', REACTOR1_THETA2(Q)
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF
! Convert angles from degrees to radians
               REACTOR1_THETA1(Q) = REACTOR1_THETA1(Q)/180.0D0*PI
               REACTOR1_THETA2(Q) = REACTOR1_THETA2(Q)/180.0D0*PI


            CASE DEFAULT
               IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, &
                     ' HAS INCORRECT FORM: ',quadric_form(Q)
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)

         END SELECT

         IF(BC_ID_Q(Q) == UNDEFINED_I) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, &
                  ' HAS NO ASSIGNED BC ID.'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF

      ENDDO


      IF(N_QUADRIC>0) THEN


         IF(N_GROUP > DIM_GROUP) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'INPUT ERROR: INVALID VALUE OF N_GROUP =', N_GROUP
               WRITE(*,*)'MAXIMUM ACCEPTABLE VALUE IS DIM_GROUP =', DIM_GROUP
               WRITE(*,*)'CHANGE MAXIMUM VALUE IN QUADRIC_MOD.F IF NECESSARY.'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF


         DO I = 1,N_GROUP

            IF(GROUP_SIZE(I) < 1 .OR. GROUP_SIZE(I) > N_QUADRIC) THEN
              IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT SIZE:', GROUP_SIZE(I)
                  WRITE(*,*)'VALID GROUP SIZE RANGE IS:', 1, ' TO ', N_QUADRIC
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)
            ENDIF

            DO J = 1,GROUP_SIZE(I)
               IF(GROUP_Q(I,J) < 1 .OR. GROUP_Q(I,J) > N_QUADRIC) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: GROUP_Q(', I,',',J, ') HAS INCORRECT VALUE:', GROUP_Q(I,J)
                     WRITE(*,*)'VALID GROUP_Q RANGE IS:', 1, ' TO ', N_QUADRIC
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF
            ENDDO

            GR = TRIM(GROUP_RELATION(I))

            IF(GR/='OR'.AND.GR/='AND'.AND.GR/='PIECEWISE') THEN
               IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT GROUP RELATION: ', GR
                  WRITE(*,*)'VALID GROUP RELATIONS ARE ''OR'',''AND'', AND ''PIECEWISE''. '
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)
            ENDIF

         ENDDO

         DO I = 2,N_GROUP

            GR = TRIM(RELATION_WITH_PREVIOUS(I))

            IF(GR/='OR'.AND.GR/='AND') THEN
               IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT RELATION WITH PREVIOUS: ', GR
                  WRITE(*,*)'VALID GROUP RELATIONS ARE ''OR'', AND ''AND''. '
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)
            ENDIF

         ENDDO

      ENDIF


      IF(TOL_SNAP(1)<ZERO.OR.TOL_SNAP(1)>HALF) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN X-DIRECTION =', TOL_SNAP(1)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SNAP(2)==UNDEFINED) TOL_SNAP(2)=TOL_SNAP(1)

      IF(TOL_SNAP(2)<ZERO.OR.TOL_SNAP(2)>HALF) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN Y-DIRECTION =', TOL_SNAP(2)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SNAP(3)==UNDEFINED) TOL_SNAP(3)=TOL_SNAP(1)

      IF(TOL_SNAP(3)<ZERO.OR.TOL_SNAP(3)>HALF) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN Z-DIRECTION =', TOL_SNAP(3)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


      IF(TOL_DELH<ZERO.OR.TOL_DELH>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_DELH =', TOL_DELH
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SMALL_CELL<ZERO.OR.TOL_SMALL_CELL>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SMALL_CELL =', TOL_SMALL_CELL
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SMALL_AREA<ZERO.OR.TOL_SMALL_AREA>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SMALL_AREA =', TOL_SMALL_AREA
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(ALPHA_MAX<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF ALPHA_MAX =', ALPHA_MAX
            WRITE(*,*)'ACCEPTABLE VALUES ARE POSITIVE NUMBERS (E.G. 1.0).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


      IF(TOL_F<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF TOL_F =', TOL_F
            WRITE(*,*)'ACCEPTABLE VALUES ARE SMALL POSITIVE NUMBERS (E.G. 1.0E-9).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_POLY<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF TOL_POLY =', TOL_POLY
            WRITE(*,*)'ACCEPTABLE VALUES ARE SMALL POSITIVE NUMBERS (E.G. 1.0E-9).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(ITERMAX_INT<0) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF ITERMAX_INT =', ITERMAX_INT
            WRITE(*,*)'ACCEPTABLE VALUES ARE LARGE POSITIVE INTEGERS (E.G. 10000).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(FAC_DIM_MAX_CUT_CELL<0.05.OR.FAC_DIM_MAX_CUT_CELL>5.0) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF FAC_DIM_MAX_CUT_CELL =', FAC_DIM_MAX_CUT_CELL
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.05 AND 5.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


      IF((CG_SAFE_MODE(1)==1).AND.(PG_OPTION/=0)) THEN
         PG_OPTION = 0
         IF(MyPE == PE_IO) WRITE(*,*)'WARNING: SAFE_MODE ACTIVATED FOR GAS PRESSURE, REVERTING TO PG_OPTION = 0'
      ENDIF

      IF(PG_OPTION <0 .OR. PG_OPTION>2) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF PG_OPTION =', PG_OPTION
            WRITE(*,*)'ACCEPTABLE VALUES ARE 0,1,AND 2.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(CG_UR_FAC(2)<ZERO.OR.CG_UR_FAC(2)>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF CG_UR_FAC(2) =', CG_UR_FAC(2)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(BAR_WIDTH<10.OR.BAR_WIDTH>80) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF BAR_WIDTH =', BAR_WIDTH
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 10 AND 80.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(BAR_RESOLUTION<ONE.OR.BAR_RESOLUTION>100.0) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF BAR_RESOLUTION =', BAR_RESOLUTION
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 100.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(F_DASHBOARD<1) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF F_DASHBOARD =', F_DASHBOARD
            WRITE(*,*)'ACCEPTABLE VALUES ARE INTEGERS >= 1.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF



!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER

      DTMIN =  LARGE_NUMBER
      DTMAX = -LARGE_NUMBER

      NIT_MIN = MAX_NIT
      NIT_MAX = 0

      N_DASHBOARD = 0



      CG_MI_CONVERTED_TO_PS = .FALSE.

      DO BCV = 1, DIMENSION_BC

         IF (BC_TYPE_ENUM(BCV) == CG_MI) THEN
            BC_TYPE_ENUM(BCV) = CG_NSW
            CG_MI_CONVERTED_TO_PS(BCV) = .TRUE.
            if(MyPE==0) print*,'From check_data_cartesian: Converted CG_MI to CG_NSW for BC#',BCV
         ENDIF

      ENDDO


      IF(RE_INDEXING) THEN

         IF(MyPE==0) THEN
            WRITE(*,*)' From check_data_cartesian: RE_INDEXING is turned on.'
            WRITE(*,*)' The preconditionner will be turned off for all equations'
            WRITE(*,*)' regardless of the mfix.dat setting.'
            LEQ_PC = 'NONE'
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE CHECK_DATA_CARTESIAN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_BC_FLAGS                                         C
!  Purpose: check the boundary conditions flags                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_BC_FLAGS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE cutcell
      USE dashboard
      USE fldvar
      USE functions
      USE funits
      USE indices
      USE leqsol
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE polygon
      USE quadric
      USE run
      USE scalars
      USE toleranc
      USE vtk
      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: IJK,IJKW,IJKS,IJKB,M,NN
      INTEGER :: IJKWW,IJKSS,IJKBB
      INTEGER :: BCV,BCV_U,BCV_V,BCV_W
!-----------------------------------------------
      DOUBLE PRECISION SUM, SUM_EP
!-----------------------------------------------
!======================================================================
! Boundary conditions
!======================================================================

      CALL INIT_ERR_MSG("CHECK_BC_FLAGS")

      DO BCV = 1, DIMENSION_BC
         IF(CG_MI_CONVERTED_TO_PS(BCV)) THEN
            IF(BC_MASSFLOW_g(BCV)==UNDEFINED) THEN
               WRITE(ERR_MSG, 1710) BCV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            DO M = 1, MMAX
               IF(BC_MASSFLOW_s(BCV,M)==UNDEFINED) THEN
                  WRITE(ERR_MSG, 1711) BCV, M
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

1710 FORMAT('Error 1110: BC :',I3,'. When using CG_MI, the gas mass flow rate',/1X, &
         'must be specified, including when it is zero.',/1X, &
         ' Please correct the mfix.dat file.')

1711 FORMAT('Error 1111: BC :',I3,'. When using CG_MI, the solids mass flow rate',/1X, &
         'for M=',I4,' must be specified, including when it is zero.',/1X, &
         ' Please correct the mfix.dat file.')

      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN

            IF(BC_TYPE_ENUM(BCV)  == CG_MI) THEN

               print*,'CG_MI at', IJK  ! This should not be printed on the screen anymore after conversion to point source.

!               FLAG(IJK) = 20
!               FLAG_E(IJK) = UNDEFINED_I
!               FLAG_N(IJK) = UNDEFINED_I
!               FLAG_T(IJK) = UNDEFINED_I

            ELSEIF(BC_TYPE_ENUM(BCV)  == CG_PO) THEN

               FLAG(IJK) = 11
               FLAG_E(IJK) = UNDEFINED_I
               FLAG_N(IJK) = UNDEFINED_I
               FLAG_T(IJK) = UNDEFINED_I

               IJKW = WEST_OF(IJK)
               BCV_U = BC_U_ID(IJKW)
               IF(BCV_U>0) THEN
                  IF(BC_TYPE_ENUM(BCV_U)  == CG_PO) THEN
                     FLAG(IJKW) = 11
                     FLAG_E(IJKW) = UNDEFINED_I
                     FLAG_N(IJKW) = UNDEFINED_I
                     FLAG_T(IJKW) = UNDEFINED_I
                  ENDIF
               ENDIF

               IJKS = SOUTH_OF(IJK)
               BCV_V = BC_V_ID(IJKS)
               IF(BCV_V>0) THEN
                  IF(BC_TYPE_ENUM(BCV_V)  == CG_PO) THEN
                     FLAG(IJKS) = 11
                     FLAG_E(IJKS) = UNDEFINED_I
                     FLAG_N(IJKS) = UNDEFINED_I
                     FLAG_T(IJKS) = UNDEFINED_I
                  ENDIF
               ENDIF

               IF(DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  BCV_W = BC_W_ID(IJKB)
                  IF(BCV_W>0) THEN
                     IF(BC_TYPE_ENUM(BCV_W)  == CG_PO) THEN
                        FLAG(IJKB) = 11
                        FLAG_E(IJKB) = UNDEFINED_I
                        FLAG_N(IJKB) = UNDEFINED_I
                        FLAG_T(IJKB) = UNDEFINED_I
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
         ENDIF
      ENDDO


      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN
            IF(BC_TYPE_ENUM(BCV)  == CG_MI) THEN

!               IJKW = WEST_OF(IJK)
!               IF(FLUID_AT(IJKW)) THEN
!                  FLAG_E(IJKW) = 2020
!               ENDIF

!               IJKS = SOUTH_OF(IJK)
!               IF(FLUID_AT(IJKS)) THEN
!                  FLAG_N(IJKS) = 2020
!               ENDIF

!               IJKB = BOTTOM_OF(IJK)
!               IF(FLUID_AT(IJKB)) THEN
!                  FLAG_T(IJKB) = 2020
!               ENDIF

               IF (BC_U_G(BCV) == UNDEFINED) THEN
                   IF (NO_I) THEN
                       BC_U_G(BCV) = ZERO
                   ELSEIF(BC_VOLFLOW_g(BCV)==UNDEFINED.AND. &
                          BC_MASSFLOW_g(BCV)==UNDEFINED.AND.&
                          BC_VELMAG_g(BCV)==UNDEFINED) THEN
                       IF(DMP_LOG)WRITE (UNIT_LOG, 900) 'BC_U_g', BCV
                       call mfix_exit(myPE)
                   ENDIF
               ENDIF
               IF (BC_V_G(BCV) == UNDEFINED) THEN
                   IF (NO_J) THEN
                       BC_V_G(BCV) = ZERO
                   ELSEIF(BC_VOLFLOW_g(BCV)==UNDEFINED.AND. &
                          BC_MASSFLOW_g(BCV)==UNDEFINED.AND.&
                          BC_VELMAG_g(BCV)==UNDEFINED) THEN
                       IF(DMP_LOG)WRITE (UNIT_LOG, 900) 'BC_V_g', BCV
                       call mfix_exit(myPE)
                   ENDIF
               ENDIF
               IF (BC_W_G(BCV) == UNDEFINED) THEN
                   IF (NO_K) THEN
                       BC_W_G(BCV) = ZERO
                   ELSEIF(BC_VOLFLOW_g(BCV)==UNDEFINED.AND. &
                          BC_MASSFLOW_g(BCV)==UNDEFINED.AND.&
                          BC_VELMAG_g(BCV)==UNDEFINED) THEN
                       IF(DMP_LOG)WRITE (UNIT_LOG, 900) 'BC_W_g', BCV
                       call mfix_exit(myPE)
                   ENDIF
               ENDIF
               IF (K_Epsilon .AND. BC_K_Turb_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_K_Turb_G', BCV
                   call mfix_exit(myPE)
               ENDIF
               IF (K_Epsilon .AND. BC_E_Turb_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_E_Turb_G', BCV
                   call mfix_exit(myPE)
               ENDIF

!               Check whether the bc velocity components have the correct sign
!               SELECT CASE (BC_PLANE(BCV))
!               CASE ('W')
!                   IF (BC_U_G(BCV) > ZERO) THEN
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '<'
!                       CALL MFIX_EXIT(myPE)
!                   ENDIF
!               CASE ('E')
!                   IF (BC_U_G(BCV) < ZERO) THEN
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '>'
!                       CALL MFIX_EXIT(myPE)
!                   ENDIF
!               CASE ('S')
!                   IF (BC_V_G(BCV) > ZERO) THEN
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '<'
!                       CALL MFIX_EXIT(myPE)
!                   ENDIF
!               CASE ('N')
!                   IF (BC_V_G(BCV) < ZERO) THEN
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '>'
!                       CALL MFIX_EXIT(myPE)
!                   ENDIF
!               CASE ('B')
!                   IF (BC_W_G(BCV) > ZERO) THEN
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '<'
!                       CALL MFIX_EXIT(myPE)
!                   ENDIF
!               CASE ('T')
!                   IF (BC_W_G(BCV) < ZERO) THEN
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '>'
!                       CALL MFIX_EXIT(myPE)
!                   ENDIF
!               END SELECT

               SUM_EP = BC_EP_G(BCV)
               DO M = 1, MMAX
                  IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN
                     IF (BC_EP_G(BCV) == ONE) THEN
                        BC_ROP_S(BCV,M) = ZERO
                     ELSEIF (MMAX == 1) THEN
                         BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S0(M)
                     ELSE
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_ROP_s', BCV, M
                         call mfix_exit(myPE)
                     ENDIF
                  ENDIF

                  SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S0(M)
                  IF (SPECIES_EQ(M)) THEN
                     SUM = ZERO
                        DO NN = 1, NMAX(M)
                           IF(BC_X_S(BCV,M,NN)/=UNDEFINED)SUM=SUM+BC_X_S(BCV,M,NN)
                        ENDDO

                     IF (BC_ROP_S(BCV,M)==ZERO .AND. SUM==ZERO) THEN
                        BC_X_S(BCV,M,1) = ONE
                        SUM = ONE
                     ENDIF

                     DO NN = 1, NMAX(M)
                        IF (BC_X_S(BCV,M,NN) == UNDEFINED) THEN
                           IF(.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)WRITE (UNIT_LOG,1110)BCV,M,NN
                              BC_X_S(BCV,M,NN) = ZERO
                        ENDIF
                     ENDDO

                     IF (.NOT.COMPARE(ONE,SUM)) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV, M
                           call mfix_exit(myPE)
                     ENDIF
                  ENDIF

                  IF (BC_U_S(BCV,M) == UNDEFINED) THEN
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_I) THEN
                        BC_U_S(BCV,M) = ZERO
                   ELSEIF(BC_VOLFLOW_s(BCV,M)==UNDEFINED.AND. &
                          BC_MASSFLOW_s(BCV,M)==UNDEFINED.AND.&
                          BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 910) 'BC_U_s', BCV, M
                            call mfix_exit(myPE)
                     ENDIF
                  ENDIF

                  IF (BC_V_S(BCV,M) == UNDEFINED) THEN
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_J) THEN
                        BC_V_S(BCV,M) = ZERO
                   ELSEIF(BC_VOLFLOW_s(BCV,M)==UNDEFINED.AND. &
                          BC_MASSFLOW_s(BCV,M)==UNDEFINED.AND.&
                          BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 910) 'BC_V_s', BCV, M
                            call mfix_exit(myPE)
                     ENDIF
                  ENDIF

                  IF (BC_W_S(BCV,M) == UNDEFINED) THEN
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_K) THEN
                        BC_W_S(BCV,M) = ZERO
                   ELSEIF(BC_VOLFLOW_s(BCV,M)==UNDEFINED.AND. &
                          BC_MASSFLOW_s(BCV,M)==UNDEFINED.AND.&
                          BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 910) 'BC_W_s', BCV, M
                           call mfix_exit(myPE)
                     ENDIF
                  ENDIF

                  IF (ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN
                     IF (BC_ROP_S(BCV,M) == ZERO) THEN
                        BC_T_S(BCV,M) = BC_T_G(BCV)
                     ELSE
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_T_s', BCV, M
                           call mfix_exit(myPE)
                     ENDIF
                  ENDIF

                  IF (GRANULAR_ENERGY .AND. BC_THETA_M(BCV,M)==UNDEFINED) THEN
                     IF (BC_ROP_S(BCV,M) == ZERO) THEN
                        BC_THETA_M(BCV,M) = ZERO
                     ELSE
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Theta_m', BCV, M
                          call mfix_exit(myPE)
                     ENDIF
                  ENDIF

!                   Check whether the bc velocity components have the correct sign
!                    SELECT CASE (TRIM(BC_PLANE(BCV)))
!                    CASE ('W')
!                        IF (BC_U_S(BCV,M) > ZERO) THEN
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '<'
!                            CALL MFIX_EXIT(myPE)
!                        ENDIF
!                    CASE ('E')
!                        IF (BC_U_S(BCV,M) < ZERO) THEN
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '>'
!                            CALL MFIX_EXIT(myPE)
!                        ENDIF
!                    CASE ('S')
!                        IF (BC_V_S(BCV,M) > ZERO) THEN
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '<'
!                            CALL MFIX_EXIT(myPE)
!                        ENDIF
!                    CASE ('N')
!                        IF (BC_V_S(BCV,M) < ZERO) THEN
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '>'
!                            CALL MFIX_EXIT(myPE)
!                        ENDIF
!                    CASE ('B')
!                        IF (BC_W_S(BCV,M) > ZERO) THEN
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '<'
!                            CALL MFIX_EXIT(myPE)
!                        ENDIF
!                    CASE ('T')
!                        IF (BC_W_S(BCV,M) < ZERO) THEN
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '>'
!                            CALL MFIX_EXIT(myPE)
!                        ENDIF
!                    END SELECT


               ENDDO

               IF (.NOT.COMPARE(ONE,SUM_EP)) THEN
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1125) BCV
                     call mfix_exit(myPE)
               ENDIF

               DO NN = 1, NScalar
                  IF (BC_Scalar(BCV,NN) == UNDEFINED) THEN
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Scalar', BCV, NN
                        CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDDO


            ELSEIF(BC_TYPE_ENUM(BCV)  == CG_PO) THEN

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN
                  FLAG_E(IJKW) = 2011
               ENDIF

               BCV_U = BC_U_ID(IJKW)
               IF(BCV_U>0) THEN
                  IF(BC_TYPE_ENUM(BCV_U)  == CG_PO) THEN
                    IJKWW = WEST_OF(IJKW)
                    IF(FLUID_AT(IJKWW)) THEN
                       FLAG_E(IJKWW) = 2011
                    ENDIF
                  ENDIF
               ENDIF

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  FLAG_N(IJKS) = 2011
               ENDIF

               BCV_V = BC_V_ID(IJKS)
               IF(BCV_V>0) THEN
                  IF(BC_TYPE_ENUM(BCV_V)  == CG_PO) THEN
                    IJKSS = SOUTH_OF(IJKS)
                    IF(FLUID_AT(IJKSS)) THEN
                       FLAG_N(IJKSS) = 2011
                    ENDIF
                  ENDIF
               ENDIF


               IF(DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  IF(FLUID_AT(IJKB)) THEN
                     FLAG_T(IJKB) = 2011
                  ENDIF

                  BCV_W = BC_W_ID(IJKB)
                  IF(BCV_W>0) THEN
                     IF(BC_TYPE_ENUM(BCV_W)  == CG_PO) THEN
                       IJKBB = BOTTOM_OF(IJKB)
                       IF(FLUID_AT(IJKBB)) THEN
                          FLAG_T(IJKBB) = 2011
                       ENDIF
                     ENDIF
                  ENDIF

               ENDIF

               IF (BC_P_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV
                   call mfix_exit(myPE)
               ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV)
                   call mfix_exit(myPE)
               ENDIF

            ENDIF

         ENDIF

      ENDDO

      RETURN


 900 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') not specified',/1X,'One of the following must be specified:',/1X,&
         'BC_VOLFLOW_g, BC_MASSFLOW_g or BC_VELMAG_g',/1X,70('*')/)

 910 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,'One of the following must be specified:',/1X,&
         'BC_VOLFLOW_g, BC_MASSFLOW_g or BC_VELMAG_g',/1X,70('*')/)

 1000 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/)
 1001 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: Illegal BC_TYPE for BC # = ',I2,/'   BC_TYPE = ',A,/&
         '  Valid BC_TYPE are: ')
 1002 FORMAT(5X,A16)
 1003 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/)
 1004 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ') not specified',/1X,70('*')/)
 1005 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ') value is unphysical',/1X,70('*')/)
 1010 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_P_g( ',I2,&
         ') = ',G12.5,/&
         ' Pressure should be greater than zero for compressible flow',/1X,70(&
         '*')/)
 1050 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - ',A,' should be ',A,' zero.',/1X,70('*')/)
 1060 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/)
 1065 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/)
 1100 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/)
 1103 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') value is unphysical',/1X,70('*')/)
 1104 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') not specified',/1X,70('*')/)
 1105 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') value is unphysical',/1X,70('*')/)
 1110 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/)
 1120 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/)
 1125 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/)
 1150 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - ',A,I1,' should be ',A,' zero.',/1X,70('*')/)
 1160 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: Boundary condition no', &
         I2,' is a second outflow condition.',/1X,&
         '  Only one outflow is allowed.  Consider using P_OUTFLOW.',/1X, 70('*')/)
 1200 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined BC location',/1X,70('*')/)
 1300 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined BC location',/1X,70('*')/)
 1400 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K')
 1410 FORMAT(I5,3X,I5,3X,I5)
 1420 FORMAT(/1X,70('*')/)

 1500 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K')


      END SUBROUTINE CHECK_BC_FLAGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALL BUILD_CONE_FOR_C2C                                C
!  Purpose: Define cone parameters for Cylider to Cylinder junction.   C
!           The C2C quadric ID must be between the two cylinders ID    C
!           (e.g., if Quadric 4 is a C2C, then Quadrics 3 and 5        C
!            must be cylinders). The two cylinders must be aligned     C
!           in the same direction and be clipped to define the extent  C
!           of the conical junction.                                   C
!           This method is currentl available for internal flow only.  C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE BUILD_CONE_FOR_C2C(Q)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: Q,QM1,QP1
      DOUBLE PRECISION :: x1,x2,y1,y2,z1,z2,R1,R2
      DOUBLE PRECISION :: tan_half_angle
      LOGICAL :: aligned
!-----------------------------------------------
!

      QM1 = Q-1
      QP1 = Q+1

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)' INFO FOR QUADRIC', Q
         WRITE(*,*)' Defining Cone for Cylinder to Cylinder junction'
         WRITE(*,*)' Between Quadrics ',QM1,' AND ', QP1
      ENDIF


      IF((TRIM(QUADRIC_FORM(QM1))=='X_CYL_INT').AND.  &
         (TRIM(QUADRIC_FORM(QP1))=='X_CYL_INT')) THEN       !Internal flow x-direction

         QUADRIC_FORM(Q) = 'X_CONE'

         aligned = (t_y(QM1)==t_y(QP1)).AND.(t_z(QM1)==t_z(QP1))
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         x1 = piece_xmax(QM1)
         x2 = piece_xmin(QP1)
         IF(x2<=x1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(x2-x1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = -ONE
         lambda_y(Q) = ONE/(tan_half_angle)**2
         lambda_z(Q) = ONE/(tan_half_angle)**2
         dquadric(Q) = ZERO

         piece_xmin(Q) = x1
         piece_xmax(Q) = x2

         t_x(Q) = x1 - R1/tan_half_angle
         t_y(Q) = t_y(QM1)
         t_z(Q) = t_z(QM1)

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF


      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='X_CYL_EXT').AND.  &
         (TRIM(QUADRIC_FORM(QP1))=='X_CYL_EXT')) THEN     !External flow x-direction

         QUADRIC_FORM(Q) = 'X_CONE'

         aligned = (t_y(QM1)==t_y(QP1)).AND.(t_z(QM1)==t_z(QP1))
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         x1 = piece_xmax(QM1)
         x2 = piece_xmin(QP1)
         IF(x2<=x1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(x2-x1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = ONE
         lambda_y(Q) = -ONE/(tan_half_angle)**2
         lambda_z(Q) = -ONE/(tan_half_angle)**2
         dquadric(Q) = ZERO

         piece_xmin(Q) = x1
         piece_xmax(Q) = x2

         t_x(Q) = x1 - R1/tan_half_angle
         t_y(Q) = t_y(QM1)
         t_z(Q) = t_z(QM1)

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF



      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='Y_CYL_INT').AND.  &
             (TRIM(QUADRIC_FORM(QP1))=='Y_CYL_INT')) THEN     !Internal flow y-direction

         QUADRIC_FORM(Q) = 'Y_CONE'

         aligned = (t_x(QM1)==t_x(QP1)).AND.(t_z(QM1)==t_z(QP1))
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         y1 = piece_ymax(QM1)
         y2 = piece_ymin(QP1)
         IF(y2<=y1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(y2-y1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = ONE/(tan_half_angle)**2
         lambda_y(Q) = -ONE
         lambda_z(Q) = ONE/(tan_half_angle)**2
         dquadric(Q) = ZERO

         piece_ymin(Q) = y1
         piece_ymax(Q) = y2

         t_x(Q) = t_x(QM1)
         t_y(Q) = y1 - R1/tan_half_angle
         t_z(Q) = t_z(QM1)

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF

      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='Y_CYL_EXT').AND.  &
             (TRIM(QUADRIC_FORM(QP1))=='Y_CYL_EXT')) THEN     !External flow y-direction

         QUADRIC_FORM(Q) = 'Y_CONE'

         aligned = (t_x(QM1)==t_x(QP1)).AND.(t_z(QM1)==t_z(QP1))
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         y1 = piece_ymax(QM1)
         y2 = piece_ymin(QP1)
         IF(y2<=y1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(y2-y1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = -ONE/(tan_half_angle)**2
         lambda_y(Q) = ONE
         lambda_z(Q) = -ONE/(tan_half_angle)**2
         dquadric(Q) = ZERO

         piece_ymin(Q) = y1
         piece_ymax(Q) = y2

         t_x(Q) = t_x(QM1)
         t_y(Q) = y1 - R1/tan_half_angle
         t_z(Q) = t_z(QM1)

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF



      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='Z_CYL_INT').AND.  &
             (TRIM(QUADRIC_FORM(QP1))=='Z_CYL_INT')) THEN     !Internal flow z-direction

         QUADRIC_FORM(Q) = 'Z_CONE'

         aligned = (t_x(QM1)==t_x(QP1)).AND.(t_y(QM1)==t_y(QP1))
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         z1 = piece_zmax(QM1)
         z2 = piece_zmin(QP1)
         IF(z2<=z1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(z2-z1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = ONE/(tan_half_angle)**2
         lambda_y(Q) = ONE/(tan_half_angle)**2
         lambda_z(Q) = -ONE
         dquadric(Q) = ZERO

         piece_zmin(Q) = z1
         piece_zmax(Q) = z2

         t_x(Q) = t_x(QM1)
         t_y(Q) = t_y(QM1)
         t_z(Q) = z1 - R1/tan_half_angle

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF

      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='Z_CYL_EXT').AND.  &
             (TRIM(QUADRIC_FORM(QP1))=='Z_CYL_EXT')) THEN     !External flow z-direction

         QUADRIC_FORM(Q) = 'Z_CONE'

         aligned = (t_x(QM1)==t_x(QP1)).AND.(t_y(QM1)==t_y(QP1))
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         z1 = piece_zmax(QM1)
         z2 = piece_zmin(QP1)
         IF(z2<=z1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(z2-z1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = -ONE/(tan_half_angle)**2
         lambda_y(Q) = -ONE/(tan_half_angle)**2
         lambda_z(Q) = ONE
         dquadric(Q) = ZERO

         piece_zmin(Q) = z1
         piece_zmax(Q) = z2

         t_x(Q) = t_x(QM1)
         t_y(Q) = t_y(QM1)
         t_z(Q) = z1 - R1/tan_half_angle

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF


      ELSE
         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' ERROR: C2C MUST BE DEFINED BETWEEN 2 CYLINDERS'
            WRITE(*,*) ' QUADRIC:',QM1, ' IS ',  TRIM(QUADRIC_FORM(QM1))
            WRITE(*,*) ' QUADRIC:',QP1, ' IS ',  TRIM(QUADRIC_FORM(QP1))
         ENDIF
         call mfix_exit(myPE)

      ENDIF

      RETURN
    END SUBROUTINE BUILD_CONE_FOR_C2C

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_DXYZ_FROM_CONTROL_POINTS                           C
!  Purpose: Define DX, DY, and DZ using control points                 C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_DXYZ_FROM_CONTROL_POINTS
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NN,NX,NY,NZ
      INTEGER :: I,I1,I2,J,J1,J2,K,K1,K2
      DOUBLE PRECISION :: L,CELL_RATIO

      LOGICAL,DIMENSION(MAX_CP) :: INDEPENDENT_SEGMENT

!-----------------------------------------------
!

!======================================================================
! X-DIRECTION
!======================================================================

! Step 1.  Input verification
!      1.1 Shift control points arrays such that the user only needs to enter
!          CPX(1) and above, and CPX(0) is automatically set to zero.

      DO NN = MAX_CP,1,-1
         CPX(nn) = CPX(NN-1)
      ENDDO

      CPX(0) = ZERO

!      1.2. Last control point must match domain length.

      NX = 0
      DO NN = 1,MAX_CP
         IF(CPX(nn)>ZERO) NX = NX + 1
      ENDDO

      IF(NX>0) THEN
         IF(MyPE==0)  WRITE(*,*)' INFO: DEFINING GRID SPACING IN X-DIRECTION... '
         IF(MyPE==0)  WRITE(*,*)' INFO: NUMBER OF CONTROL POINTS IN X-DIRECTION = ',NX
         IF(CPX(NX)/=XLENGTH) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: LAST CONTROL POINT MUST BE EQUAL TO XLENGTH.'
            IF(MyPE==0)  WRITE(*,*)' XLENGTH = ',XLENGTH
            IF(MyPE==0)  WRITE(*,*)' LAST CONTROL POINT = ',CPX(NX)
            call mfix_exit(myPE)
         ENDIF
      ENDIF

!      1.3. Check for acceptable values, and identify independent segments. If
!           the first or last cell dimension is given, it is converted into an
!           expansion ratio.

      INDEPENDENT_SEGMENT = .TRUE.

      DO NN = 1,NX   ! For each segment

         IF(CPX(nn) <= CPX(NN-1)) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: CONTROL POINTS ALONG X MUST BE SORTED IN ASCENDING ORDER.'
            IF(MyPE==0)  WRITE(*,*)' CPX = ',CPX(0:NX)
            call mfix_exit(myPE)
         ENDIF

         IF(NCX(nn) <= 1) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: NUMBER OF CELLS MUST BE LARGER THAN 1 IN X-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' NCX = ',NCX(nn)
            call mfix_exit(myPE)
         ENDIF

         IF(ERX(nn) <= ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: EXPANSION RATIO MUST BE POSITIVE IN X-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' ERX = ',ERX(nn)
            call mfix_exit(myPE)
         ENDIF

      ENDDO

      DO NN = 1,NX   ! For each segment

         IF(FIRST_DX(nn)/=ZERO.AND.LAST_DX(nn)/=ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST AND LAST DX ARE DEFINED, WHICH IS NOT ALLOWED IN X-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' FIRST DX = ',FIRST_DX(nn)
            IF(MyPE==0)  WRITE(*,*)' LAST  DX = ',LAST_DX(nn)
            call mfix_exit(myPE)
         ELSEIF(FIRST_DX(nn)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DX DEFINED IN X-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' FIRST DX = ',FIRST_DX(nn)
            L = CPX(nn) - CPX(NN-1)  ! Size of the current segment
            IF(L<=FIRST_DX(nn)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DX IS NOT SMALLER THAN SEGMENT LENGTH IN X-SEGMENT :',NN
               IF(MyPE==0)  WRITE(*,*)' FIRST DX = ',FIRST_DX(nn)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('FIRST',FIRST_DX(nn),L,NCX(nn),CELL_RATIO)
            ERX(nn) = CELL_RATIO**(NCX(nn)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERX(nn)
         ELSEIF(LAST_DX(nn)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: LAST DX DEFINED IN X-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' LAST DX = ',LAST_DX(nn)
            L = CPX(nn) - CPX(NN-1)  ! Size of the current segment
            IF(L<=LAST_DX(nn)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DX IS NOT SMALLER THAN SEGMENT LENGTH IN X-SEGMENT :',NN
               IF(MyPE==0)  WRITE(*,*)' LAST DX = ',LAST_DX(nn)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('LAST ',LAST_DX(nn),L,NCX(NN),CELL_RATIO)
            ERX(nn) = CELL_RATIO**(NCX(nn)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERX(nn)
         ELSEIF(FIRST_DX(nn)<ZERO) THEN
            IF(NN==1) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DX CANNOT MATCH PREVIOUS DX FOR FIRST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DX WILL ATTEMPT TO MATCH PREVIOUS DX FOR SEGMENT :',NN
               INDEPENDENT_SEGMENT(nn) = .FALSE.
            ENDIF
         ELSEIF(LAST_DX(nn)<ZERO) THEN
            IF(NN==NX) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DX CANNOT MATCH NEXT DX FOR LAST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: LAST DX WILL ATTEMPT TO MATCH NEXT DX FOR SEGMENT :',NN
               INDEPENDENT_SEGMENT(nn) = .FALSE.
            ENDIF
         ENDIF

      ENDDO

! Step 3.  Computation of cell sizes.

!      3.1 First pass: Set-up all independent segments


      I1 = 0  ! First index of segment
      I2 = 0  ! Last index of segment

      DO NN = 1,NX   ! For each segment

         I2 = I1 + NCX(nn) - 1

         IF(INDEPENDENT_SEGMENT(nn)) THEN

            L = CPX(nn) - CPX(NN-1)  ! Size of the current segment

            IF(ERX(nn)/=ONE) THEN
               CELL_RATIO = ERX(nn)**(ONE/DBLE(NCX(nn)-1))                     ! Ratio between two consecutive cells
               DX(I1) = L * (ONE - CELL_RATIO) / (ONE - CELL_RATIO**NCX(nn))     ! First cell size

               DO I = I1+1,I2                                                   ! All other cell sizes, geometric series
                 DX(I) = DX(I-1) * CELL_RATIO
               ENDDO

            ELSE
               DX(I1:I2) = L / NCX(nn)                                           ! Uniform size if expansion ratio is unity.
            ENDIF

         ENDIF

         I1 = I2 + 1                                                            ! Prepare First index for next segment

      ENDDO

!      3.2 Second pass: Set-up all dependent segments


      I1 = 0  ! First index of segment
      I2 = 0  ! Last index of segment

      DO NN = 1,NX   ! For each segment

         I2 = I1 + NCX(nn) - 1

         IF(.NOT.INDEPENDENT_SEGMENT(nn)) THEN

            L = CPX(nn) - CPX(NN-1)  ! Size of the current segment

            IF(FIRST_DX(nn)<ZERO) THEN
               DX(I1) = DX(I1-1)                                                ! First cell size
               CALL FIND_CELL_RATIO('FIRST',DX(I1),L,NCX(NN),CELL_RATIO)
               DO I = I1+1,I2                                                   ! All other cell sizes, geometric series
                 DX(I) = DX(I-1) * CELL_RATIO
               ENDDO
            ELSEIF(LAST_DX(nn)<ZERO) THEN
               DX(I2) = DX(I2+1)                                                ! Last cell size
               CALL FIND_CELL_RATIO('LAST ',DX(I2),L,NCX(nn),CELL_RATIO)
               DO I = I2-1,I1,-1                                                ! All other cell sizes, geometric series
                 DX(I) = DX(I+1) / CELL_RATIO
               ENDDO
            ENDIF

         ENDIF

         I1 = I2 + 1                                                  ! Prepare First index for next segment

      ENDDO


! Step 4. Verify that the sum of cells among all segment matches the total number of cells

      IF(I1>0.AND.I1/=IMAX) THEN
         IF(MyPE==0)  WRITE(*,*)' ERROR: IMAX MUST BE EQUAL TO THE SUM OF NCX.'
         IF(MyPE==0)  WRITE(*,*)' IMAX = ', IMAX
         IF(MyPE==0)  WRITE(*,*)' SUM OF NCX = ', I1
         call mfix_exit(myPE)
      ENDIF


!======================================================================
! Y-DIRECTION
!======================================================================

! Step 1.  Input verification
!      1.1 Shift control points arrays such that the user only needs to enter
!          CPY(1) and above, and CPY(0) is automatically set to zero.

      DO NN = MAX_CP,1,-1
         CPY(nn) = CPY(NN-1)
      ENDDO

      CPY(0) = ZERO

!      1.2. Last control point must match domain length.

      NY = 0
      DO NN = 1,MAX_CP
         IF(CPY(nn)>ZERO) NY = NY + 1
      ENDDO

      IF(NY>0) THEN
         IF(MyPE==0)  WRITE(*,*)' INFO: DEFINING GRID SPACING IN Y-DIRECTION... '
         IF(MyPE==0)  WRITE(*,*)' INFO: NUMBER OF CONTROL POINTS IN Y-DIRECTION = ',NY
         IF(CPY(NY)/=YLENGTH) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: LAST CONTROL POINT MUST BE EQUAL TO YLENGTH.'
            IF(MyPE==0)  WRITE(*,*)' YLENGTH = ',YLENGTH
            IF(MyPE==0)  WRITE(*,*)' LAST CONTROL POINT = ',CPY(NY)
            call mfix_exit(myPE)
         ENDIF
      ENDIF

!      1.3. Check for acceptable values, and identify independent segments. If
!           the first or last cell dimension is given, it is converted into an
!           expansion ratio.

      INDEPENDENT_SEGMENT = .TRUE.

      DO NN = 1,NY   ! For each segment

         IF(CPY(nn) <= CPY(NN-1)) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: CONTROL POINTS ALONG Y MUST BE SORTED IN ASCENDING ORDER.'
            IF(MyPE==0)  WRITE(*,*)' CPY = ',CPY(0:NY)
            call mfix_exit(myPE)
         ENDIF

         IF(NCY(nn) <= 1) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: NUMBER OF CELLS MUST BE LARGER THAN 1 IN Y-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' NCY = ',NCY(nn)
            call mfix_exit(myPE)
         ENDIF

         IF(ERY(nn) <= ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: EXPANSION RATIO MUST BE POSITIVE IN Y-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' ERY = ',ERY(nn)
            call mfix_exit(myPE)
         ENDIF

      ENDDO

      DO NN = 1,NY   ! For each segment

         IF(FIRST_DY(nn)/=ZERO.AND.LAST_DY(nn)/=ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST AND LAST DY ARE DEFINED, WHICH IS NOT ALLOWED IN Y-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' FIRST DY = ',FIRST_DY(nn)
            IF(MyPE==0)  WRITE(*,*)' LAST  DY = ',LAST_DY(nn)
            call mfix_exit(myPE)
         ELSEIF(FIRST_DY(nn)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DY DEFINED IN Y-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' FIRST DY = ',FIRST_DY(nn)
            L = CPY(nn) - CPY(NN-1)  ! Size of the current segment
            IF(L<=FIRST_DY(nn)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DY IS NOT SMALLER THAN SEGMENT LENGTH IN Y-SEGMENT :',NN
               IF(MyPE==0)  WRITE(*,*)' FIRST DY = ',FIRST_DY(nn)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('FIRST',FIRST_DY(nn),L,NCY(nn),CELL_RATIO)
            ERY(nn) = CELL_RATIO**(NCY(nn)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERY(nn)
         ELSEIF(LAST_DY(nn)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: LAST DY DEFINED IN Y-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' LAST DY = ',LAST_DY(nn)
            L = CPY(nn) - CPY(NN-1)  ! Size of the current segment
            IF(L<=LAST_DY(nn)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DY IS NOT SMALLER THAN SEGMENT LENGTH IN Y-SEGMENT :',NN
               IF(MyPE==0)  WRITE(*,*)' LAST DY = ',LAST_DY(nn)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('LAST ',LAST_DY(nn),L,NCY(nn),CELL_RATIO)
            ERY(nn) = CELL_RATIO**(NCY(nn)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERY(nn)
         ELSEIF(FIRST_DY(nn)<ZERO) THEN
            IF(NN==1) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DY CANNOT MATCH PREVIOUS DY FOR FIRST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DY WILL ATTEMPT TO MATCH PREVIOUS DY FOR SEGMENT :',NN
               INDEPENDENT_SEGMENT(nn) = .FALSE.
            ENDIF
         ELSEIF(LAST_DY(nn)<ZERO) THEN
            IF(NN==NY) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DY CANNOT MATCH NEXT DY FOR LAST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: LAST DY WILL ATTEMPT TO MATCH NEXT DY FOR SEGMENT :',NN
               INDEPENDENT_SEGMENT(nn) = .FALSE.
            ENDIF
         ENDIF

      ENDDO

! Step 3.  Computation of cell sizes.

!      3.1 First pass: Set-up all independent segments


      J1 = 0  ! First index of segment
      J2 = 0  ! Last index of segment

      DO NN = 1,NY   ! For each segment

         J2 = J1 + NCY(nn) - 1

         IF(INDEPENDENT_SEGMENT(nn)) THEN

            L = CPY(nn) - CPY(NN-1)  ! Size of the current segment

            IF(ERY(nn)/=ONE) THEN
               CELL_RATIO = ERY(nn)**(ONE/DBLE(NCY(nn)-1))                     ! Ratio between two consecutive cells
               DY(J1) = L * (ONE - CELL_RATIO) / (ONE - CELL_RATIO**NCY(nn))     ! First cell size

               DO J = J1+1,J2                                                   ! All other cell sizes, geometric series
                 DY(J) = DY(J-1) * CELL_RATIO
               ENDDO

            ELSE
               DY(J1:J2) = L / NCY(nn)                                           ! Uniform size if expansion ratio is unity.
            ENDIF

         ENDIF

         J1 = J2 + 1                                                            ! Prepare First index for next segment

      ENDDO

!      3.2 Second pass: Set-up all dependent segments


      J1 = 0  ! First index of segment
      J2 = 0  ! Last index of segment

      DO NN = 1,NY   ! For each segment

         J2 = J1 + NCY(nn) - 1

         IF(.NOT.INDEPENDENT_SEGMENT(nn)) THEN

            L = CPY(nn) - CPY(NN-1)  ! Size of the current segment

            IF(FIRST_DY(nn)<ZERO) THEN
               DY(J1) = DY(J1-1)                                                ! First cell size
               CALL FIND_CELL_RATIO('FIRST',DY(J1),L,NCY(nn),CELL_RATIO)
               DO J = J1+1,J2                                                   ! All other cell sizes, geometric series
                 DY(J) = DY(J-1) * CELL_RATIO
               ENDDO
            ELSEIF(LAST_DY(nn)<ZERO) THEN
               DY(J2) = DY(J2+1)                                                ! Last cell size
               CALL FIND_CELL_RATIO('LAST ',DY(J2),L,NCY(nn),CELL_RATIO)
               DO J = J2-1,J1,-1                                                ! All other cell sizes, geometric series
                 DY(J) = DY(J+1) / CELL_RATIO
               ENDDO
            ENDIF

         ENDIF

         J1 = J2 + 1                                                  ! Prepare First index for next segment

      ENDDO


! Step 4. Verify that the sum of cells among all segment matches the total number of cells

      IF(J1>0.AND.J1/=JMAX) THEN
         IF(MyPE==0)  WRITE(*,*)' ERROR: JMAX MUST BE EQUAL TO THE SUM OF NCY.'
         IF(MyPE==0)  WRITE(*,*)' JMAX = ', JMAX
         IF(MyPE==0)  WRITE(*,*)' SUM OF NCY = ', J1
         call mfix_exit(myPE)
      ENDIF


!======================================================================
! Z-DIRECTION
!======================================================================

      IF(NO_K) RETURN

! Step 1.  Input verification
!      1.1 Shift control points arrays such that the user only needs to enter
!          CPZ(1) and above, and CPZ(0) is automatically set to zero.

      DO NN = MAX_CP,1,-1
         CPZ(nn) = CPZ(NN-1)
      ENDDO

      CPZ(0) = ZERO

!      1.2. Last control point must match domain length.

      NZ = 0
      DO NN = 1,MAX_CP
         IF(CPZ(nn)>ZERO) NZ = NZ + 1
      ENDDO

      IF(NZ>0) THEN
         IF(MyPE==0)  WRITE(*,*)' INFO: DEFINING GRID SPACING IN Z-DIRECTION... '
         IF(MyPE==0)  WRITE(*,*)' INFO: NUMBER OF CONTROL POINTS IN Z-DIRECTION = ',NZ
         IF(CPZ(NZ)/=ZLENGTH) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: LAST CONTROL POINT MUST BE EQUAL TO ZLENGTH.'
            IF(MyPE==0)  WRITE(*,*)' ZLENGTH = ',ZLENGTH
            IF(MyPE==0)  WRITE(*,*)' LAST CONTROL POINT = ',CPZ(NZ)
            call mfix_exit(myPE)
         ENDIF
      ENDIF

!      1.3. Check for acceptable values, and identify independent segments. If
!           the first or last cell dimension is given, it is converted into an
!           expansion ratio.

      INDEPENDENT_SEGMENT = .TRUE.

      DO NN = 1,NZ   ! For each segment

         IF(CPZ(nn) <= CPZ(NN-1)) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: CONTROL POINTS ALONG Z MUST BE SORTED IN ASCENDING ORDER.'
            IF(MyPE==0)  WRITE(*,*)' CPZ = ',CPZ(0:NZ)
            call mfix_exit(myPE)
         ENDIF

         IF(NCZ(nn) <= 1) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: NUMBER OF CELLS MUST BE LARGER THAN 1 IN Z-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' NCZ = ',NCZ(nn)
            call mfix_exit(myPE)
         ENDIF

         IF(ERZ(nn) <= ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: EXPANSION RATIO MUST BE POSITIVE IN Z-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' ERZ = ',ERZ(nn)
            call mfix_exit(myPE)
         ENDIF

      ENDDO

      DO NN = 1,NZ   ! For each segment

         IF(FIRST_DZ(nn)/=ZERO.AND.LAST_DZ(nn)/=ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST AND LAST DZ ARE DEFINED, WHICH IS NOT ALLOWED IN Z-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' FIRST DZ = ',FIRST_DZ(nn)
            IF(MyPE==0)  WRITE(*,*)' LAST  DZ = ',LAST_DZ(nn)
            call mfix_exit(myPE)
         ELSEIF(FIRST_DZ(nn)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DZ DEFINED IN Z-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' FIRST DZ = ',FIRST_DZ(nn)
            L = CPZ(nn) - CPZ(NN-1)  ! Size of the current segment
            IF(L<=FIRST_DZ(nn)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DZ IS NOT SMALLER THAN SEGMENT LENGTH IN Z-SEGMENT :',NN
               IF(MyPE==0)  WRITE(*,*)' FIRST DZ = ',FIRST_DZ(nn)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('FIRST',FIRST_DZ(nn),L,NCZ(nn),CELL_RATIO)
            ERZ(nn) = CELL_RATIO**(NCZ(nn)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERZ(nn)
         ELSEIF(LAST_DZ(nn)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: LAST DZ DEFINED IN Z-SEGMENT :',NN
            IF(MyPE==0)  WRITE(*,*)' LAST DZ = ',LAST_DZ(nn)
            L = CPZ(nn) - CPZ(NN-1)  ! Size of the current segment
            IF(L<=LAST_DZ(nn)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DZ IS NOT SMALLER THAN SEGMENT LENGTH IN Z-SEGMENT :',NN
               IF(MyPE==0)  WRITE(*,*)' LAST DZ = ',LAST_DZ(nn)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('LAST ',LAST_DZ(nn),L,NCZ(nn),CELL_RATIO)
            ERZ(nn) = CELL_RATIO**(NCZ(nn)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERZ(nn)
         ELSEIF(FIRST_DZ(nn)<ZERO) THEN
            IF(NN==1) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DZ CANNOT MATCH PREVIOUS DZ FOR FIRST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DZ WILL ATTEMPT TO MATCH PREVIOUS DZ FOR SEGMENT :',NN
               INDEPENDENT_SEGMENT(nn) = .FALSE.
            ENDIF
         ELSEIF(LAST_DZ(nn)<ZERO) THEN
            IF(NN==NZ) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DZ CANNOT MATCH NEXT DZ FOR LAST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: LAST DZ WILL ATTEMPT TO MATCH NEXT DZ FOR SEGMENT :',NN
               INDEPENDENT_SEGMENT(nn) = .FALSE.
            ENDIF
         ENDIF

      ENDDO

! Step 3.  Computation of cell sizes.

!      3.1 First pass: Set-up all independent segments


      K1 = 0  ! First index of segment
      K2 = 0  ! Last index of segment

      DO NN = 1,NZ   ! For each segment

         K2 = K1 + NCZ(nn) - 1

         IF(INDEPENDENT_SEGMENT(nn)) THEN

            L = CPZ(nn) - CPZ(NN-1)  ! Size of the current segment

            IF(ERZ(nn)/=ONE) THEN
               CELL_RATIO = ERZ(nn)**(ONE/DBLE(NCZ(nn)-1))                     ! Ratio between two consecutive cells
               DZ(K1) = L * (ONE - CELL_RATIO) / (ONE - CELL_RATIO**NCZ(nn))     ! First cell size

               DO K = K1+1,K2                                                   ! All other cell sizes, geometric series
                 DZ(K) = DZ(K-1) * CELL_RATIO
               ENDDO

            ELSE
               DZ(K1:K2) = L / NCZ(nn)                                           ! Uniform size if expansion ratio is unity.
            ENDIF

         ENDIF

         K1 = K2 + 1                                                            ! Prepare First index for next segment

      ENDDO

!      3.2 Second pass: Set-up all dependent segments


      K1 = 0  ! First index of segment
      K2 = 0  ! Last index of segment

      DO NN = 1,NZ   ! For each segment

         K2 = K1 + NCZ(nn) - 1

         IF(.NOT.INDEPENDENT_SEGMENT(nn)) THEN

            L = CPZ(nn) - CPZ(NN-1)  ! Size of the current segment

            IF(FIRST_DZ(nn)<ZERO) THEN
               DZ(K1) = DZ(K1-1)                                                ! First cell size
               CALL FIND_CELL_RATIO('FIRST',DZ(K1),L,NCZ(nn),CELL_RATIO)
               DO K = K1+1,K2                                                   ! All other cell sizes, geometric series
                 DZ(K) = DZ(K-1) * CELL_RATIO
               ENDDO
            ELSEIF(LAST_DZ(nn)<ZERO) THEN
               DZ(K2) = DZ(K2+1)                                                ! Last cell size
               CALL FIND_CELL_RATIO('LAST ',DZ(K2),L,NCZ(nn),CELL_RATIO)
               DO K = K2-1,K1,-1                                                ! All other cell sizes, geometric series
                 DZ(K) = DZ(K+1) / CELL_RATIO
               ENDDO
            ENDIF

         ENDIF

         K1 = K2 + 1                                                  ! Prepare First index for next segment

      ENDDO


! Step 4. Verify that the sum of cells among all segment matches the total number of cells

      IF(K1>0.AND.K1/=KMAX) THEN
         IF(MyPE==0)  WRITE(*,*)' ERROR: KMAX MUST BE EQUAL TO THE SUM OF NCZ.'
         IF(MyPE==0)  WRITE(*,*)' KMAX = ', KMAX
         IF(MyPE==0)  WRITE(*,*)' SUM OF NCZ = ', K1
         call mfix_exit(myPE)
      ENDIF



      RETURN

      END SUBROUTINE GET_DXYZ_FROM_CONTROL_POINTS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FIND_CELL_RATIO                                        C
!  Purpose: Given the interval length L, number of cells N, and the    C
!           target value of D_target, find the cell ratio alpha3       C
!           such that D(POS) matches D_Target. POS can be either       C
!           FIRST or LAST cell in the segment.                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE FIND_CELL_RATIO(POS,D_Target,L,NN,ALPHA3)

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

      IMPLICIT NONE
      LOGICAL :: SOLUTION_FOUND

      DOUBLE PRECISION :: f1,f2,f3
      DOUBLE PRECISION :: ALPHA1,ALPHA2,ALPHA3,D_Target,L,DU
      DOUBLE PRECISION, PARAMETER :: ALPHAMAX = 10.0D0  ! maximum  value of cell ratio
      INTEGER :: NN,niter
      CHARACTER (LEN=5) :: POS

      DU = L / DBLE(NN)                  ! Cell size if uniform distribution

      IF(DU==D_TARGET) THEN
         ALPHA3 = 1.0
         SOLUTION_FOUND = .TRUE.
         RETURN
      ELSE

         IF(TRIM(POS)=='FIRST') THEN     ! Determine two initial guesses
            IF(D_TARGET<DU) THEN
               ALPHA1 = ONE
               ALPHA2 = ALPHAMAX
            ELSE
               ALPHA1 = ONE/ALPHAMAX
               ALPHA2 = ONE
            ENDIF
         ELSEIF(TRIM(POS)=='LAST') THEN
            IF(D_TARGET>DU) THEN
               ALPHA1 = ONE
               ALPHA2 = ALPHAMAX
            ELSE
               ALPHA1 = ONE/ALPHAMAX
               ALPHA2 = ONE
            ENDIF
         ELSE
            IF(MyPE==0) WRITE(*,*)' ERROR, IN FUNCTION F: POS MUST BE FIRST OR LAST.'
            call mfix_exit(myPE)
         ENDIF

      ENDIF

      f1 = F(POS,ALPHA1,D_Target,L,NN)
      f2 = F(POS,ALPHA2,D_Target,L,NN)

!======================================================================
!  The cell ratio is solution of F(alpha) = zero. The root is found by
!  the secant method, based on two inital guesses.
!======================================================================

      niter = 1
      SOLUTION_FOUND = .FALSE.

        if(DABS(f1)<TOL_F) then         ! First guess is solution
           SOLUTION_FOUND = .TRUE.
           ALPHA3 = ALPHA1
        elseif(DABS(f2)<TOL_F) then    ! Second guess is solution
           SOLUTION_FOUND = .TRUE.
           ALPHA3 = ALPHA2
        elseif(f1*f2 < ZERO) then       ! Solution is between two guesses
          niter = 0
          f3 = 2.0d0*TOL_F
          do while (   (abs(f3) > TOL_F)   .AND.   (niter<ITERMAX_INT)       )

            ALPHA3 = ALPHA1 - f1*(ALPHA2-ALPHA1)/(f2-f1)  ! secant point

            f3 = F(POS,ALPHA3,D_Target,L,NN)

            if(f1*f3<0) then            ! Reduce size of interval
              ALPHA2 = ALPHA3
              f2 = f3
            else
              ALPHA1 = ALPHA3
              f1 = f3
            endif
            niter = niter + 1

          end do
          if (niter < ITERMAX_INT) then
            SOLUTION_FOUND = .TRUE.
          else
             WRITE(*,*)   'Unable to find a solution'
             WRITE(*,1000)'between ALPHA1 = ', ALPHA1
             WRITE(*,1000)'   and  ALPHA2 = ', ALPHA2
             WRITE(*,1000)'Current value of ALPHA3 = ', ALPHA3
             WRITE(*,1000)'Current value of abs(f) = ', DABS(f3)
             WRITE(*,1000)'Tolerance = ', TOL_F
             WRITE(*,*)'Maximum number of iterations = ', ITERMAX_INT
             WRITE(*,*)   'Please increase the intersection tolerance, '
             WRITE(*,*)   'or the maximum number of iterations, and try again.'
             WRITE(*,*)   'MFiX will exit now.'
             CALL MFIX_EXIT(myPE)
             SOLUTION_FOUND = .FALSE.
          endif
        else
          WRITE(*,*)   'Unable to find a solution'
          WRITE(*,*)   'MFiX will exit now.'
          CALL MFIX_EXIT(myPE)
          SOLUTION_FOUND = .FALSE.
        endif

 1000 FORMAT(A,3(2X,G12.5))

      RETURN

    contains

      DOUBLE PRECISION Function F(POS,ALPHAC,D_Target,L,N)
        use param1, only: one
        USE constant
        USE mpi_utility

        IMPLICIT NONE
        DOUBLE PRECISION:: ALPHAC,D,D_Target,DU,L
        INTEGER:: N
        CHARACTER (LEN=5) :: POS

        DU = L / DBLE(nn)    ! Cell size if uniform distribution

        IF(ALPHAC==ONE) THEN
           D = DU
        ELSE
           IF(TRIM(POS)=='FIRST') THEN
              D = L * (ONE - ALPHAC) / (ONE -ALPHAC**N)
           ELSEIF(TRIM(POS)=='LAST') THEN
              D = L * (ONE - ALPHAC) / (ONE -ALPHAC**N) * ALPHAC**(N-1)
           ELSE
              IF(MyPE==0) WRITE(*,*)' ERROR, IN FUNCTION F: POS MUST BE FIRST OR LAST.'
              call mfix_exit(myPE)
           ENDIF
        ENDIF

        F = D - D_Target

        RETURN

      END FUNCTION F

    END SUBROUTINE FIND_CELL_RATIO

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_IJK_SIZE                                        C
!  Purpose: Adjust domain size of each processor, based on an          C
!           estimate on the total number of useful cells.              C
!           The objective is to reduce load imbalance.                 C
!           This is the parallel version                               C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ADJUST_IJK_SIZE
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE cut_cell_preproc, only: eval_f, print_cg_header
      USE cutcell
      USE dashboard
      USE discretelement
      USE funits
      USE gridmap
      USE indices
      USE leqsol
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE polygon
      USE quadric
      USE run
      USE scalars
      USE stl
      USE vtk

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LC,I,J,K,Q_ID,TOTAL_NUC,IDEAL_NCPP
      DOUBLE PRECISION :: X_COPY,Y_COPY,Z_COPY,F_COPY
      LOGICAL :: SHIFT,CLIP_FLAG,PRESENT
      DOUBLE PRECISION, DIMENSION(0:DIM_I) :: DXT ,XCC
      DOUBLE PRECISION, DIMENSION(0:DIM_J) :: DYT ,YCC
      DOUBLE PRECISION, DIMENSION(0:DIM_K) :: DZT, ZCC

      INTEGER, DIMENSION(0:DIM_I) :: NUC_I,GLOBAL_NUC_I
      INTEGER, DIMENSION(0:DIM_J) :: NUC_J,GLOBAL_NUC_J
      INTEGER, DIMENSION(0:DIM_K) :: NUC_K,GLOBAL_NUC_K

      INTEGER :: IPROC,PSUM

        INTEGER, DIMENSION(0:numPEs-1) :: NCPP_OLD,NCPP,NCPP_WITH_GHOST

        INTEGER, DIMENSION(0:NODESJ-1) :: JSIZE_OLD

      INTEGER :: JSIZE, JREMAIN
      INTEGER :: MAXVAL_NCPP_OLD,MINVAL_NCPP_OLD,MAXVAL_NCPP,MINVAL_NCPP
      INTEGER :: AVG_NCPP_OLD,AVG_NCPP
      DOUBLE PRECISION :: LIP_OLD,LIP,MAXSPEEDUP_OLD,P

      INTEGER :: I_OFFSET,J_OFFSET,K_OFFSET,IERR

      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount

      LOGICAL :: PRINT_STATISTICS


!-----------------------------------------------
!

      IF(.NOT.CARTESIAN_GRID) RETURN            ! Perform adjustement only when both CG
      IF(.NOT.ADJUST_PROC_DOMAIN_SIZE) RETURN   ! and domain adjustment
      IF(NODESI*NODESJ*NODESK==1) RETURN         ! and parallel run are active


      INQUIRE(FILE='gridmap.dat',EXIST=PRESENT)
      IF(PRESENT) THEN
         IF (myPE == PE_IO) WRITE(*,*)'gridmap was assigned from grimap.dat. Skipping the adjustment.'
         RETURN

      ENDIF

      IF (myPE == PE_IO) CALL PRINT_CG_HEADER

      allocate( NCPP_UNIFORM(0:NumPEs-1))

      SHORT_GRIDMAP_INIT = .TRUE.
      CALL gridmap_init
      SHORT_GRIDMAP_INIT = .FALSE.

      allocate( ISIZE_ALL(0:NODESI-1))
      allocate( JSIZE_ALL(0:NODESJ-1))
      allocate( KSIZE_ALL(0:NODESK-1))



      NCPP_UNIFORM = ijksize3_all

!      print*,'------> MyPE,NCPP_UNIFORM=',MyPE,NCPP_UNIFORM

      DIMENSION_I   = IMAX3
      DIMENSION_J   = JMAX3
      DIMENSION_K   = KMAX3

      PARTIAL_CHECK_03 = .TRUE.
!      CALL CHECK_DATA_03(SHIFT)
      PARTIAL_CHECK_03 = .FALSE.
      CALL CHECK_DATA_CARTESIAN                                ! Make sure CG data is valid

      CALL DEFINE_QUADRICS

      SHIFT = .TRUE.

      IF(SHIFT) THEN                                           ! Shift DX,DY,DZ and store it into temporary DXT,DYT,DZT

         IF (DO_I) THEN
            DXT(IMAX3) = DX(IMAX-1)
            DXT(IMAX2) = DX(IMAX-1)
            DO LC = IMAX1, IMIN1, -1
                 DXT(LC) = DX(LC-2)
            ENDDO
            DXT(IMIN2) = DX(IMIN1)
            DXT(IMIN3) =DX(IMIN2)

            XCC(IMIN1) = HALF*DXT(IMIN1)
            DO I=IMIN1+1,IMAX1
               XCC(I) = XCC(I-1) + HALF*(DXT(I-1) + DXT(I))
            ENDDO

         ENDIF
   !
         IF (DO_J) THEN
            DYT(JMAX3) = DY(JMAX-1)
            DYT(JMAX2) = DY(JMAX-1)
            DO LC = JMAX1, JMIN1, -1
               DYT(LC) = DY(LC-2)
            ENDDO
            DYT(JMIN2) = DY(JMIN1)
            DYT(JMIN3) =DY(JMIN2)

            YCC(JMIN1) = HALF*DYT(JMIN1)
            DO J=JMIN1+1,JMAX1
               YCC(J) = YCC(J-1) + HALF*(DYT(J-1) + DYT(J))
            ENDDO


         ENDIF
   !
         IF (DO_K) THEN
            DZT(KMAX3) = DZ(KMAX-1)
            DZT(KMAX2) = DZ(KMAX-1)
            DO LC = KMAX1, KMIN1, -1
               DZT(LC) = DZ(LC-2)
            ENDDO
            DZT(KMIN2) = DZ(KMIN1)
            DZT(KMIN3) =DZ(KMIN2)

            ZCC(KMIN1) = HALF*DZT(KMIN1)
            DO K=KMIN1+1,KMAX1
               ZCC(K) = ZCC(K-1) + HALF*(DZT(K-1) + DZT(K))
            ENDDO

         ENDIF

      ENDIF  ! SHIFT



      IF(NODESI>1) THEN                                ! DOMAIN DECOMPOSITION IN I-DIRECTION

         IF(myPE == 0.AND.NODESJ*NODESK/=1) THEN
            WRITE(*,*)'ERROR IN SUBROUTINE ADJUST_IJK_SIZE.'
            WRITE(*,*)'ADJUSTMENT POSSIBLE ONLY FOR DOMAIN DECOMPOSITION In ONE DIRECTION.'
            WRITE(*,*)'NODESI,NODESJ,NODESK =',NODESI,NODESJ,NODESK
            WRITE(*,*)'MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         ENDIF


         JSIZE_ALL(0:NODESJ-1) = jmax1-jmin1+1         ! Assign size in J and K-direction
         KSIZE_ALL(0:NODESK-1) = kmax1-kmin1+1

! Assign size in I-direction
         IF (myPE == 0) THEN
            WRITE(*,1000)'INFO FROM ADJUST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO ADJUST DOMAIN SIZE IN I-DIRECTION ...'
            WRITE(*,1000)'THIS IS BASED ON AN ESTIMATED (NOT EXACT) NUMBER OF USEFUL CELLS,'
            WRITE(*,1000)'AND IT INCLUDES GHOST LAYERS.'

         ENDIF

         DO I = ISTART1,IEND1

            NUC_I(I) = 0                     ! NUC : Number of Useful Cells

            DO J = JSTART1,JEND1
               DO K = KSTART1,KEND1
                  Q_ID = 1
                  CLIP_FLAG = .FALSE.

                  X_COPY = XCC(I)
                  Y_COPY = YCC(J)
                  Z_COPY = ZCC(K)


                  CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

   !               CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

   !               CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)

   !                CALL EVAL_STL_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)


                  IF (F_COPY < TOL_F ) THEN      ! Interior point, counted as useful
                     NUC_I(I) = NUC_I(I) + 1
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

! Gather NUC onto the head node

         CALL allgather_1i (IEND1-ISTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            I_OFFSET = 0
         ELSE
            I_OFFSET = 0
            DO iproc=0,myPE-1
               I_OFFSET = I_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (I_OFFSET,disp,IERR)

         call gatherv_1i( NUC_I(ISTART1:IEND1), IEND1-ISTART1+1, GLOBAL_NUC_I(IMIN1:IMAX1), rcount, disp, PE_IO, ierr )



      ELSEIF(NODESJ>1) THEN                            ! DOMAIN DECOMPOSITION IN J-DIRECTION


         IF(myPE == 0.AND.NODESI*NODESK/=1) THEN
            WRITE(*,*)'ERROR IN SUBROUTINE ADJUST_IJK_SIZE.'
            WRITE(*,*)'ADJUSTMENT POSSIBLE ONLY FOR DOMAIN DECOMPOSITION In ONE DIRECTION.'
            WRITE(*,*)'NODESI,NODESJ,NODESK =',NODESI,NODESJ,NODESK
            WRITE(*,*)'MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         ENDIF


         ISIZE_ALL(0:NODESI-1) = imax1-imin1+1         ! Assign size in I and K-direction
         KSIZE_ALL(0:NODESK-1) = kmax1-kmin1+1



!         WRITE(*,*)'Attempting to read gridmap from grimap.dat...',MyPE
!         INQUIRE(FILE='gridmap.dat',EXIST=PRESENT)
!         IF(PRESENT) THEN
!          WRITE(*,*)'Reading gridmap from grimap.dat...'
!            OPEN(CONVERT='BIG_ENDIAN',UNIT=777, FILE='gridmap.dat', STATUS='OLD')
!            DO IPROC = 0,NumPEs-1
!                  READ(777,*) jsize_all(IPROC)
!            ENDDO
!            CLOSE(777)
!            GOTO 9999
!         ENDIF


! Assign size in J-direction
         IF (myPE == 0) THEN
            WRITE(*,1000)'INFO FROM ADJUST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO ADJUST DOMAIN SIZE IN J-DIRECTION ...'
            WRITE(*,1000)'THIS IS BASED ON AN ESTIMATED (NOT EXACT) NUMBER OF USEFUL CELLS,'
            WRITE(*,1000)'AND IT INCLUDES GHOST LAYERS.'

         ENDIF

         DO J = JSTART1,JEND1

            NUC_J(J) = 0                     ! NUC : Number of Useful Cells

            DO I = ISTART1,IEND1
               DO K = KSTART1,KEND1
                  Q_ID = 1
                  CLIP_FLAG = .FALSE.

                  X_COPY = XCC(I)
                  Y_COPY = YCC(J)
                  Z_COPY = ZCC(K)


                  CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

   !               CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

   !               CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)

   !                CALL EVAL_STL_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)


                  IF (F_COPY < TOL_F ) THEN      ! Interior point, counted as useful
                     NUC_J(J) = NUC_J(J) + 1
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

! Gather NUC onto the head node

         CALL allgather_1i (JEND1-JSTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            J_OFFSET = 0
         ELSE
            J_OFFSET = 0
            DO iproc=0,myPE-1
               J_OFFSET = J_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (J_OFFSET,disp,IERR)

         call gatherv_1i( NUC_J(JSTART1:JEND1), JEND1-JSTART1+1, GLOBAL_NUC_J(JMIN1:JMAX1), rcount, disp, PE_IO, ierr )


      ELSEIF(NODESK>1) THEN                            ! DOMAIN DECOMPOSITION IN K-DIRECTION

         IF(myPE == 0.AND.NODESI*NODESJ/=1) THEN
            WRITE(*,*)'ERROR IN SUBROUTINE ADJUST_IJK_SIZE.'
            WRITE(*,*)'ADJUSTMENT POSSIBLE ONLY FOR DOMAIN DECOMPOSITION In ONE DIRECTION.'
            WRITE(*,*)'NODESI,NODESJ,NODESK =',NODESI,NODESJ,NODESK
            WRITE(*,*)'MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         ENDIF


         ISIZE_ALL(0:NODESI-1) = imax1-imin1+1         ! Assign size in I and J-direction
         JSIZE_ALL(0:NODESJ-1) = jmax1-jmin1+1

! Assign size in K-direction
         IF (myPE == 0) THEN
            WRITE(*,1000)'INFO FROM ADJUST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO ADJUST DOMAIN SIZE IN K-DIRECTION ...'
            WRITE(*,1000)'THIS IS BASED ON AN ESTIMATED (NOT EXACT) NUMBER OF USEFUL CELLS,'
            WRITE(*,1000)'AND IT INCLUDES GHOST LAYERS.'
         ENDIF

         DO K = KSTART1,KEND1

            NUC_K(K) = 0                     ! NUC : Number of Useful Cells

            DO I = ISTART1,IEND1
               DO J = JSTART1,JEND1
                  Q_ID = 1
                  CLIP_FLAG = .FALSE.

                  X_COPY = XCC(I)
                  Y_COPY = YCC(J)
                  Z_COPY = ZCC(K)


                  CALL EVAL_F('QUADRIC',X_COPY,Y_COPY,Z_COPY,Q_ID,F_COPY,CLIP_FLAG)

   !               CALL EVAL_F('POLYGON',X_COPY,Y_COPY,Z_COPY,N_POLYGON,F_COPY,CLIP_FLAG)

   !               CALL EVAL_F('USR_DEF',X_COPY,Y_COPY,Z_COPY,N_USR_DEF,F_COPY,CLIP_FLAG)

   !                CALL EVAL_STL_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)


                  IF (F_COPY < TOL_F ) THEN      ! Interior point, counted as useful
                     NUC_K(K) = NUC_K(K) + 1
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

! Gather NUC onto the head node

         CALL allgather_1i (KEND1-KSTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            K_OFFSET = 0
         ELSE
            K_OFFSET = 0
            DO iproc=0,myPE-1
               K_OFFSET = K_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (K_OFFSET,disp,IERR)

         call gatherv_1i( NUC_K(KSTART1:KEND1), KEND1-KSTART1+1, GLOBAL_NUC_K(KMIN1:KMAX1), rcount, disp, PE_IO, ierr )

      ENDIF !(NODESI,NODESJ, OR NODESK)

      IF (myPE == PE_IO) THEN                              ! DETERMINE BEST DOMAIN DECOMPOSITION FROM HEAD NODE

         IF(NODESI>1) THEN                                 ! DOMAIN DECOMPOSITION IN I-DIRECTION

         ELSEIF(NODESJ>1) THEN                            ! DOMAIN DECOMPOSITION IN J-DIRECTION

   ! For comparison with old grid size, determine the size in j direction (without adjustment) and add the remainder sequentially

            JSIZE = (JMAX1-JMIN1+1)/NODESJ
            JSIZE_OLD(0:NODESJ-1) = JSIZE

            JREMAIN = (JMAX1-JMIN1+1) - NODESJ*JSIZE
            IF (JREMAIN.GE.1) THEN
               JSIZE_OLD( 0:(JREMAIN-1) ) = JSIZE + 1
            ENDIF


!            DO IPROC = 0 ,numPEs-1
!               NCPP_UNIFORM(IPROC) = (imax3-imin3+1)*(JSIZE_OLD(IPROC)+4)
!               IF(DO_K) NCPP_UNIFORM(IPROC) = NCPP_UNIFORM(IPROC)*(kmax3-kmin3+1)
!            ENDDO



           JSIZE_ALL = JSIZE_OLD

           CALL MINIMIZE_LOAD_IMBALANCE(NODESJ,GLOBAL_NUC_J(JMIN1:JMAX1),JMIN1,JMAX1,JSIZE_ALL,NCPP,NCPP_WITH_GHOST)


            TOTAL_NUC  = SUM(NCPP_WITH_GHOST(0:NumPEs-1))
            IDEAL_NCPP = TOTAL_NUC / NumPEs

            WRITE (*, 1000) 'AFTER OPTIMIZATION:'
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/PROC.  = ',IDEAL_NCPP
            WRITE (*, 1000) 'ACTUALL CELL COUNT:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,numPEs-1
               WRITE (*, 1020) IPROC,JSIZE_ALL(IPROC),NCPP_WITH_GHOST(IPROC), &
                    DBLE(NCPP_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO

            DO IPROC = 0 ,numPEs-1
               NCPP_OLD(IPROC) = (imax3-imin3+1)*(JSIZE_ALL(IPROC)+4)
               IF(DO_K) NCPP_OLD(IPROC) = NCPP_OLD(IPROC)*(kmax3-kmin3+1)
            ENDDO


   ! Verify that the sum of all JSIZE_ALL matches JMAX

            PSUM = 0
            DO IPROC = 0,numPEs-1
               PSUM = PSUM + JSIZE_ALL(IPROC)
               IF(JSIZE_ALL(IPROC)<5) THEN
                  WRITE (*, 1010) 'ERROR: J-SIZE TOO SMALL FOR PROCESSOR:',IPROC
                  WRITE (*, 1010) 'J-SIZE = ',JSIZE_ALL(IPROC)
                  CALL MFIX_EXIT(myPE)
               ENDIF

            ENDDO

            IF(PSUM/=JMAX) THEN
               WRITE (*, 1000) 'ERROR IN ADJUST_IJK_SIZE: UNABLE TO ASSIGN JSIZE TO PROCESSORS.'
               WRITE (*, 1000) 'SUM OF JSIZE_ALL DOES NOT MATCH JMAX:'
               WRITE (*, 1010) 'SUM OF JSIZE_ALL = ',PSUM
               WRITE (*, 1010) 'JMAX1 = ',JMAX
               CALL MFIX_EXIT(myPE)
            ENDIF

         ELSEIF(NODESK>1) THEN                            ! DOMAIN DECOMPOSITION IN K-DIRECTION

         ENDIF !(NODESI,NODESJ, OR NODESK)

         PRINT_STATISTICS = .TRUE.

         IF(PRINT_STATISTICS) THEN

!            MAXVAL_NCPP_OLD = MAXVAL(NCPP_OLD)
!            MINVAL_NCPP_OLD = MINVAL(NCPP_OLD)
!            AVG_NCPP_OLD    = SUM(NCPP_OLD)/NUMPES

!            LIP_OLD = DBLE(MAXVAL_NCPP_OLD-AVG_NCPP_OLD)/DBLE(AVG_NCPP_OLD)*100.0D0

!            P = DBLE(MAXVAL_NCPP_OLD)/DBLE(AVG_NCPP_OLD)

!            MAXSPEEDUP_OLD = DBLE(NumPes)*(ONE-LIP_OLD/100.0D0)



            MAXVAL_NCPP_OLD = MAXVAL(NCPP_UNIFORM)
            MINVAL_NCPP_OLD = MINVAL(NCPP_UNIFORM)
            AVG_NCPP_OLD    = SUM(NCPP_UNIFORM)/NUMPES

            LIP_OLD = DBLE(MAXVAL_NCPP_OLD-AVG_NCPP_OLD)/DBLE(AVG_NCPP_OLD)*100.0D0

            P = DBLE(MAXVAL_NCPP_OLD)/DBLE(AVG_NCPP_OLD)

            MAXSPEEDUP_OLD = DBLE(NumPes)*(ONE-LIP_OLD/100.0D0)




            MAXVAL_NCPP = MAXVAL(NCPP_WITH_GHOST)
            MINVAL_NCPP = MINVAL(NCPP_WITH_GHOST)
            AVG_NCPP    = SUM(NCPP_WITH_GHOST(0:NumPEs-1))/NUMPES

            LIP = DBLE(MAXVAL_NCPP-AVG_NCPP)/DBLE(AVG_NCPP)*100.0D0

            P = DBLE(MAXVAL_NCPP)/DBLE(AVG_NCPP)

!            MAXSPEEDUP = DBLE(NumPes)*(ONE-LIP/100.0D0)


            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) 'ESTIMATED PARALLEL LOAD BALANCING STATISTICS'
            WRITE (*, 1000) 'COMPARISION BETWEEN UNIFORM SIZE (OLD)'
            WRITE (*, 1000) 'AND ADJUSTED SIZE (NEW)'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '                               OLD       NEW'
            WRITE (*, 1010) 'MAX CELL COUNT        : ',MAXVAL_NCPP_OLD,MAXVAL_NCPP
            WRITE (*, 1010) 'AT PROCESSOR          : ',MAXLOC(NCPP_OLD)-1,MAXLOC(NCPP_WITH_GHOST)-1
            WRITE (*, 1010) 'MIN CELL COUNT        : ',MINVAL_NCPP_OLD,MINVAL_NCPP
            WRITE (*, 1010) 'AT PROCESSOR          : ',MINLOC(NCPP_OLD)-1,MINLOC(NCPP_WITH_GHOST)-1
            WRITE (*, 1010) 'AVG CELL COUNT        : ',AVG_NCPP_OLD,AVG_NCPP
            WRITE (*, 1000) ''
            WRITE (*, 1030) 'LOAD IMBALANCE (%)    : ',LIP_OLD,LIP
            WRITE (*, 1000) ''
!            WRITE (*, 1030) 'IDEAL SPEEDUP         : ',DBLE(NumPEs),DBLE(NumPEs)
!            WRITE (*, 1030) 'MAX SPEEDUP           : ',MAXSPEEDUP_OLD,MAXSPEEDUP
!            WRITE (*, 1030) 'MAX EFFICIENCY (%)    : ',100.0D0 - LIP_OLD,100.0D0 - LIP

            WRITE (*, 1000) '================================================='

            WRITE (*, 1000) 'NOTE: ACTUAL LOAD BALANCING WILL BE COMPUTED AFTER PRE_PROCESSING.'


        ENDIF !(PRINT_STATISTICS)

     ENDIF  ! (MyPE==PE_IO)

9999  CONTINUE

! Broadcast Domain sizes to all processors

      CALL BCAST(ISIZE_ALL)
      CALL BCAST(JSIZE_ALL)
      CALL BCAST(KSIZE_ALL)

      DOMAIN_SIZE_ADJUSTED = .TRUE.


1000  FORMAT(1x,A)
1010  FORMAT(1x,A,I10,I10)
1020  FORMAT(1X,I8,2(I12),F12.2)
1030  FORMAT(1X,A,2(F10.1))
1040  FORMAT(F10.1)
1050  FORMAT(1X,3(A))


      RETURN

      END SUBROUTINE ADJUST_IJK_SIZE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REPORT_BEST_PROCESSOR_SIZE                             C
!  Purpose: Adjust domain size of each processor, based on an          C
!           estimate on the total number of useful cells.              C
!           The objective is to reduce load imbalance.                 C
!           This subroutine can be called in serial and parallel mode  C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE REPORT_BEST_PROCESSOR_SIZE
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE parallel

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!

      IS_SERIAL=(numPEs==1)

      IF(IS_SERIAL) THEN   ! Temporarily mimick a parallel run

         NODESI = NODESI_REPORT
         NODESJ = NODESJ_REPORT
         NODESK = NODESK_REPORT
         numPEs = NODESI*NODESJ*NODESK

         IF(numPEs>1.AND.myPE==0) THEN
            WRITE(*,1000)'TEMPORARILY SETTING:'
            WRITE(*,1010)'NODESI = ',NODESI
            WRITE(*,1010)'NODESJ = ',NODESJ
            WRITE(*,1010)'NODESK = ',NODESK
            WRITE(*,1000)'TO REPORT BEST DOMAIN SIZE FOR PARALLEL RUN'
         ENDIF


      ENDIF



      IF(numPEs>1) CALL REPORT_BEST_IJK_SIZE


      IF(IS_SERIAL) THEN   ! Revert to serial values
                            ! These values were changed to allow reporting
                            ! optimized sizes even for a serial run

         NODESI = 1
         NODESJ = 1
         NODESK = 1
         numPEs = 1

      ENDIF

1000  FORMAT(1x,A)
1010  FORMAT(1x,A,I10)

      RETURN

      END SUBROUTINE REPORT_BEST_PROCESSOR_SIZE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REPORT_BEST_IJK_SIZE                                   C
!  Purpose: Adjust domain size of each processor, based on an          C
!           estimate on the total number of useful cells.              C
!           The objective is to reduce load imbalance.                 C
!           This is the parallel version                               C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE REPORT_BEST_IJK_SIZE
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE parallel

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I,J,K,TOTAL_NUC,IDEAL_NCPP

      INTEGER, DIMENSION(0:DIM_I) :: NUC_I
      INTEGER, DIMENSION(0:DIM_J) :: NUC_J
      INTEGER, DIMENSION(0:DIM_K) :: NUC_K


      INTEGER :: ilistsize,jlistsize,klistsize             ! size of list of cells
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_NUC_I      ! Number of Useful Cells at I for all processors
                                                           ! (I will repeat if decomposing in J or K direction)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_LIST_I     ! List of I for all processors
      INTEGER, ALLOCATABLE, DIMENSION(:) :: GLOBAL_NUC_I   ! Number of Useful Cells at Global I

      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_NUC_J      ! Number of Useful Cells at J for all processors
                                                           ! (J will repeat if decomposing in I or K direction)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_LIST_J     ! List of J for all processors
      INTEGER, ALLOCATABLE, DIMENSION(:) :: GLOBAL_NUC_J   ! Number of Useful Cells at Global J

      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_NUC_K      ! Number of Useful Cells at K for all processors
                                                           ! (K will repeat if decomposing in I or J direction)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_LIST_K     ! List of K for all processors
      INTEGER, ALLOCATABLE, DIMENSION(:) :: GLOBAL_NUC_K   ! Number of Useful Cells at Global K


      INTEGER :: IPROC,PSUM

        INTEGER, DIMENSION(0:numPEs-1) :: NCPP_OLD,NCPP,NCPP_OLD_WITH_GHOST,NCPP_WITH_GHOST

        INTEGER, DIMENSION(0:NODESI-1) :: ISIZE_OLD
        INTEGER, DIMENSION(0:NODESJ-1) :: JSIZE_OLD
        INTEGER, DIMENSION(0:NODESK-1) :: KSIZE_OLD

      INTEGER :: JSIZE, IREMAIN,ISIZE, JREMAIN,KSIZE, KREMAIN
      DOUBLE PRECISION :: LIP_OLD

      INTEGER :: I_OFFSET,J_OFFSET,K_OFFSET,IERR

      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount

      INTEGER :: IPROC_OF_MAX_OLD,IPROC_OF_MIN_OLD

!-----------------------------------------------
!
      IF(.NOT.REPORT_BEST_DOMAIN_SIZE)   RETURN

      IF(.NOT.CARTESIAN_GRID) RETURN            ! Perform adjustement only when both CG
!      IF(.NOT.ADJUST_PROC_DOMAIN_SIZE) RETURN   ! and domain adjustment
      IF(NODESI*NODESJ*NODESK==1) RETURN         ! and parallel run are active

      IF(.not.allocated(ISIZE_ALL)) allocate( ISIZE_ALL(0:NODESI-1))
      IF(.not.allocated(JSIZE_ALL)) allocate( JSIZE_ALL(0:NODESJ-1))
      IF(.not.allocated(KSIZE_ALL)) allocate( KSIZE_ALL(0:NODESK-1))

      ISIZE_ALL(0:NODESI-1) = imax1-imin1+1  ! Assign default sizes in I, J and K-direction
      JSIZE_ALL(0:NODESJ-1) = jmax1-jmin1+1
      KSIZE_ALL(0:NODESK-1) = kmax1-kmin1+1


      IF(NODESI>1) THEN                                ! DOMAIN DECOMPOSITION IN I-DIRECTION

! Assign size in I-direction

         DO I = ISTART1,IEND1
            NUC_I(I) = 0                     ! NUC : Number of Useful Cells
            DO J = JSTART3,JEND3
               DO K = KSTART3,KEND3
                  IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN   ! Count number of useful cells
                     NUC_I(I) = NUC_I(I) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! Gather NUC onto the head node
! rcount is the size for each processor
! disp is the cumulative sum for each processor

         CALL allgather_1i (IEND1-ISTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            I_OFFSET = 0
         ELSE
            I_OFFSET = 0
            DO iproc=0,myPE-1
               I_OFFSET = I_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (I_OFFSET,disp,IERR)

         ilistsize=SUM(rcount)

         allocate( ALL_NUC_I(ilistsize))
         allocate( ALL_LIST_I(ilistsize))
         allocate( GLOBAL_NUC_I(IMIN1:IMAX1))

! Gather list of I and NUC, each processor has its own list
         call gatherv_1i( (/(I,I=ISTART1,IEND1)/), IEND1-ISTART1+1, ALL_LIST_I(:), rcount, disp, PE_IO, ierr )
         call gatherv_1i( NUC_I(ISTART1:IEND1), IEND1-ISTART1+1, ALL_NUC_I(:), rcount, disp, PE_IO, ierr )

! Get the glocal NUC for each unique value of I
         IF (myPE == 0) THEN
            GLOBAL_NUC_I = 0
            DO I=1,ilistsize
               GLOBAL_NUC_I(ALL_LIST_I(I)) = GLOBAL_NUC_I(ALL_LIST_I(I)) + ALL_NUC_I(I)
            ENDDO
         ENDIF

      ENDIF ! NODESI

      IF(NODESJ>1) THEN                            ! DOMAIN DECOMPOSITION IN J-DIRECTION

! Assign size in J-direction

         DO J = JSTART1,JEND1
            NUC_J(J) = 0                     ! NUC : Number of Useful Cells
            DO I = ISTART3,IEND3
               DO K = KSTART3,KEND3
                  IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN   ! Count number of useful cells
                     NUC_J(J) = NUC_J(J) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! Gather NUC onto the head node
! rcount is the size for each processor
! disp is the cumulative sum for each processor

         CALL allgather_1i (JEND1-JSTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            J_OFFSET = 0
         ELSE
            J_OFFSET = 0
            DO iproc=0,myPE-1
               J_OFFSET = J_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (J_OFFSET,disp,IERR)

         Jlistsize=SUM(rcount)

         allocate( ALL_NUC_J(Jlistsize))
         allocate( ALL_LIST_J(Jlistsize))
         allocate( GLOBAL_NUC_J(JMIN1:JMAX1))

! Gather list of J and NUC, each processor has its own list
         call gatherv_1i( (/(J,J=JSTART1,JEND1)/), JEND1-JSTART1+1, ALL_LIST_J(:), rcount, disp, PE_IO, ierr )
         call gatherv_1i( NUC_J(JSTART1:JEND1), JEND1-JSTART1+1, ALL_NUC_J(:), rcount, disp, PE_IO, ierr )

! Get the glocal NUC for each unique value of J
         IF (myPE == 0) THEN
            GLOBAL_NUC_J = 0
            DO J=1,Jlistsize
               GLOBAL_NUC_J(ALL_LIST_J(J)) = GLOBAL_NUC_J(ALL_LIST_J(J)) + ALL_NUC_J(J)
            ENDDO
         ENDIF

      ENDIF ! NODESJ

      IF(NODESK>1) THEN                            ! DOMAIN DECOMPOSITION IN K-DIRECTION

! Assign size in K-direction

         DO K = KSTART1,KEND1
            NUC_K(K) = 0                     ! NUC : Number of Useful Cells
            DO I = ISTART3,IEND3
               DO J = JSTART3,JEND3
                  IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN   ! Count number of useful cells
                     NUC_K(K) = NUC_K(K) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! Gather NUC onto the head node
! rcount is the size for each processor
! disp is the cumulative sum for each processor

         CALL allgather_1i (KEND1-KSTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            K_OFFSET = 0
         ELSE
            K_OFFSET = 0
            DO iproc=0,myPE-1
               K_OFFSET = K_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (K_OFFSET,disp,IERR)

         Klistsize=SUM(rcount)

         allocate( ALL_NUC_K(Klistsize))
         allocate( ALL_LIST_K(Klistsize))
         allocate( GLOBAL_NUC_K(KMIN1:KMAX1))

! Gather list of K and NUC, each processor has its own list
         call gatherv_1i( (/(K,K=KSTART1,KEND1)/), KEND1-KSTART1+1, ALL_LIST_K(:), rcount, disp, PE_IO, ierr )
         call gatherv_1i( NUC_K(KSTART1:KEND1), KEND1-KSTART1+1, ALL_NUC_K(:), rcount, disp, PE_IO, ierr )

! Get the glocal NUC for each unique value of K
         IF (myPE == 0) THEN
            GLOBAL_NUC_K = 0
            DO K=1,Klistsize
               GLOBAL_NUC_K(ALL_LIST_K(K)) = GLOBAL_NUC_K(ALL_LIST_K(K)) + ALL_NUC_K(K)
            ENDDO
         ENDIF

      ENDIF ! NODESK



! DETERMINE BEST DOMAIN DECOMPOSITION FROM HEAD NODE

      IF (myPE == PE_IO) THEN

         IF(NODESI>1) THEN      ! DOMAIN DECOMPOSITION IN I-DIRECTION

! For comparison with old grid size, determine the size in i direction (without adjustment) and add the remainder sequentially

            ISIZE = (IMAX1-IMIN1+1)/NODESI
            ISIZE_OLD(0:NODESI-1) = ISIZE

            IREMAIN = (IMAX1-IMIN1+1) - NODESI*ISIZE
            IF (IREMAIN.GE.1) THEN
               ISIZE_OLD( 0:(IREMAIN-1) ) = ISIZE + 1
            ENDIF

            ISIZE_ALL = ISIZE_OLD

! Get load imbalance before optimization
            CALL GET_LIP_WITH_GHOST_LAYERS(NODESI,GLOBAL_NUC_I(IMIN1:IMAX1),IMIN1, &
                                           IMAX1,ISIZE_ALL,NCPP_OLD,NCPP_OLD_WITH_GHOST, &
                                           LIP_OLD,IPROC_OF_MAX_OLD,IPROC_OF_MIN_OLD)

            TOTAL_NUC  = SUM(NCPP_OLD_WITH_GHOST(0:NODESI-1))
            IDEAL_NCPP = TOTAL_NUC / NODESI
            WRITE(*,1000)'INFO FROM REPORT_BEST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO REPORT BEST DOMAIN SIZE IN I-DIRECTION ...'
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/I-NODE = ',IDEAL_NCPP
            WRITE (*, 1000) 'BEFORE OPTIMIZATION:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   I-NODE       I-SIZE   CELLS/NODE    DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,NODESI-1
               WRITE (*, 1020) IPROC,ISIZE_ALL(IPROC),NCPP_OLD_WITH_GHOST(IPROC), &
                    DBLE(NCPP_OLD_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO
            WRITE (*, 1000) '================================================='

! Optimize load imbalance
            CALL MINIMIZE_LOAD_IMBALANCE(NODESI,GLOBAL_NUC_I(IMIN1:IMAX1),IMIN1,IMAX1,ISIZE_ALL,NCPP,NCPP_WITH_GHOST)

            TOTAL_NUC  = SUM(NCPP_WITH_GHOST(0:NODESI-1))
            IDEAL_NCPP = TOTAL_NUC / NODESI
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/I-NODE = ',IDEAL_NCPP
            WRITE (*, 1010) 'SINCE GHOST CELLS ARE INCLUDED, THE TOTALS'
            WRITE (*, 1010) 'BEFORE AND AFTER OPTIMIZATION MAY NOT MATCH.'
            WRITE (*, 1000) 'AFTER OPTIMIZATION:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   I-NODE       I-SIZE   CELLS/NODE    DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,NODESI-1
               WRITE (*, 1020) IPROC,ISIZE_ALL(IPROC),NCPP_WITH_GHOST(IPROC), &
                    DBLE(NCPP_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO
            WRITE (*, 1000) '================================================='

! Verify that the sum of all ISIZE_ALL matches IMAX
            PSUM = 0
            DO IPROC = 0,NODESI-1
               PSUM = PSUM + ISIZE_ALL(IPROC)
               IF(ISIZE_ALL(IPROC)<5) THEN
                  WRITE (*, 1010) 'ERROR: I-SIZE TOO SMALL FOR I-NODE:',IPROC
                  WRITE (*, 1010) 'I-SIZE = ',ISIZE_ALL(IPROC)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO

            IF(PSUM/=IMAX) THEN
               WRITE (*, 1000) 'ERROR IN ADJUST_IJK_SIZE: UNABLE TO ASSIGN ISIZE TO PROCESSORS.'
               WRITE (*, 1000) 'SUM OF ISIZE_ALL DOES NOT MATCH IMAX:'
               WRITE (*, 1010) 'SUM OF ISIZE_ALL = ',PSUM
               WRITE (*, 1010) 'IMAX = ',IMAX
               CALL MFIX_EXIT(myPE)
            ENDIF

         ENDIF                  ! DOMAIN DECOMPOSITION IN I-DIRECTION


         IF(NODESJ>1) THEN      ! DOMAIN DECOMPOSITION IN J-DIRECTION
! For comparison with old grid size, determine the size in i direction (without adjustment) and add the remainder sequentially

            JSIZE = (JMAX1-JMIN1+1)/NODESJ
            JSIZE_OLD(0:NODESJ-1) = JSIZE

            JREMAIN = (JMAX1-JMIN1+1) - NODESJ*JSIZE
            IF (JREMAIN.GE.1) THEN
               JSIZE_OLD( 0:(JREMAIN-1) ) = JSIZE + 1
            ENDIF

            JSIZE_ALL = JSIZE_OLD

! Get load imbalance before optimization
            CALL GET_LIP_WITH_GHOST_LAYERS(NODESJ,GLOBAL_NUC_J(JMIN1:JMAX1),JMIN1,JMAX1,JSIZE_ALL, &
                 NCPP_OLD,NCPP_OLD_WITH_GHOST,LIP_OLD,IPROC_OF_MAX_OLD,IPROC_OF_MIN_OLD)

            TOTAL_NUC  = SUM(NCPP_OLD_WITH_GHOST(0:NODESJ-1))
            IDEAL_NCPP = TOTAL_NUC / NODESJ
            WRITE(*,1000)'INFO FROM REPORT_BEST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO REPORT BEST DOMAIN SIZE IN J-DIRECTION ...'
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/J-NODE = ',IDEAL_NCPP
            WRITE (*, 1000) 'BEFORE OPTIMIZATION:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   J-NODE       J-SIZE   CELLS/NODE    DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,NODESJ-1
               WRITE (*, 1020) IPROC,JSIZE_ALL(IPROC),NCPP_OLD_WITH_GHOST(IPROC), &
                    DBLE(NCPP_OLD_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO
            WRITE (*, 1000) '================================================='

! Optimize load imbalance
            CALL MINIMIZE_LOAD_IMBALANCE(NODESJ,GLOBAL_NUC_J(JMIN1:JMAX1),JMIN1,JMAX1,JSIZE_ALL,NCPP,NCPP_WITH_GHOST)

            TOTAL_NUC  = SUM(NCPP_WITH_GHOST(0:NODESJ-1))
            IDEAL_NCPP = TOTAL_NUC / NODESJ
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/J-NODE = ',IDEAL_NCPP
            WRITE (*, 1010) 'SINCE GHOST CELLS ARE INCLUDED, THE TOTALS'
            WRITE (*, 1010) 'BEFORE AND AFTER OPTIMIZATION MAY NOT MATCH.'
            WRITE (*, 1000) 'AFTER OPTIMIZATION:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   J-NODE       J-SIZE   CELLS/NODE    DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,NODESJ-1
               WRITE (*, 1020) IPROC,JSIZE_ALL(IPROC),NCPP_WITH_GHOST(IPROC), &
                    DBLE(NCPP_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO
            WRITE (*, 1000) '================================================='

! Verify that the sum of all JSIZE_ALL matches JMAX
            PSUM = 0
            DO IPROC = 0,NODESJ-1
               PSUM = PSUM + JSIZE_ALL(IPROC)
               IF(JSIZE_ALL(IPROC)<5) THEN
                  WRITE (*, 1010) 'ERROR: J-SIZE TOO SMALL FOR J-NODE:',IPROC
                  WRITE (*, 1010) 'J-SIZE = ',JSIZE_ALL(IPROC)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO

            IF(PSUM/=JMAX) THEN
               WRITE (*, 1000) 'ERROR IN ADJUST_IJK_SIZE: UNABLE TO ASSIGN JSIZE TO PROCESSORS.'
               WRITE (*, 1000) 'SUM OF JSIZE_ALL DOES NOT MATCH JMAX:'
               WRITE (*, 1010) 'SUM OF JSIZE_ALL = ',PSUM
               WRITE (*, 1010) 'JMAX = ',JMAX
               CALL MFIX_EXIT(myPE)
            ENDIF

         ENDIF                  ! DOMAIN DECOMPOSITION IN J-DIRECTION

         IF(NODESK>1) THEN      ! DOMAIN DECOMPOSITION IN K-DIRECTION
! For comparison with old grid size, determine the size in i direction (without adjustment) and add the remainder sequentially

            KSIZE = (KMAX1-KMIN1+1)/NODESK
            KSIZE_OLD(0:NODESK-1) = KSIZE

            KREMAIN = (KMAX1-KMIN1+1) - NODESK*KSIZE
            IF (KREMAIN.GE.1) THEN
               KSIZE_OLD( 0:(KREMAIN-1) ) = KSIZE + 1
            ENDIF

            KSIZE_ALL = KSIZE_OLD

! Get load imbalance before optimization
            CALL GET_LIP_WITH_GHOST_LAYERS(NODESK,GLOBAL_NUC_K(KMIN1:KMAX1),KMIN1,KMAX1,KSIZE_ALL, &
                 NCPP_OLD,NCPP_OLD_WITH_GHOST,LIP_OLD,IPROC_OF_MAX_OLD,IPROC_OF_MIN_OLD)

            TOTAL_NUC  = SUM(NCPP_OLD_WITH_GHOST(0:NODESK-1))
            IDEAL_NCPP = TOTAL_NUC / NODESK
            WRITE(*,1000)'INFO FROM REPORT_BEST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO REPORT BEST DOMAIN SIZE IN K-DIRECTION ...'
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/K_NODE = ',IDEAL_NCPP
            WRITE (*, 1000) 'BEFORE OPTIMIZATION:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   K-NODE       K-SIZE   CELLS/NODE    DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,NODESK-1
               WRITE (*, 1020) IPROC,KSIZE_ALL(IPROC),NCPP_OLD_WITH_GHOST(IPROC), &
                    DBLE(NCPP_OLD_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO
            WRITE (*, 1000) '================================================='

! Optimize load imbalance
            CALL MINIMIZE_LOAD_IMBALANCE(NODESK,GLOBAL_NUC_K(KMIN1:KMAX1),KMIN1,KMAX1,KSIZE_ALL,NCPP,NCPP_WITH_GHOST)

            TOTAL_NUC  = SUM(NCPP_WITH_GHOST(0:NODESK-1))
            IDEAL_NCPP = TOTAL_NUC / NODESK
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/K_NODE = ',IDEAL_NCPP
            WRITE (*, 1010) 'SINCE GHOST CELLS ARE INCLUDED, THE TOTALS'
            WRITE (*, 1010) 'BEFORE AND AFTER OPTIMIZATION MAY NOT MATCH.'
            WRITE (*, 1000) 'AFTER OPTIMIZATION:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   K-NODE       K-SIZE   CELLS/NODE    DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,NODESK-1
               WRITE (*, 1020) IPROC,KSIZE_ALL(IPROC),NCPP_WITH_GHOST(IPROC), &
                    DBLE(NCPP_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO
            WRITE (*, 1000) '================================================='

! Verify that the sum of all KSIZE_ALL matches KMAX
            PSUM = 0
            DO IPROC = 0,NODESK-1
               PSUM = PSUM + KSIZE_ALL(IPROC)
               IF(KSIZE_ALL(IPROC)<5) THEN
                  WRITE (*, 1010) 'ERROR: K-SIZE TOO SMALL FOR K-NODE:',IPROC
                  WRITE (*, 1010) 'K-SIZE = ',KSIZE_ALL(IPROC)
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO

            IF(PSUM/=KMAX) THEN
               WRITE (*, 1000) 'ERROR IN ADJUST_IJK_SIZE: UNABLE TO ASSIGN KSIZE TO PROCESSORS.'
               WRITE (*, 1000) 'SUM OF KSIZE_ALL DOES NOT MATCH KMAX:'
               WRITE (*, 1010) 'SUM OF KSIZE_ALL = ',PSUM
               WRITE (*, 1010) 'KMAX = ',KMAX
               CALL MFIX_EXIT(myPE)
            ENDIF


         ENDIF                  ! DOMAIN DECOMPOSITION IN K-DIRECTION


         OPEN(CONVERT='BIG_ENDIAN',UNIT=777, FILE='suggested_gridmap.dat')
         WRITE (777, 1005) NODESI,NODESJ,NODESK, '     ! NODESI, NODESJ, NODESK'
         DO IPROC = 0,NODESI-1
               WRITE(777,1060) IPROC,Isize_all(IPROC)
         ENDDO
         DO IPROC = 0,NODESJ-1
               WRITE(777,1060) IPROC,Jsize_all(IPROC)
         ENDDO
         DO IPROC = 0,NODESK-1
               WRITE(777,1060) IPROC,Ksize_all(IPROC)
         ENDDO

         CLOSE(777)
         WRITE (*, 1000) '================================================='
         WRITE (*, 1000) 'GRID PARTITION SAVED IN FILE: suggested_gridmap.dat'
         WRITE (*, 1000) 'TO USE THIS DISTRIBUTION, RENAME THE FILE AS: gridmap.dat'
         WRITE (*, 1000) 'AND RUN MFIX AGAIN.'
!         WRITE (*, 1000) 'MFIX WILL STOP NOW.'
         WRITE (*, 1000) '================================================='


      ENDIF  ! (MyPE==PE_IO)

! Finalize and terminate MPI
      call parallel_fin

      STOP

! Broadcast Domain sizes to all processors

!      CALL BCAST(ISIZE_ALL)
!      CALL BCAST(JSIZE_ALL)
!      CALL BCAST(KSIZE_ALL)

!      DOMAIN_SIZE_ADJUSTED = .TRUE.

1000  FORMAT(1x,A)
1005  FORMAT(1x,I10,I10,I10,A)
1010  FORMAT(1x,A,I10,I10)
1020  FORMAT(1X,I8,2(I12),F12.2)
1030  FORMAT(1X,A,2(F10.1))
1040  FORMAT(F10.1)
1050  FORMAT(1X,3(A))
1060  FORMAT(1x,I10,I10)

      RETURN

      END SUBROUTINE REPORT_BEST_IJK_SIZE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_LIP_WITH_GHOST_LAYERS                              C
!  Purpose: Compute Load Imbalance percentage                          C
!           by including size of ghost layers                          C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_LIP_WITH_GHOST_LAYERS(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST,LIP,IPROC_OF_MAX,IPROC_OF_MIN)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NODESL,L,LMIN1,LMAX1,TOTAL_NUC,TOTAL_NUC_WITH_GHOST,IPROC_OF_MAX,IPROC_OF_MIN

      INTEGER :: LCOUNT1,LCOUNT2,MINVAL_NCPP,MAXVAL_NCPP,IDEAL_NCPP
      INTEGER, DIMENSION(LMIN1:LMAX1) :: NUC_L

      INTEGER :: IPROC

        INTEGER, DIMENSION(0:NODESL-1) :: NCPP,NCPP_WITH_GHOST,L_SIZE,L1,L2


      DOUBLE PRECISION :: LIP

!-----------------------------------------------

      LCOUNT1 = LMAX1 - LMIN1 + 1
      LCOUNT2 = SUM(L_SIZE(0:NODESL-1))

      IF(LCOUNT1/=LCOUNT2) THEN
         WRITE(*,*)' ERROR: SUM OF CELLS DO NOT MATCH:',LCOUNT1,LCOUNT2
         CALL MFIX_EXIT(myPE)
      ENDIF

      L1(0) = LMIN1
      L2(0) = L1(0) + L_SIZE(0) - 1

      DO IPROC = 1,NODESL-1
         L1(IPROC) = L2(IPROC-1) + 1
         L2(IPROC) = L1(IPROC) + L_SIZE(IPROC) - 1
      ENDDO

      DO IPROC = 0,NODESL-1
         NCPP(IPROC) = SUM(NUC_L(L1(IPROC):L2(IPROC)))
!         print*,'NUC=',NUC_L(L1(IPROC):L2(IPROC))
!         print*,'L1,L2=',IPROC,L1(IPROC),L2(IPROC),NCPP(IPROC)
      ENDDO

      TOTAL_NUC = 0

      DO L=LMIN1,LMAX1
         TOTAL_NUC = TOTAL_NUC + NUC_L(L)
      ENDDO

      NCPP_WITH_GHOST(0) = NCPP(0) + 2*NUC_L(L1(0)) + NUC_L(L1(1)) + NUC_L(L1(1)+1)

      DO IPROC = 1,NODESL-2
         NCPP_WITH_GHOST(IPROC) =   NCPP(IPROC)  &
                                  + NUC_L(L2(IPROC-1)) + NUC_L(L2(IPROC-1)-1) &
                                  + NUC_L(L1(IPROC+1)) + NUC_L(L1(IPROC+1)+1)
      ENDDO

      NCPP_WITH_GHOST(NODESL-1) = NCPP(NODESL-1) + 2*NUC_L(L2(NODESL-1)) + NUC_L(L2(NODESL-2)) + NUC_L(L2(NODESL-2)-1)

      TOTAL_NUC_WITH_GHOST = 0
      DO IPROC = 0,NODESL-1
!         print*,'NCPP_WITH_GHOST=',IPROC,L_SIZE(IPROC),NCPP(IPROC),NCPP_WITH_GHOST(IPROC)
         TOTAL_NUC_WITH_GHOST = TOTAL_NUC_WITH_GHOST + NCPP_WITH_GHOST(IPROC)
      ENDDO


      IDEAL_NCPP = TOTAL_NUC_WITH_GHOST / NumPEs

      MAXVAL_NCPP = MAXVAL(NCPP_WITH_GHOST)
      MINVAL_NCPP = MINVAL(NCPP_WITH_GHOST)

      LIP = DBLE(MAXVAL_NCPP-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
!      LIP = DBLE(MAXVAL_NCPP-MINVAL_NCPP)/DBLE(MINVAL_NCPP)*100.0D0


      IPROC_OF_MAX = MAXLOC(NCPP_WITH_GHOST,1)-1
      IPROC_OF_MIN = MINLOC(NCPP_WITH_GHOST,1)-1


!      print*,'IPROC_OF_MIN=',IPROC_OF_MIN,MINVAL_NCPP
!      print*,'IPROC_OF_MAX=',IPROC_OF_MAX,MAXVAL_NCPP

      RETURN
      END SUBROUTINE GET_LIP_WITH_GHOST_LAYERS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MINIMIZE_LOAD_IMBALANCE                                C
!  Purpose: Rearrange L_SIZE to minimize load imbalance                C
!           by including size of ghost layers                          C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MINIMIZE_LOAD_IMBALANCE(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NODESL,LMIN1,LMAX1,IPROC_OF_MAX,IPROC_OF_MIN

      INTEGER, DIMENSION(LMIN1:LMAX1) :: NUC_L

      INTEGER :: NN,NOIMPROVEMENT

      INTEGER,PARAMETER :: NAMAX=10000  ! maximum number of adjustments, increase if optimized load is not reached

      INTEGER, DIMENSION(0:numPEs-1) :: NCPP,NCPP_WITH_GHOST,L_SIZE,BEST_L_SIZE,BEST_NCPP,BEST_NCPP_WITH_GHOST

      DOUBLE PRECISION :: LIP,BEST_LIP

!-----------------------------------------------


!      NA = NAMAX   ! Number of adjustments

! Initial estimate of LIP

      CALL GET_LIP_WITH_GHOST_LAYERS(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST,LIP,IPROC_OF_MAX,IPROC_OF_MIN)

      BEST_LIP = LIP
      BEST_L_SIZE = L_SIZE

!      print*,'INITIAL ESTIMATE OF LIP:',LIP
!         WRITE (*, 1000) '================================================='
!         WRITE (*, 1010) 'AFTER STEP:',N
!         WRITE (*, 1010) 'NOIMPROVEMENT=',NOIMPROVEMENT
!         WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!         WRITE (*, 1000) '================================================='
!         DO IPROC = 0,numPEs-1
!            WRITE (*, 1020) IPROC,L_SIZE(IPROC),NCPP_WITH_GHOST(IPROC)
!         ENDDO
!         WRITE (*, 1000) '================================================='

!      print*,'MIN ',IPROC_OF_MIN,NCPP_WITH_GHOST(IPROC_OF_MIN)
!      print*,'MAX ',IPROC_OF_MAX,NCPP_WITH_GHOST(IPROC_OF_MAX)

      print*,'ATTEMPTING TO OPTIMIZE LOAD BALANCE...'

      NOIMPROVEMENT=0

      DO NN = 1,NAMAX

         L_SIZE(IPROC_OF_MAX) = L_SIZE(IPROC_OF_MAX) - 1
         L_SIZE(IPROC_OF_MIN) = L_SIZE(IPROC_OF_MIN) + 1

         CALL GET_LIP_WITH_GHOST_LAYERS(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST,LIP,IPROC_OF_MAX,IPROC_OF_MIN)

!      print*,'After adjustment of LIP:',N, LIP
!      print*,'MIN ',IPROC_OF_MIN,NCPP_WITH_GHOST(IPROC_OF_MIN)
!      print*,'MAX ',IPROC_OF_MAX,NCPP_WITH_GHOST(IPROC_OF_MAX)

         IF(LIP<BEST_LIP) THEN
            BEST_LIP    = LIP
            BEST_L_SIZE = L_SIZE
            BEST_NCPP   = NCPP
            BEST_NCPP_WITH_GHOST   = NCPP_WITH_GHOST
            NOIMPROVEMENT=0
         ELSE
            NOIMPROVEMENT = NOIMPROVEMENT + 1
         ENDIF

!         WRITE (*, 1000) '================================================='
!         WRITE (*, 1010) 'AFTER STEP:',N
!         WRITE (*, 1010) 'NOIMPROVEMENT=',NOIMPROVEMENT
!         WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!         WRITE (*, 1000) '================================================='
!         DO IPROC = 0,numPEs-1
!            WRITE (*, 1020) IPROC,L_SIZE(IPROC),NCPP_WITH_GHOST(IPROC)
!         ENDDO

         IF(NOIMPROVEMENT==10) THEN
            WRITE (*, 1000) 'OPTIMIZED LOAD BALANCE REACHED.'
            EXIT
         ENDIF

      ENDDO

!      print*,'Best LIP = ',BEST_LIP
      L_SIZE = BEST_L_SIZE
      NCPP   = BEST_NCPP
      NCPP_WITH_GHOST = BEST_NCPP_WITH_GHOST


!      WRITE (*, 1000) '================================================='
!      WRITE (*, 1000) 'AFTER OPTIMIZATION'
!      WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!      WRITE (*, 1000) '================================================='
!      DO IPROC = 0,numPEs-1
!         WRITE (*, 1020) IPROC,L_SIZE(IPROC),NCPP_WITH_GHOST(IPROC)
!      ENDDO

1000  FORMAT(1x,A)
1010  FORMAT(1x,A,I10,I10)
1020  FORMAT(1X,I8,2(I12),F12.2)

      RETURN
      END SUBROUTINE MINIMIZE_LOAD_IMBALANCE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REPORT_BEST_IJK_SIZE0                                   C
!  Purpose: Adjust domain size of each processor, based on an          C
!           estimate on the total number of useful cells.              C
!           The objective is to reduce load imbalance.                 C
!           This is the parallel version                               C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE REPORT_BEST_IJK_SIZE0
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE parallel

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I,J,K,TOTAL_NUC,IDEAL_NCPP

      INTEGER, DIMENSION(0:DIM_I) :: NUC_I,GLOBAL_NUC_I
      INTEGER, DIMENSION(0:DIM_J) :: NUC_J,GLOBAL_NUC_J
      INTEGER, DIMENSION(0:DIM_K) :: NUC_K,GLOBAL_NUC_K

      INTEGER :: IPROC,PSUM

      INTEGER, DIMENSION(0:numPEs-1) :: NCPP_OLD,NCPP,NCPP_OLD_WITH_GHOST,NCPP_WITH_GHOST

      INTEGER, DIMENSION(0:NODESJ-1) :: JSIZE_OLD

      INTEGER :: JSIZE, JREMAIN
      INTEGER :: MAXVAL_NCPP_OLD,MINVAL_NCPP_OLD,MAXVAL_NCPP,MINVAL_NCPP
      INTEGER :: AVG_NCPP_OLD,AVG_NCPP
      DOUBLE PRECISION :: LIP_OLD,LIP,MAXSPEEDUP_OLD,MAXSPEEDUP,P

      INTEGER :: I_OFFSET,J_OFFSET,K_OFFSET,IERR

      INTEGER, DIMENSION(0:numPEs-1) :: disp,rcount

      INTEGER :: IPROC_OF_MAX_OLD,IPROC_OF_MIN_OLD

!-----------------------------------------------
!
      RETURN

      IF(.NOT.CARTESIAN_GRID) RETURN            ! Perform adjustement only when both CG
!      IF(.NOT.ADJUST_PROC_DOMAIN_SIZE) RETURN   ! and domain adjustment
      IF(NODESI*NODESJ*NODESK==1) RETURN         ! and parallel run are active



      IF(.not.allocated(ISIZE_ALL)) allocate( ISIZE_ALL(0:NODESI-1))
      IF(.not.allocated(JSIZE_ALL)) allocate( JSIZE_ALL(0:NODESJ-1))
      IF(.not.allocated(KSIZE_ALL)) allocate( KSIZE_ALL(0:NODESK-1))

!      allocate( ISIZE_ALL(0:NODESI-1))
!      allocate( JSIZE_ALL(0:NODESJ-1))
!      allocate( KSIZE_ALL(0:NODESK-1))




      IF(NODESI>1) THEN                                ! DOMAIN DECOMPOSITION IN I-DIRECTION

         IF(myPE == 0.AND.NODESJ*NODESK/=1) THEN
            WRITE(*,*)'ERROR IN SUBROUTINE REPORT_BEST_IJK_SIZE.'
            WRITE(*,*)'ADJUSTMENT POSSIBLE ONLY FOR DOMAIN DECOMPOSITION In ONE DIRECTION.'
            WRITE(*,*)'NODESI,NODESJ,NODESK =',NODESI,NODESJ,NODESK
            WRITE(*,*)'MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         ENDIF


         JSIZE_ALL(0:NODESJ-1) = jmax1-jmin1+1         ! Assign size in J and K-direction
         KSIZE_ALL(0:NODESK-1) = kmax1-kmin1+1

! Assign size in I-direction
         IF (myPE == 0) THEN
            WRITE(*,1000)'INFO FROM REPORT_BEST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO REPORT BEST DOMAIN SIZE IN I-DIRECTION ...'
         ENDIF

         DO I = ISTART1,IEND1

            NUC_I(I) = 0                     ! NUC : Number of Useful Cells

            DO J = JSTART3,JEND3
               DO K = KSTART3,KEND3

                  IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN   ! Count number of useful cells
                     NUC_I(I) = NUC_I(I) + 1
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

! Gather NUC onto the head node

         CALL allgather_1i (IEND1-ISTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            I_OFFSET = 0
         ELSE
            I_OFFSET = 0
            DO iproc=0,myPE-1
               I_OFFSET = I_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (I_OFFSET,disp,IERR)

         call gatherv_1i( NUC_I(ISTART1:IEND1), IEND1-ISTART1+1, GLOBAL_NUC_I(IMIN1:IMAX1), rcount, disp, PE_IO, ierr )



      ELSEIF(NODESJ>1) THEN                            ! DOMAIN DECOMPOSITION IN J-DIRECTION


         IF(myPE == 0.AND.NODESI*NODESK/=1) THEN
            WRITE(*,*)'ERROR IN SUBROUTINE REPORT_BEST_IJK_SIZE.'
            WRITE(*,*)'ADJUSTMENT POSSIBLE ONLY FOR DOMAIN DECOMPOSITION In ONE DIRECTION.'
            WRITE(*,*)'NODESI,NODESJ,NODESK =',NODESI,NODESJ,NODESK
            WRITE(*,*)'MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         ENDIF


         ISIZE_ALL(0:NODESI-1) = imax1-imin1+1         ! Assign size in I and K-direction
         KSIZE_ALL(0:NODESK-1) = kmax1-kmin1+1

! Assign size in J-direction
         IF (myPE == 0) THEN
            WRITE(*,1000)'INFO FROM REPORT_BEST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO REPORT BEST DOMAIN SIZE IN J-DIRECTION ...'
         ENDIF

         DO J = JSTART1,JEND1

            NUC_J(J) = 0                     ! NUC : Number of Useful Cells

            DO I = ISTART3,IEND3
               DO K = KSTART3,KEND3

                  IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN   ! Count number of useful cells
                     NUC_J(J) = NUC_J(J) + 1
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

! Gather NUC onto the head node

         IF(IS_SERIAL) THEN

            DO J=JMIN1,JMAX1
               GLOBAL_NUC_J(J) = NUC_J(J)
            ENDDO


         ELSE

            CALL allgather_1i (JEND1-JSTART1+1,rcount,IERR)

            IF (myPE == 0) THEN
               J_OFFSET = 0
            ELSE
               J_OFFSET = 0
               DO iproc=0,myPE-1
                  J_OFFSET = J_OFFSET + rcount(iproc)
               ENDDO
            ENDIF

            CALL allgather_1i (J_OFFSET,disp,IERR)

            call gatherv_1i( NUC_J(JSTART1:JEND1), JEND1-JSTART1+1, GLOBAL_NUC_J(JMIN1:JMAX1), rcount, disp, PE_IO, ierr )

         ENDIF

!         IF (myPE == 0) THEN
!            DO J=JMIN1,JMAX1
!               print*,'J,NUC=',J,GLOBAL_NUC_J(J)
!            ENDDO
!         ENDIF



      ELSEIF(NODESK>1) THEN                            ! DOMAIN DECOMPOSITION IN K-DIRECTION

         IF(myPE == 0.AND.NODESI*NODESJ/=1) THEN
            WRITE(*,*)'ERROR IN SUBROUTINE REPORT_BEST_IJK_SIZE.'
            WRITE(*,*)'ADJUSTMENT POSSIBLE ONLY FOR DOMAIN DECOMPOSITION In ONE DIRECTION.'
            WRITE(*,*)'NODESI,NODESJ,NODESK =',NODESI,NODESJ,NODESK
            WRITE(*,*)'MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         ENDIF


         ISIZE_ALL(0:NODESI-1) = imax1-imin1+1         ! Assign size in I and J-direction
         JSIZE_ALL(0:NODESJ-1) = jmax1-jmin1+1

! Assign size in K-direction
         IF (myPE == 0) THEN
            WRITE(*,1000)'INFO FROM REPORT_BEST_IJK_SIZE:'
            WRITE(*,1000)'ATTEMPTING TO REPORT BEST DOMAIN SIZE IN K-DIRECTION ...'
         ENDIF

         DO K = KSTART1,KEND1

            NUC_K(K) = 0                     ! NUC : Number of Useful Cells

            DO I = ISTART1,IEND1
               DO J = JSTART1,JEND1

                  IF( .NOT.DEAD_CELL_AT(I,J,K)) THEN   ! Count number of useful cells
                     NUC_K(K) = NUC_K(K) + 1
                  ENDIF

               ENDDO
            ENDDO

         ENDDO

! Gather NUC onto the head node

         CALL allgather_1i (KEND1-KSTART1+1,rcount,IERR)

         IF (myPE == 0) THEN
            K_OFFSET = 0
         ELSE
            K_OFFSET = 0
            DO iproc=0,myPE-1
               K_OFFSET = K_OFFSET + rcount(iproc)
            ENDDO
         ENDIF

         CALL allgather_1i (K_OFFSET,disp,IERR)

         call gatherv_1i( NUC_K(KSTART1:KEND1), KEND1-KSTART1+1, GLOBAL_NUC_K(KMIN1:KMAX1), rcount, disp, PE_IO, ierr )



      ENDIF !(NODESI,NODESJ, OR NODESK)




      IF (myPE == PE_IO) THEN                              ! DETERMINE BEST DOMAIN DECOMPOSITION FROM HEAD NODE

         IF(NODESI>1) THEN                                 ! DOMAIN DECOMPOSITION IN I-DIRECTION



         ELSEIF(NODESJ>1) THEN                            ! DOMAIN DECOMPOSITION IN J-DIRECTION

   ! For comparison with old grid size, determine the size in j direction (without adjustment) and add the remainder sequentially

            JSIZE = (JMAX1-JMIN1+1)/NODESJ
            JSIZE_OLD(0:NODESJ-1) = JSIZE

            JREMAIN = (JMAX1-JMIN1+1) - NODESJ*JSIZE
            IF (JREMAIN.GE.1) THEN
               JSIZE_OLD( 0:(JREMAIN-1) ) = JSIZE + 1
            ENDIF



           JSIZE_ALL = JSIZE_OLD

      CALL GET_LIP_WITH_GHOST_LAYERS(NODESJ,GLOBAL_NUC_J(JMIN1:JMAX1),JMIN1,JMAX1,JSIZE_ALL, &
           NCPP_OLD,NCPP_OLD_WITH_GHOST,LIP_OLD,IPROC_OF_MAX_OLD,IPROC_OF_MIN_OLD)


!      print*,'INITIAL ESTIMATE OF LIP, before minimizing LIP:',LIP_OLD
!         WRITE (*, 1000) '================================================='
!         WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!         WRITE (*, 1000) '================================================='
!         DO IPROC = 0,numPEs-1
!            WRITE (*, 1020) IPROC,JSIZE_ALL(IPROC),NCPP_OLD_WITH_GHOST(IPROC)
!         ENDDO
!         WRITE (*, 1000) '================================================='

           CALL MINIMIZE_LOAD_IMBALANCE(nodesj,GLOBAL_NUC_J(JMIN1:JMAX1),JMIN1,JMAX1,JSIZE_ALL,NCPP,NCPP_WITH_GHOST)


            TOTAL_NUC  = SUM(NCPP_WITH_GHOST(0:NumPEs-1))
            IDEAL_NCPP = TOTAL_NUC / NumPEs

            WRITE (*, 1000) 'AFTER OPTIMIZATION:'
            WRITE (*, 1010) 'TOTAL NUMBER OF USEFUL CELLS = ',TOTAL_NUC
            WRITE (*, 1010) 'IDEAL NUMBER OF CELLS/PROC.  = ',IDEAL_NCPP
            WRITE (*, 1000) 'ACTUALL CELL COUNT:'
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   DIFF. (%)'
            WRITE (*, 1000) '================================================='
            DO IPROC = 0,numPEs-1
               WRITE (*, 1020) IPROC,JSIZE_ALL(IPROC),NCPP_WITH_GHOST(IPROC), &
                    DBLE(NCPP_WITH_GHOST(IPROC)-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
            ENDDO


!            DO IPROC = 0 ,numPEs-1
!               NCPP_OLD(IPROC) = (imax3-imin3+1)*(JSIZE_ALL(IPROC)+4)
!               IF(DO_K) NCPP_OLD(IPROC) = NCPP_OLD(IPROC)*(kmax3-kmin3+1)
!            ENDDO


   ! Verify that the sum of all JSIZE_ALL matches JMAX

            PSUM = 0
            DO IPROC = 0,numPEs-1
               PSUM = PSUM + JSIZE_ALL(IPROC)
               IF(JSIZE_ALL(IPROC)<5) THEN
                  WRITE (*, 1010) 'ERROR: J-SIZE TOO SMALL FOR PROCESSOR:',IPROC
                  WRITE (*, 1010) 'J-SIZE = ',JSIZE_ALL(IPROC)
                  CALL MFIX_EXIT(myPE)
               ENDIF

            ENDDO

            IF(PSUM/=JMAX) THEN
               WRITE (*, 1000) 'ERROR IN ADJUST_IJK_SIZE: UNABLE TO ASSIGN JSIZE TO PROCESSORS.'
               WRITE (*, 1000) 'SUM OF JSIZE_ALL DOES NOT MATCH JMAX:'
               WRITE (*, 1010) 'SUM OF JSIZE_ALL = ',PSUM
               WRITE (*, 1010) 'JMAX1 = ',JMAX
               CALL MFIX_EXIT(myPE)
            ENDIF


            OPEN(CONVERT='BIG_ENDIAN',UNIT=777, FILE='suggested_gridmap.dat')
            WRITE (777, 1000) 'J-SIZE DISTRIBUTION'
            WRITE (777, 1010) 'NUMBER OF PROCESSORS = ',NumPEs
            WRITE (777, 1000) '================================================='
            WRITE (777, 1000) '   PROCESSOR    J-SIZE'
            WRITE (777, 1000) '================================================='

            DO IPROC = 0,NumPEs-1
                  WRITE(777,1060) IPROC,jsize_all(IPROC)
            ENDDO
            CLOSE(777)
            WRITE (*, 1000) '================================================='
            WRITE (*, 1000) 'J-SIZE DISTRIBUTION SAVED IN FILE: suggested_gridmap.dat'
            WRITE (*, 1000) 'TO USE THIS DISTRIBUTION, RENAME THE FILE AS: gridmap.dat'
            WRITE (*, 1000) 'AND RUN MFIX AGAIN.'
            WRITE (*, 1000) '================================================='


         ELSEIF(NODESK>1) THEN                            ! DOMAIN DECOMPOSITION IN K-DIRECTION




         ENDIF !(NODESI,NODESJ, OR NODESK)


         MAXVAL_NCPP_OLD = MAXVAL(NCPP_OLD_WITH_GHOST)
         MINVAL_NCPP_OLD = MINVAL(NCPP_OLD_WITH_GHOST)
         AVG_NCPP_OLD    = SUM(NCPP_OLD_WITH_GHOST)/NUMPES

!         LIP_OLD = DBLE(MAXVAL_NCPP_OLD-AVG_NCPP_OLD)/DBLE(AVG_NCPP_OLD)*100.0D0
         LIP_OLD = DBLE(MAXVAL_NCPP_OLD-MINVAL_NCPP_OLD)/DBLE(MINVAL_NCPP_OLD)*100.0D0

!         P = DBLE(MAXVAL_NCPP_OLD)/DBLE(AVG_NCPP_OLD)
         P = DBLE(MINVAL_NCPP_OLD)/DBLE(MAXVAL_NCPP_OLD)

!         MAXSPEEDUP_OLD = DBLE(NumPes)*(ONE-LIP_OLD/100.0D0)
         MAXSPEEDUP_OLD = ONE / ((ONE-P) + P/NumPes)

         MAXVAL_NCPP = MAXVAL(NCPP_WITH_GHOST)
         MINVAL_NCPP = MINVAL(NCPP_WITH_GHOST)
         AVG_NCPP    = SUM(NCPP_WITH_GHOST(0:NumPEs-1))/NUMPES

!         LIP = DBLE(MAXVAL_NCPP-AVG_NCPP)/DBLE(AVG_NCPP)*100.0D0
         LIP = DBLE(MAXVAL_NCPP-MINVAL_NCPP)/DBLE(MINVAL_NCPP)*100.0D0

!         P = DBLE(MAXVAL_NCPP)/DBLE(AVG_NCPP)
         P = DBLE(MINVAL_NCPP)/DBLE(MAXVAL_NCPP)

!         MAXSPEEDUP = DBLE(NumPes)*(ONE-LIP/100.0D0)
         MAXSPEEDUP = ONE / ((ONE-P) + P/NumPes)

         WRITE (*, 1000) '================================================='
         WRITE (*, 1000) 'ESTIMATED PARALLEL LOAD BALANCING STATISTICS'
         WRITE (*, 1000) 'COMPARISION BETWEEN UNIFORM SIZE (OLD)'
         WRITE (*, 1000) 'AND SUGGESTED SIZE (NEW)'
         WRITE (*, 1000) '================================================='
         WRITE (*, 1000) '                               OLD       NEW'
         WRITE (*, 1010) 'MAX CELL COUNT        : ',MAXVAL_NCPP_OLD,MAXVAL_NCPP
         WRITE (*, 1010) 'AT PROCESSOR          : ',MAXLOC(NCPP_OLD_WITH_GHOST)-1,MAXLOC(NCPP_WITH_GHOST)-1
         WRITE (*, 1010) 'MIN CELL COUNT        : ',MINVAL_NCPP_OLD,MINVAL_NCPP
         WRITE (*, 1010) 'AT PROCESSOR          : ',MINLOC(NCPP_OLD_WITH_GHOST)-1,MINLOC(NCPP_WITH_GHOST)-1
         WRITE (*, 1010) 'AVG CELL COUNT        : ',AVG_NCPP_OLD,AVG_NCPP
         WRITE (*, 1000) ''
         WRITE (*, 1030) 'LOAD IMBALANCE (%)    : ',LIP_OLD,LIP
         WRITE (*, 1000) ''
!         WRITE (*, 1030) 'IDEAL SPEEDUP         : ',DBLE(NumPEs),DBLE(NumPEs)
!         WRITE (*, 1030) 'MAX SPEEDUP           : ',MAXSPEEDUP_OLD,MAXSPEEDUP
!         WRITE (*, 1030) 'MAX EFFICIENCY (%)    : ',MAXSPEEDUP_OLD/DBLE(NumPEs)*100.0,MAXSPEEDUP/DBLE(NumPEs)*100.0

         WRITE (*, 1000) '================================================='



     ENDIF  ! (MyPE==PE_IO)



! Broadcast Domain sizes to all processors

!      CALL BCAST(ISIZE_ALL)
!      CALL BCAST(JSIZE_ALL)
!      CALL BCAST(KSIZE_ALL)

!      DOMAIN_SIZE_ADJUSTED = .TRUE.




1000  FORMAT(1x,A)
1010  FORMAT(1x,A,I10,I10)
1020  FORMAT(1X,I8,2(I12),F12.2)
1030  FORMAT(1X,A,2(F10.1))
1040  FORMAT(F10.1)
1050  FORMAT(1X,3(A))
1060  FORMAT(1x,I8,I12)

      RETURN

      END SUBROUTINE REPORT_BEST_IJK_SIZE0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_LIP_WITH_GHOST_LAYERS0                             C
!  Purpose: Compute Load Imbalance percentage                          C
!           by including size of ghost layers                          C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_LIP_WITH_GHOST_LAYERS0(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST,LIP,IPROC_OF_MAX,IPROC_OF_MIN)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NODESL,L,LMIN1,LMAX1,TOTAL_NUC,TOTAL_NUC_WITH_GHOST,IPROC_OF_MAX,IPROC_OF_MIN

      INTEGER :: LCOUNT1,LCOUNT2,MINVAL_NCPP,MAXVAL_NCPP,IDEAL_NCPP
      INTEGER, DIMENSION(LMIN1:LMAX1) :: NUC_L

      INTEGER :: IPROC

        INTEGER, DIMENSION(0:NODESL-1) :: NCPP,NCPP_WITH_GHOST,L_SIZE,L1,L2


      DOUBLE PRECISION :: LIP

!-----------------------------------------------

      LCOUNT1 = LMAX1 - LMIN1 + 1
      LCOUNT2 = SUM(L_SIZE(0:NODESL-1))

      IF(LCOUNT1/=LCOUNT2) THEN
         WRITE(*,*)' ERROR: SUM OF CELLS DO NOT MATCH:',LCOUNT1,LCOUNT2
         CALL MFIX_EXIT(myPE)
      ENDIF

      L1(0) = LMIN1
      L2(0) = L1(0) + L_SIZE(0) - 1

      DO IPROC = 1,NODESL-1
         L1(IPROC) = L2(IPROC-1) + 1
         L2(IPROC) = L1(IPROC) + L_SIZE(IPROC) - 1
      ENDDO

      DO IPROC = 0,NODESL-1
         NCPP(IPROC) = SUM(NUC_L(L1(IPROC):L2(IPROC)))
!         print*,'NUC=',NUC_L(L1(IPROC):L2(IPROC))
!         print*,'L1,L2=',IPROC,L1(IPROC),L2(IPROC),NCPP(IPROC)
      ENDDO

      TOTAL_NUC = 0

      DO L=LMIN1,LMAX1
         TOTAL_NUC = TOTAL_NUC + NUC_L(L)
      ENDDO

      NCPP_WITH_GHOST(0) = NCPP(0) + 2*NUC_L(L1(0)) + NUC_L(L1(1)) + NUC_L(L1(1)+1)

      DO IPROC = 1,NODESL-2
         NCPP_WITH_GHOST(IPROC) =   NCPP(IPROC)  &
                                  + NUC_L(L2(IPROC-1)) + NUC_L(L2(IPROC-1)-1) &
                                  + NUC_L(L1(IPROC+1)) + NUC_L(L1(IPROC+1)+1)
      ENDDO

      NCPP_WITH_GHOST(NODESL-1) = NCPP(NODESL-1) + 2*NUC_L(L2(NODESL-1)) + NUC_L(L2(NODESL-2)) + NUC_L(L2(NODESL-2)-1)

      TOTAL_NUC_WITH_GHOST = 0
      DO IPROC = 0,NODESL-1
!         print*,'NCPP_WITH_GHOST=',IPROC,L_SIZE(IPROC),NCPP(IPROC),NCPP_WITH_GHOST(IPROC)
         TOTAL_NUC_WITH_GHOST = TOTAL_NUC_WITH_GHOST + NCPP_WITH_GHOST(IPROC)
      ENDDO


      IDEAL_NCPP = TOTAL_NUC_WITH_GHOST / NumPEs

      MAXVAL_NCPP = MAXVAL(NCPP_WITH_GHOST)
      MINVAL_NCPP = MINVAL(NCPP_WITH_GHOST)

!      LIP = DBLE(MAXVAL_NCPP-IDEAL_NCPP)/DBLE(IDEAL_NCPP)*100.0D0
      LIP = DBLE(MAXVAL_NCPP-MINVAL_NCPP)/DBLE(MINVAL_NCPP)*100.0D0


      IPROC_OF_MAX = MAXLOC(NCPP_WITH_GHOST,1)-1
      IPROC_OF_MIN = MINLOC(NCPP_WITH_GHOST,1)-1


!      print*,'IPROC_OF_MIN=',IPROC_OF_MIN,MINVAL_NCPP
!      print*,'IPROC_OF_MAX=',IPROC_OF_MAX,MAXVAL_NCPP

      RETURN
      END SUBROUTINE GET_LIP_WITH_GHOST_LAYERS0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MINIMIZE_LOAD_IMBALANCE0                               C
!  Purpose: Rearrange L_SIZE to minimize load imbalance                C
!           by including size of ghost layers                          C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MINIMIZE_LOAD_IMBALANCE0(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE gridmap
      USE param
      USE param1
      USE constant
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar
      USE mpi_utility
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NODESL,LMIN1,LMAX1,IPROC_OF_MAX,IPROC_OF_MIN

      INTEGER, DIMENSION(LMIN1:LMAX1) :: NUC_L

      INTEGER :: NN,NOIMPROVEMENT

      INTEGER,PARAMETER :: NAMAX=10000  ! maximum number of adjustments, increase if optimized load is not reached

        INTEGER, DIMENSION(0:numPEs-1) :: NCPP,NCPP_WITH_GHOST,L_SIZE,BEST_L_SIZE,BEST_NCPP,BEST_NCPP_WITH_GHOST


      DOUBLE PRECISION :: LIP,BEST_LIP

!-----------------------------------------------


!      NA = NAMAX   ! Number of adjustments

! Initial estimate of LIP

      CALL GET_LIP_WITH_GHOST_LAYERS0(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST,LIP,IPROC_OF_MAX,IPROC_OF_MIN)

      BEST_LIP = LIP
      BEST_L_SIZE = L_SIZE

!      print*,'INITIAL ESTIMATE OF LIP:',LIP
!         WRITE (*, 1000) '================================================='
!         WRITE (*, 1010) 'AFTER STEP:',N
!         WRITE (*, 1010) 'NOIMPROVEMENT=',NOIMPROVEMENT
!         WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!         WRITE (*, 1000) '================================================='
!         DO IPROC = 0,numPEs-1
!            WRITE (*, 1020) IPROC,L_SIZE(IPROC),NCPP_WITH_GHOST(IPROC)
!         ENDDO
!         WRITE (*, 1000) '================================================='

!      print*,'MIN ',IPROC_OF_MIN,NCPP_WITH_GHOST(IPROC_OF_MIN)
!      print*,'MAX ',IPROC_OF_MAX,NCPP_WITH_GHOST(IPROC_OF_MAX)

      print*,'ATTEMPTING TO OPTIMIZE LOAD BALANCE...'

      NOIMPROVEMENT=0

      DO NN = 1,NAMAX

         L_SIZE(IPROC_OF_MAX) = L_SIZE(IPROC_OF_MAX) - 1
         L_SIZE(IPROC_OF_MIN) = L_SIZE(IPROC_OF_MIN) + 1

         CALL GET_LIP_WITH_GHOST_LAYERS(NODESL,NUC_L,LMIN1,LMAX1,L_SIZE,NCPP,NCPP_WITH_GHOST,LIP,IPROC_OF_MAX,IPROC_OF_MIN)

!      print*,'After adjustment of LIP:',N, LIP
!      print*,'MIN ',IPROC_OF_MIN,NCPP_WITH_GHOST(IPROC_OF_MIN)
!      print*,'MAX ',IPROC_OF_MAX,NCPP_WITH_GHOST(IPROC_OF_MAX)

         IF(LIP<BEST_LIP) THEN
            BEST_LIP    = LIP
            BEST_L_SIZE = L_SIZE
            BEST_NCPP   = NCPP
            BEST_NCPP_WITH_GHOST   = NCPP_WITH_GHOST
            NOIMPROVEMENT=0
         ELSE
            NOIMPROVEMENT = NOIMPROVEMENT + 1
         ENDIF

!         WRITE (*, 1000) '================================================='
!         WRITE (*, 1010) 'AFTER STEP:',N
!         WRITE (*, 1010) 'NOIMPROVEMENT=',NOIMPROVEMENT
!         WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!         WRITE (*, 1000) '================================================='
!         DO IPROC = 0,numPEs-1
!            WRITE (*, 1020) IPROC,L_SIZE(IPROC),NCPP_WITH_GHOST(IPROC)
!         ENDDO

         IF(NOIMPROVEMENT==10) THEN
            WRITE (*, 1000) 'OPTIMIZED LOAD BALANCE REACHED.'
            EXIT
         ENDIF

      ENDDO

!      print*,'Best LIP = ',BEST_LIP
      L_SIZE = BEST_L_SIZE
      NCPP   = BEST_NCPP
      NCPP_WITH_GHOST = BEST_NCPP_WITH_GHOST


!      WRITE (*, 1000) '================================================='
!      WRITE (*, 1000) 'AFTER OPTIMIZATION'
!      WRITE (*, 1000) '   PROCESSOR    J-SIZE   CELLS/PROC.   ERROR (%)'
!      WRITE (*, 1000) '================================================='
!      DO IPROC = 0,numPEs-1
!         WRITE (*, 1020) IPROC,L_SIZE(IPROC),NCPP_WITH_GHOST(IPROC)
!      ENDDO

1000  FORMAT(1x,A)
1010  FORMAT(1x,A,I10,I10)
1020  FORMAT(1X,I8,2(I12),F12.2)

      RETURN
      END SUBROUTINE MINIMIZE_LOAD_IMBALANCE0
END MODULE CHECK_DATA_CG
