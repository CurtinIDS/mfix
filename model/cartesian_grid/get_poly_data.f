
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_POLY_DATA                                          C
!  Purpose: reads polygon(s) coordinates from poly.dat (2D only)       C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!  Reviewer:                                          Date: **-***-**  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_POLY_DATA

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------

      USE compar
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE polygon
      USE progress_bar
      USE run
      USE rxns
      USE scalars
      IMPLICIT NONE

      INTEGER :: POLY,V,NN,NSKIP
      LOGICAL :: PRESENT

      WRITE(*,2000) 'READING polygon geometry from poly.dat...'

      INQUIRE(FILE='poly.dat',EXIST=PRESENT)
      IF(.NOT.PRESENT) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,"('(PE ',I3,'): input data file, ',A11,' is missing: run aborted')") &
            myPE,'poly.dat'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF
!
!
!     OPEN poly.dat ASCII FILE
!
      OPEN(CONVERT='BIG_ENDIAN',UNIT=333, FILE='poly.dat', STATUS='OLD', ERR=910)

      NSKIP = 13
      DO NN=1,NSKIP
         READ(333,*,ERR=920,END=930)   ! Skip first NSKIP comment lines
      ENDDO

      READ(333,*,ERR=920,END=930) N_POLYGON

      DO POLY = 1, N_POLYGON
         READ(333,*,ERR=920,END=930) N_VERTEX(POLY), POLY_SIGN(POLY)
         IF(DABS(POLY_SIGN(POLY))/=ONE) THEN
            WRITE(*,*)'ERROR WHILE READIND poly.dat:'
            WRITE(*,*)'POLYGON SIGN MUST BE +1.0 or -1.0.'
            CALL MFIX_EXIT(myPE)
         ENDIF
         DO V = 1,N_VERTEX(POLY)
            READ(333,*,ERR=920,END=930) X_VERTEX(POLY,V),Y_VERTEX(POLY,V),BC_ID_P(POLY,V)
         ENDDO
         X_VERTEX(POLY,N_VERTEX(POLY) + 1) = X_VERTEX(POLY, 1)  ! Copy first vertex as last vertex to close the polygon
         Y_VERTEX(POLY,N_VERTEX(POLY) + 1) = Y_VERTEX(POLY, 1)
      ENDDO

      CLOSE(333)

      WRITE(*,2010) 'Polygon geometry successfully read for ', N_POLYGON, ' polygon(s).'

      RETURN

!======================================================================
!     HERE IF AN ERROR OCCURED OPENNING/READING THE FILE
!======================================================================
!
 910  CONTINUE
      WRITE (*, 1500)
      CALL MFIX_EXIT(myPE)
 920  CONTINUE
      WRITE (*, 1600)
      CALL MFIX_EXIT(myPE)
 930  CONTINUE
      WRITE (*, 1700)
      CALL MFIX_EXIT(myPE)
!
 1500 FORMAT(/1X,70('*')//' From: GET_POLY_DATA',/' Message: ',&
      'Unable to open poly.dat file',/1X,70('*')/)
 1600 FORMAT(/1X,70('*')//' From: GET_POLY_DATA',/' Message: ',&
      'Error while reading poly.dat file',/1X,70('*')/)
 1700 FORMAT(/1X,70('*')//' From: GET_POLY_DATA',/' Message: ',&
      'End of file reached while reading poly.dat file',/1X,70('*')/)
 2000 FORMAT(1X,A)
 2010 FORMAT(1X,A,I4,A)

      END SUBROUTINE GET_POLY_DATA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EVAL_POLY_FCT                                          C
!  Purpose: Evaluates a user-defined function                          C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE EVAL_POLY_FCT(x1,x2,x3,Q,f_pol,CLIP_FLAG,BCID)

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
      USE fldvar
      USE quadric
      USE cutcell
      USE polygon

      IMPLICIT NONE
      DOUBLE PRECISION x1,x2,x3
      DOUBLE PRECISION f_pol
      DOUBLE PRECISION, DIMENSION(DIM_POLYGON) :: F_POLY
      INTEGER :: Q
      LOGICAL :: CLIP_X,CLIP_Y,CLIP_FLAG
      LOGICAL :: test1,test2

      INTEGER :: POLY,N_EDGE,COUNTER,EDGE
      DOUBLE PRECISION :: X_VERTEX_MIN,X_VERTEX_MAX,Y_VERTEX_MIN,Y_VERTEX_MAX
      DOUBLE PRECISION :: XV1,XV2,YV1,YV2
      DOUBLE PRECISION :: x_west,x_east,y_north,y_south,slope,x_star,D1,D2
      DOUBLE PRECISION :: max_slope = 200.0D0
      DOUBLE PRECISION :: TOL_xstar
      INTEGER :: BCID

!======================================================================
! This subroutine checks whether a point P(x1,x2,x3) lies inside any
! of the N_POLGON polygon(s), using the Shimrat's algorithm:
! 1) Draw a horizontal line from point P to -infinity (here taken as
!    the minimum x-location of the polygon vertices)
! 2) Count the number of intersection(s) with all edges (COUNTER)
! 3) If COUNTER is even, the point is outside
!    If COUNTER is odd,  the point is inside
! Once the point's location is determined, a value is assigned to f_pol:
!    fpol = 0 if point lies on top of one of the polygon's edges
!                (within tolerance TOL_POLY)
!    fpol = POLY_SIGN(POLY) if point lies inside the polygon
!    fpol = - POLY_SIGN(POLY) if point lies outside the polygon
!    where
!           POLY_SIGN(POLY) = -1.0 means the interior of the polygon
!                                  is part of the computational domain
!           POLY_SIGN(POLY) = +1.0 means the interior of the polygon
!                                  is excluded from the computational domain
!======================================================================

      IF(N_POLYGON < 1) RETURN

      BCID = -1

      DO POLY = 1,N_POLYGON

         N_EDGE   = N_VERTEX(POLY)

         X_VERTEX_MIN = MINVAL(X_VERTEX(POLY,:)) - 2.0D0 * TOL_POLY
         X_VERTEX_MAX = MAXVAL(X_VERTEX(POLY,:)) + 2.0D0 * TOL_POLY

         Y_VERTEX_MIN = MINVAL(Y_VERTEX(POLY,:)) - 2.0D0 * TOL_POLY
         Y_VERTEX_MAX = MAXVAL(Y_VERTEX(POLY,:)) + 2.0D0 * TOL_POLY

         CLIP_X = (X_VERTEX_MIN <= x1).AND.( x1 <= X_VERTEX_MAX)
         CLIP_Y = (Y_VERTEX_MIN <= x2).AND.( x2 <= Y_VERTEX_MAX)


         IF(CLIP_X.AND.CLIP_Y) THEN

            COUNTER = 0

            DO EDGE = 1, N_EDGE

               XV1 = X_VERTEX(POLY,EDGE)
               XV2 = X_VERTEX(POLY,EDGE + 1)

               YV1 = Y_VERTEX(POLY,EDGE)
               YV2 = Y_VERTEX(POLY,EDGE + 1)


               D1 = DSQRT((x1 - XV1)**2 + (x2 - YV1)**2)
               D2 = DSQRT((x1 - XV2)**2 + (x2 - YV2)**2)
               IF(DMIN1(D1,D2)<TOL_POLY) THEN
                  F_POL = ZERO
                  BCID = BC_ID_P(POLY,EDGE)
                  RETURN
               ENDIF

               IF(DABS(YV2 - YV1) < TOL_POLY) THEN
                  test1 = (DABS(x2 - YV1) < TOL_POLY)
                  x_west = DMIN1(XV1,XV2)
                  x_east = DMAX1(XV1,XV2)
                  test2 = (x_west <= x1).AND.( x1 <= x_east)
                  IF (test1.and.test2) THEN
                     F_POL = ZERO
                     BCID = BC_ID_P(POLY,EDGE)
                     RETURN
                  ENDIF
               ELSE
                  y_south = DMIN1(YV1,YV2)
                  y_north = DMAX1(YV1,YV2)
                  test1 = (y_south <= x2).AND.( x2 <= y_north)
                  IF(test1) THEN
                     slope = (XV2 - XV1) / (YV2 - YV1)
                     x_star = XV1 + slope * (x2 - YV1)

                     IF(DABS(slope)>max_slope) THEN
                        TOL_xstar = max_slope * TOL_POLY
                     ELSE
                        TOL_xstar = TOL_POLY
                     ENDIF

                     IF (x_star < x1 - TOL_xstar) THEN
                        COUNTER =COUNTER + 1
                     ELSEIF(DABS(x_star - x1) <= TOL_xstar) THEN
                        F_POL = ZERO
                        BCID = BC_ID_P(POLY,EDGE)
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            IF((COUNTER/2) * 2 == COUNTER) THEN
               F_POLY(POLY) = -POLY_SIGN(POLY)    ! even ==> outside
            ELSE
               F_POLY(POLY) = POLY_SIGN(POLY)     ! odd  ==> inside
            ENDIF

         ELSE
            F_POLY(POLY) = -POLY_SIGN(POLY)       ! outside
         ENDIF

      ENDDO

      f_pol = MAXVAL(F_POLY(1:N_POLYGON))

      CLIP_FLAG = .TRUE.

      RETURN

      END SUBROUTINE EVAL_POLY_FCT
