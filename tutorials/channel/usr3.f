!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is
!           user-definable.  The user may insert code in this routine
!           or call appropriate user defined subroutines.
!           This routine is not called from an IJK loop, hence
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR3
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      Use usr
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!

      CALL WRITE_VELOCITY_PROFILES

      RETURN
      END SUBROUTINE USR3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VELOCITY_PROFILES                                C
!  Purpose: Save Velocity data at exact velocity node locations        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VELOCITY_PROFILES

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
      USE functions

      IMPLICIT NONE
      INTEGER :: I,J,K,IJK
      INTEGER I1,I2,J1,J2,N,NP
      CHARACTER (LEN=32) :: U_datafile,V_datafile
      DOUBLE PRECISION :: left_boundary,half_width,x_norm,u_norm,v_norm
      DOUBLE PRECISION :: error_u,error_v,u_exact,v_exact
      DOUBLE PRECISION :: angle,UREF,LREF
      DOUBLE PRECISION :: xu,yu,xv,yv
      INTEGER, PARAMETER :: ARRAY_SIZE = 10000
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: x_normalized
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: u_normalized
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: v_normalized
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: ue_normalized
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: ve_normalized
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: u_error
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE) :: v_error
      DOUBLE PRECISION, DIMENSION(ARRAY_SIZE,10) :: Temp_Array
      CHARACTER(LEN=3) :: BC

      OPEN(UNIT = 777, FILE= 'extract_velocity.inp')
      READ(777,*)BC
      READ(777,*)angle,UREF,LREF
      READ(777,*)I1,I2
      READ(777,*)J1,J2
      READ(777,*)U_datafile
      READ(777,*)V_datafile
      CLOSE(777)


      angle = HALF*DACOS(-1.0D0) - angle/180.0D0*DACOS(-1.0D0)

      half_width = LREF / DSIN(angle)

!======================================================================
!  Composite U_g profile
!======================================================================

      K = 1

      N = 0

      DO J = J1,J2
         DO I = I1,I2

            IJK = FUNIJK(I,J,K)

            XU = X_U(IJK)
            YU = Y_U(IJK)

            left_boundary = YU / DTAN(angle)

            x_norm = ( XU - left_boundary - half_width ) / (half_width)

            IF((x_norm>-ONE).AND.(x_norm<ONE)) THEN
               N = N + 1
               x_normalized(N) = x_norm
               u_normalized(N) = U_g(IJK) / UREF
            ENDIF
         ENDDO
      ENDDO

      NP = N

      Temp_Array(1:NP,1) = x_normalized(1:NP)
      Temp_Array(1:NP,2) = u_normalized(1:NP)
      CALL SORT(Temp_Array, NP,1,2)
      x_normalized(1:NP) = Temp_Array(1:NP,1)
      u_normalized(1:NP) = Temp_Array(1:NP,2)


      DO N = 1,NP
         IF(BC == 'NSW') THEN
            ue_normalized(N) = 1.5d0 * UREF * (1.0d0 - x_normalized(N)**2) * DCOS(angle)
         ELSEIF(BC == 'FSW') THEN
            ue_normalized(N) = UREF * DCOS(angle)
         ENDIF
         u_error(N) = dabs(u_normalized(N) - ue_normalized(N))
      ENDDO


      OPEN(UNIT = 777, FILE= U_datafile)
      WRITE(777,*)'       X         U_g        u_exact        Error'

      DO N = 1,NP
         WRITE(777,1000) x_normalized(N),u_normalized(N),ue_normalized(N),u_error(N)
      ENDDO

      CLOSE(777)

!======================================================================
!  Composite V_g profile
!======================================================================

      N = 0

      DO J = J1,J2
         DO I = I1,I2

            IJK = FUNIJK(I,J,K)

            XV = X_V(IJK)
            YV = Y_V(IJK)

            left_boundary = YV / DTAN(angle)
            half_width = LREF / DSIN(angle)

            x_norm = ( XV - left_boundary  - half_width ) / (half_width)

             IF((x_norm>-ONE).AND.(x_norm<ONE)) THEN
               N = N + 1
               x_normalized(N) = x_norm
               v_normalized(N) = V_g(IJK) / UREF
            ENDIF
         ENDDO
      ENDDO

      NP = N

      Temp_Array(1:NP,1) = x_normalized(1:NP)
      Temp_Array(1:NP,2) = v_normalized(1:NP)
      CALL SORT(Temp_Array, NP,1,2)
      x_normalized(1:NP) = Temp_Array(1:NP,1)
      v_normalized(1:NP) = Temp_Array(1:NP,2)


      DO N = 1,NP
         IF(BC == 'NSW') THEN
            ve_normalized(N) = 1.5d0 * UREF * (1.0d0 - x_normalized(N)**2) * DSIN(angle)
         ELSEIF(BC == 'FSW') THEN
            ve_normalized(N) = UREF * DSIN(angle)
         ENDIF
         v_error(N) = dabs(v_normalized(N) - ve_normalized(N))
      ENDDO


      OPEN(UNIT = 777, FILE= V_datafile)
      WRITE(777,*)'       X         V_g        v_exact        Error'

      DO N = 1,NP
         WRITE(777,1000) x_normalized(N),v_normalized(N),ve_normalized(N),v_error(N)
      ENDDO

      CLOSE(777)


1000  FORMAT(5(F12.5,2X))


      RETURN


      END SUBROUTINE WRITE_VELOCITY_PROFILES


