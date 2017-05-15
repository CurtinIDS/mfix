!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_MU_s                                               C
!  Purpose: calculate MU_s from the RES file                           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-MAR-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJKMAX2, EP_s, RUN_NAME, TIME           C
!  Variables modified: M, IJK                                          C
!                                                                      C
!  Local variables: ARRAY                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_MU_s

      USE compar
      USE constant
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE rdf
      USE run, only: run_name, time
      USE visc_s

      IMPLICIT NONE
      INCLUDE 'xforms.inc'

      DOUBLE PRECISION    ARRAY(DIMENSION_3,DIMENSION_M), K_1m
      LOGICAL             INTER
      INTEGER             IER, M, IJK

      INTER = .FALSE.
!
      IF (DO_XFORMS) THEN
         C_E = E_PASS
      ELSE
         IF (C_E.EQ.UNDEFINED) THEN
            WRITE (*,*) ' Enter Coefficient of restitution (e) value'
            READ  (*,*) C_e
         END IF
      END IF
!
      CLOSE(UNIT_OUT)
      IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(TEMP_FILE)
      OPEN (UNIT=UNIT_OUT,FILE=TEMP_FILE,STATUS='UNKNOWN',convert='big_endian')
      CALL READ_RES1
      WRITE (UNIT_OUT,*) RUN_NAME
      WRITE (UNIT_OUT,*) ' '
      WRITE (UNIT_OUT,*) ' data from restart file'
      WRITE (UNIT_OUT,*) ' TIME = ' , TIME
!
      DO M = 1, MMAX
        CALL CALC_MU_s(M, IER)
      ENDDO
!
      DO 200 M = 1,MMAX
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) RETURN
         END IF
         DO 100 IJK = 1,IJKMAX2
            IF (EP_s(IJK,M) .NE. 0.0) THEN
               K_1m = 2.D0 * (ONE + C_e) * RO_s(IJK,M) * G_0(IJK, M, M)
               ARRAY(IJK,M) = P_s(IJK,M) / K_1m / EP_s(IJK,M)**2
            ELSE
               ARRAY(IJK,M) = 0.0
            END IF
100      CONTINUE
200   CONTINUE
!
      DO 300 M = 1,MMAX
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) RETURN
         END IF
         WRITE (UNIT_OUT,*) ' M    = ' , M
         WRITE (UNIT_OUT,*) ' '
         CALL OUT_ARRAY(ARRAY(1,M),'THETA')
         CALL OUT_ARRAY(MU_s(1,M),'MU_s')
300   CONTINUE
!
      CLOSE (UNIT=UNIT_OUT)
      RETURN
      END
