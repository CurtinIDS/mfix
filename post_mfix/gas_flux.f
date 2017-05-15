!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GAS_FLUX                                               C
!                                                                      C
!  Purpose: Calculate the time averaged gas flux (2 methods)           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 05-JUN-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: V_g, EP_g, IMIN1, IMAX1, JMIN1, JMAX1         C
!                        KMIN1, KMAX1, IJKMAX2                         C
!  Variables modified: PLOT_TYPE, VAR_INDEX, LOC_X, LOC_Y, LOC_Z       C
!                      I, J, K, IJK                                    C
!                                                                      C
!  Local variables: FILE_NAME, FILE1_INDEX, FILE2_INDEX, NX, NY, NZ    C
!                   L, NT                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GAS_FLUX
!
!
      Use param
      Use param1
      Use fldvar
      Use run
      Use geometry
      Use indices
      Use post3d
      Use physprop
      Use compar
      Use functions

      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      DOUBLE PRECISION  TAVG(DIMENSION_3,3)
      REAL              TIME_REAL(N_SPX)
      REAL              TIME_FOUND, TIME_NOW
      INTEGER           NSTEP_1
      INTEGER           NX , NY , NZ
      INTEGER           REC_POINTER(N_SPX) , L , NT
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
      INTEGER           I, J, K, IJK
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A)',ADVANCE='NO') &
                 ' Enter time to start and end time averaging > '
         READ  (*,*) TIME_START, TIME_END
         CALL GET_FILE_NAME(TEMP_FILE)
      END IF
!
      OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN',convert='big_endian')
!
      DO IJK = 1,IJKMAX2
         TAVG(IJK,1) = ZERO
         TAVG(IJK,2) = ZERO
         TAVG(IJK,3) = ZERO
      END DO
!
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
      READ_SPX(3) = .TRUE.         !  V_g
      READ_SPX(1) = .TRUE.         !  EP_g
      CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
      NT = 0
!
100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                          TIME_NOW, TIME_REAL, NSTEP_1)
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_END) GOTO 200
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      NT = NT + 1
      DO IJK = 1,IJKMAX2
         TAVG(IJK,1) = TAVG(IJK,1) + V_g(IJK)
         TAVG(IJK,2) = TAVG(IJK,2) + EP_g(IJK)
         TAVG(IJK,3) = TAVG(IJK,3) + EP_g(IJK) * V_g(IJK)
      END DO
      GOTO 100
!
200   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF
!
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE (40,*)' Avg(EP_g)*Avg(V_g) and Avg(EP_g*V_g) '
      WRITE (40,'(1X, 2(A, 2X, G12.5))')&
         'Start time = ', TIME_start,'End time = ',TIME_end
      WRITE (40,'(1X, 3(2X, A, I4))')&
        'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               TAVG(IJK,1) = TAVG(IJK,1) / REAL(NT)
               TAVG(IJK,2) = TAVG(IJK,2) / REAL(NT)
               TAVG(IJK,3) = TAVG(IJK,3) / REAL(NT)
               WRITE (40,*) TAVG(IJK,1)*TAVG(IJK,2) , TAVG(IJK,3)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
!
      WRITE(*,*)' Number of data points used = ', NT
      WRITE(*,*)' Avg(EP_g)*Avg(V_g) and Avg(EP_g*V_g) written'
!
      RETURN
      END
