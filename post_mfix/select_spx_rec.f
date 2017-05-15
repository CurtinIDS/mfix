!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SELECT_SPX_REC                                         C
!  Purpose: Select records from SPx files interactively and write to a C
!           new SPx file                                               C
!           NOTE : USES N_SPX                                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-APR-94   C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2                                       C
!  Variables modified: IJK                                             C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SELECT_SPX_REC
!
!
      Use param
      Use param1
      Use geometry
      Use indices
      Use run
      Use machine
      Use funits
      Use post3d
      Use physprop
      Use fldvar
      IMPLICIT NONE
      INCLUDE 'xforms.inc'

      LOGICAL  AT_EOF(N_SPX), READ_SPX(N_SPX),SELECT
      INTEGER  REC_POINTER(N_SPX)
      INTEGER  NSTEP_1
!      CHARACTER(LEN=1)  IANS
      REAL     TIME_REAL(N_SPX) , tstart
      LOGICAL  ERROR
!!      double precision :: uavg(20000) , vavg(20000) , wavg(20000)
!!      double precision :: uavg2(20000) , vavg2(20000) , wavg2(20000)
!
      INTEGER L, L_SPX
!!      integer :: tcount
      integer :: unit_add = 10
!
!!      uavg(:) = 0.0
!!      vavg(:) = 0.0
!!      wavg(:) = 0.0
!!      tcount  = 0
!
      ERROR  = .FALSE.
      SELECT = .TRUE.
      L_SPX  = SPX_NUM
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A)',ADVANCE='NO') 'Enter the number of Spx file > '
         READ  (*,*) L_SPX
10       WRITE (*,'(A)',ADVANCE='NO') 'Enter the name of new file > '
         READ  (*,'(A60)') TEMP_FILE
         unit_add = 10
         OPEN (UNIT=UNIT_SPX+L_SPX+unit_add,FILE=TEMP_FILE, &
            STATUS='NEW',&
            RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=10,convert='big_endian')
         write (*,*) ' enter starting time to write output'
         read  (*,*) tstart
      END IF
!
!      IF (DO_XFORMS.AND.GET_TIMES) THEN
!         OPEN (UNIT=UNIT_SPX+L_SPX,FILE=TEMP_FILE,STATUS='UNKNOWN',
!     &      RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=10)
!      END IF
!
      CALL WRITE_SPX0(L_SPX,unit_add)
!
      DO 20 L = 1, N_SPX
        READ_SPX(L) = .FALSE.
        REC_POINTER(L) = 4
        AT_EOF(L) = .FALSE.
20    CONTINUE
      READ_SPX(L_SPX) = .TRUE.
!
      L = 0
100   continue
      CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL,NSTEP_1)
      IF (.NOT.AT_EOF(L_SPX)) THEN
         L = L + 1
         IF (.NOT.DO_XFORMS) THEN

!            WRITE(*,'(A,G12.5,A)',ADVANCE='NO')&
!                 'Write time ', TIME_REAL(L_SPX), '(Y/N) [N]'
!            READ(*,'(A1)')IANS
!            IF (IANS .EQ. 'y' .OR. IANS .EQ. 'Y') THEN

             if (time_real(l_spx).ge.tstart) then
!!                tcount = tcount + 1
                TIME = DBLE(TIME_REAL(L_SPX))

!!               uavg(:) = uavg(:) + u_s(:,1)
!!               uavg2(:) = uavg(:) / real(tcount)

                CALL WRITE_SPX1(L_SPX,unit_add)
                write (*,*) 'write data ... time = ' , time_real(l_spx)

            END IF
         ENDIF
         GOTO 100
      ENDIF
      CLOSE(UNIT_SPX+L_SPX+unit_add)
      WRITE(*,*)
      WRITE(*,'(A,I1,A)')' A New SP',L_SPX, ' file written '
      RETURN
      END
