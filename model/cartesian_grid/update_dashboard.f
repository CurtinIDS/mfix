!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_DASHBOARD                                       C
!  Purpose: Updates and writes dashboard file                          C
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
      SUBROUTINE UPDATE_DASHBOARD(NIT,TLEFT,TUNIT)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE dashboard
      USE leqsol
      USE machine
      USE parallel
      USE run, ONLY: get_tunit, description, dt, dt_dir, run_name, time, tstop
      USE sendrecv
      USE time_cpu
      USE vtk
      Use residual

      IMPLICIT NONE

      DOUBLE PRECISION :: Smass,TLEFT
      DOUBLE PRECISION :: CPU_NOW
      INTEGER :: NIT
      CHARACTER(LEN=4)  :: TUNIT,CPU_TUNIT
!     temporary array to hold time data
      INTEGER DAT(8)
      CHARACTER(LEN=10) DATE, TIM, ZONE
      LOGICAL :: Sm_flag

! Intel Linux compiler supports this function thru it's portability library
      CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
      ID_YEAR   = DAT(1)
      ID_MONTH  = DAT(2)
      ID_DAY    = DAT(3)
      ID_HOUR   = DAT(5)
      ID_MINUTE = DAT(6)
      ID_SECOND = DAT(7)

      CALL CPU_TIME (CPU_NOW)
      CALL GET_TUNIT(CPU_NOW,CPU_TUNIT)

      IF (CPU_NOW > 60.0d0 .AND. CPU_TUNIT == 's') THEN
         CPU_NOW = CPU_NOW/60.0d0
         CPU_TUNIT = 'min'
      ENDIF

      IF (TLEFT > 60.0d0 .AND. TUNIT == 's') THEN
         TLEFT = TLEFT/60.0d0
         TUNIT = 'min'
      ENDIF

      Sm_flag = (TRIM(RUN_STATUS)=='In Progress...'.OR.TRIM(RUN_STATUS)=='Complete.')

      IF(Sm_flag) THEN
         CALL GET_SMASS (SMASS)

         SMMIN = DMIN1(SMASS,SMMIN)
         SMMAX = DMAX1(SMASS,SMMAX)
      ENDIF

      DTMIN = DMIN1(DT,DTMIN)
      DTMAX = DMAX1(DT,DTMAX)

!      IF(NIT==0) NIT=NIT_MIN

      NIT_MIN = MIN0(NIT,NIT_MIN)
      NIT_MAX = MAX0(NIT,NIT_MAX)

      IF(myPE /= PE_IO) RETURN



      OPEN(CONVERT='BIG_ENDIAN',UNIT     =  111           , &
           FILE     = 'DASHBOARD.TXT', &
           FORM     = 'FORMATTED'    , &
           ACCESS   = 'SEQUENTIAL'   , &
           STATUS   = 'REPLACE'      , &
           ACTION   = 'WRITE')


!      OPEN(CONVERT='BIG_ENDIAN',UNIT=111,FILE='DASHBOARD.TXT',STATUS='REPLACE')
      WRITE(111,30) '  _____________________________________________________________________________ '
      WRITE(111,30) ' |                                                                             |'
      WRITE(111,30) ' |                                MFIX DASHBOARD                               |'
      WRITE(111,30) ' |_____________________________________________________________________________|'
      WRITE(111,30) ' |                                                                             |'
      WRITE (111,10)' | RUN_NAME             = ',RUN_NAME,                                         '|'
      WRITE (111,10)' | Description          = ',DESCRIPTION,                                      '|'
      WRITE (111,10)' | Run Status           = ',RUN_STATUS,                                       '|'
      IF(RUN_STATUS/='Complete.') THEN
      WRITE (111,15)' | CPU time elapsed     = ',CPU_NOW,CPU_TUNIT,                                '|'
      WRITE (111,15)' | CPU time left        = ',TLEFT,TUNIT,                                      '|'
      ELSE
      WRITE (111,15)' | CPU time used        = ',TLEFT,TUNIT,                                      '|'
      ENDIF
      IF(WRITE_VTK_FILES) WRITE(111,10) ' | Latest vtu file      = ',VTU_FILENAME,                 '|'
      IF(IS_SERIAL) THEN
         WRITE (111,30)' | Serial run                                                                  |'
      ELSE
         WRITE (111,25)' | Parallel run, numPEs = ',numPEs,                                        '|'
      ENDIF


      WRITE(111,30) ' |_____________________________________________________________________________|'
      WRITE(111,30) ' |         |         |         |         |         |                           |'
      WRITE(111,30) ' |  Name   |  Value  |   Min   |   Max   | % of max|0%       Progress      100%|'
      WRITE(111,30) ' |_________|_________|_________|_________|_________|___________________________|'
      WRITE(111,30) ' |         |         |         |         |         |                           |'
      WRITE(111,40,ADVANCE='NO')' Time    ',Time,Init_Time,Tstop
      CALL WRITE_SIMPLE_PROGRESS_BAR(Time,TStop-Init_Time)
      IF(DT_DIR>0) THEN
         WRITE(111,40,ADVANCE='NO')' DT (+)  ',DT,DTMIN,DTMAX
!         WRITE(111,40,ADVANCE='NO')' DT      ',DT,DTMIN,DTMAX
      ELSE
         WRITE(111,40,ADVANCE='NO')' DT (-)  ',DT,DTMIN,DTMAX
      ENDIF
      CALL WRITE_SIMPLE_PROGRESS_BAR(DT,DTMAX)
      IF(Sm_flag) THEN
         WRITE(111,40,ADVANCE='NO')' Sm      ',SMASS,SMMIN,SMMAX
         CALL WRITE_SIMPLE_PROGRESS_BAR(SMASS,SMMAX)
      ELSE
      WRITE(111,30) ' | Sm      |         |         |         |         |                           |'
      ENDIF
      WRITE(111,50,ADVANCE='NO')' NIT     ',NIT,NIT_MIN,NIT_MAX
      CALL WRITE_SIMPLE_PROGRESS_BAR(dble(NIT),dble(NIT_MAX))
      IF (RESID_INDEX(8,1) == UNDEFINED_I) THEN
         WRITE (111,55) ' Max res ',RESID_STRING(8)
      ENDIF
      WRITE(111,30)' |_________|_________|_________|_________|_________|___________________________|'
      WRITE(111,60)'  Last updated at: ',ID_HOUR,ID_MINUTE,ID_SECOND,' on: ',ID_MONTH,ID_DAY,ID_YEAR
10    FORMAT(A,A,T80,A)
15    FORMAT(A, F6.1, 1X, A,T80,A)
25    FORMAT(A,I6,T80,A)
30    FORMAT(A)
40    FORMAT(' |',A,'|',3(E9.2,'|'))
50    FORMAT(' |',A,'|',3(I9,'|'))
55    FORMAT(' |',A,'|',A7,'  |         |         |         |                           |')
60    FORMAT(A,I2.2,':',I2.2,':',I2.2,A,I2.2,'/',I2.2,'/',I4)
      CLOSE(111)

      RETURN
      END SUBROUTINE UPDATE_DASHBOARD

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SIMPLE_PROGRESS_BAR                              C
!  Purpose: Displays a progress bar on the screen                      C
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
      SUBROUTINE WRITE_SIMPLE_PROGRESS_BAR(x,x_MAX)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param1, only: zero


      IMPLICIT NONE

      INTEGER            :: BAR_WIDTH
      CHARACTER (LEN=1)  :: BAR_CHAR
      DOUBLE PRECISION   :: BAR_RESOLUTION

      INTEGER :: PROGRESS
      INTEGER :: P
      CHARACTER (LEN=9) :: TEXT
      CHARACTER (LEN=69) :: PROGRESSBAR
      DOUBLE PRECISION :: PERCENT,PTEST
      DOUBLE PRECISION :: x,x_max


      IF(X_MAX==ZERO) THEN
         WRITE(111,5)
5        FORMAT(9X,'|',27X,'|')
         RETURN
      ENDIF

      BAR_WIDTH = 27
      BAR_CHAR = '='
      BAR_RESOLUTION = 1.0

      PERCENT  = x/x_MAX * 100.0
      PROGRESS = INT(PERCENT * BAR_WIDTH)

      WRITE(TEXT,10) PERCENT
10    FORMAT(' ',F5.1,' % ')

      DO P = 1, BAR_WIDTH
         PTEST = FLOAT(P)/FLOAT(BAR_WIDTH) * 100.0
         IF(PERCENT<PTEST-BAR_RESOLUTION) THEN
            PROGRESSBAR(P:P)= ' '
         ELSE
            PROGRESSBAR(P:P)= BAR_CHAR
         ENDIF
      ENDDO

      WRITE(111,15)TEXT,'|',TRIM(PROGRESSBAR),'|'

15    FORMAT(A,A,A27,A)

      RETURN
      END SUBROUTINE WRITE_SIMPLE_PROGRESS_BAR



