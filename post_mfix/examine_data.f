!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EXAMINE_DATA                                           C
!  Purpose: Examine/print selected data                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 03-NOV-93  C
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

      ! format_one
      subroutine format_one(spec,nPrec,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec,w,sLen

      if (nPrec.lt.0) then
         spec = '(a,g12.5)'
         sLen = 12
      else
         w = nPrec + 8
         spec = '(a,exx.xxE4)'
         write (spec(5:6),'(i2.2)') w
         write (spec(8:9),'(i2.2)') nPrec
         sLen = w + 5
      end if
      return
      end subroutine format_one

      ! format_oneB
      subroutine format_oneB(spec,nPrec,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec,w,sLen

      if (nPrec.lt.0) then
         spec = '(1x,g12.5)'
         sLen = 13
      else
         w = nPrec + 8
         spec = '(1x,exx.xxE4)'
         write (spec(6:7),'(i2.2)') w
         write (spec(9:10),'(i2.2)') nPrec
         sLen = w + 1
      end if
      return
      end subroutine format_oneB


      ! format_oneC
      subroutine format_oneC(spec,nPrec,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec,w,sLen

      if (nPrec.lt.0) then
         spec = '(1x,1a8,a,g12.5)'
         sLen = 21
      else
         w = nPrec + 8
         spec = '(1x,1a8,a,exx.xxE4)'
         write (spec(12:13),'(i2.2)') w
         write (spec(15:16),'(i2.2)') nPrec
         sLen = w + 9
      end if
      return
      end subroutine format_oneC

      ! format_two
      subroutine format_two(spec,nPrec1,nPrec2,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec1,nPrec2,w1,w2,sLen

      if (nPrec1.lt.0) then
         spec = '(a,g12.5,a,g12.5)'
         sLen = 24
      else
         w1 = nPrec1 + 8
         w2 = nPrec1 + 8
         spec = '(a,exx.xxE4,a,exx.xxE4)'
         write (spec(5:6),'(i2.2)') w1
         write (spec(8:9),'(i2.2)') nPrec1
         write (spec(16:17),'(i2.2)') w2
         write (spec(19:20),'(i2.2)') nPrec2
         sLen = w1 + w2 + 10
      end if
      return
      end subroutine format_two

      ! format_twoB
      subroutine format_twoB(spec,nPrec1,nPrec2,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec1,nPrec2,w1,w2,sLen

      if (nPrec1.lt.0) then
         spec = '(1X,A,1A8,A,G12.5,A,G12.5)'
         sLen = 24
      else
         w1 = nPrec1 + 8
         w2 = nPrec2 + 8
         spec = '(1x,a,1A8,A,exx.xxE4,a,exx.xxE4)'
         !       1234567890123456789012345678901
         write (spec(14:15),'(i2.2)') w1
         write (spec(17:18),'(i2.2)') nPrec1
         write (spec(25:26),'(i2.2)') w2
         write (spec(28:29),'(i2.2)') nPrec2
         sLen = w1 + w2 + 10
      end if
      return
      end subroutine format_twoB

       ! format_twoC
      subroutine format_twoC(spec,nPrec1,nPrec2,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec1,nPrec2,w1,w2,sLen

      if (nPrec1.lt.0) then
         spec = '(1X,G12.5,2X,G12.5)'
         sLen = 27
      else
         w1 = nPrec1 + 8
         w2 = nPrec2 + 8
         spec = '(1x,exx.xxE4,2x,exx.xxE4)'
         !       123456789012345678901234567890
         write (spec(6:7),'(i2.2)') w1
         write (spec(9:10),'(i2.2)') nPrec1
         write (spec(18:19),'(i2.2)') w2
         write (spec(21:22),'(i2.2)') nPrec2
         sLen = w1 + w2 + 3
      end if
      return
      end subroutine format_twoC


      ! format_four
      subroutine format_four(spec,nPrec1,nPrec2,nPrec3,nPrec4,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec1,nPrec2,nPrec3,nPrec4,w1,w2,w3,w4,sLen

      if (nPrec1.lt.0) then
         spec = '(1X,4(G12.5,2X))'
         sLen = 57 ! = 1 + 4*14
      else
         w1 = nPrec1 + 8
         w2 = nPrec2 + 8
         w3 = nPrec3 + 8
         w4 = nPrec4 + 8
         spec = &
         '(1x,exx.xxE4,2x,exx.xxE4,2x,exx.xxE4,2x,exx.xxE4)'
         !123  67 90123456 89 12 4567  01 3456789  23 56
         write (spec(6:7),'(i2.2)')   w1
         write (spec(9:10),'(i2.2)')   nPrec1
         write (spec(18:19),'(i2.2)') w2
         write (spec(21:22),'(i2.2)') nPrec2
         write (spec(30:31),'(i2.2)') w3
         write (spec(33:34),'(i2.2)') nPrec3
         write (spec(42:43),'(i2.2)') w4
         write (spec(45:46),'(i2.2)') nPrec4
         sLen = w1 + w2 + w3 + w4 + 7
      end if
      return
      end subroutine format_four



      ! format_five
      subroutine format_five(spec,nPrec1,nPrec2,nPrec3,nPrec4,nPrec5,sLen)
      implicit none
      character(len=*) :: spec
      integer       :: nPrec1,nPrec2,nPrec3,nPrec4,nPrec5,w1,w2,sLen
      integer       :: w3,w4,w5

      if (nPrec1.lt.0) then
         spec = '(1X,5(G12.5,2X))'
         sLen = 71 ! = 1 + 5*14
      else
         w1 = nPrec1 + 8
         w2 = nPrec2 + 8
         w3 = nPrec3 + 8
         w4 = nPrec4 + 8
         w5 = nPrec5 + 8
         spec = &
         '(1x,exx.xxE4,2x,exx.xxE4,2x,exx.xxE4,2x,exx.xxE4,2x,exx.xxE4)'
         !123 56 89012345 78 01234567 90 23456789 12 45678901 34 67890
         write (spec(6:7),'(i2.2)')   w1
         write (spec(9:10),'(i2.2)')   nPrec1
         write (spec(18:19),'(i2.2)') w2
         write (spec(21:22),'(i2.2)') nPrec2
         write (spec(30:31),'(i2.2)') w3
         write (spec(33:34),'(i2.2)') nPrec3
         write (spec(42:43),'(i2.2)') w4
         write (spec(45:46),'(i2.2)') nPrec4
         write (spec(54:55),'(i2.2)') w5
         write (spec(57:58),'(i2.2)') nPrec5
         sLen = w1 + w2 + w3 + w4 + w5 + 9
      end if
      return
      end subroutine format_five


!        WRITE(LINE,'(1X,5(G12.5,2X))')TIME_NOW, XTMP, YTMP, ZTMP, &
!                                      VALUE_TMP


      SUBROUTINE EXAMINE_DATA
!
      Use param, only: dimension_3
      Use param1
      Use constant
      Use physprop
      Use fldvar
      Use indices
      Use run, only: any_solve_ros, k_epsilon, time, run_name
      Use geometry
      Use post3d, only: xdist_vec, xdist_sc, ydist_vec, ydist_sc, zdist_vec, zdist_sc
      Use rxns
      Use scalars
      Use compar
      use post_precision
      Use functions

      IMPLICIT NONE
      INTEGER  N_VAR
      PARAMETER (N_VAR=52)
      INCLUDE 'xforms.inc'

      CHARACTER(LEN=80)  :: LINE
      CHARACTER(LEN=120)  :: STRING, SUBSTR
      CHARACTER(LEN=8)   :: VAR, VAR_DAT(N_VAR)
      CHARACTER(LEN=120) :: FILE_NAME
      INTEGER      L, L3, L4, LMAX, IANS, NSTEP_1
      REAL         DIST(DIMENSION_3), VALUE(DIMENSION_3)
      REAL         TIME_IN_RES
      INTEGER      DISPLAY, DIRECTION, NT
      INTEGER      I2d, J2d, K2d
      LOGICAL      FILE_EXIST,END_AVERAGE
      LOGICAL      SUM
      LOGICAL      STRCMP,INTER
      REAL         XTMP, YTMP, ZTMP, VALUE_TMP
      INTEGER      IJK1
      INTEGER      M_LOCAL, mIJK, lIJK, IER
      INTEGER      I, J, K, IJK, M, N
      REAL         DELm, DELl, FAC1, FAC2
!
      REAL              TIME_REAL(N_SPX), TIME_FOUND, TIME_NOW , TIME_OLD
      INTEGER           REC_POINTER(N_SPX)
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
!
!  Function subroutines
!
      REAL XFLOW_gx, XFLOW_gy, XFLOW_gz
      REAL VFLOW_gx, VFLOW_gy, VFLOW_gz, MFLOW_gx, MFLOW_gy, MFLOW_gz
      REAL XFLOW_sx, XFLOW_sy, XFLOW_sz
      REAL VFLOW_sx, VFLOW_sy, VFLOW_sz, MFLOW_sx, MFLOW_sy, MFLOW_sz
      REAL FLUX_gx, FLUX_gy, FLUX_gz
      REAL FLUX_sx, FLUX_sy, FLUX_sz
      REAL CALC_RO_g
      INTEGER, EXTERNAL :: FUNIJK_LOC
!
!                   1       2      3         4      5      6
      DATA VAR_DAT/'EP_g', 'P_g', 'P_star', 'U_g', 'V_g', 'W_g', &

!                   7      8      9      10       11     12
                   'U_s', 'V_s', 'W_s', 'ROP_s', 'T_g', 'T_s', &

!                   13      14     15     16          17
                   'T_s2', 'X_g', 'X_s', 'XFLOW_gx', 'XFLOW_gy', &

!                   18          19          20          21
                   'XFLOW_gz', 'XFLOW_sx', 'XFLOW_sy', 'XFLOW_sz', &

!                   22          23          24          25
                   'MFLOW_gx', 'MFLOW_gy', 'MFLOW_gz', 'MFLOW_sx', &

!                   26          27         28           29
                   'MFLOW_sy', 'MFLOW_sz', 'VFLOW_gx', 'VFLOW_gy', &

!                   30          31          32          33
                   'VFLOW_gz', 'VFLOW_sx', 'VFLOW_sy', 'VFLOW_sz', &

!                   34        35        36         37
                   'MASS_g', 'MASS_s', 'FLUX_gx', 'FLUX_gy', &

!                   38         39         40         41
                   'FLUX_gz', 'FLUX_sx', 'FLUX_sy', 'FLUX_sz' ,&

!                   42      43     44    45     46     47
                   'KE_g', 'KE_s','P_s','PE_g','PE_s','BERN_s', &

!                   48          49          50          51
                   'Theta_m', 'Scalar' , 'RRates' , 'K_Turb_G', &

!                   52
                   'E_Turb_G'/

      integer :: gas_species_index , solid_species_index , solid_index
      logical :: bRead_all

      common /fast_sp7/ gas_species_index , solid_species_index , &
                         solid_index , bRead_all

      solid_species_index = 0
      solid_index    = 0
      gas_species_index = 0
      bRead_all = .true.


      CALL READ_RES1
      TIME_IN_RES = TIME
!
      SUM  = .FALSE.
      LMAX = LEN(STRING)
      INTER = .FALSE.
!
      IF (.NOT.DO_XFORMS) THEN
         TIME_START = 0.
         TIME_END   = 0.
         TIME_AVERAGE = .FALSE.
         I1         = 1
         I2         = 1
         I_AVERAGE  = .FALSE.
         J1         = 1
         J2         = 1
         J_AVERAGE  = .FALSE.
         K1         = 1
         K2         = 1
         K_AVERAGE  = .FALSE.
         M = 1
         N = 1
         VAR     = 'EP_g'
      ELSE
         VAR     = VAR_DAT(VAR_NO)
      END IF
!
!
      FILE_NAME  = '*'
      IF (.NOT.DO_XFORMS) THEN
         RES_P_g    = .FALSE.
         RES_T_g    = .FALSE.
         RES_X_g    = .FALSE.
         MINMAX     = -1
      END IF
      DO L = 1, N_SPX
         REC_POINTER(L) = 4
      END DO
!
!
      IF (DO_XFORMS) THEN
        M          = M_USE
        N          = N_USE
        GOTO 5500
      ENDIF
!
!
      WRITE(*,*)
      WRITE(*,*)&
      ' Interactive data retrieval program. Type ? any time for help,'
      WRITE(*,*)&
      ' or press RETURN to select default values shown in parenthesis.'

      WRITE(*,*)

9     write (*,'(A)',ADVANCE='NO') &
         ' Write output using user-supplied precision? (T/F) '
!      read  (*,*) bPrecision
        READ(*,'(1A60)',ERR=9) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(9)
          GOTO 9
        ENDIF
        L3 = 1
        CALL GET_SUBSTR(STRING, L3, SUBSTR)
        IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=9)bPrecision

      nPrec_location = -1  ! default value if not using precision
      nPrec_time     = -1  ! default value if not using precision
      nPrec_variable =  5  ! default value if not using precision

      if (bPrecision) then
         write (*,'(A)',ADVANCE='NO') ' Enter precision for location values: '
         read  (*,*) nPrec_location
         write (*,'(A)',ADVANCE='NO') ' Enter precision for time values: '
         read  (*,*) nPrec_time
      end if
!
!  Read time
!
10    continue

      solid_species_index = 0
      solid_index         = 0
      gas_species_index   = 0
      bRead_all           = .true.

      WRITE(*,*)
      WRITE(*,'(A,F7.3,A,F7.3,A)',ADVANCE='NO')&
       ' Time: (',TIME_START,',',TIME_END,') > '
      READ(*,'(1A60)',ERR=10) STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(10)
        GOTO 10
      ELSEIF(STRING(1:1) .EQ. 'e' .OR. STRING(1:1) .EQ. 'E' .OR. &
             STRING(1:1) .EQ. 'q' .OR. STRING(1:1) .EQ. 'Q') THEN
        IF(FILE_NAME(1:1) .NE. '*') CLOSE(40)
        bRead_all = .true.
        RETURN
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=10)TIME_START
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=10)TIME_END
      IF(TIME_START .LT. ZERO .OR. TIME_END .LT. ZERO) THEN
        IF(FILE_NAME(1:1) .NE. '*') CLOSE(40)
        RETURN
      ENDIF
      IF(TIME_START .GE. TIME_IN_RES) THEN
        TIME_START = TIME_IN_RES
        TIME_END   = TIME_IN_RES
      ENDIF
      IF(TIME_END .LT. TIME_START)GOTO 10
      IF(TIME_START .NE. TIME_END) THEN
        IF(TIME_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
11      WRITE(*, '(A,1A1,A)',ADVANCE='NO')' Time average ? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=11) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(11)
          GOTO 11
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          TIME_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          TIME_AVERAGE = .FALSE.
        ENDIF
      ENDIF
!
!  Read variable name
!
20    CONTINUE
      IF(MINMAX .EQ. 1) THEN
        SUBSTR(1:1) = '1'
        SUBSTR(2:9) = VAR
      ELSEIF(MINMAX .EQ. 0) THEN
        SUBSTR(1:1) = '0'
        SUBSTR(2:9) = VAR
      ELSE
        SUBSTR(1:8) = VAR
        SUBSTR(9:9) = ' '
      ENDIF
      WRITE(*,'(A,1A9,A)',ADVANCE='NO')&
       ' Variable: (', SUBSTR(1:9), ') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(20)
        GOTO 20
      ENDIF
      L3 = 1
!
      IF(STRING(1:1) .EQ. '1') THEN
        MINMAX = 1
        L3 =2
      ELSEIF(STRING(1:1) .EQ. '0')THEN
        MINMAX = 0
        L3 = 2
      ELSEIF(STRING(1:1) .NE. ' ')THEN
        MINMAX = -1
      ENDIF
!
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,'(1A8)',ERR=20) VAR
!
!     Identify variable number
!
      DO 22 L = 1, N_VAR
        IF (STRCMP(VAR,VAR_DAT(L)))THEN
          VAR_NO = L
          GOTO 23
        ENDIF
22    CONTINUE
      WRITE(*,'(A,1A8,A)')' Variable ', VAR, ' not found'
      GOTO 20
23    CONTINUE

!
      IF((VAR_NO .GE.  7 .AND. VAR_NO .LE. 10) .OR. &
         (VAR_NO .EQ. 12                     ) .OR. &
         (VAR_NO .EQ. 15                     ) .OR. &
         (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR. &
         (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. &
         (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR. &
         (VAR_NO .EQ. 35                     ) .OR. &
         (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR. &
         (VAR_NO .EQ. 43                     ) .OR. &
         (VAR_NO .EQ. 44                     ) .OR. &
         (VAR_NO .EQ. 46                     ) .OR. &
         (VAR_NO .EQ. 47                     ) .OR. &
         (VAR_NO .EQ. 48                     ) &
                                                   )THEN
        IF(MMAX .GT. 1) THEN
24        WRITE(*,'(A,I2,A)',ADVANCE='NO') ' Solids phase: (', M, ') > '
          READ(*,'(1A60)',ERR=24) STRING
          IF(STRING(1:1) .EQ. '?') THEN
            CALL HELP(24)
            GOTO 24
          ENDIF
          L3 = 1
          CALL GET_SUBSTR(STRING, L3, SUBSTR)
          IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=24)M
          IF(M .GT. MMAX) THEN
            WRITE(*,*)' Value should not exceed ', MMAX
            M = MMAX
            GOTO 24
          ENDIF
        ELSE
          M = 1
        ENDIF
        solid_index = m
      ENDIF
!
      IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21)THEN
25      WRITE(*,'(A,I2,A)',ADVANCE='NO') ' Species: (', N, ') > '
        READ(*,'(1A60)',ERR=25) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(25)
          GOTO 25
        ENDIF
        L3 = 1
        CALL GET_SUBSTR(STRING, L3, SUBSTR)
        IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=25)N
        IF((VAR_NO .EQ. 14 .OR. VAR_NO .EQ. 16 .OR. VAR_NO .EQ. 17 .OR.&
           VAR_NO .EQ. 18 ) .AND. N .GT. NMAX(0)) THEN
          WRITE(*,*)' Value should not exceed ', NMAX(0)
          N = NMAX(0)
          GOTO 25
        ELSEIF((VAR_NO .EQ. 15 .OR. VAR_NO .EQ. 19 .OR. VAR_NO .EQ. 20&
                .OR. VAR_NO .EQ. 21) .AND. N .GT. NMAX(M)) THEN
          WRITE(*,*)' Value should not exceed ', NMAX(M)
          N = NMAX(M)
          GOTO 25
        ENDIF
        if (var_no .eq. 14) gas_species_index   = N
        if (var_no .eq. 15) solid_species_index = N
      ENDIF
!
      IF(VAR_NO .EQ. 49)THEN
        IF(NScalar .LE. 0) THEN
          write(*,'(i40)')Nscalar
          WRITE(*,'(A)')' No user-defined Scalar found'
          GOTO 20
        ENDIF
        IF(N .GT. NScalar)N = 1
29      WRITE(*,'(A,I2,A)',ADVANCE='NO') ' Scalar: (', N, ') > '
        READ(*,'(1A60)',ERR=29) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(29)
          GOTO 29
        ENDIF
        L3 = 1
        CALL GET_SUBSTR(STRING, L3, SUBSTR)
        IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=29)N
        IF(N .GT. NScalar) THEN
          WRITE(*,*)' Value should not exceed ', NScalar
          N = NScalar
          GOTO 29
        ENDIF
      ENDIF
!
      IF(VAR_NO .EQ. 50)THEN
        IF(nRR .LE. 0) THEN
          write(*,'(i40)')nRR
          WRITE(*,'(A)')' No Reaction Rate data found'
          GOTO 20
        ENDIF
        IF(N .GT. nRR)N = 1
229     WRITE(*,'(A,I2,A)',ADVANCE='NO') ' Reaction Rate: (', N, ') > '
        READ(*,'(1A60)',ERR=229) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(32)
          GOTO 229
        ENDIF
        L3 = 1
        CALL GET_SUBSTR(STRING, L3, SUBSTR)
        IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=229)N
        IF(N .GT. nRR) THEN
          WRITE(*,*)' Value should not exceed ', nRR
          N = nRR
          GOTO 229
        ENDIF
      ENDIF
!
      IF(VAR_NO .EQ. 51)THEN
        IF( .NOT. K_Epsilon) THEN
          WRITE(*,'(A)')' K_Turb_G not found'
          GOTO 20
        ENDIF
      ENDIF
!
      IF(VAR_NO .EQ. 52)THEN
        IF( .NOT. K_Epsilon) THEN
          WRITE(*,'(A)')' E_Turb_G not found'
          GOTO 20
        ENDIF
      ENDIF
!
!
 5500 CONTINUE
!
!
      IF(VAR_NO .EQ.  4 .OR. VAR_NO .EQ.  7 .OR. VAR_NO .EQ. 16 .OR.&
         VAR_NO .EQ. 19 .OR. VAR_NO .EQ. 22 .OR. VAR_NO .EQ. 25 .OR.&
         VAR_NO .EQ. 28 .OR. VAR_NO .EQ. 31 .OR. VAR_NO .EQ. 36 .OR.&
         VAR_NO .EQ. 39)THEN
        DIRECTION = 1
      ELSEIF &
        (VAR_NO .EQ.  5 .OR. VAR_NO .EQ.  8 .OR. VAR_NO .EQ. 17 .OR.&
         VAR_NO .EQ. 20 .OR. VAR_NO .EQ. 23 .OR. VAR_NO .EQ. 26 .OR.&
         VAR_NO .EQ. 29 .OR. VAR_NO .EQ. 32 .OR. VAR_NO .EQ. 37 .OR.&
         VAR_NO .EQ. 40)THEN
        DIRECTION = 2
      ELSEIF&
        (VAR_NO .EQ.  6 .OR. VAR_NO .EQ.  9 .OR. VAR_NO .EQ. 18 .OR.&
         VAR_NO .EQ. 21 .OR. VAR_NO .EQ. 24 .OR. VAR_NO .EQ. 27 .OR.&
         VAR_NO .EQ. 30 .OR. VAR_NO .EQ. 33 .OR. VAR_NO .EQ. 38 .OR.&
         VAR_NO .EQ. 41)THEN
        DIRECTION = 3
      ELSE
         DIRECTION = 0
      ENDIF
!
      IF(VAR_NO .GE. 16 .AND. VAR_NO .LE. 41) THEN
        SUM = .TRUE.
      ELSE
        SUM = .FALSE.
      ENDIF
!
!  Enable the required SPX file
!
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         AT_EOF(L)      = .FALSE.
      END DO

      if (var_no.eq.14 .or. var_no.eq.15) then
          bRead_all = .false.
      else
          bRead_all = .true.
      end if
!
!
      IF(VAR_NO .EQ. 1 .OR. (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.&
         (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR. &
         (VAR_NO .GE. 28 .AND. VAR_NO .LE. 30) .OR.&
          VAR_NO .EQ. 34                       .OR.&
         (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.&
         (VAR_NO .GE. 42 .AND. VAR_NO .LE. 47)    )            THEN
        READ_SPX(1) = .TRUE.    ! EP_g
      ENDIF
      IF(VAR_NO .EQ. 2 .OR. VAR_NO .EQ. 3  .OR.&
         VAR_NO .EQ. 44) THEN
        READ_SPX(2) = .TRUE.    ! P_g, P_star
      ENDIF
      IF((VAR_NO .GE.  4 .AND. VAR_NO .LE.  6) .OR.&
         (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.&
         (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR. &
         (VAR_NO .GE. 28 .AND. VAR_NO .LE. 30) .OR.&
         (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.&
         (VAR_NO .EQ. 42                     ) &
                                                 ) THEN
        READ_SPX(3) = .TRUE.    ! U_g, V_g, W_g
      ENDIF
      IF((VAR_NO .GE.  7 .AND. VAR_NO .LE.  9) .OR.&
         (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.&
         (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. &
         (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.&
         (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR.&
         (VAR_NO .EQ. 43                     ) .OR.&
         (VAR_NO .EQ. 44                     ) .OR.&
         (VAR_NO .EQ. 47                     ) &
                                                   ) THEN
        READ_SPX(4) = .TRUE.    ! U_s, V_s, W_s
      ENDIF
      IF (VAR_NO .EQ. 44) THEN
        READ_SPX(5) = .TRUE.
      END IF
      IF (VAR_NO .EQ. 10 .OR.&
         (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.&
         (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. &
         (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.&
          VAR_NO .EQ. 35                       .OR.&
         (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41)    ) THEN
        IF(MMAX .EQ. 1) THEN
          READ_SPX(1) = .TRUE.    ! EP_g
        ELSE
          READ_SPX(5) = .TRUE.    ! ROP_s
        ENDIF
      ENDIF
      IF(VAR_NO .GE. 11 .AND. VAR_NO .LE. 13) THEN
        READ_SPX(6) = .TRUE.    ! T_g, T_s, T_s2
      ENDIF
!     When the solids density is not constant, it is computed
!     from solids species mass fractions
!     This is needed for volumetric flow rates VFLOW_sx, VFLOW_sy or VFLOW_sz
!     because EP_S is needed and EP_S = ROP_S / RO_S
!     ROP_S is read from the SP5 file, but RO_S is not saved anywhere
      IF((VAR_NO .GE. 14 .AND. VAR_NO .LE. 21) .OR. &
         (ANY_SOLVE_ROs.AND.(VAR_NO .GE. 31 .AND. VAR_NO .LE. 33))) THEN
        READ_SPX(7) = .TRUE.    ! X_g, X_s
      ENDIF

      IF(VAR_NO .EQ. 48 ) THEN
        READ_SPX(8) = .TRUE.    ! Theta_m
      ENDIF

      IF(VAR_NO .EQ. 49 ) THEN
        READ_SPX(9) = .TRUE.    ! Scalar
      ENDIF

      IF(VAR_NO .EQ. 50 ) THEN
        READ_SPX(10) = .TRUE.    ! Reaction Rates
      ENDIF

      IF(VAR_NO .EQ. 51 ) THEN
        READ_SPX(11) = .TRUE.    ! K_Turb_G
      ENDIF

      IF(VAR_NO .EQ. 52 ) THEN
        READ_SPX(11) = .TRUE.    ! E_Turb_G
      ENDIF
!
!  Open P_g, T_g, and X_g files, if gas density needs to be determined
!
      IF (DO_XFORMS) GOTO 1125
      IF( RO_g0 .EQ. UNDEFINED .AND. &
         ( (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.&
           (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR.&
           (VAR_NO .EQ. 34                     ) .OR.&
           (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.&
           (VAR_NO .EQ. 42                     ) .OR.&
           (VAR_NO .EQ. 45                     ) &
                                                    ) ) THEN
          IF (.NOT.DO_XFORMS) THEN
             WRITE(*,*)&
           ' To calculate gas density P_g, T_g, and X_g are needed'
!
           IF(RES_P_g)THEN
             SUBSTR(1:1) = 'Y'
           ELSE
             SUBSTR(1:1) = 'N'
           ENDIF
26         WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
            ' P_g from RES file? (',SUBSTR(1:1),') > '
           READ(*,'(1A60)',ERR=26) STRING
           IF(STRING(1:1) .EQ. '?') THEN
             CALL HELP(26)
             GOTO 26
           ENDIF
           IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
             RES_P_g = .TRUE.
           ELSEIF(STRING(1:1) .NE. ' ')THEN
             RES_P_g = .FALSE.
           ENDIF
!
           IF(RES_T_g)THEN
             SUBSTR(1:1) = 'Y'
           ELSE
             SUBSTR(1:1) = 'N'
           ENDIF
27         WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
             ' T_g from RES file? (',SUBSTR(1:1),') > '
           READ(*,'(1A60)',ERR=27) STRING
           IF(STRING(1:1) .EQ. '?') THEN
             CALL HELP(27)
             GOTO 27
           ENDIF
           IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
             RES_T_g = .TRUE.
           ELSEIF(STRING(1:1) .NE. ' ')THEN
             RES_T_g = .FALSE.
           ENDIF
        ELSE
           CONTINUE
        END IF
!
        RES_X_G = .FALSE.
        IF(MW_avg .EQ. UNDEFINED) THEN
          IF (.NOT.DO_XFORMS) THEN
             IF(RES_X_g)THEN
               SUBSTR(1:1) = 'Y'
             ELSE
               SUBSTR(1:1) = 'N'
             ENDIF
28           WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
               ' X_g from RES file? (',SUBSTR(1:1),') > '
             READ(*,'(1A60)',ERR=28) STRING
             IF(STRING(1:1) .EQ. '?') THEN
               CALL HELP(28)
               GOTO 28
             ENDIF
             IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
               RES_X_g = .TRUE.
             ELSEIF(STRING(1:1) .NE. ' ')THEN
               RES_X_g = .FALSE.
             ENDIF
          ELSE
             CONTINUE
          END IF
        ENDIF
!
!
        IF(.NOT.RES_P_g) READ_SPX(2) = .TRUE.    ! P_g, P_star
        IF(.NOT.RES_X_g) READ_SPX(7) = .TRUE.    ! X_g, X_s
        IF(.NOT.RES_T_g) READ_SPX(6) = .TRUE.    ! T_g, T_s, T_s2
        IF(RES_P_g .OR. RES_T_g .OR. RES_X_g) CALL READ_RES1
!
!
      ENDIF
!
1125  IF (DO_XFORMS) GOTO 5501
!
! if doing user-specifiied precision output, get value for this variable
!
      if (bPrecision) then
          write (*,'(a,i2,a)',ADVANCE='NO') ' Enter precision ( ' , nPrec_variable , ') >'
          read (*,'(1a60)') string
          L3 = 1
          CALL GET_SUBSTR(STRING, L3, SUBSTR)
          IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) nPrec_variable
      end if
!
!  Read I range
!
30    WRITE(*,'(A,I3,A,I3,A)',ADVANCE='NO')&
       ' I range: (', I1, ',', I2, ') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(30)
        GOTO 30
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) I1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) I2
!
!     Check bounds
!
      IF(I2 .LT. I1) I2 = I1
      IF(I1 .LT. 1 .OR. I1 .GT. IMAX2 .OR.&
         I2 .LT. 1 .OR. I2 .GT. IMAX2     ) THEN
        WRITE(*,'(A,I3)')&
          ' I1 and I2 should be in the range 1 to ', IMAX2
        GOTO 30
      ENDIF
!
      IF(I1 .NE. I2) THEN
        IF(I_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
31      WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
        ' Average or sum over I? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=31) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(31)
          GOTO 31
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          I_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          I_AVERAGE = .FALSE.
        ENDIF
      ENDIF
!
!  Read J range
!
40    WRITE(*,'(A,I3,A,I3,A)',ADVANCE='NO')&
       ' J range: (', J1,',',J2, ') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(40)
        GOTO 40
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=40) J1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=40) J2
!
!     Check bounds
!
      IF(J2 .LT. J1) J2 = J1
      IF(J1 .LT. 1 .OR. J1 .GT. JMAX2 .OR.&
         J2 .LT. 1 .OR. J2 .GT. JMAX2     ) THEN
        WRITE(*,'(A,I3)')&
          ' J1 and J2 should be in the range 1 to ', JMAX2
        GOTO 40
      ENDIF
!
      IF(J1 .NE. J2) THEN
        IF(J_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
41      WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
        ' Average or sum over J? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=41) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(41)
          GOTO 41
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          J_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          J_AVERAGE = .FALSE.
        ENDIF
      ENDIF
!
!  Read K range
!
50    WRITE(*,'(A,I3,A,I3,A)',ADVANCE='NO')&
       ' K range: (', K1,',',K2,') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(50)
        GOTO 50
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=50) K1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=50) K2
!
!     Check bounds
!
      IF(K2 .LT. K1) K2 = K1
      IF(K1 .LT. 1 .OR. K1 .GT. KMAX2 .OR.&
         K2 .LT. 1 .OR. K2 .GT. KMAX2     ) THEN
        WRITE(*,'(A,I3)')&
          ' K1 and K2 should be in the range 1 to ', KMAX2
        GOTO 50
      ENDIF
!
      IF(K1 .NE. K2) THEN
        IF(K_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
51      WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
        ' Average or sum over K? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=51) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(51)
          GOTO 51
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          K_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          K_AVERAGE = .FALSE.
        ENDIF
      ENDIF
!
!
 5501 CONTINUE
!
!
!  Read file name
!
!      IF (DO_XFORMS) THEN
!         L3 = INDEX(TEMP_FILE,'*')
!         IF (L3.EQ.0) THEN
!            IF (APPEND_MODE) THEN
!               OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN', &
!                             POSITION='APPEND',CONVERT='BIG_ENDIAN')
!            ELSE
!               OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
!            END IF
!            FILE_NAME(1:1) = 'A'
!         ELSE
!            FILE_NAME(1:1) = '*'
!         END IF
!         GOTO 5502
!      END IF
!
!
70    WRITE(*,'(A,1A30,A)',ADVANCE='NO') ' File: (', FILE_NAME,') >'
      READ(*,'(1A60)') STRING
      IF (STRING(1:1) .EQ. '?') THEN
         CALL HELP(70)
         GOTO 70
      ELSE IF (STRING(1:1) .EQ. '!') THEN
         GOTO 10
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF (SUBSTR(1:1) .NE. ' ')THEN
         IF (.NOT.STRCMP(STRING(1:30),FILE_NAME)) THEN
            IF (FILE_NAME(1:1) .NE. '*') CLOSE(40)
            CALL STREQS(FILE_NAME,STRING(1:30))
            IF (FILE_NAME(1:1) .NE. '*') THEN
               INQUIRE (FILE=FILE_NAME,EXIST=FILE_EXIST)
               IF (FILE_EXIST .AND. .NOT.DO_XFORMS) THEN
                  WRITE(*,'(A)',ADVANCE='NO')' File exists.  Over write? (1=Yes) >'
                  READ(*,*)IANS
                  IF(IANS .NE. 1)GOTO 70
               ENDIF
               OPEN (UNIT=40,FILE=FILE_NAME,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
            ENDIF
         ENDIF
      ENDIF
!
!
 5502 CONTINUE
!
      IF (TIME_START .LT. TIME_IN_RES) THEN
         CALL SEEK_TIME(READ_SPX,TIME_START,REC_POINTER,TIME_FOUND)
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) THEN
               RETURN
            END IF
         END IF
        IF(TIME_FOUND .LT. ZERO) THEN
          WRITE(*,*)' Could not find record for TIME_START'
          GOTO 10
        ENDIF
      ENDIF
!
!  write initial data
!
      CALL WRITE_LINE(FILE_NAME,' ',1)
      CALL WRITE_LINE(FILE_NAME,' ',1)
      DISPLAY = 15
      IF(TIME_START .EQ. TIME_END) THEN
        DISPLAY = DISPLAY - 8
      ENDIF
      IF(I1 .EQ. I2) THEN
        IF(DIRECTION .EQ. 1)THEN
          XTMP = XDIST_VEC(I1)
        ELSE
          XTMP = XDIST_SC(I1)
        ENDIF
        DISPLAY = DISPLAY - 4
        call format_one(spec,nPrec_location,nPrec_length)
        WRITE(LINE,spec)' X = ', XTMP
        CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length + 5)
      ELSEIF(I_AVERAGE) THEN
        call format_two(spec,nPrec_location,nPrec_location,nPrec_length)
        IF(DIRECTION .EQ. 1)THEN
          XTMP = XDIST_VEC(I1)
          IF(SUM) THEN
            WRITE(LINE,spec)&
            ' Sum of values for X = ', XTMP, ' to ',XDIST_VEC(I2)
          ELSE
            WRITE(LINE,spec)&
            ' Average value for X = ', XTMP, ' to ',XDIST_VEC(I2)
          ENDIF
        ELSE
          XTMP = XDIST_SC(I1)
          IF(SUM) THEN
            WRITE(LINE,spec)&
            ' Sum of values for X = ', XTMP, ' to ',XDIST_SC(I2)
          ELSE
            WRITE(LINE,spec)&
            ' Average value for X = ', XTMP, ' to ',XDIST_SC(I2)
          ENDIF
        ENDIF
        DISPLAY = DISPLAY - 4
        CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+27)
      ENDIF
!
      IF(J1 .EQ. J2) THEN
        IF(DIRECTION .EQ. 2)THEN
          YTMP = YDIST_VEC(J1)
        ELSE
          YTMP = YDIST_SC(J1)
        ENDIF
        DISPLAY = DISPLAY - 2
        call format_one(spec,nPrec_location,nPrec_length)
        WRITE(LINE,spec)' Y = ', YTMP
        CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+5)
      ELSEIF(J_AVERAGE) THEN
        call format_two(spec,nPrec_location,nPrec_location,nPrec_length)
        IF(DIRECTION .EQ. 2)THEN
          YTMP = YDIST_VEC(J1)
          IF(SUM) THEN
            WRITE(LINE,spec)&
            ' Sum of values for Y = ', YTMP, ' to ',YDIST_VEC(J2)
          ELSE
            WRITE(LINE,spec)&
            ' Average value for Y = ', YTMP, ' to ',YDIST_VEC(J2)
          ENDIF
        ELSE
          YTMP = YDIST_SC(J1)
          IF(SUM) THEN
            WRITE(LINE,spec)&
            ' Sum of values for Y = ', YTMP, ' to ',YDIST_SC(J2)
          ELSE
            WRITE(LINE,spec)&
            ' Average value for Y = ', YTMP, ' to ',YDIST_SC(J2)
          ENDIF
        ENDIF
        DISPLAY = DISPLAY - 2
        CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+27)
      ENDIF
!
      IF(K1 .EQ. K2) THEN
        IF(DIRECTION .EQ. 3)THEN
          ZTMP = ZDIST_VEC(K1)
        ELSE
          ZTMP = ZDIST_SC(K1)
        ENDIF
        DISPLAY = DISPLAY - 1
        call format_one(spec,nPrec_location,nPrec_length)
        WRITE(LINE,spec)' Z = ', ZTMP
        CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+5)
      ELSEIF(K_AVERAGE) THEN
        call format_two(spec,nPrec_location,nPrec_location,nPrec_length)
        IF(DIRECTION .EQ. 3)THEN
          ZTMP = ZDIST_VEC(K1)
          IF(SUM) THEN
            WRITE(LINE,spec)&
            ' Sum of values for Z = ', ZTMP, ' to ',ZDIST_VEC(K2)
          ELSE
            WRITE(LINE,spec)&
            ' Average value for Z = ', ZTMP, ' to ',ZDIST_VEC(K2)
          ENDIF
        ELSE
          ZTMP = ZDIST_SC(K1)
          IF(SUM) THEN
            WRITE(LINE,spec)&
            ' Sum of values for Z = ', ZTMP, ' to ',ZDIST_SC(K2)
          ELSE
            WRITE(LINE,spec)&
            ' Average value for Z = ', ZTMP, ' to ',ZDIST_SC(K2)
          ENDIF
        ENDIF
        DISPLAY = DISPLAY - 1
        CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+27)
      ENDIF
!
      IF(DISPLAY .EQ. 8 .AND. .NOT.TIME_AVERAGE)THEN
        WRITE(LINE,'(5X,A,10X,1A8)')'Time', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 27)
      ENDIF
      IF((VAR_NO .GE.  7 .AND. VAR_NO .LE. 10) .OR.&
         (VAR_NO .EQ. 12                     ) .OR.&
         (VAR_NO .EQ. 13                     ) .OR.&
         (VAR_NO .EQ. 15                     ) .OR.&
         (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.&
         (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR.&
         (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.&
         (VAR_NO .EQ. 35                     ) .OR.&
         (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR.&
         (VAR_NO .EQ. 43                     ) .OR.&
         (VAR_NO .EQ. 44                     ) .OR.&
         (VAR_NO .EQ. 46                     ) .OR.&
         (VAR_NO .EQ. 47                     ) .OR.&
         (VAR_NO .EQ. 48                     ) &
                                                   )THEN
        WRITE(LINE,'(A,I2)') ' Solids phase = ', M
        CALL WRITE_LINE(FILE_NAME, LINE, 18)
      ENDIF
      IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21) THEN
        WRITE(LINE,'(A,I2)') ' Species = ', N
        CALL WRITE_LINE(FILE_NAME, LINE, 13)
      ENDIF
      IF(VAR_NO .EQ. 50) THEN
        WRITE(LINE,'(A,I2)') ' Rrates = ', N
        CALL WRITE_LINE(FILE_NAME, LINE, 12)
      ENDIF
!
!  Read data
!
      END_AVERAGE = .FALSE.
      TIME_OLD = -1.
      NT = 0
      IF(TIME_AVERAGE) THEN
        DO K = K1, K2
           DO J = J1, J2
              DO I = I1, I2
                 IJK = FUNIJK_LOC(I,J,K)
                 VALUE(IJK) = ZERO
              ENDDO
           ENDDO
        ENDDO
      ENDIF
!
      IF(MINMAX .EQ. 0) THEN
        WRITE(LINE,'(6X,A, 10X, A,13X,A,13X,A,5X,A,1A8)')&
            'Time', 'X', 'Y', 'Z', 'Minimum ', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 70)
      ELSEIF(MINMAX .EQ. 1) THEN
        WRITE(LINE,'(6X,A, 10X, A,13X,A,13X,A,5X,A,1A8)')&
            'Time', 'X', 'Y', 'Z', 'Maximum ', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 70)
      ENDIF
!
100   CONTINUE
      IF (TIME_START .LT. TIME_IN_RES) THEN
         CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                          TIME_NOW, TIME_REAL, NSTEP_1)
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) THEN
               RETURN
            END IF
         END IF
      ELSE
         CALL READ_RES1
         TIME_NOW = TIME_IN_RES
      ENDIF
!
      IF (TIME_NOW .LT. ZERO) THEN
        IF(.NOT.TIME_AVERAGE) THEN
           IF (.NOT.DO_XFORMS) THEN
              GOTO 10
           ELSE
              IF(FILE_NAME(1:1) .NE. '*') CLOSE(40)
              RETURN
           END IF
        END IF
        TIME_NOW = TIME_OLD
        END_AVERAGE = .TRUE.
        GOTO 106
      ENDIF
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      TIME_OLD = TIME_NOW
!
      IF(MMAX .EQ. 1.AND.(.NOT.ANY_SOLVE_ROs)) THEN
        IF (VAR_NO .EQ. 10 .OR.&
         (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.&
         (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. &
         (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.&
         (VAR_NO .EQ. 35                     ) .OR.&
         (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR. &
         (VAR_NO .EQ. 43                     ) .OR.&
         (VAR_NO .EQ. 44                     ) .OR.&
         (VAR_NO .EQ. 46                     ) .OR.&
         (VAR_NO .EQ. 47                     ) &
                                                  ) THEN
          ! loop over the entire domain because mass flux calculations
          ! need ROP_s outside i,j,k limits specified by the user
          DO K = KMIN1, KMAX2
             DO J = JMIN1, JMAX2
                DO I = IMIN1, IMAX2
                   IJK = FUNIJK_LOC(I, J, K)
                   ROP_s(IJK, 1) = (ONE - EP_g(IJK)) * RO_s0(1)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
      ENDIF
!
!     FIND THETA IF CALCULATING SOLIDS PRESSURE, P_s
!
      IF(VAR_NO .EQ. 44 .OR. VAR_NO .EQ. 47) THEN
        M_LOCAL=M
        DO M = 1, MMAX
          CALL CALC_MU_s(M, IER)
        ENDDO
        M=M_LOCAL
      ENDIF
!
      NT = NT + 1
      DO K = K1, K2
         DO J = J1, J2
            DO I = I1, I2
               IJK = FUNIJK_LOC(I, J, K)
               IF(VAR_NO .EQ. 1) THEN
                  VALUE_TMP = EP_g(IJK)
               ELSEIF(VAR_NO .EQ. 2)THEN
                  VALUE_TMP = P_g(IJK)
               ELSEIF(VAR_NO .EQ. 3)THEN
                  VALUE_TMP = P_star(IJK)
               ELSEIF(VAR_NO .EQ. 4)THEN
                  VALUE_TMP = U_g(IJK)
               ELSEIF(VAR_NO .EQ. 5)THEN
                  VALUE_TMP = V_g(IJK)
               ELSEIF(VAR_NO .EQ. 6)THEN
                  VALUE_TMP = W_g(IJK)
               ELSEIF(VAR_NO .EQ. 7)THEN
                  VALUE_TMP = U_s(IJK, M)
               ELSEIF(VAR_NO .EQ. 8)THEN
                  VALUE_TMP = V_s(IJK, M)
               ELSEIF(VAR_NO .EQ. 9)THEN
                  VALUE_TMP = W_s(IJK, M)
               ELSEIF(VAR_NO .EQ. 10)THEN
                  VALUE_TMP = ROP_s(IJK, M)
               ELSEIF(VAR_NO .EQ. 11)THEN
                  VALUE_TMP = T_g(IJK)
               ELSEIF(VAR_NO .EQ. 12)THEN
                  VALUE_TMP = T_s(IJK, M)
               ELSEIF(VAR_NO .EQ. 13)THEN
                  VALUE_TMP = T_s(IJK, 2)
               ELSEIF(VAR_NO .EQ. 14)THEN
                  VALUE_TMP = X_g(IJK, N)
               ELSEIF(VAR_NO .EQ. 15)THEN
                  VALUE_TMP = X_s(IJK, M, N)
               ELSEIF(VAR_NO .EQ. 16)THEN
                  VALUE_TMP = XFLOW_gx(I, J, K, IJK, N)
               ELSEIF(VAR_NO .EQ. 17)THEN
                  VALUE_TMP = XFLOW_gy(I, J, K, IJK, N)
               ELSEIF(VAR_NO .EQ. 18)THEN
                  VALUE_TMP = XFLOW_gz(I, J, K, IJK, N)
               ELSEIF(VAR_NO .EQ. 19)THEN
                  VALUE_TMP = XFLOW_sx(I, J, K, IJK, M, N)
               ELSEIF(VAR_NO .EQ. 20)THEN
                  VALUE_TMP = XFLOW_sy(I, J, K, IJK, M, N)
               ELSEIF(VAR_NO .EQ. 21)THEN
                  VALUE_TMP = XFLOW_sz(I, J, K, IJK, M, N)
               ELSEIF(VAR_NO .EQ. 22)THEN
                  VALUE_TMP = MFLOW_gx(I, J, K, IJK)
               ELSEIF(VAR_NO .EQ. 23)THEN
                  VALUE_TMP = MFLOW_gy(I, J, K, IJK)
               ELSEIF(VAR_NO .EQ. 24)THEN
                  VALUE_TMP = MFLOW_gz(I, J, K, IJK)
               ELSEIF(VAR_NO .EQ. 25)THEN
                  VALUE_TMP = MFLOW_sx(I, J, K, IJK, M)
               ELSEIF(VAR_NO .EQ. 26)THEN
                  VALUE_TMP = MFLOW_sy(I, J, K, IJK, M)
               ELSEIF(VAR_NO .EQ. 27)THEN
                  VALUE_TMP = MFLOW_sz(I, J, K, IJK, M)
               ELSEIF(VAR_NO .EQ. 28)THEN
                  VALUE_TMP = VFLOW_gx(I, J, K, IJK)
               ELSEIF(VAR_NO .EQ. 29)THEN
                  VALUE_TMP = VFLOW_gy(I, J, K, IJK)
               ELSEIF(VAR_NO .EQ. 30)THEN
                  VALUE_TMP = VFLOW_gz(I, J, K, IJK)
               ELSEIF(VAR_NO .EQ. 31)THEN
                  VALUE_TMP = VFLOW_sx(I, J, K, IJK, M)
               ELSEIF(VAR_NO .EQ. 32)THEN
                  VALUE_TMP = VFLOW_sy(I, J, K, IJK, M)
               ELSEIF(VAR_NO .EQ. 33)THEN
                  VALUE_TMP = VFLOW_sz(I, J, K, IJK, M)
               ELSEIF(VAR_NO .EQ. 34)THEN
                  VALUE_TMP = EP_g(IJK) * CALC_RO_g(IJK) * VOL(IJK)
               ELSEIF(VAR_NO .EQ. 35)THEN
                  VALUE_TMP = ROP_s(IJK, M) * VOL(IJK)
               ELSEIF(VAR_NO .EQ. 36)THEN
                  VALUE_TMP = FLUX_gx(IJK)
               ELSEIF(VAR_NO .EQ. 37)THEN
                  VALUE_TMP = FLUX_gy(IJK)
               ELSEIF(VAR_NO .EQ. 38)THEN
                  VALUE_TMP = FLUX_gz(IJK)
               ELSEIF(VAR_NO .EQ. 39)THEN
                  VALUE_TMP = FLUX_sx(IJK, M)
               ELSEIF(VAR_NO .EQ. 40)THEN
                  VALUE_TMP = FLUX_sy(IJK, M)
               ELSEIF(VAR_NO .EQ. 41)THEN
                  VALUE_TMP = FLUX_sz(IJK, M)
               ELSEIF(VAR_NO .EQ. 42)THEN
                  VALUE_TMP = 0.5*CALC_RO_g(IJK)*EP_g(IJK)&
                       *(U_g(IJK)**2 + V_g(IJK)**2 + W_g(IJK)**2)
               ELSEIF(VAR_NO .EQ. 43)THEN
                  VALUE_TMP = 0.5*ROP_s(IJK,M)&
                       *(U_s(IJK,M)**2 + V_s(IJK,M)**2 + W_s(IJK,M)**2)
               ELSEIF(VAR_NO .EQ. 44)THEN
                  VALUE_TMP = P_s(IJK,M)
                  IF(EP_g(IJK) .LT. EP_star) THEN
                     VALUE_TMP = P_star(IJK)
                  ENDIF
               ELSEIF(VAR_NO .EQ. 45)THEN
                  VALUE_TMP = GRAVITY*YDIST_SC(J)*CALC_RO_g(IJK)*EP_g(IJK)
               ELSEIF(VAR_NO .EQ. 46)THEN
                  VALUE_TMP = GRAVITY*YDIST_SC(J)*ROP_s(IJK,M)
               ELSEIF(VAR_NO .EQ. 47)THEN
                  mIJK= FUNIJK_LOC(I,J+1,K)
                  lIJK= FUNIJK_LOC(I,J-1,K)
                  DELm=YDIST_SC(J+1)-YDIST_SC(J)
                  DELl=YDIST_SC(J)-YDIST_SC(J-1)
                  FAC1= (ROP_s(mIJK,M)-ROP_s(IJK,M))/DELm &
                       +(ROP_s(IJK,M)-ROP_s(lIJK,M))/DELl
                  FAC2=P_s(IJK,M)/(ROP_s(IJK,M)**2)
                  !          IF(EP_g(IJK) .LT. EP_star) THEN
                  !             FAC2=P_star(IJK)/(ROP_s(IJK,M)**2)
                  !          ENDIF
                  VALUE_TMP = DY(J)*FAC1*FAC2/2
               ELSEIF(VAR_NO .EQ. 48)THEN
                  VALUE_TMP = Theta_m(IJK, M)
               ELSEIF(VAR_NO .EQ. 49)THEN
                  VALUE_TMP = Scalar(IJK, N)
               ELSEIF(VAR_NO .EQ. 50)THEN
                  VALUE_TMP = ReactionRates(IJK, N)
               ELSEIF(VAR_NO .EQ. 51)THEN
                  VALUE_TMP = K_Turb_G(IJK)
               ELSEIF(VAR_NO .EQ. 52)THEN
                  VALUE_TMP = E_Turb_G(IJK)
               ENDIF

               IF(TIME_AVERAGE)THEN
                  VALUE(IJK) = VALUE(IJK) + VALUE_TMP
               ELSE
                  VALUE(IJK) = VALUE_TMP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

106   IF(TIME_AVERAGE)THEN
        IF(TIME_NOW .GE. TIME_END .OR. END_AVERAGE)THEN
          IF(NT .EQ. 0) THEN
            WRITE(*,*)' Could not do time averaging'
            GOTO 10
          ENDIF
          DO K = K1, K2
             DO J = J1, J2
                DO I = I1, I2
                   IJK = FUNIJK_LOC(I,J,K)
                   VALUE(IJK) = VALUE(IJK) / REAL(NT)
                ENDDO
             ENDDO
          ENDDO
          call format_twoB(spec,nPrec_time,nPrec_time,nPrec_length)
          WRITE(LINE,spec) &
           'Time average of ',VAR, ' from Time = ',TIME_START, &
           ' to ', TIME_NOW
          CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+33)
        ELSE
          GOTO 100
        ENDIF
      ENDIF
!
!  DO I, J, or K averaging
!
      IF(K_AVERAGE)THEN
        K = K1
        DO J = J1, J2
           DO I = I1, I2
              IJK = FUNIJK_LOC(I, J, K)
              IF(WALL_AT(IJK))THEN
                 VALUE(IJK) = ZERO
                 DIST(IJK)  = ZERO
              ELSEIF(.NOT. SUM) THEN
                 IF(DIRECTION .EQ. 3)THEN
                    VALUE(IJK) = VALUE(IJK) * DZ_T(K)
                    DIST(IJK)  = DZ_T(K)
                 ELSE
                    VALUE(IJK) = VALUE(IJK) * DZ(K)
                    DIST(IJK)  = DZ(K)
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
!
        DO K = K1+1, K2
           DO J = J1, J2
              DO I = I1, I2
                 IJK = FUNIJK_LOC(I, J, K)
                 IJK1 = FUNIJK_LOC(I, J, K1)
                 IF(.NOT.WALL_AT(IJK)) THEN
                    IF(SUM) THEN
                       VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
                    ELSE
                       IF(DIRECTION .EQ. 3)THEN
                          VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DZ_T(K)
                          DIST(IJK1)  = DIST(IJK1)  + DZ_T(K)
                       ELSE
                          VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DZ(K)
                          DIST(IJK1)  = DIST(IJK1)  + DZ(K)
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        IF(.NOT. SUM) THEN
          K = K1
          DO J = J1, J2
             DO I = I1, I2
                IJK = FUNIJK_LOC(I, J, K)
                IF(DIST(IJK) .NE. ZERO) THEN
                   VALUE(IJK) = VALUE(IJK) / DIST(IJK)
                ELSEIF(VALUE(IJK) .NE. ZERO) THEN
                   WRITE(*,*)' Error in K-averaging'
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       K2d = K1
    ELSE
       K2d = K2
    ENDIF

    IF(J_AVERAGE)THEN
       J = J1
       DO K = K1, K2d
          DO I = I1, I2
             IJK = FUNIJK_LOC(I, J, K)
             IF(WALL_AT(IJK) .AND. .NOT.K_AVERAGE)THEN
                VALUE(IJK) = ZERO
                DIST(IJK)  = ZERO
             ELSEIF(.NOT. SUM) THEN
                IF(DIRECTION .EQ. 2)THEN
                   VALUE(IJK) = VALUE(IJK) * DY_N(J)
                   DIST(IJK)  = DY_N(J)
                ELSE
                   VALUE(IJK) = VALUE(IJK) * DY(J)
                   DIST(IJK)  = DY(J)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       DO K = K1, K2d
          DO J = J1+1, J2
             DO I = I1, I2
                IJK = FUNIJK_LOC(I, J, K)
                IJK1 = FUNIJK_LOC(I, J1, K)
                IF(.NOT.WALL_AT(IJK) .OR. K_AVERAGE) THEN
                   IF(SUM) THEN
                      VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
                   ELSE
                      IF(DIRECTION .EQ. 2)THEN
                         VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DY_N(J)
                         DIST(IJK1)  = DIST(IJK1)  + DY_N(J)
                      ELSE
                         VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DY(J)
                         DIST(IJK1)  = DIST(IJK1)  + DY(J)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       IF(.NOT.SUM) THEN
          J = J1
          DO K = K1, K2d
             DO I = I1, I2
                IJK = FUNIJK_LOC(I, J, K)
                IF(DIST(IJK) .NE. ZERO) THEN
                   VALUE(IJK) = VALUE(IJK) / DIST(IJK)
                ELSEIF(VALUE(IJK) .NE. ZERO) THEN
                   WRITE(*,*)' Error in J-averaging'
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       J2d = J1
    ELSE
       J2d = J2
    ENDIF

    IF(I_AVERAGE)THEN
       I = I1
       DO K = K1, K2d
          DO J = J1, J2d
             IJK = FUNIJK_LOC(I, J, K)
             IF(WALL_AT(IJK) .AND. .NOT.J_AVERAGE .AND. .NOT.K_AVERAGE)THEN
                VALUE(IJK) = ZERO
                DIST(IJK)  = ZERO
             ELSEIF(.NOT. SUM) THEN
                IF(DIRECTION .EQ. 1)THEN
                   VALUE(IJK) = VALUE(IJK) * DX_E(I) * X_E(I)
                   DIST(IJK)  = DX_E(I) * X_E(I)
                ELSE
                   VALUE(IJK) = VALUE(IJK) * DX(I) * X(I)
                   DIST(IJK)  = DX(I) * X(I)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       DO K = K1, K2d
          DO J = J1, J2d
             DO I = I1+1, I2
                IJK = FUNIJK_LOC(I, J, K)
                IJK1 = FUNIJK_LOC(I1, J, K)
                IF(.NOT.WALL_AT(IJK) .OR. J_AVERAGE .OR. K_AVERAGE) THEN
                   IF(SUM) THEN
                      VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
                   ELSE
                      IF(DIRECTION .EQ. 1)THEN
                         VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DX_E(I)* X_E(I)
                         DIST(IJK1)  = DIST(IJK1)  + DX_E(I) * X_E(I)
                      ELSE
                         VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DX(I) * X(I)
                         DIST(IJK1)  = DIST(IJK1)  + DX(I) * X(I)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       IF(.NOT.SUM) THEN
          I = I1
          DO K = K1, K2d
             DO J = J1, J2d
                IJK = FUNIJK_LOC(I, J, K)
                IF(DIST(IJK) .NE. ZERO) THEN
                   VALUE(IJK) = VALUE(IJK) / DIST(IJK)
                ELSEIF(VALUE(IJK) .NE. ZERO) THEN
                   WRITE(*,*)' Error in I-averaging'
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       I2d = I1
    ELSE
       I2d = I2
    ENDIF
!
!  Display data or write data to file
!
    IF(MINMAX .GE. 0) THEN
       IF(MINMAX .EQ. 0)THEN
          VALUE_TMP = 1E32
          DO K = K1, K2d
             DO J = J1, J2d
                DO I = I1, I2d
                   IJK = FUNIJK_LOC(I,J,K)
                   IF(WALL_AT(IJK) .AND. &
                        .NOT. (I_AVERAGE .OR. J_AVERAGE .OR. K_AVERAGE) ) CYCLE
                   IF(VALUE(IJK) .GE. VALUE_TMP) CYCLE
                   XTMP = XDIST_SC(I)
                   YTMP = YDIST_SC(J)
                   ZTMP = ZDIST_SC(K)
                   IF(DIRECTION .EQ. 1) THEN
                      XTMP = XDIST_VEC(I)
                   ELSEIF(DIRECTION .EQ. 2) THEN
                      YTMP = YDIST_VEC(J)
                   ELSEIF(DIRECTION .EQ. 3) THEN
                      ZTMP = ZDIST_VEC(K)
                   ENDIF
                   VALUE_TMP = VALUE(IJK)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF(MINMAX .EQ. 1)THEN
          VALUE_TMP = -1E32
          DO K = K1, K2d
             DO J = J1, J2d
                DO I = I1, I2d
                   IJK = FUNIJK_LOC(I,J,K)
                   IF(WALL_AT(IJK) .AND. &
                        .NOT. (I_AVERAGE .OR. J_AVERAGE .OR. K_AVERAGE) ) CYCLE
                   IF(VALUE(IJK) .LE. VALUE_TMP) CYCLE
                   XTMP = XDIST_SC(I)
                   YTMP = YDIST_SC(J)
                   ZTMP = ZDIST_SC(K)
                   IF(DIRECTION .EQ. 1) THEN
                      XTMP = XDIST_VEC(I)
                   ELSEIF(DIRECTION .EQ. 2) THEN
                      YTMP = YDIST_VEC(J)
                   ELSEIF(DIRECTION .EQ. 3) THEN
                      ZTMP = ZDIST_VEC(K)
                   ENDIF
                   VALUE_TMP = VALUE(IJK)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       call format_five(spec,nPrec_time,nPrec_location,nPrec_location, &
            nPrec_location,nPrec_variable,nPrec_length)
       WRITE(LINE,spec)TIME_NOW, XTMP, YTMP, ZTMP, &
            VALUE_TMP
       CALL WRITE_LINE(FILE_NAME, LINE, nPrec_length)
    ELSEIF(DISPLAY .EQ. 8) THEN
       IJK = FUNIJK_LOC(I1, J1, K1)
       IF(TIME_AVERAGE) THEN
          call format_oneB(spec,nPrec_variable,nPrec_length)
          WRITE(LINE,spec)VALUE(IJK)
          CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length)
       ELSE
          call format_twoC(spec,nPrec_time,nPrec_variable,nPrec_length)
          WRITE(LINE,spec)TIME_NOW, VALUE(IJK)
          CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length)
       ENDIF
    ELSEIF(DISPLAY .EQ. 4 .OR. DISPLAY .EQ. 12) THEN
       call format_one(spec,nPrec_time,nPrec_length)
       WRITE(LINE,spec)' Time = ',TIME_NOW
       CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length + 8)
       WRITE(LINE,'(6X,A,13X,1A8)')'X', VAR
       CALL WRITE_LINE(FILE_NAME, LINE, 28)
       DO I = I1, I2d
          IJK = FUNIJK_LOC(I, J1, K1)
          call format_twoC(spec,nPrec_location,nPrec_variable,nPrec_length)
          IF(DIRECTION .EQ. 1)THEN
             WRITE(LINE,spec)XDIST_VEC(I), VALUE(IJK)
          ELSE
             WRITE(LINE,spec)XDIST_SC(I), VALUE(IJK)
          ENDIF
          CALL WRITE_LINE(FILE_NAME, LINE, nPrec_length)
       ENDDO
    ELSEIF(DISPLAY .EQ. 2 .OR. DISPLAY .EQ. 10) THEN
       call format_one(spec,nPrec_time,nPrec_length)
       WRITE(LINE,spec)' Time = ',TIME_NOW
       CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+8)
       WRITE(LINE,'(6X,A,13X,1A8)')'Y', VAR
       CALL WRITE_LINE(FILE_NAME, LINE, 28)
       DO J = J1, J2d
          IJK = FUNIJK_LOC(I1, J, K1)
          call format_twoC(spec,nPrec_location,nPrec_variable,nPrec_length)
          IF(DIRECTION .EQ. 2)THEN
             WRITE(LINE,spec)YDIST_VEC(J), VALUE(IJK)
          ELSE
             WRITE(LINE,spec)YDIST_SC(J), VALUE(IJK)
          ENDIF
          CALL WRITE_LINE(FILE_NAME, LINE, nPrec_length)
       ENDDO
    ELSEIF(DISPLAY .EQ. 1 .OR. DISPLAY .EQ. 9) THEN
       call format_one(spec,nPrec_time,nPrec_length)
       WRITE(LINE,spec)' Time = ',TIME_NOW
       CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+8)
       WRITE(LINE,'(6X,A,13X,1A8)')'Z', VAR
       CALL WRITE_LINE(FILE_NAME, LINE, 28)
       DO K = K1, K2d
          IJK = FUNIJK_LOC(I1, J1, K)
          call format_twoC(spec,nPrec_location,nPrec_variable,nPrec_length)
          IF(DIRECTION .EQ. 3)THEN
             WRITE(LINE,spec)ZDIST_VEC(K), VALUE(IJK)
          ELSE
             WRITE(LINE,spec)ZDIST_SC(K), VALUE(IJK)
          ENDIF
          CALL WRITE_LINE(FILE_NAME, LINE, nPrec_length)
       ENDDO
    ELSEIF(DISPLAY .EQ. 0) THEN
       call format_one(spec,nPrec_time,nPrec_length)
       WRITE(LINE,spec)' Time = ',TIME_NOW
       CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+8)
       IJK = FUNIJK_LOC(I1, J1, K1)
       call format_oneC(spec,nPrec_variable,nPrec_length)
       WRITE(LINE,spec)VAR,' = ',VALUE(IJK)
       CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length + 3)
    ELSE
       call format_one(spec,nPrec_time,nPrec_length)
       WRITE(LINE,spec)' Time = ',TIME_NOW
       CALL WRITE_LINE(FILE_NAME,LINE,nPrec_length+8)
       WRITE(LINE,'(6X,A,13X,A,13X,A,13X,1A8)')'X','Y','Z',VAR
       CALL WRITE_LINE(FILE_NAME, LINE, 56)
       DO K = K1, K2d
          DO J = J1, J2d
             DO I = I1, I2d
                IJK = FUNIJK_LOC(I,J,K)
                XTMP = XDIST_SC(I)
                YTMP = YDIST_SC(J)
                ZTMP = ZDIST_SC(K)
                IF(DIRECTION .EQ. 1) THEN
                   XTMP = XDIST_VEC(I)
                ELSEIF(DIRECTION .EQ. 2) THEN
                   YTMP = YDIST_VEC(J)
                ELSEIF(DIRECTION .EQ. 3) THEN
                   ZTMP = ZDIST_VEC(K)
                ENDIF
                call format_four(spec,nPrec_location,nPrec_location, &
                     nPrec_location, &
                     nPrec_variable,nPrec_length)

                WRITE(LINE,spec)XTMP, YTMP, ZTMP, VALUE(IJK)
                CALL WRITE_LINE(FILE_NAME, LINE, nPrec_length)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF (TIME_AVERAGE .AND. END_AVERAGE) THEN
       IF (DO_XFORMS) THEN
          IF (FILE_NAME(1:1) .NE. '*') CLOSE(40)
          RETURN
       ELSE
          IF (FILE_NAME(1:1) .NE. '*') then
             close (40)
             open (unit=40,file=file_name,position='append',convert='big_endian')
          end if
          GOTO 10
       END IF
    END IF

    IF(TIME_NOW .GE. TIME_END) THEN
       IF (DO_XFORMS) THEN
          IF (FILE_NAME(1:1) .NE. '*') CLOSE(UNIT=40)
          RETURN
       ELSE
          IF (FILE_NAME(1:1) .NE. '*') then
             close (40)
             open (unit=40,file=file_name,position='append',convert='big_endian')
          end if
          GOTO 10
       END IF
    END IF

    GOTO 100

    CONTAINS

      DOUBLE PRECISION FUNCTION DZ_T(K)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: K
        DZ_T = HALF * (DZ(K) + DZ(KP1(K)))
      END FUNCTION DZ_T

      DOUBLE PRECISION FUNCTION DY_N(J)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: J
        DY_N = HALF * (DY(J) + DY(JP1(J)))
      END FUNCTION DY_N

      DOUBLE PRECISION FUNCTION DX_E(I)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: I
        DX_E = HALF * (DX(I) + DX(IP1(I)))
      END FUNCTION DX_E

      END

      ! This routine appears to need its own version of the funijk routine
      INTEGER FUNCTION FUNIJK_LOC(LI, LJ, LK)
        USE geometry
	IMPLICIT NONE
        INTEGER, INTENT(IN) :: LI, LJ, LK
        FUNIJK_LOC = LI + (LJ-1)*IMAX2 + (LK-1)*IJMAX2
      END FUNCTION FUNIJK_LOC
!
      SUBROUTINE WRITE_LINE(FILE_NAME,LINE,NCHARS)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: FILE_NAME, LINE
      CHARACTER(LEN=81):: LINE2
      INTEGER       NCHARS
!
      INCLUDE 'xforms.inc'
!
      IF (DO_XFORMS) THEN
         LINE2 = LINE
         LINE2(NCHARS+1:NCHARS+1) = CHAR(0)
         CALL ADD_TO_RESULTS_BROWSER(LINE2(1:NCHARS+1))
      END IF
      IF(FILE_NAME(1:1) .EQ. '*')THEN
        IF (.NOT.DO_XFORMS) WRITE(*,*) LINE(1:NCHARS)
      ELSE
        WRITE(40,*) LINE(1:NCHARS)
      ENDIF
      RETURN
      END
!
      SUBROUTINE HELP(N)
      IMPLICIT NONE
      INTEGER N
!
      WRITE(*,*)
!
      IF(N .EQ. 9)THEN
        WRITE(*,*) &
        ' If true, then user is asked for the number of digits of precision'
        WRITE(*,*) &
        ' for the location, time and variable.'
      ELSEIF(N .EQ. 10)THEN
        WRITE(*,*) &
        ' Enter start time and end time for data retrieval.  If the'
        WRITE(*,*) &
        ' time entered is negative, control returns to the main menu.'
      ELSEIF(N .EQ. 11)THEN
        WRITE(*,*) &
        ' Enter Y or y for time averaging.'
      ELSEIF(N .EQ. 20)THEN
        WRITE(*,'(42(A,/))') &
        ' Valid variable names:', &
        '   EP_g      - Void fraction', &
        '   P_g       - Gas pressure, dyne/cm^2', &
        '   P_star    - Solids pressure (frictional regime), dyne/cm^2', &
        '   Gas velocity, cm/s', &
        '   U_g       - X component', &
        '   V_g       - Y component', &
        '   W_g       - Z component', &
        '   Solids velocity, cm/s', &
        '   U_s       - X component', &
        '   V_s       - Y component', &
        '   W_s       - Z component', &
        '   ROP_s     - Solids density x volume fraction, g/cm^3', &
        '   T_g       - Gas temperature, K', &
        '   T_s       - Solids temperature, K', &
        '   Theta_m   - Granular temperature, cm^2/s^2', &
        '   X_g       - Gas species mass fraction', &
        '   X_s       - Solids species mass fraction', &
        '   Scalar    - User defined scalar', &
        '   RRates    - Reaction Rate data (SPA)' , &
        '   Gas species mass flow rates, g/s', &
        '   K_Turb_G    - Turbulent kinetic energy', &
        '   E_Turb_G    - Dissipation of turbulence', &
        '   XFLOW_gx  - in X direction', &
        '   XFLOW_gy  - in Y direction', &
        '   XFLOW_gz  - in Z direction', &
        '   Solids species mass flow rates, g/s', &
        '   XFLOW_sx  - in X direction', &
        '   XFLOW_sy  - in Y direction', &
        '   XFLOW_sz  - in Z direction', &
        '   Gas total mass flow rates, g/s', &
        '   MFLOW_gx  - in X direction', &
        '   MFLOW_gy  - in Y direction', &
        '   MFLOW_gz  - in Z direction', &
        '   Solids total mass flow rates, g/s', &
        '   MFLOW_sx  - in X direction', &
        '   MFLOW_sy  - in Y direction', &
        '   MFLOW_sz  - in Z direction', &
        '   Gas total volumetric flow rates, cm^3/s', &
        '   VFLOW_gx  - in X direction', &
        '   VFLOW_gy  - in Y direction', &
        '   VFLOW_gz  - in Z direction', &
        '   Solids total volumetric flow rates, cm^3/s', &
        '   VFLOW_sx  - in X direction', &
        '   VFLOW_sy  - in Y direction', &
        '   VFLOW_sz  - in Z direction', &
        '   Gas total mass flux, g/(cm^2.s)',  &
        '   FLUX_gx   - in X direction', &
        '   FLUX_gy   - in Y direction', &
        '   FLUX_gz   - in Z direction', &
        '   Solids total mass flux, g/(cm^2.s)',  &
        '   FLUX_sx   - in X direction', &
        '   FLUX_sy   - in Y direction', &
        '   FLUX_sz   - in Z direction', &
        '   MASS_g    - Total mass of gas in specified volume, g', &
        '   MASS_s    - Total mass of solids in specified volume, g', &
        ' To determine the minimum or maximum preceed the variable', &
        ' name with a 0 or 1.'
      ELSEIF(N .EQ. 24)THEN
        WRITE(*,*)' Enter the solids phase number'
      ELSEIF(N .EQ. 25)THEN
        WRITE(*,*)' Enter the species number'
      ELSEIF(N .EQ. 26)THEN
        WRITE(*,*)' Enter Y or y to use P_g values from the RES file'
        WRITE(*,*)' The values may differ from current values of P_g'
      ELSEIF(N .EQ. 27)THEN
        WRITE(*,*)' Enter Y or y to use T_g values from the RES file'
        WRITE(*,*)' The values may differ from current values of T_g'
      ELSEIF(N .EQ. 28)THEN
        WRITE(*,*)' Enter Y or y to use X_g values from the RES file'
        WRITE(*,*)' The values may differ from current values of X_g'
      ELSEIF(N .EQ. 29)THEN
        WRITE(*,*)' Enter the scalar number'
      ELSEIF(N .EQ. 30)THEN
        WRITE(*,*)' Enter the starting and ending values of I'
      ELSEIF(N .EQ. 31)THEN
        WRITE(*,*)' Intensive variables, such as EP_g, are averaged'
        WRITE(*,*)' Extensive variables, such as VFLOW_gy, are summed'
      ELSEIF(N .EQ. 32)THEN
        WRITE(*,*)' Enter the Reaction Rate index'
      ELSEIF(N .EQ. 40)THEN
        WRITE(*,*)' Enter the starting and ending values of J'
      ELSEIF(N .EQ. 41)THEN
        WRITE(*,*)' Intensive variables, such as EP_g, are averaged'
        WRITE(*,*)' Extensive variables, such as VFLOW_gy, are summed'
      ELSEIF(N .EQ. 50)THEN
        WRITE(*,*)' Enter the starting and ending values of K'
      ELSEIF(N .EQ. 51)THEN
        WRITE(*,*)' Intensive variables, such as EP_g, are averaged'
        WRITE(*,*)' Extensive variables, such as VFLOW_gy, are summed'
      ELSEIF(N .EQ. 70)THEN
        WRITE(*,*)' Enter the file name.  Enter * for displaying'
        WRITE(*,*)' the values at the terminal.  Enter ! for going'
        WRITE(*,*)' to the begining of this menu.'
      ENDIF
!
      WRITE(*,*)
!
      RETURN
      END
