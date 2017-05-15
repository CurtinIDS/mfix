
!
!**********************************************************************
!*************************** init_usr_input ***************************
!**********************************************************************
!
subroutine init_usr_input

  use usr_input

  implicit none
  !
  !
  usr_t1 = 0.0
  usr_t2 = 0.0
  usr_tavg       = .false.
  usr_allow_tavg = .true.

  usr_done       = .false.

  usr_var     = 'EP_g'
  usr_var_num = 1

  usr_m = 1
  usr_n = 1

  usr_direction = 0

  usr_sum = .false.


  usr_t1_a = 0.0
  usr_t2_a = 0.0
  usr_tavg_a       = .false.
  usr_allow_tavg_a = .true.

  usr_i1         = 1
  usr_i2         = 1
  usr_i_avg      = .false.
  usr_must_i_avg = .false.

  usr_j1         = 1
  usr_j2         = 1
  usr_j_avg      = .false.
  usr_must_j_avg = .false.

  usr_k1         = 1
  usr_k2         = 1
  usr_k_avg      = .false.
  usr_must_k_avg = .false.

  usr_i1_a         = 1
  usr_i2_a         = 1
  usr_i_avg_a      = .false.
  usr_must_i_avg_a = .false.

  usr_j1_a         = 1
  usr_j2_a         = 1
  usr_j_avg_a      = .false.
  usr_must_j_avg_a = .false.

  usr_k1_a         = 1
  usr_k2_a         = 1
  usr_k_avg_a      = .false.
  usr_must_k_avg_a = .false.

  usr_status     = 0

  usr_fname        = '*'
  usr_ask_for_file = .true.

  return
end subroutine init_usr_input

!
!**********************************************************************
!*************************** get_usr_input ****************************
!**********************************************************************
!
subroutine get_usr_input(time_in_res,ask_for_times,init_read_res)

  use usr_input

  implicit none

  real    :: time_in_res
  logical :: ask_for_times , init_read_res
  !
  !
  ! get : usr_t1 , usr_t2 , usr_tavg
  !
  if (ask_for_times) then
     call get_usr_time(time_in_res)
     if (usr_done) return
  end if
  !
  ! get : usr_var , usr_var_num
  !
  call get_usr_variable

  ! get : usr_m
  !
  call get_usr_solid_phase
  !
  ! get : usr_n
  !
  call get_usr_species
  !
  ! calculate : usr_direction  and  usr_sum
  !
  call determine_usr_direction
  call determine_usr_sum
  !
  !
  ! set : read_spx(*) , at_eof(*)
  !
  call usr_enable_spx_files(init_read_res)
  !
  ! get : usr_i1 , usr_i2 , usr_i_avg
  !
  call get_usr_i_range
  !
  ! get : usr_j1 , usr_j2 , usr_j_avg
  !
  call get_usr_j_range
  !
  ! get : usr_k1 , usr_k2 , usr_k_avg
  !
  call get_usr_k_range
  !
  ! get usr_fname
  !
  call usr_get_fname
  if (usr_done) return
  !
  return
end subroutine get_usr_input
!
!**********************************************************************
!**************************** get_usr_time ****************************
!**********************************************************************
!
subroutine get_usr_time(time_in_res)

  use param1
  use usr_input

  implicit none

  CHARACTER(LEN=120) :: STRING, SUBSTR

  integer :: L

  real :: time_in_res
  !
  !  Read time
  !
10 WRITE(*,*)
  WRITE(*,'(A,F7.3,A,F7.3,A)',ADVANCE='NO')&
       ' Time: (',usr_t1,',',usr_t2,') > '
  READ(*,'(1A60)',ERR=10) STRING
  IF(STRING(1:1) .EQ. '?') THEN
     CALL HELP(10)
     GOTO 10
  ELSEIF(STRING(1:1) .EQ. 'e' .OR. STRING(1:1) .EQ. 'E' .OR. &
       STRING(1:1) .EQ. 'q' .OR. STRING(1:1) .EQ. 'Q') THEN
     usr_done = .true.
     IF (usr_fname(1:1) .NE. '*') CLOSE(40)
     RETURN
  ENDIF
  L = 1
  CALL GET_SUBSTR(STRING, L, SUBSTR)
  IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=10)usr_t1
  CALL GET_SUBSTR(STRING, L, SUBSTR)
  IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=10)usr_t2
  IF(usr_t1 .LT. ZERO .OR. usr_t2 .LT. ZERO) THEN
     usr_done = .true.
     IF (usr_fname(1:1) .NE. '*') CLOSE(40)
     RETURN
  ENDIF
  IF(usr_t1 .GE. TIME_IN_RES) THEN
     usr_t1 = TIME_IN_RES
     usr_t2   = TIME_IN_RES
  ENDIF

  IF(usr_t2 .LT. usr_t1)GOTO 10

  IF(usr_t1 .NE. usr_t2) THEN

     if (.not.usr_allow_tavg) then
        usr_tavg = .false.
        goto 100
     end if

     IF(usr_tavg)THEN
        SUBSTR(1:1) = 'Y'
     ELSE
        SUBSTR(1:1) = 'N'
     ENDIF
11   WRITE(*, '(A,1A1,A)',ADVANCE='NO')' Time average ? (',SUBSTR(1:1),') > '
     READ(*,'(1A60)',ERR=11) STRING
     IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(11)
        GOTO 11
     ENDIF
     IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
        usr_tavg = .TRUE.
     ELSEIF(STRING(1:1) .NE. ' ')THEN
        usr_tavg = .FALSE.
     ENDIF
  ENDIF
  !
100 continue
  usr_done = .false.
  !
  return
end subroutine get_usr_time
!
!**********************************************************************
!**************************** get_usr_variable ************************
!**********************************************************************
!
subroutine get_usr_variable

  use usr_input

  implicit none

  INTEGER  N_VAR
  PARAMETER (N_VAR=49)

  CHARACTER(LEN=120) :: STRING, SUBSTR
  CHARACTER(LEN=8) :: VAR      , VAR_DAT(N_VAR)

  integer :: L3 , L , var_no
  integer :: minmax    !?????

  logical :: strcmp

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

                                !                   48          49
       'Theta_m', 'Scalar' /
  !
  !  Read variable name
  !
  var = usr_var

20 CONTINUE

  SUBSTR(1:8) = var
  SUBSTR(9:9) = ' '

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
     L3 = 2
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
  DO L = 1, N_VAR
     IF (STRCMP(VAR,VAR_DAT(L)))THEN
        VAR_NO = L
        GOTO 23
     ENDIF
  ENDDO
  WRITE(*,'(A,1A8,A)')' Variable ', VAR, ' not found'
  GOTO 20
23 CONTINUE

  usr_var = var
  usr_var_num = var_no
  !
  return
end subroutine get_usr_variable
!
!**********************************************************************
!************************ get_usr_solid_phase *************************
!**********************************************************************
!
subroutine get_usr_solid_phase

  use physprop
  use usr_input

  implicit none

  integer   :: var_no , m , L3
  CHARACTER(LEN=120) :: STRING, SUBSTR
  !
  if (mmax .eq. 1) then
     usr_m = 1
     return
  end if
  !
  var_no = usr_var_num
  m      = usr_m
  !
  IF((VAR_NO .GE.  7 .AND. VAR_NO .LE. 10) .OR. &
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
24      WRITE(*,'(A,I2,A)',ADVANCE='NO') ' Solids phase: (', M, ') > '
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
  ENDIF




  return
end subroutine get_usr_solid_phase
!
!**********************************************************************
!*************************** get_usr_species **************************
!**********************************************************************
!
subroutine get_usr_species

  use physprop
  use usr_input

  implicit none
  !
  CHARACTER(LEN=120) :: STRING, SUBSTR
  integer   :: var_no , m , n , L3
  !
  var_no = usr_var_num
  m      = usr_m
  n      = usr_n
  !
  IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21)THEN
25   WRITE(*,'(A,I2,A)',ADVANCE='NO') ' Species: (', N, ') > '
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
  ENDIF
  !
  return
end subroutine get_usr_species
!
!**********************************************************************
!************************* determine_usr_direction ********************
!**********************************************************************
!
subroutine determine_usr_direction

  use usr_input

  implicit none

  integer :: var_no

  var_no = usr_var_num

  IF(VAR_NO .EQ.  4 .OR. VAR_NO .EQ.  7 .OR. VAR_NO .EQ. 16 .OR.&
       VAR_NO .EQ. 19 .OR. VAR_NO .EQ. 22 .OR. VAR_NO .EQ. 25 .OR.&
       VAR_NO .EQ. 28 .OR. VAR_NO .EQ. 31 .OR. VAR_NO .EQ. 36 .OR.&
       VAR_NO .EQ. 39)THEN
     usr_DIRECTION = 1
  ELSEIF &
       (VAR_NO .EQ.  5 .OR. VAR_NO .EQ.  8 .OR. VAR_NO .EQ. 17 .OR.&
       VAR_NO .EQ. 20 .OR. VAR_NO .EQ. 23 .OR. VAR_NO .EQ. 26 .OR.&
       VAR_NO .EQ. 29 .OR. VAR_NO .EQ. 32 .OR. VAR_NO .EQ. 37 .OR.&
       VAR_NO .EQ. 40)THEN
     usr_DIRECTION = 2
  ELSEIF&
       (VAR_NO .EQ.  6 .OR. VAR_NO .EQ.  9 .OR. VAR_NO .EQ. 18 .OR.&
       VAR_NO .EQ. 21 .OR. VAR_NO .EQ. 24 .OR. VAR_NO .EQ. 27 .OR.&
       VAR_NO .EQ. 30 .OR. VAR_NO .EQ. 33 .OR. VAR_NO .EQ. 38 .OR.&
       VAR_NO .EQ. 41)THEN
     usr_DIRECTION = 3
  ELSE
     usr_DIRECTION = 0
  ENDIF
  !
  return
end subroutine determine_usr_direction
!
!
!**********************************************************************
!************************** determine_usr_sum *************************
!**********************************************************************
!
subroutine determine_usr_sum

  use usr_input

  implicit none

  IF(usr_var_num .GE. 16 .AND. usr_var_num .LE. 41) THEN
     usr_sum = .TRUE.
  ELSE
     usr_sum = .FALSE.
  ENDIF
  !
  return
end subroutine determine_usr_sum
!
!**********************************************************************
!************************* usr_enable_spx_files ***********************
!**********************************************************************
!
subroutine usr_enable_spx_files(init_read_res)
  !
  !  Enable the required SPX file
  !
  use param1
  use physprop
  use usr_input

  implicit none

  include 'xforms.inc'

  CHARACTER(LEN=120) :: STRING, SUBSTR
  !
  integer :: L
  logical :: init_read_res
  if (init_read_res) then
     DO L = 1, N_SPX
        READ_SPX(L)    = .FALSE.
        AT_EOF(L)      = .FALSE.
     END DO
     !
     RES_P_g    = .FALSE.
     RES_T_g    = .FALSE.
     RES_X_g    = .FALSE.
  end if
  !
  var_no = usr_var_num
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
  IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21) THEN
     READ_SPX(7) = .TRUE.    ! X_g, X_s
  ENDIF

  IF(VAR_NO .EQ. 48 ) THEN
     READ_SPX(8) = .TRUE.    ! Theta_m
  ENDIF

  IF(VAR_NO .EQ. 49 ) THEN
     READ_SPX(9) = .TRUE.    ! Scalar
  ENDIF
  !
  !  Open P_g, T_g, and X_g files, if gas density needs to be determined
  !
  IF( RO_g0 .EQ. UNDEFINED .AND. &
       ( (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.&
       (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR.&
       (VAR_NO .EQ. 34                     ) .OR.&
       (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.&
       (VAR_NO .EQ. 42                     ) .OR.&
       (VAR_NO .EQ. 45                     ) &
       ) ) THEN
     WRITE(*,*)&
          ' To calculate gas density P_g, T_g, and X_g are needed'
     !
     IF(RES_P_g)THEN
        SUBSTR(1:1) = 'Y'
     ELSE
        SUBSTR(1:1) = 'N'
     ENDIF
26   WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
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
27   WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
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
     !
     RES_X_G = .FALSE.
     IF(MW_avg .EQ. UNDEFINED) THEN
        IF(RES_X_g)THEN
           SUBSTR(1:1) = 'Y'
        ELSE
           SUBSTR(1:1) = 'N'
        ENDIF
28      WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
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



  return
end subroutine usr_enable_spx_files
!
!**********************************************************************
!**************************** get_usr_i_range *************************
!**********************************************************************
!
subroutine get_usr_i_range

  use geometry
  use usr_input

  implicit none

  CHARACTER(LEN=120) :: STRING, SUBSTR

  integer :: L3 , i1 , i2
  logical :: i_average

  i1 = usr_i1
  i2 = usr_i2
  i_average = usr_i_avg
  if (usr_must_i_avg) i_average = .true.


30 WRITE(*,'(A,I3,A,I3,A)',ADVANCE='NO')&
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
  IF(I1 .NE. I2 .and. .not.usr_must_i_avg) THEN
     IF(I_AVERAGE)THEN
        SUBSTR(1:1) = 'Y'
     ELSE
        SUBSTR(1:1) = 'N'
     ENDIF
31   WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
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

  if (i1 .eq. i2) i_average = .false.

  usr_i1    = i1
  usr_i2    = i2
  usr_i_avg = i_average

  return
end subroutine get_usr_i_range

!
!**********************************************************************
!**************************** get_usr_j_range *************************
!**********************************************************************
!
subroutine get_usr_j_range

  use geometry
  use usr_input

  implicit none

  CHARACTER(LEN=120) :: STRING, SUBSTR

  integer :: L3 , j1 , j2
  logical :: j_average

  j1 = usr_j1
  j2 = usr_j2
  j_average = usr_j_avg
  if (usr_must_j_avg) j_average = .true.


30 WRITE(*,'(A,I3,A,I3,A)',ADVANCE='NO')&
       ' J range: (', j1, ',', j2, ') >'
  READ(*,'(1A60)') STRING
  IF(STRING(1:1) .EQ. '?') THEN
     CALL HELP(40)
     GOTO 30
  ENDIF
  L3 = 1
  CALL GET_SUBSTR(STRING, L3, SUBSTR)
  IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) j1
  CALL GET_SUBSTR(STRING, L3, SUBSTR)
  IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) j2
  !
  !     Check bounds
  !
  IF(j2 .LT. j1) j2 = j1
  IF(j1 .LT. 1 .OR. j1 .GT. JMAX2 .OR.&
       j2 .LT. 1 .OR. j2 .GT. JMAX2     ) THEN
     WRITE(*,'(A,I3)')&
          ' j1 and j2 should be in the range 1 to ', JMAX2
     GOTO 30
  ENDIF
  !
  IF(j1 .NE. j2 .and. .not.usr_must_j_avg) THEN
     IF(J_AVERAGE)THEN
        SUBSTR(1:1) = 'Y'
     ELSE
        SUBSTR(1:1) = 'N'
     ENDIF
31   WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
          ' Average or sum over J? (',SUBSTR(1:1),') > '
     READ(*,'(1A60)',ERR=31) STRING
     IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(41)
        GOTO 31
     ENDIF
     IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
        J_AVERAGE = .TRUE.
     ELSEIF(STRING(1:1) .NE. ' ')THEN
        J_AVERAGE = .FALSE.
     ENDIF
  ENDIF

  if (j1 .eq. j2) j_average = .false.

  usr_j1    = j1
  usr_j2    = j2
  usr_j_avg = j_average

  return
end subroutine get_usr_j_range

!
!**********************************************************************
!**************************** get_usr_k_range *************************
!**********************************************************************
!
subroutine get_usr_k_range

  use geometry
  use usr_input

  implicit none

  CHARACTER(LEN=120) :: STRING, SUBSTR

  integer :: L3 , k1 , k2
  logical :: k_average

  k1 = usr_k1
  k2 = usr_k2
  k_average = usr_k_avg
  if (usr_must_k_avg) k_average = .true.


30 WRITE(*,'(A,I3,A,I3,A)',ADVANCE='NO')&
       ' K range: (', k1, ',', k2, ') >'
  READ(*,'(1A60)') STRING
  IF(STRING(1:1) .EQ. '?') THEN
     CALL HELP(50)
     GOTO 30
  ENDIF
  L3 = 1
  CALL GET_SUBSTR(STRING, L3, SUBSTR)
  IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) k1
  CALL GET_SUBSTR(STRING, L3, SUBSTR)
  IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) k2
  !
  !     Check bounds
  !
  IF(k2 .LT. k1) k2 = k1
  IF(k1 .LT. 1 .OR. k1 .GT. KMAX2 .OR.&
       k2 .LT. 1 .OR. k2 .GT. KMAX2     ) THEN
     WRITE(*,'(A,I3)')&
          ' k1 and k2 should be in the range 1 to ', KMAX2
     GOTO 30
  ENDIF
  !
  IF(k1 .NE. k2 .and. .not.usr_must_k_avg) THEN
     IF(K_AVERAGE)THEN
        SUBSTR(1:1) = 'Y'
     ELSE
        SUBSTR(1:1) = 'N'
     ENDIF
31   WRITE(*, '(A,1A1,A)',ADVANCE='NO')&
          ' Average or sum over K? (',SUBSTR(1:1),') > '
     READ(*,'(1A60)',ERR=31) STRING
     IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(51)
        GOTO 31
     ENDIF
     IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
        K_AVERAGE = .TRUE.
     ELSEIF(STRING(1:1) .NE. ' ')THEN
        K_AVERAGE = .FALSE.
     ENDIF
  ENDIF

  if (k1 .eq. k2) k_average = .false.

  usr_k1    = k1
  usr_k2    = k2
  usr_k_avg = k_average

  return
end subroutine get_usr_k_range
!
! ******************************************************************
! **************************** usr_set_array ***********************
! ******************************************************************
!
subroutine usr_set_array(nsize,arr,var_no,m,n)

  use geometry
  use fldvar
  use indices
  use compar
  use constant
  use post3d, only: ydist_sc
  use functions

  implicit none

  integer :: var_no , nsize , m , n , mIJK , lIJK
  real    :: arr(nsize)
  real    :: DELm , DELl , FAC1 , FAC2

  integer :: i , j , k , ijk

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
  !
  if (var_no .eq.  1) arr(:) = EP_g(:)
  if (var_no .eq.  2) arr(:) = P_g(:)
  if (var_no .eq.  3) arr(:) = P_star(:)
  if (var_no .eq.  4) arr(:) = U_g(:)
  if (var_no .eq.  5) arr(:) = V_g(:)
  if (var_no .eq.  6) arr(:) = W_g(:)
  if (var_no .eq.  7) arr(:) = U_s(:,m)
  if (var_no .eq.  8) arr(:) = V_s(:,m)
  if (var_no .eq.  9) arr(:) = W_s(:,m)
  if (var_no .eq. 10) arr(:) = ROP_s(:,m)
  if (var_no .eq. 11) arr(:) = T_g(:)
  if (var_no .eq. 12) arr(:) = T_s(:,m)
  if (var_no .eq. 13) arr(:) = T_s(:,2)
  if (var_no .eq. 14) arr(:) = X_g(:,n)
  if (var_no .eq. 15) arr(:) = X_s(:,m,n)

  if (var_no .eq. 48) arr(:) = Theta_m(:,m)
  if (var_no .eq. 49) arr(:) = Scalar(:,n)   ! need a usr_scalar

  if (var_no.ge.16 .and. var_no.le.47) then

     do ijk = 1,nsize
        i = i_of(ijk)
        j = j_of(ijk)
        k = k_of(ijk)
        if (var_no .eq. 16) arr(ijk) = XFLOW_gx(I, J, K, IJK, N)
        if (var_no .eq. 17) arr(ijk) = XFLOW_gy(I, J, K, IJK, N)
        if (var_no .eq. 18) arr(ijk) = XFLOW_gz(I, J, K, IJK, N)
        if (var_no .eq. 19) arr(ijk) = XFLOW_sx(I, J, K, IJK, M , N)
        if (var_no .eq. 20) arr(ijk) = XFLOW_sy(I, J, K, IJK, M , N)
        if (var_no .eq. 21) arr(ijk) = XFLOW_sz(I, J, K, IJK, M , N)
        if (var_no .eq. 22) arr(ijk) = MFLOW_gx(I, J, K, IJK)
        if (var_no .eq. 23) arr(ijk) = MFLOW_gy(I, J, K, IJK)
        if (var_no .eq. 24) arr(ijk) = MFLOW_gz(I, J, K, IJK)
        if (var_no .eq. 25) arr(ijk) = MFLOW_sx(I, J, K, IJK, M)
        if (var_no .eq. 26) arr(ijk) = MFLOW_sy(I, J, K, IJK, M)
        if (var_no .eq. 27) arr(ijk) = MFLOW_sz(I, J, K, IJK, M)
        if (var_no .eq. 28) arr(ijk) = VFLOW_gx(I, J, K, IJK)
        if (var_no .eq. 29) arr(ijk) = VFLOW_gy(I, J, K, IJK)
        if (var_no .eq. 30) arr(ijk) = VFLOW_gz(I, J, K, IJK)
        if (var_no .eq. 31) arr(ijk) = VFLOW_sx(I, J, K, IJK, M)
        if (var_no .eq. 32) arr(ijk) = VFLOW_sy(I, J, K, IJK, M)
        if (var_no .eq. 33) arr(ijk) = VFLOW_sz(I, J, K, IJK, M)
        if (var_no .eq. 34) arr(ijk) = EP_g(IJK) * CALC_RO_g(IJK) * VOL(IJK)
        if (var_no .eq. 35) arr(ijk) = ROP_s(IJK, M) * VOL(IJK)
        if (var_no .eq. 36) arr(ijk) = FLUX_gx(IJK)
        if (var_no .eq. 37) arr(ijk) = FLUX_gy(IJK)
        if (var_no .eq. 38) arr(ijk) = FLUX_gz(IJK)
        if (var_no .eq. 39) arr(ijk) = FLUX_sx(IJK, M)
        if (var_no .eq. 40) arr(ijk) = FLUX_sy(IJK, M)
        if (var_no .eq. 41) arr(ijk) = FLUX_sz(IJK, M)
        if (var_no .eq. 42) arr(ijk) = 0.5*CALC_RO_g(IJK)*EP_g(IJK) &
             *(U_g(IJK)**2 + V_g(IJK)**2 + W_g(IJK)**2)
        if (var_no .eq. 43) arr(ijk) = 0.5*ROP_s(IJK,M) &
             *(U_s(IJK,M)**2 + V_s(IJK,M)**2 + W_s(IJK,M)**2)
        if (var_no .eq. 44) then
           arr(ijk) = P_s(IJK,M)
           if (EP_g(IJK) .LT. EP_star) arr(ijk) = P_star(IJK)
        end if
        if (var_no .eq. 45) arr(ijk) = GRAVITY * YDIST_SC(J) &
             * CALC_RO_g(IJK) * EP_g(IJK)
        if (var_no .eq. 46) arr(ijk) = GRAVITY * YDIST_SC(J) &
             * ROP_s(IJK,M)
        if (var_no .eq. 47) then
           mIJK= FUNIJK(I,J+1,K)
           lIJK= FUNIJK(I,J-1,K)
           DELm = YDIST_SC(J+1) - YDIST_SC(J)
           DELl = YDIST_SC(J)   - YDIST_SC(J-1)
           FAC1= (ROP_s(mIJK,M)-ROP_s(IJK,M))/DELm &
                +(ROP_s(IJK,M)-ROP_s(lIJK,M))/DELl
           FAC2=P_s(IJK,M)/(ROP_s(IJK,M)**2)
           arr(ijk) = DY(J)*FAC1*FAC2/2
        end if

     end do


     !        ELSEIF(VAR_NO .EQ. 47)THEN
     !         mIJK= FUNIJK(I,J+1,K)
     !         lIJK= FUNIJK(I,J-1,K)
     !         DELm=YDIST_SC(J+1)-YDIST_SC(J)
     !         DELl=YDIST_SC(J)-YDIST_SC(J-1)
     !         FAC1= (ROP_s(mIJK,M)-ROP_s(IJK,M))/DELm &
     !              +(ROP_s(IJK,M)-ROP_s(lIJK,M))/DELl
     !        FAC2=P_s(IJK,M)/(ROP_s(IJK,M)**2)
     !         VALUE_TMP = DY(J)*FAC1*FAC2/2


  end if



  return
end subroutine usr_set_array
!
! ***************************************************************
! ************************* spatial_averaging *******************
! ***************************************************************
!
subroutine spatial_averaging(value,i1,i2,j1,j2,k1,k2, &
     usr_ii_avg,usr_jj_avg,usr_kk_avg)
  !
  use geometry
  use indices
  use compar
  use usr_input
  use functions

  implicit none
  !
  logical :: usr_ii_avg , usr_jj_avg , usr_kk_avg
  integer :: i , j , k , ijk , ijk1 , i1 , i2 , j1 , j2
  integer :: k1 , k2 , k2d , j2d , i2d
  real    :: value(*)
  real, allocatable :: dist(:)

  allocate (dist(ijkmax2))

  !
  !  DO I, J, or K averaging
  !
  IF(usr_kk_avg)THEN
     K = K1
     DO J = J1, J2
        DO I = I1, I2
           IJK = FUNIJK(I, J, K)
           IF(WALL_AT(IJK))THEN
              VALUE(IJK) = ZERO
              DIST(IJK)  = ZERO
           ELSEIF(.NOT. usr_SUM) THEN
              IF(usr_DIRECTION .EQ. 3)THEN
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
              IJK = FUNIJK(I, J, K)
              IJK1 = FUNIJK(I, J, K1)
              IF(.NOT.WALL_AT(IJK)) THEN
                 IF(usr_SUM) THEN
                    VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
                 ELSE
                    IF(usr_DIRECTION .EQ. 3)THEN
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
     !
     IF(.NOT. usr_sum) THEN
        K = K1
        DO J = J1, J2
           DO I = I1, I2
              IJK = FUNIJK(I, J, K)
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
  !
  IF(usr_jj_avg)THEN
     J = J1
     DO K = K1, K2d
        DO I = I1, I2
           IJK = FUNIJK(I, J, K)
           IF(WALL_AT(IJK) .AND. .NOT.usr_k_avg)THEN
              VALUE(IJK) = ZERO
              DIST(IJK)  = ZERO
           ELSEIF(.NOT. usr_sum) THEN
              IF(usr_DIRECTION .EQ. 2)THEN
                 VALUE(IJK) = VALUE(IJK) * DY_N(J)
                 DIST(IJK)  = DY_N(J)
              ELSE
                 VALUE(IJK) = VALUE(IJK) * DY(J)
                 DIST(IJK)  = DY(J)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !
     DO K = K1, K2d
        DO J = J1+1, J2
           DO I = I1, I2
              IJK = FUNIJK(I, J, K)
              IJK1 = FUNIJK(I, J1, K)
              IF(.NOT.WALL_AT(IJK) .OR. usr_k_avg) THEN
                 IF(usr_SUM) THEN
                    VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
                 ELSE
                    IF(usr_DIRECTION .EQ. 2)THEN
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
     !
     IF(.NOT.usr_sum) THEN
        J = J1
        DO K = K1, K2d
           DO I = I1, I2
              IJK = FUNIJK(I, J, K)
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
  !
  IF(usr_ii_avg)THEN
     I = I1
     DO K = K1, K2d
        DO J = J1, J2d
           IJK = FUNIJK(I, J, K)
           IF(WALL_AT(IJK) .AND. .NOT.usr_j_avg .AND. .NOT.usr_k_avg)THEN
              VALUE(IJK) = ZERO
              DIST(IJK)  = ZERO
           ELSEIF(.NOT. usr_sum) THEN
              IF(usr_DIRECTION .EQ. 1)THEN
                 VALUE(IJK) = VALUE(IJK) * DX_E(I) * X_E(I)
                 DIST(IJK)  = DX_E(I) * X_E(I)
              ELSE
                 VALUE(IJK) = VALUE(IJK) * DX(I) * X(I)
                 DIST(IJK)  = DX(I) * X(I)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !
     DO K = K1, K2d
        DO J = J1, J2d
           DO I = I1+1, I2
              IJK = FUNIJK(I, J, K)
              IJK1 = FUNIJK(I1, J, K)
              IF(.NOT.WALL_AT(IJK) .OR. usr_j_avg .OR. usr_k_avg) THEN
                 IF(usr_SUM) THEN
                    VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
                 ELSE
                    IF(usr_DIRECTION .EQ. 1)THEN
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

     IF(.NOT.usr_SUM) THEN
        I = I1
        DO K = K1, K2d
           DO J = J1, J2d
              IJK = FUNIJK(I, J, K)
              IF (DIST(IJK) .NE. ZERO) THEN
                 VALUE(IJK) = VALUE(IJK) / DIST(IJK)
              ELSE IF (VALUE(IJK) .NE. ZERO) THEN
                 WRITE(*,*)' Error in I-averaging'
              ENDIF
           end do
        end do
     ENDIF
     I2d = I1
  ELSE
     I2d = I2
  ENDIF
  !
  deallocate(dist)

  return

  contains

    double precision function DZ_T(K)
      implicit none
      integer, intent(in) :: k
      DZ_T = HALF * (DZ(K) + DZ(Kp1(K)))
    end function DZ_T

    double precision function DY_N(J)
      implicit none
      integer, intent(in) :: j
      DY_N = HALF * (DY(J) + DY(Jp1(J)))
    end function DY_N

    double precision function DX_E(I)
      implicit none
      integer, intent(in) :: i
      DX_E = HALF * (DX(I) + DX(Ip1(I)))
    end function DX_E

end subroutine spatial_averaging
!
! *********************************************************************
! ************************** usr_start_time ***************************
! *********************************************************************
!
subroutine usr_start_time(time_in_res)

  use usr_input

  implicit none

  real :: time_in_res , time_found

  IF (usr_t1 .LT. TIME_IN_RES) THEN
     CALL SEEK_TIME(READ_SPX,usr_t1,REC_POINTER,TIME_FOUND)
     IF (TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        usr_status = -1
        return
     ENDIF
  ENDIF
  !
  usr_status = 0
  !
  return
end subroutine usr_start_time
!
! *********************************************************************
! *************************** usr_get_fname ***************************
! *********************************************************************
!
subroutine usr_get_fname

  use usr_input

  implicit none

  integer   :: L3         , ians
  logical   :: strcmp     , file_exist
  CHARACTER(LEN=120) :: STRING, SUBSTR

  if (.not.usr_ask_for_file) return

  usr_done = .false.

70 WRITE(*,'(A,1A30,A)',ADVANCE='NO') ' File: (', usr_fname,') >'
  READ(*,'(1A60)') STRING
  IF (STRING(1:1) .EQ. '?') THEN
     CALL HELP(70)
     GOTO 70
  ELSE IF (STRING(1:1) .EQ. '!') THEN
     usr_done = .true.
     return
  ENDIF
  L3 = 1
  CALL GET_SUBSTR(STRING, L3, SUBSTR)
  IF (SUBSTR(1:1) .NE. ' ')THEN
     IF (.NOT.STRCMP(STRING(1:30),usr_fname)) THEN
        IF (usr_fname(1:1) .NE. '*') CLOSE(40)
        CALL STREQS(usr_fname,STRING(1:30))
        IF (usr_fname(1:1) .NE. '*') THEN
           INQUIRE (FILE=usr_fname,EXIST=FILE_EXIST)
           IF (FILE_EXIST) THEN
              WRITE(*,'(A)',ADVANCE='NO')' File exists.  Over write? (1=Yes) >'
              READ(*,*)IANS
              IF(IANS .NE. 1)GOTO 70
           ENDIF
           OPEN (UNIT=40,FILE=usr_fname,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
        ENDIF
     ENDIF
  ENDIF
  !
  return
end subroutine usr_get_fname
!
!**********************************************************************
!************************** usr_write_input ***************************
!**********************************************************************
!
subroutine usr_write_input(code)

  use usr_input

  implicit none

  integer :: code
  !
  if (code .eq. 1) then
     call usr_write_input2(usr_t1,usr_t2,usr_tavg, &
          usr_var,usr_m,usr_n, &
          usr_i1,usr_i2,usr_i_avg, &
          usr_j1,usr_j2,usr_j_avg, &
          usr_k1,usr_k2,usr_k_avg, &
          usr_fname(1:1))
  end if
  !
  if (code .eq. 2) then
     call usr_write_input2(usr_t1,usr_t2,usr_tavg, &
          usr_var_a,usr_m_a,usr_n_a, &
          usr_i1_a,usr_i2_a,usr_i_avg_a, &
          usr_j1_a,usr_j2_a,usr_j_avg_a, &
          usr_k1_a,usr_k2_a,usr_k_avg_a, &
          usr_fname(1:1))
  end if

  return
end subroutine usr_write_input

!
!**********************************************************************
!************************** usr_write_input2 **************************
!**********************************************************************
!
subroutine usr_write_input2(usr_t1_a,usr_t2_a,usr_tavg_a, &
     usr_var_a,usr_m_a,usr_n_a, &
     usr_i1_a,usr_i2_a,usr_i_avg_a, &
     usr_j1_a,usr_j2_a,usr_j_avg_a, &
     usr_k1_a,usr_k2_a,usr_k_avg_a, &
     output_flag)
  implicit none
  !

  real :: usr_t1_a
  real :: usr_t2_a

  logical :: usr_tavg_a

  character(len=8) :: usr_var_a
  character(len=1) :: output_flag

  integer   :: usr_m_a
  integer   :: usr_n_a

  integer   :: usr_i1_a
  integer   :: usr_i2_a
  logical   :: usr_i_avg_a

  integer   :: usr_j1_a
  integer   :: usr_j2_a
  logical   :: usr_j_avg_a

  integer   :: usr_k1_a
  integer   :: usr_k2_a
  logical   :: usr_k_avg_a
  !
  if (output_flag .eq. '*') then

     write (*,*) ' '
     write (*,*) ' time        : ' , usr_t1_a , ' to ' , usr_t2_a
     if (usr_tavg_a) write (*,*) '                     time averaged'

     write (*,*) ' variable    : ' , usr_var_a

     write (*,*) ' solid phase : ' , usr_m_a
     write (*,*) ' species     : ' , usr_n_a

     write (*,*) ' I-range     : ' , usr_i1_a , ' to ' , usr_i2_a
     if (usr_i_avg_a) write (*,*) '                   i-spatially averaged'

     write (*,*) ' J-range     : ' , usr_j1_a , ' to ' , usr_j2_a
     if (usr_j_avg_a) write (*,*) '                   j-spatially averaged'

     write (*,*) ' K-range     : ' , usr_k1_a , ' to ' , usr_k2_a
     if (usr_k_avg_a) write (*,*) '                   k-spatially averaged'

     write (*,*) ' '

  else


     write (40,*) ' '
     write (40,*) ' time        : ' , usr_t1_a , ' to ' , usr_t2_a
     if (usr_tavg_a) write (40,*) '                     time averaged'

     write (40,*) ' '
     write (40,*) ' variable    : ' , usr_var_a

     write (40,*) ' solid phase : ' , usr_m_a
     write (40,*) ' species     : ' , usr_n_a

     write (40,*) ' I-range     : ' , usr_i1_a , ' to ' , usr_i2_a
     if (usr_i_avg_a) write (40,*) '                   i-spatially averaged'

     write (40,*) ' J-range     : ' , usr_j1_a , ' to ' , usr_j2_a
     if (usr_j_avg_a) write (40,*) '                   j-spatially averaged'

     write (40,*) ' K-range     : ' , usr_k1_a , ' to ' , usr_k2_a
     if (usr_k_avg_a) write (40,*) '                   k-spatially averaged'

     write (40,*) ' '

  endif

  return
end subroutine usr_write_input2

