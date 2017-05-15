!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: output                                                      !
!  Author: M. Syamlal                                 Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Contain data for output control.                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE output

      use param, only: DIMENSION_USR
      use param1, only: N_SPX

! Interval at which restart (.RES) file data is updated.
      DOUBLE PRECISION :: RES_DT
! Interval to create a backup copy of the RES file.
      DOUBLE PRECISION :: RES_BACKUP_DT
! Number of RES file copies to retain.
      INTEGER :: RES_BACKUPS
! Interval at which REAL restart (.SPx) files data are updated.
      DOUBLE PRECISION :: SPX_DT(N_SPX)
! Interval at which standard output (.OUT) file data is updated.
      DOUBLE PRECISION :: OUT_DT
! Interval at which user-defined output files are updated.
      DOUBLE PRECISION :: USR_DT (DIMENSION_USR)
! Interval to check and report the mass balance
      DOUBLE PRECISION :: REPORT_MASS_BALANCE_DT
! Interval in number of time steps at which LOG file is written
      INTEGER :: NLOG
! Flag to display messages and residuals on the screen
      LOGICAL :: FULL_LOG
! Flag to enable all ranks to write private LOG files.
      LOGICAL :: ENABLE_DMP_LOG
! Flag to print the index layout for  ijk<=>i,j,k  debugging tasks
      LOGICAL :: DBGPRN_LAYOUT
! Generate log files when negative gas density is detected.
      LOGICAL :: REPORT_NEG_DENSITY
! Generate log files when negative specific heat is detected.
      LOGICAL :: REPORT_NEG_SPECIFICHEAT


! Time at which special output is to be written
      DOUBLE PRECISION :: USR_TIME(DIMENSION_USR)
! Time at which restart file is to be written
      DOUBLE PRECISION :: RES_TIME
! Time at which restart backup file is to be written
      DOUBLE PRECISION :: RES_BACKUP_TIME
! Time at which REAL restart file is to be written
      DOUBLE PRECISION :: SPX_TIME(N_SPX)
! Time at which standard output is to be written
      DOUBLE PRECISION :: OUT_TIME

! The approximate amount (in MB) of space needed to write
! one time step of data into the indexed SPX file.
      DOUBLE PRECISION :: DISK(N_SPX) = 0.0d0
! The approximated total disk space (in MB)
      DOUBLE PRECISION :: DISK_TOT = 0.0d0
! One megabite (MB)
      DOUBLE PRECISION, PARAMETER :: ONEMEG = 1048576

! The following have no direct usage in the code. They are generic
! hooks wereby a user can specify some information. It may be useful
! to have the USR_I/J/K variables calculated in a check routine
! the same way that coordinates for BCs, ICs, PSs, are calculated.
!--------------------------------------------------------------------//
! X coordinate of the west face of user output region
      DOUBLE PRECISION :: USR_X_w (DIMENSION_USR)
! X coordinate of the east face of user output region
      DOUBLE PRECISION :: USR_X_e (DIMENSION_USR)
! Y coordinate of the south face of user output region
      DOUBLE PRECISION :: USR_Y_s (DIMENSION_USR)
! Y coordinate of the north face of user output region
      DOUBLE PRECISION :: USR_Y_n (DIMENSION_USR)
! Z coordinate of the bottom face of user output region
      DOUBLE PRECISION :: USR_Z_b (DIMENSION_USR)
! Z coordinate of the top face of user output region
      DOUBLE PRECISION :: USR_Z_t (DIMENSION_USR)
! I index of the west face of user output region
      INTEGER :: USR_I_w (DIMENSION_USR)
! I index of the east face of user output region
      INTEGER :: USR_I_e (DIMENSION_USR)
! J index of the south face of user output region
      INTEGER :: USR_J_s (DIMENSION_USR)
! J index of the north face of user output region
      INTEGER :: USR_J_n (DIMENSION_USR)
! K index of the bottom face of user output region
      INTEGER :: USR_K_b (DIMENSION_USR)
! K index of the top face of user output region
      INTEGER :: USR_K_t (DIMENSION_USR)
! Type of user-defined output: BINARY or ASCII.
      CHARACTER(LEN=6) :: USR_TYPE (DIMENSION_USR)
! Variables to be written in the user-defined output file.
      CHARACTER(LEN=60) :: USR_VAR (DIMENSION_USR)
! Format for writing user-defined (ASCII) output file.
      CHARACTER(LEN=60) :: USR_FORMAT (DIMENSION_USR)
! Extension for the user-defined output file.
      CHARACTER(LEN=16) :: USR_EXT (DIMENSION_USR)


      END MODULE output
