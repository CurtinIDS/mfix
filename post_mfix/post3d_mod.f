      MODULE post3d

      Use param
      Use param1

      INTEGER       N_PLOT_MENU
      PARAMETER     (N_PLOT_MENU = 7)
!
      INTEGER       N_VAR_MENU
      PARAMETER     (N_VAR_MENU = 13)
!
!                   Flag to tell whether an SPX file is open
      LOGICAL       SPX_OPEN(N_SPX)
!
!                   Last record in an SPX file, number of files
      INTEGER       LAST_REC(N_SPX), NUM_REC(N_SPX)
!
!                   index of variable within the file to plot
      INTEGER       VAR_INDEX
      INTEGER       LOC_X,LOC_Y,LOC_Z
      INTEGER       PLOT_TYPE,PLOT_TYPE2
      INTEGER       SELECTION
      INTEGER       COUNT
      CHARACTER(LEN=60) ::  PLOT_MENU(N_PLOT_MENU)
      CHARACTER(LEN=60) ::  VAR_MENU(N_VAR_MENU)
      REAL          VERSION_NUMBER
      REAL          XDIST_SC(DIM_I) , XDIST_VEC(DIM_I)
      REAL          YDIST_SC(DIM_J) , YDIST_VEC(DIM_J)
      REAL          ZDIST_SC(DIM_K) , ZDIST_VEC(DIM_K)

      END MODULE post3d
