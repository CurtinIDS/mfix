!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: DEPRECATED_OR_UNKNOWN                               !
!     Author: J.Musser                                Date:  5-SEPT-14 !
!                                                                      !
!     Purpose: This routine is called when a keyword was not matched   !
!     to any of the keywords in the namelist files. This routine       !
!     reports if the keyword was deprecated or incorrect.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEPRECATED_OR_UNKNOWN(LINE_NO, INPUT)

      use param
      use param1
      use compar, only: myPE
      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT

      CHARACTER(len=256) :: STRING
      INTEGER :: IOS

! Old keyword for solids density :: Replaced by RO_s0
      DOUBLE PRECISION :: RO_s
      LOGICAL :: BC_APPLY_TO_MPPIC(DIMENSION_BC)
      INTEGER :: COHESION_DEBUG, DES_MMAX,DES_NMAX_s(DIM_M),DIMN
      CHARACTER(LEN=16) :: DES_BC_TYPE(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_MASSFLOW_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION :: DES_BC_ROP_s (DIMENSION_BC, DIM_M)
      DOUBLE PRECISION :: DES_BC_T_s (DIMENSION_IC, DIM_M)
      DOUBLE PRECISION :: DES_BC_VOLFLOW_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION :: DES_BC_X_e(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_X_w(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Y_n(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Y_s(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Z_b(DIMENSION_BC)
      DOUBLE PRECISION :: DES_BC_Z_t(DIMENSION_BC)
      DOUBLE PRECISION :: DES_C_ps0(DIM_M), DES_K_s0(DIM_M)
      LOGICAL :: DES_CALC_BEDHEIGHT, DES_CONV_EQ, DES_ENERGY_EQ
      LOGICAL :: DES_COND_EQ,DES_COND_EQ_PFP,DES_COND_EQ_PP, DES_RADI_EQ

      DOUBLE PRECISION :: DES_D_P0 (DIM_M), DES_F, DES_GAMMA
      DOUBLE PRECISION :: DES_EPS_XSTART,DES_EPS_YSTART,DES_EPS_ZSTART

      DOUBLE PRECISION :: DES_IC_X_e(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_X_w(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Y_n(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Y_s(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Z_b(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_Z_t(DIMENSION_IC)
      DOUBLE PRECISION :: DES_IC_T_s (DIMENSION_IC, DIM_M)
      DOUBLE PRECISION :: DES_MW_s(DIM_M, DIM_N_s), DES_RO_s (DIM_M)
      DOUBLE PRECISION :: DTSOLID_FACTOR, LID_VEL
      DOUBLE PRECISION :: MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH
      CHARACTER(len=18) DES_SPECIES_s(DIM_M, DIM_N_s)
      CHARACTER(len=32)  DES_SPECIES_ALIAS_s(DIM_M, DIM_N_s)
      LOGICAL :: DES_PERIODIC_WALLS, DES_PERIODIC_WALLS_X,DES_SPECIES_EQ
      LOGICAL :: DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z,TSUJI_DRAG

      INTEGER :: MAX_DES_BC_CELL, NPC_PIC(DIM_M)
      INTEGER :: QLM, QLN, INIT_QUAD_COUNT
      LOGICAL :: MPPIC_CONSTANTNPC, MPPIC_CONSTANTWT, SQUARE_WELL
      LOGICAL :: WALLDTSPLIT,WALLREFLECT, CALL_DI, CALL_GROW, CALL_ISAT
      DOUBLE PRECISION :: MQUAD_FACTOR,STATWT_PIC(DIM_M), ISATdt
      DOUBLE PRECISION :: RADIUS_RATIO,WALL_RADIUS_RATIO
      DOUBLE PRECISION :: pvel_mean, PVEL_StDev,VOL_FRAC(DIM_M)
      CHARACTER(LEN=64) :: REACTION_MODEL

      LOGICAL :: DISCRETE_ELEMENT, MPPIC, DES_CONTINUUM_HYBRID

      DOUBLE PRECISION :: PARTICLES_FACTOR, DES_RES_DT, DES_SPX_DT
      INTEGER :: MAX_PIS, MAX_FACETS_PER_CELL_DES
      LOGICAL :: USE_STL_DES, DES_CONTINUUM_COUPLED, &
         DES_CONVERT_BOX_TO_FACETS


! 2014-1 Deprecated list:
!-----------------------------------------------------------------------
      NAMELIST / DEP_2014_1 / RO_s, BC_APPLY_TO_MPPIC, COHESION_DEBUG, &
         DES_BC_MASSFLOW_s, DES_BC_ROP_s,DES_BC_T_s,DES_BC_TYPE,       &
         DES_BC_VOLFLOW_s,DES_BC_X_e,DES_BC_X_w, DES_BC_Y_n,DES_BC_Y_s,&
         DES_BC_Z_b,DES_BC_Z_t,DES_CALC_BEDHEIGHT,QLM,QLN,DES_COND_EQ, &
         DES_COND_EQ_PFP,DES_COND_EQ_PP,DES_CONV_EQ,DES_C_ps0,DES_D_p0,&
         DES_ENERGY_EQ,DES_EPS_XSTART,DES_F,DES_EPS_YSTART,DIMN,       &
         DES_EPS_ZSTART,DES_GAMMA,DES_IC_T_s,DES_IC_X_e,DES_IC_X_w,    &
         DES_IC_Y_n,DES_IC_Y_s,DES_IC_Z_b,DES_IC_Z_t,DES_K_s0,DES_MMAX,&
         DES_MW_s,DES_NMAX_s,DES_PERIODIC_WALLS,DES_SPECIES_s,LID_VEL, &
         DES_PERIODIC_WALLS_X,DES_PERIODIC_WALLS_Y,WALLREFLECT,        &
         DES_PERIODIC_WALLS_Z,DES_RADI_EQ,DES_RO_s,DES_SPECIES_ALIAS_s,&
         DES_SPECIES_EQ,DTSOLID_FACTOR,INIT_QUAD_COUNT,MAX_DES_BC_CELL,&
         MASTER_WALL_WELL_DEPTH,MASTER_WELL_DEPTH,MPPIC_CONSTANTNPC,   &
         MPPIC_CONSTANTWT,MQUAD_FACTOR,NPC_PIC,pvel_mean,pvel_stdev,   &
         RADIUS_RATIO,REACTION_MODEL,SQUARE_WELL,STATWT_PIC,STATWT_PIC,&
         TSUJI_DRAG,VOL_FRAC,WALLDTSPLIT,WALL_RADIUS_RATIO, CALL_DI,   &
         CALL_GROW, CALL_ISAT, ISATdt


! 2015-1 Deprecated list:
!-----------------------------------------------------------------------
      NAMELIST / DEP_2015_1 / MAX_PIS, PARTICLES_FACTOR, USE_STL_DES,  &
         DISCRETE_ELEMENT, MPPIC, DES_CONTINUUM_HYBRID

! 2015-2 Deprecated list:
!-----------------------------------------------------------------------
      NAMELIST / DEP_2015_2 / DES_CONTINUUM_COUPLED, DES_RES_DT,       &
         DES_SPX_DT

! 2016-1 Deprecated list:
!-----------------------------------------------------------------------
      NAMELIST / DEP_2016_1 / MAX_FACETS_PER_CELL_DES,                 &
          DES_CONVERT_BOX_TO_FACETS

! 2014-1 Release Deprecated keywords.
      STRING=''; STRING = '&DEP_2014_1 '//trim(adjustl(INPUT))//'/'
      READ(STRING,NML=DEP_2014_1,IOSTAT=IOS)
      IF(IOS == 0) CALL DEPRECATED(LINE_NO, INPUT, '2014-1')


! 2015-1 Release Deprecated keywords.
      STRING=''; STRING = '&DEP_2015_1 '//trim(adjustl(INPUT))//'/'
      READ(STRING,NML=DEP_2015_1,IOSTAT=IOS)
      IF(IOS == 0) CALL DEPRECATED(LINE_NO, INPUT, '2015-1')

! 2015-2 Release Deprecated keywords.
      STRING=''; STRING = '&DEP_2015_2 '//trim(adjustl(INPUT))//'/'
      READ(STRING,NML=DEP_2015_1,IOSTAT=IOS)
      IF(IOS == 0) CALL DEPRECATED(LINE_NO, INPUT, '2015-2')

! Everything else...  This should be the last call in this routine.
      CALL UNKNOWN_KEYWORD(LINE_NO, INPUT)


      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: DEPRECATED                                          !
!     Author: J.Musser                                Date:  5-SEPT-14 !
!                                                                      !
!     Purpose: Write the error message for deprecated keywords.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEPRECATED(LINE_NO, INPUT, RELEASE)

      INTEGER, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT
      CHARACTER(len=*), INTENT(IN) :: RELEASE

      IF(myPE == 0) &
         WRITE(*,1000) trim(iVAL(LINE_NO)), RELEASE, trim(INPUT)

      CALL MFIX_EXIT(myPE)

 1000 FORMAT(//1X,70('*')/' From DEPRECATED',/' Error 1000:',          &
         ' A keyword pair on line ',A,' of the mfix.dat file was',/    &
         ' identified as being deprecated as of the ',A,' Release.',// &
         3x,A,//' Please see the user documentation and update the ',  &
         'mfix.dat file.',/1X,70('*')//)

      END SUBROUTINE DEPRECATED


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: UNKNOWN_KEYWORD                                     !
!     Author: J.Musser                                Date:  5-SEPT-14 !
!                                                                      !
!     Purpose: Write the error message for deprecated keywords.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE UNKNOWN_KEYWORD(LINE_NO, INPUT)

      INTEGER, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT

      IF(myPE == 0) WRITE(*,2000) trim(iVAL(LINE_NO)), trim(INPUT)

      CALL MFIX_EXIT(myPE)

 2000 FORMAT(//1X,70('*')/' From: UNKNOWN_KEYWORD',/' Error 2000: ',   &
         'Unable to process line ',A,' of the mfix.dat file.',2/3x,    &
         A,2/1x,'Possible causes are',/3x,'* Incorrect or illegal ',   &
         'keyword format',/3x,'* Unknown or mistyped name',/3x,'* ',   &
         'The mensioned item is too small (array overflow).', 2/1x,    &
         'Please see the user documentation and update the mfix.dat ', &
         'file. ',/1X,70('*')//)


      END SUBROUTINE UNKNOWN_KEYWORD


      END SUBROUTINE DEPRECATED_OR_UNKNOWN
