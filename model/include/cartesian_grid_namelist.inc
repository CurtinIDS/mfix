!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CARTESIAN_GRID_NAMELIST.INC                            C
!  Purpose: Variables needed for CARTESIAN_GRID input                  C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 26-Aug-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
       NAMELIST / CARTESIAN_GRID_INPUT_DATA /&

      CARTESIAN_GRID , N_QUADRIC,USE_STL,USE_MSH,USE_POLYGON,N_USR_DEF,&
      quadric_form, quadric_scale, &
      lambda_x, lambda_y, lambda_z, dquadric, &
      theta_x, theta_y, theta_z, &
      Radius, Half_angle,n_x,n_y,n_z, &
      Torus_R1,Torus_R2, &
      UCOIL_R1, UCOIL_R2, UCOIL_Y1, UCOIL_Y2, &
      BEND_R1, BEND_R2, BEND_THETA1, BEND_THETA2, &
      C2C_R1, C2C_R2, C2C_Y1,C2C_Y2, &
      REACTOR1_R1,REACTOR1_R2,REACTOR1_Y1,REACTOR1_Y2, &
      REACTOR1_YR1,REACTOR1_YR2, &
      REACTOR1_RR1, REACTOR1_RR2,REACTOR1_THETA1,REACTOR1_THETA2, &
      t_x, t_y, t_z, &
      clip_xmin, clip_xmax, clip_ymin, clip_ymax, clip_zmin, clip_zmax, &
      piece_xmin, piece_xmax, piece_ymin, piece_ymax, piece_zmin, piece_zmax, &
      FLUID_IN_CLIPPED_REGION, &
      BC_ID_Q,&
      N_GROUP, GROUP_SIZE,GROUP_Q,GROUP_RELATION,RELATION_WITH_PREVIOUS, &
      TOL_SNAP,TOL_DELH, TOL_SMALL_CELL,TOL_MERGE,&
      TOL_SMALL_AREA, ALPHA_MAX, &
      TOL_F, TOL_POLY,ITERMAX_INT, &
      TOL_STL,STL_SMALL_ANGLE,TOL_STL_DP,RAY_DIR,DIM_FACETS_PER_CELL,&
      OUT_STL_VALUE,STL_BC_ID,CAD_PROPAGATE_ORDER,&
      TX_STL,TY_STL,TZ_STL,SCALE_STL,&
      TOL_MSH,MSH_SMALL_ANGLE,OUT_MSH_VALUE,&
      TX_MSH,TY_MSH,TZ_MSH,SCALE_MSH,&
      SET_CORNER_CELLS,&
      FAC_DIM_MAX_CUT_CELL, &
      WRITE_VTK_FILES,TIME_DEPENDENT_FILENAME, VTK_DT, VTK_DBG_FILE, &
      VTK_VAR,FRAME,VTU_DIR,&
      VTK_X_w, VTK_X_e, VTK_Y_s, VTK_Y_n, VTK_Z_b, VTK_Z_t,VTK_VARLIST, &
      VTK_EP_g, VTK_P_g, VTK_P_star, VTK_VEL_g, VTK_VEL_S, &
      VTK_U_G, VTK_V_G, VTK_W_G, VTK_U_S, VTK_V_S, VTK_W_S, &
      VTK_ROP_s,  VTK_T_g,  VTK_T_s, VTK_X_g, VTK_X_s, &
      VTK_Theta_m, VTK_Scalar, VTK_RRate, VTK_K_Turb_G, VTK_E_Turb_G, &
      VTK_VORTICITY, VTK_LAMBDA_2, VTK_PARTITION, VTK_BC_ID, VTK_DWALL, &
      VTK_FACET_COUNT_DES, VTK_NB_FACET_DES, VTK_IJK, VTK_NORMAL, VTK_DEBUG, &
      VTK_PART_DIAMETER, VTK_PART_VEL,VTK_PART_ANGULAR_VEL, VTK_PART_ORIENTATION, &
      VTK_PART_USR_VAR, VTK_PART_TEMP, VTK_PART_X_s, VTK_PART_COHESION, &
      VTK_FILEBASE,VTK_DATA, VTK_CUTCELL_ONLY, &
      VTK_NXS, VTK_NYS, VTK_NZS,VTK_SLICE_TOL, VTK_SELECT_MODE, &
      PG_OPTION,&
      CG_SAFE_MODE,PRINT_WARNINGS,&
      CG_UR_FAC,&
      PRINT_PROGRESS_BAR,BAR_WIDTH, BAR_CHAR, BAR_RESOLUTION,  &
      WRITE_DASHBOARD, F_DASHBOARD, &
      CPX,CPY,CPZ,NCX,NCY,NCZ,ERX,ERY,ERZ, &
      FIRST_DX,LAST_DX,FIRST_DY,LAST_DY,FIRST_DZ,LAST_DZ, &
      RE_INDEXING, ADJUST_PROC_DOMAIN_SIZE, REPORT_BEST_DOMAIN_SIZE,&
      NODESI_REPORT,NODESJ_REPORT,NODESK_REPORT,MINIMIZE_SEND_RECV, &
      DWALL_BRUTE_FORCE
