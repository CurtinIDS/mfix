!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: DEallocate_ARRAYS
!  Purpose: deallocate arrays
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98
!  Reviewer:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE Deallocate_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE ambm
      USE cont
      USE drag
      USE energy
      USE fldvar
      USE geometry
      USE ghdtheory
      USE indices
      USE kintheory
      USE mflux
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE residual
      USE run
      USE rxns
      USE scalars
      USE tau_g
      USE tau_s
      USE trace
      USE visc_g
      USE visc_s
      USE vshear
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER M
!-----------------------------------------------

!ambm
      Deallocate( A_m )
      Deallocate( B_m )

!cont
      Deallocate( DO_CONT )

!drag
      Deallocate(  F_gs )
      Deallocate(  F_ss )

!energy
      Deallocate(  HOR_g  )
      Deallocate(  HOR_s  )
      Deallocate(  GAMA_gs  )
      Deallocate(  GAMA_Rg  )
      Deallocate(  GAMA_Rs  )
      Deallocate(  T_Rg  )
      Deallocate(  T_Rs  )

!fldvar
      Deallocate(  EP_g  )
      Deallocate(  epg_jfac  )
      Deallocate(  epg_ifac  )
      Deallocate(  eps_ifac  )
      Deallocate(  EP_go  )
      Deallocate(  P_g  )
      Deallocate(  P_go  )
      Deallocate(  RO_g  )
      Deallocate(  RO_go  )
      Deallocate(  ROP_g  )
      Deallocate(  ROP_go  )
      Deallocate(  ROP_s  )
      Deallocate(  ROP_so  )

!      Deallocate(  RO_SV )
!      Deallocate(  RO_SVo )
      Deallocate(  EP_SS )
      Deallocate(  ERR_ARRAY )

      Deallocate(  T_g  )
      Deallocate(  T_s  )
      Deallocate(  T_go  )
      Deallocate(  T_so  )
      Deallocate(  X_g  )
      Deallocate(  X_s  )
      Deallocate(  X_go  )
      Deallocate(  X_so  )
      Deallocate(  U_g  )
      Deallocate(  U_go  )
      Deallocate(  U_s  )
      Deallocate(  U_so  )
      Deallocate(  V_g  )
      Deallocate(  V_go  )
      Deallocate(  V_s  )
      Deallocate(  V_so  )
      Deallocate(  W_g  )
      Deallocate(  W_go  )
      Deallocate(  W_s  )
      Deallocate(  W_so  )
      Deallocate(  P_s  )
      Deallocate(  P_s_c  )
      Deallocate(  P_s_v  )
      Deallocate(  P_s_f  )
      Deallocate(  P_s_p  )
      Deallocate(  P_star  )
      Deallocate(  P_staro  )
      Deallocate(  THETA_m  )
      Deallocate(  THETA_mo  )

      IF(DIMENSION_Scalar /= 0)then
        Deallocate(  Scalar  )
        Deallocate(  Scalaro  )
      ENDIF

      IF(K_Epsilon)then
        Deallocate(  K_Turb_G  )
        Deallocate(  E_Turb_G  )
      ENDIF

!geometry
      Deallocate(  FLAG  )
      Deallocate(  FLAG_E  )
      Deallocate(  FLAG_N  )
      Deallocate(  FLAG_T  )
      Deallocate(  ICBC_FLAG  )
      Deallocate(  oDX  )
      Deallocate(  oDY  )
      Deallocate(  oDZ  )
      Deallocate(  oDX_E  )
      Deallocate(  oDY_N  )
      Deallocate(  oDZ_T  )
      Deallocate(  X  )
      Deallocate(  X_E  )
      Deallocate(  oX  )
      Deallocate(  oX_E  )
      Deallocate(  Z  )
      Deallocate(  Z_T  )
      Deallocate(  FX  )
      Deallocate(  FX_bar  )
      Deallocate(  FX_E  )
      Deallocate(  FX_E_bar  )
      Deallocate(  FY_N  )
      Deallocate(  FY_N_bar  )
      Deallocate(  FZ_T  )
      Deallocate(  FZ_T_bar  )
      Deallocate(  AYZ  )
      Deallocate(  AXZ  )
      Deallocate(  AXY  )
      Deallocate(  VOL  )
      Deallocate(  AYZ_U  )
      Deallocate(  AXZ_U  )
      Deallocate(  AXY_U  )
      Deallocate(  VOL_U  )
      Deallocate(  AYZ_V  )
      Deallocate(  AXZ_V  )
      Deallocate(  AXY_V  )
      Deallocate(  VOL_V  )
      Deallocate(  AYZ_W  )
      Deallocate(  AXZ_W  )
      Deallocate(  AXY_W  )
      Deallocate(  VOL_W  )

      Deallocate(  IJK_ARRAY_OF, DEAD_CELL_AT  )
      Deallocate(  ro_s, ro_so )


!indices
      Deallocate(  STORE_LM  )
      Deallocate(  CELL_CLASS  )
      Deallocate(  I_OF  )
      Deallocate(  J_OF  )
      Deallocate(  K_OF  )
      Deallocate(  Im1  )
      Deallocate(  Ip1  )
      Deallocate(  Jm1  )
      Deallocate(  Jp1  )
      Deallocate(  Km1  )
      Deallocate(  Kp1  )

!pgcor
      Deallocate(  d_e )
      Deallocate(  d_n )
      Deallocate(  d_t )
      Deallocate(  Pp_g )
      Deallocate(  PHASE_4_P_g )

!physprop
      Deallocate(  MU_g  )
      Deallocate(  C_pg  )
      Deallocate(  C_ps  )
      Deallocate(  K_g  )
      Deallocate(  K_s  )
      Deallocate(  Kth_s  )
      Deallocate(  Kphi_s  )
      Deallocate(  DIF_g  )
      Deallocate(  DIF_s  )
      Deallocate(  MW_MIX_g  )

!pscor
      Deallocate(  e_e )
      Deallocate(  e_n )
      Deallocate(  e_t )
      Deallocate(  K_cp )
      Deallocate(  EPp )
      Deallocate(  PHASE_4_P_s )

!residual
      Deallocate( RESID )
      Deallocate( MAX_RESID )
      Deallocate( IJK_RESID )
      deallocate( NUM_RESID )   ! added 17-mar-2008
      deallocate( den_resid )   ! added 17-mar-2008
      deallocate( resid_pack )  ! added 17-mar-2008

!rxns
      if (nRR .gt. 0) Deallocate( ReactionRates )
      Deallocate(  R_gp  )
      Deallocate(  R_sp  )
      Deallocate(  RoX_gc  )
      Deallocate(  RoX_sc  )
      Deallocate(  SUM_R_g  )
      Deallocate(  SUM_R_s  )
      Deallocate(  R_phase  )

!scalars
      IF(DIMENSION_Scalar /= 0)then
        Deallocate(  Scalar_c  )
        Deallocate(  Scalar_p  )
        Deallocate(  Dif_Scalar  )
      ENDIF

!dqmom
      deallocate(  D_p  )
      deallocate(  D_po )
      deallocate(  Source_a)
      deallocate(  S_bar)
      deallocate(  Matrix_a)
      deallocate(  Matrix_b)
      deallocate(  Matrix_c)
      deallocate(  Inv_a)
      deallocate(  A)
      deallocate(  omega)
      deaLLocate(  beta_a)
      deaLLocate(  ystart)

!tau_g
      Deallocate(  TAU_U_g )
      Deallocate(  TAU_V_g )
      Deallocate(  TAU_W_g )
      Deallocate(  DF_gu )
      Deallocate(  DF_gv )
      Deallocate(  DF_gw )
      Deallocate(  CTAU_U_G )
      Deallocate(  CTAU_V_G )
      Deallocate(  CTAU_W_G )

!tau_s
      Deallocate(  TAU_U_s )
      Deallocate(  TAU_V_s )
      Deallocate(  TAU_W_s )

!trace
      Deallocate(  trD_s_C  )
      Deallocate(  trD_s2  )
      Deallocate(  trD_s_Co  )
      Deallocate(  trD_s_Co2  )

!visc_g
      Deallocate(  trD_g )
      Deallocate(  MU_gt  )
      Deallocate(  EPMU_gt  )
      Deallocate(  LAMBDA_gt  )
      Deallocate(  EPLAMBDA_gt  )
      Deallocate(  L_scale  )

!visc_s
      Deallocate(  trD_s )
      Deallocate(  MU_s  )
      Deallocate(  EPMU_s  )
      Deallocate(  LAMBDA_s  )
      Deallocate(  EPLAMBDA_s  )
      Deallocate(  ALPHA_s  )
      Deallocate(  MU_s_c  )
      Deallocate(  LAMBDA_s_c  )
      Deallocate(  LAMBDA_s_v )
      Deallocate(  LAMBDA_s_f )
      Deallocate(  LAMBDA_s_p )
      Deallocate(  MU_s_v )
      Deallocate(  MU_s_f )
      Deallocate(  MU_s_p )
      Deallocate(  MU_b_v )
      Deallocate(  EP_star_array )
      Deallocate(  EP_g_blend_start )
      Deallocate(  EP_g_blend_end )
      Deallocate(  I2_devD_s )
      Deallocate(  TrM_s )
      Deallocate(  TrDM_s )

!vshear
      Deallocate(  VSH )
      Deallocate(  VSHE )

!kintheory
!additional general kinetic theory terms
      Deallocate(  KTMOM_U_s)
      Deallocate(  KTMOM_V_s)
      Deallocate(  KTMOM_W_s)

!kintheory
!kinetic theory terms Iddir & Arastoopour (2005)
      IF (TRIM(KT_TYPE) == 'IA_NONEP') THEN
        Deallocate(  trD_s2_ip)
        Deallocate(  MU_sM_ip)
        Deallocate(  MU_sL_ip)
        Deallocate(  XI_sM_ip)
        Deallocate(  XI_sL_ip)
        Deallocate(  Fnu_s_ip)
        Deallocate(  FT_sM_ip)
        Deallocate(  FT_sL_ip)
        Deallocate(  Kth_sL_ip)
        Deallocate(  Knu_sM_ip)
        Deallocate(  Knu_sL_ip)
        Deallocate(  Kvel_s_ip)
        Deallocate(  EDvel_sL_ip)
        Deallocate(  ED_ss_ip)
      ENDIF
!kinetic theory terms Iddir & Arastoopour (2005) or Garzo Dufty (1999)
      IF (TRIM(KT_TYPE) == 'IA_NONEP' .OR. TRIM(KT_TYPE) == 'GD99') THEN
        Deallocate(  EDvel_sM_ip)
        Deallocate(  EDT_s_ip)
      ENDIF

! ghdtheory
! GHD theory
      IF (TRIM(KT_TYPE) == 'GHD') THEN
         Deallocate(  Flux_nE)
         Deallocate(  Flux_nN)
         Deallocate(  Flux_nT)
         Deallocate(  Zeta0)
         Deallocate(  ZetaU)
         Deallocate(  DiT)
         Deallocate(  DijF)
         Deallocate(  Lij)
         Deallocate(  Dij)
         Deallocate(  DijQ)
         Deallocate(  JoiX)
         Deallocate(  JoiY)
         Deallocate(  JoiZ)
         Deallocate(  FiX)
         Deallocate(  FiY)
         Deallocate(  FiZ)
      ENDIF

! array allocation of add on packages, such as linear equation solvers

! array allocation for higher order implementation
      Deallocate(  FLAG3 )
      Deallocate(  CELL_CLASS3 )
      Deallocate(  I3_OF )
      Deallocate(  J3_OF )
      Deallocate(  K3_OF )
      Deallocate(  Im1_3 )
      Deallocate(  Ip1_3 )
      Deallocate(  Jm1_3 )
      Deallocate(  Jp1_3 )
      Deallocate(  Km1_3 )
      Deallocate(  Kp1_3 )

!mflux
      Deallocate(  Flux_gE)
      Deallocate(  Flux_sE)
      Deallocate(  Flux_gN)
      Deallocate(  Flux_sN)
      Deallocate(  Flux_gT)
      Deallocate(  Flux_sT)
      Deallocate(  ROP_gE)
      Deallocate(  ROP_sE)
      Deallocate(  ROP_gN)
      Deallocate(  ROP_sN)
      Deallocate(  ROP_gT)
      Deallocate(  ROP_sT)

!spill over from interp_res.f
      deallocate( ijksize3_all )
      deallocate( ijkstart3_all )
      deallocate( ijkend3_all )
      deallocate( ijksize4_all )
      deallocate( ijkstart4_all )
      deallocate( ijkend4_all )

      deallocate( istart_all)
      deallocate( jstart_all )
      deallocate( kstart_all )
      deallocate( istart1_all )
      deallocate( jstart1_all )
      deallocate( kstart1_all )
      deallocate( istart2_all )
      deallocate( jstart2_all )
      deallocate( kstart2_all )
      deallocate( istart3_all )
      deallocate( jstart3_all )
      deallocate( kstart3_all )
      deallocate( istart4_all )
      deallocate( jstart4_all )
      deallocate( kstart4_all )

      deallocate( iend_all )
      deallocate( jend_all )
      deallocate( kend_all )
      deallocate( iend1_all )
      deallocate( jend1_all )
      deallocate( kend1_all )
      deallocate( iend2_all )
      deallocate( jend2_all )
      deallocate( kend2_all )
      deallocate( iend3_all )
      deallocate( jend3_all )
      deallocate( kend3_all )
      deallocate( iend4_all )
      deallocate( jend4_all )
      deallocate( kend4_all )

      deallocate( displs )

      deallocate( imap)
      deallocate( jmap )
      deallocate( kmap )
      deallocate( imap_c )
      deallocate( jmap_c )
      deallocate( kmap_c )


      RETURN
      END SUBROUTINE Deallocate_ARRAYS

