!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Module name: ps_mod.f                                               C
!                                                                      C
!  Purpose: Common block containing point source data.                 C
!                                                                      C
!  Author: J. Musser                                  Date: 10-Jun-13  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE ps

      use param, only: dimension_ps, dim_m, dim_n_g, dim_n_s

! Run-time logical indicating that point sources are present.
      LOGICAL :: POINT_SOURCE

      LOGICAL :: PS_DEFINED(DIMENSION_PS)

! Physical location of point sources.
      DOUBLE PRECISION :: PS_X_w(DIMENSION_PS)  ! West
      DOUBLE PRECISION :: PS_X_e(DIMENSION_PS)  ! East
      DOUBLE PRECISION :: PS_Y_s(DIMENSION_PS)  ! South
      DOUBLE PRECISION :: PS_Y_n(DIMENSION_PS)  ! North
      DOUBLE PRECISION :: PS_Z_b(DIMENSION_PS)  ! Bottom
      DOUBLE PRECISION :: PS_Z_t(DIMENSION_PS)  ! Top

! Cell indices delineating point source region:
      INTEGER :: PS_I_w(DIMENSION_PS)  ! West
      INTEGER :: PS_I_e(DIMENSION_PS)  ! East
      INTEGER :: PS_J_s(DIMENSION_PS)  ! South
      INTEGER :: PS_J_n(DIMENSION_PS)  ! North
      INTEGER :: PS_K_b(DIMENSION_PS)  ! Bottom
      INTEGER :: PS_K_t(DIMENSION_PS)  ! Top

! Gas mass flow rate through the point source:
      DOUBLE PRECISION PS_MASSFLOW_g (DIMENSION_PS)

! Velocity vector for gas point source: (normalized)
      DOUBLE PRECISION :: PS_U_g(DIMENSION_PS) ! X-axis
      DOUBLE PRECISION :: PS_V_g(DIMENSION_PS) ! Y-axis
      DOUBLE PRECISION :: PS_W_g(DIMENSION_PS) ! Z-axis

! Gas phase velocity magnitude: (calculated)
      DOUBLE PRECISION :: PS_VEL_MAG_G(DIMENSION_PS)

! Gas phase species mass fractions
      DOUBLE PRECISION :: PS_X_g(DIMENSION_PS, DIM_N_g)

! Gas phase temperature.
      DOUBLE PRECISION :: PS_T_g(DIMENSION_PS)
      DOUBLE PRECISION :: PS_CpxMFLOW_g(DIMENSION_PS)

! Solids mass flow rate through the point source:
      DOUBLE PRECISION PS_MASSFLOW_s (DIMENSION_PS, DIM_M)

! Velocity vector for solids point sources: (normalized)
      DOUBLE PRECISION :: PS_U_s(DIMENSION_PS, DIM_M) ! X-axis
      DOUBLE PRECISION :: PS_V_s(DIMENSION_PS, DIM_M) ! Y-axis
      DOUBLE PRECISION :: PS_W_s(DIMENSION_PS, DIM_M) ! Z-axis

! Solids phase velocity magnitude: (calculated)
      DOUBLE PRECISION :: PS_VEL_MAG_S(DIMENSION_PS, DIM_M)

! Solids phase species mass fractions
      DOUBLE PRECISION :: PS_X_s(DIMENSION_PS, DIM_M, DIM_N_s)

! Solids phase temperature.
      DOUBLE PRECISION :: PS_T_s(DIMENSION_PS, DIM_M)
      DOUBLE PRECISION :: PS_CpxMFLOW_s(DIMENSION_PS, DIM_M)

! Total volume of cells comprising point source cells (calculated)
      DOUBLE PRECISION :: PS_VOLUME(DIMENSION_PS)

! Legacy variable... to be deleated
      INTEGER, DIMENSION(:), ALLOCATABLE :: POINT_SOURCES


      END MODULE ps
