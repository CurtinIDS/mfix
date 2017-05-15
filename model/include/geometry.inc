!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: GEOMETRY_LOCKED                                           !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / GEOMETRY_LOCKED /&

      COORDINATES, &
      IMAX, DX, XLENGTH, NO_I, CYCLIC_X, CYCLIC_X_PD, XMIN, &
      JMAX, DY, YLENGTH, NO_J, CYCLIC_Y, CYCLIC_Y_PD,       &
      KMAX, DZ, ZLENGTH, NO_K, CYCLIC_Z, CYCLIC_Z_PD,       &
      CYLINDRICAL_2D, I_CYL_NUM, I_CYL_TRANSITION,          &
      SHEAR, V_SH



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: GEOMETRY_UNLOCKED                                         !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / GEOMETRY_UNLOCKED /&

      FLUX_G, DELP_X, DELP_Y, DELP_Z
