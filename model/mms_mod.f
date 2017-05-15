!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MMS_SRC                                                !
!  Purpose: Global storage containers for MMS variables.               !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE mms

! By default the MMS functions are unavailable.
      LOGICAL :: USE_MMS = .FALSE.

!! Method of Manufactured Solutions (MMS) and Tecplot variables :

! Gas volume fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Ep_g

! Gas pressure
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_P_g

! Gas velocity components
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_g

! Gas temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_g

! Solid bulk density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_s

! Solids velocity components
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_s

! Solids temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_s

! Granular temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Theta_m

! Gas continuity MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_g_Src

! Gas Momentum source terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_g_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_g_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_g_Src

! Gas energy equation MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_g_Src

! Solid continuity MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_s_Src

! Solid momentum MMS source terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_s_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_s_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_s_Src

! Solid energy equation MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_s_Src

! Granular energy MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Theta_m_Src

! Temporary variable for pressure shifting while plotting and
! discretization error norm calculation
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_G_Sh

! Index for pressure shifting
      INTEGER :: IJK_Sh

!  x, y, z locations of top-right corner of a cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  xtr
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ytr
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ztr

      contains



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_MMS_VARS                                      !
!  Purpose: Allocate memory for allocatable variables defined inside   !
!  MMS module.                                                         !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_MMS_VARS
      use param, only   : DIMENSION_3

      IMPLICIT NONE

! Note: both fluid and solid phase variables can be of size DIMENSION_3
! since MMS tests do not work for more than one solid phase.

        allocate(MMS_Ep_g(DIMENSION_3))
        allocate(MMS_P_g(DIMENSION_3))
        allocate(MMS_U_g(DIMENSION_3))
        allocate(MMS_V_g(DIMENSION_3))
        allocate(MMS_W_g(DIMENSION_3))
        allocate(MMS_T_g(DIMENSION_3))

        allocate(MMS_ROP_s(DIMENSION_3))
        allocate(MMS_U_s(DIMENSION_3))
        allocate(MMS_V_s(DIMENSION_3))
        allocate(MMS_W_s(DIMENSION_3))
        allocate(MMS_T_s(DIMENSION_3))
        allocate(MMS_Theta_m(DIMENSION_3))

        allocate(MMS_ROP_g_Src(DIMENSION_3))
        allocate(MMS_U_g_Src(DIMENSION_3))
        allocate(MMS_V_g_Src(DIMENSION_3))
        allocate(MMS_W_g_Src(DIMENSION_3))
        allocate(MMS_T_g_Src(DIMENSION_3))

        allocate(MMS_ROP_s_Src(DIMENSION_3))
        allocate(MMS_U_s_Src(DIMENSION_3))
        allocate(MMS_V_s_Src(DIMENSION_3))
        allocate(MMS_W_s_Src(DIMENSION_3))
        allocate(MMS_T_s_Src(DIMENSION_3))
        allocate(MMS_Theta_m_Src(DIMENSION_3))

        allocate(P_g_Sh(DIMENSION_3))

        allocate(xtr(DIMENSION_3))
        allocate(ytr(DIMENSION_3))
        allocate(ztr(DIMENSION_3))


      RETURN
      END SUBROUTINE ALLOCATE_MMS_VARS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_MMS_VARS                                    !
!  Purpose: Deallocate memory for allocatable variables defined inside !
!  MMS module.                                                         !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_MMS_VARS

      IMPLICIT NONE


        deallocate(MMS_Ep_g)
        deallocate(MMS_P_g)
        deallocate(MMS_U_g)
        deallocate(MMS_V_g)
        deallocate(MMS_W_g)
        deallocate(MMS_T_g)

        deallocate(MMS_ROP_s)
        deallocate(MMS_U_s)
        deallocate(MMS_V_s)
        deallocate(MMS_W_s)
        deallocate(MMS_T_s)
        deallocate(MMS_Theta_m)

        deallocate(MMS_ROP_g_Src)
        deallocate(MMS_U_g_Src)
        deallocate(MMS_V_g_Src)
        deallocate(MMS_W_g_Src)
        deallocate(MMS_T_g_Src)

        deallocate(MMS_ROP_s_Src)
        deallocate(MMS_U_s_Src)
        deallocate(MMS_V_s_Src)
        deallocate(MMS_W_s_Src)
        deallocate(MMS_T_s_Src)
        deallocate(MMS_Theta_m_Src)

        deallocate(P_g_Sh)

        deallocate(xtr)
        deallocate(ytr)
        deallocate(ztr)

      RETURN
      END SUBROUTINE DEALLOCATE_MMS_VARS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:CALCULATE_MMS                                           !
!  Purpose: Generic stub/placeholder for MMS analytical solutions.     !
!                                                                      !
!  Author: J.Musser                                   Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALCULATE_MMS

      IMPLICIT NONE

      RETURN
      END SUBROUTINE CALCULATE_MMS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:CALCULATE_MMS_SOURCE                                    !
!  Purpose: Generic stub/placeholder for MMS source terms.             !
!                                                                      !
!  Author: J.Musser                                   Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALCULATE_MMS_SOURCE

      IMPLICIT NONE

      RETURN
      End SUBROUTINE CALCULATE_MMS_SOURCE

      END MODULE mms
