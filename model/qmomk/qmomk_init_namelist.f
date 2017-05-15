!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_INIT_NAMELIST                                    C
!  Purpose: Initializie namelist for QMOM                              C
!                                                                      C
!                                                                      C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!  Comments:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

SUBROUTINE qmomk_init_namelist

  USE param1
  USE qmom_kinetic_equation
  Implicit none
  INCLUDE 'qmomknamelist.inc'

  QMOMK = .FALSE.
  QMOMK_TYPE = 'QMOM'
  QMOMK_WALL_BC_TYPE = 'SPECULAR_REFLECTIVE'
  QMOMK_COLLISIONS = 'BGK'
  QMOMK_COUPLED = .TRUE.
  QMOMK_CFL = 0.4
  QMOMK_COLLISIONS_ORDER = 0
 RETURN

END SUBROUTINE qmomk_init_namelist
