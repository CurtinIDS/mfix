!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ADJUST_THETA                                            C
!  Purpose: Remove small negative values of theta caused by linear     C
!           solvers                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 02-APR-98  C
!                                                                      C
!  Modified: S. Benyahia                              Date: 02-AUG-06  C
!  Purpose: check for small negative numbers at walls                  C
!           (not just fluid cells)                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADJUST_THETA(M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: zero
      USE toleranc, only: zero_ep_s
      USE constant, only: pi, to_si
! granular temperature of solids phase m
      USE fldvar, only: theta_m
! material density of solids phase m
      USE fldvar, only: ro_s
! particle diameter of solids phase m
      USE fldvar, only: d_p
! number of solids phases
      USE physprop, only: smax
! kt types
      USE run, only: kt_type
      USE run, only: kt_type_enum
      USE run, only: lun_1984
      USE run, only: simonin_1996
      USE run, only: ahmadi_1995
      USE run, only: gd_1999
      USE run, only: gtsh_2012
      USE run, only: ia_2005
      USE run, only: ghd_2007
! needed for function.inc
      USE compar
      USE exit, only: mfix_exit
      USE functions
      USE geometry
      USE indices

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase
      INTEGER, INTENT(IN) :: M
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
! Solids phase index
      INTEGER :: L
! Particle mass and diameter for use with those kinetic theories
! that include mass of particle in definition of theta
      DOUBLE PRECISION :: M_PM, D_PM
! small value of theta_m
      DOUBLE PRECISION :: smallTheta
!-----------------------------------------------

      IER = 0
      smallTheta = (to_SI)**4 * ZERO_EP_S

      DO IJK = IJKSTART3, IJKEND3
        IF ( FLUID_AT(IJK) ) THEN

          SELECT CASE(KT_TYPE_ENUM)
            CASE (LUN_1984, SIMONIN_1996, AHMADI_1995, GD_1999, &
                  GTSH_2012)
              IF (THETA_M(IJK,M) < smallTheta) &
                 THETA_M(IJK,M) = smallTheta

            CASE (IA_2005)
              D_PM = D_P(IJK,M)
              M_PM = (PI/6.d0)*(D_PM**3)*RO_S(IJK,M)
              IF (THETA_M(IJK,M) < smallTheta*M_PM) &
                THETA_M(IJK,M) = smallTheta*M_PM

            CASE (GHD_2007)
              M_PM = ZERO
              DO L = 1,SMAX
                D_PM = D_P(IJK,L)
                M_PM = M_PM +(PI/6.d0)*(D_PM**3)*RO_S(IJK,L)
              ENDDO
              M_PM = M_PM/DBLE(SMAX)
              IF (THETA_M(IJK,M) < smallTheta*M_PM) &
                THETA_M(IJK,M) = smallTheta*M_PM

            CASE DEFAULT
! should never hit this
               WRITE (*, '(A)') 'ADJUST_THETA'
               WRITE (*, '(A,A)') 'Unknown KT_TYPE: ', KT_TYPE
               call mfix_exit(myPE)
          END SELECT   ! end selection of kt_type_enum
        ENDIF   ! end if (fluid_at)
      ENDDO  ! end do ijk

      RETURN
      END SUBROUTINE ADJUST_THETA

