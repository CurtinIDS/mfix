!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS1_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms have been calculated but before they are     !
!  applied. The user may insert code in this routine or call user      !
!  defined subroutines.                                                !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR1_DES

  
      USE constant, only : gravity, gravity_x, gravity_y, Pi
      USE constant, only : c
      USE param1, only : ZERO, ONE, undefined
      Use physprop, only: mmax
      USE run, only: time, dt, tstop
      USE discretelement, only : grav, dtsolid, des_mmax
      USE des_thermo_cond, only : DES_QW_Cond
      Use usr, only: QW_SAMPLE_TIME, DES_QwFlux_AVG
      Use compar, only: ijkstart3, ijkend3
      use functions, only: fluid_at
      USE vtk, only : vtk_dt, write_vtk_files
      use cutcell, only : debug_cg
      

      IMPLICIT NONE

      INCLUDE 'usrnlst.inc'

      DOUBLE PRECISION :: RPM 
      DOUBLE PRECISION :: ANGLE
      INTEGER :: IJK, M
      LOGICAL :: CHECK_TIME
      
      ! Get RPM from input file - C(1)
      RPM = C(1)
      ! Roate gravity vector for rotating drum
      ANGLE = time * RPM*2.0D0*Pi/60.0D0 !convert RPM to rad/s
      GRAVITY_X = sin(ANGLE)*GRAVITY
      GRAVITY_Y = -cos(ANGLE)*GRAVITY 
      GRAV(1) = GRAVITY_X
      GRAV(2) = GRAVITY_Y

      ! Compute running average for wall heat transfer
      DO IJK=IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK))CYCLE
         DO M = MMAX+1, DES_MMAX+MMAX
            DES_QwFlux_AVG(IJK,M)=(DES_QwFlux_AVG(IJK,M)*QW_SAMPLE_TIME&
            &   +DES_Qw_COND(IJK,M)*DTSOLID)/(QW_SAMPLE_TIME+DTSOLID)
         ENDDO
      ENDDO
      QW_SAMPLE_TIME=QW_SAMPLE_TIME+DTSOLID

      ! CHECK TO SEE IF IT IS OUTPUT TIME
      CHECK_TIME = .FALSE.
      IF(WRITE_VTK_FILES) THEN       
         IF(DT == UNDEFINED) THEN
            CHECK_TIME = .FALSE.
         ELSE
            CHECK_TIME = (TIME+0.1d0*DT>=VTK_DT(1)).OR.(TIME+0.1d0*DT>=TSTOP)
         ENDIF
      ENDIF

      IF(CHECK_TIME)THEN
         ! VTU files will be written this time step, so store average
         ! heat flux into debug array that (for output)
         ! and reset the running average
         DO M = MMAX+1, DES_MMAX+MMAX
            DEBUG_CG(:,M) =  DES_QwFlux_AVG(:,M)
            DES_QwFlux_AVG(:,M) = ZERO
            QW_SAMPLE_TIME = ZERO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE USR1_DES
