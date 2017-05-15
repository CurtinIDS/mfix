!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES

      Use discretelement, only: S_TIME, DTSOLID

      use output, only : USR_DT
      use param, only: DIMENSION_USR

      IMPLICIT NONE

! Stored times for writing user data.
      DOUBLE PRECISION, DIMENSION(DIMENSION_USR), SAVE :: USR_TIME = 0.0
      DOUBLE PRECISION :: OVERSHOOT

! Calculate if there is overshoot.
      OVERSHOOT = S_TIME + 0.1d0*DTSOLID

! Write particle mass data.
      IF(OVERSHOOT >= USR_TIME(1)) THEN
! Update the time to write the special output.
         USR_TIME(1) = (INT((OVERSHOOT)/USR_DT(1))+1)*USR_DT(1)
! Write the output.
         CALL WRITE_RXN_DATA
      ENDIF

      RETURN
      END SUBROUTINE USR2_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_DES_TP                                           !
!                                                                      !
!  THIS ROUTINE IS NOT DESIGNED TO BE CHECKED INTO CVS. IT IS STRICTLY !
!  HERE FOR PRINTING MESSAGES USED FOR DEBUGGING AND V&V WORK.         !
!                                                                      !
!  Author: J.Musser                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_RXN_DATA

      use constant, only: C
      use constant, only: RGAS => GAS_CONST_cal

      use ic, only: IC_X_s

      use des_thermo, only: DES_T_S
      use des_rxns, only: DES_X_s

      use discretelement, only: S_TIME
      use discretelement, only: PMASS
      use discretelement, only: PIJK

      use fldvar, only: T_g


      Use param1, only: UNDEFINED

      use physprop, only: MW_g, MW_s
      use run, only: RUN_NAME

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: lUNIT = 2030

! Time to freeze the analytic solution.
      DOUBLE PRECISION :: MAX_TIME

      DOUBLE PRECISION :: M_p0, Mp
      DOUBLE PRECISION :: X_B0, XB
      DOUBLE PRECISION :: X_D0, XD
      DOUBLE PRECISION :: X_I0, XI

      DOUBLE PRECISION :: ABS_ERR
      DOUBLE PRECISION T_COMP

      DOUBLE PRECISION :: rAg, rCg
      DOUBLE PRECISION :: rBs, rDs, rIs

      INTEGER :: NP
      INTEGER :: M

      include "species.inc"

! Local Functions
!---------------------------------------------------------------------//
! Local function dummy variables.
      DOUBLE PRECISION :: lt, lTs, lTg

      DOUBLE PRECISION :: M_p ! Particle Mass
      DOUBLE PRECISION :: X_B ! Solids species B
      DOUBLE PRECISION :: X_D ! Solids species D
      DOUBLE PRECISION :: X_I ! Solids species I

      M_p(lt) = (M_p0 + lt*(rBs + rDs + rIs))

      X_B(lt) = (M_p0*X_B0 + rBs*lt)/M_p(lt)
      X_D(lt) = (M_p0*X_D0 + rDs*lt)/M_p(lt)
      X_I(lt) = (M_p0*X_I0 + rIs*lt)/M_p(lt)

!---------------------------------------------------------------------//

      NP = 2
      M  = PIJK(NP,5)

      M_p0 = 0.0353429173528852d0

      X_B0 = IC_X_s(1,M,Bs)
      X_D0 = IC_X_s(1,M,Ds)
      X_I0 = IC_X_s(1,M,Is)

      rAg = -1.0d0 * MW_g(Ag) * C(1)
      rCg =  1.0d0 * MW_g(Cg) * C(1)

      rBs = -2.0d0 * MW_s(M,Bs) * C(1)
      rDs =  1.0d0 * MW_s(M,Ds) * C(1)
      rIs =  0.0d0


! Freeze the analytic solution at MAX_TIME seconds (the approximate time
! when the reaction should stop due to a lack of species B reactant).
      MAX_TIME = -(M_p0*X_B0)/rBs
      T_COMP = merge(S_TIME, MAX_TIME, S_TIME < MAX_TIME)


      FNAME = 'POST_MASS.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=lUNIT,FILE=FNAME,STATUS='NEW')

         write(lUNIT,"(3X,A)")'Initial Conditions:'
         write(lUNIT,"(15X,A,7X,A)")'Analytic','MFIX-DEM'
         write(lUNIT,"(6X,A,3X,F12.8,3X,F12.8)")'Mass',m_p0,PMASS(2)

         write(lUNIT,"(6X,'MW_Bs: ',F12.8)") MW_s(M,Bs)
         write(lUNIT,"(6X,'MW_Ds: ',F12.8)") MW_s(M,Ds)
         write(lUNIT,"(6X,'Reaction Rate: ',F12.8)") C(1)

         write(lUNIT,"(//7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','Mp','Mp-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=lUNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Calculate the absolute error
      Mp = M_p(T_COMP)
      ABS_ERR = abs(Mp - PMASS(2))
      WRITE(lUNIT,"(4(3X,F12.8))")S_TIME, Mp, PMASS(2),ABS_ERR
      CLOSE(lUNIT)


! Solids Species B:
!.....................................................................//
      FNAME = 'POST_Xs_B.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=lUNIT,FILE=FNAME,STATUS='NEW')

         write(lUNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','XB','XB-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=lUNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      XB = X_B(T_COMP)
      ABS_ERR = ABS(XB - DES_X_s(2,Bs))
      WRITE(lUNIT,"(4(3X,F12.8))")S_TIME, XB, DES_X_s(2,Bs), ABS_ERR
      CLOSE(lUNIT)



! Solids Species D:
!.....................................................................//
      FNAME = 'POST_Xs_D.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=lUNIT,FILE=FNAME,STATUS='NEW')
         write(lUNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','XD','XD-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=lUNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      XD = X_D(T_COMP)
      ABS_ERR = ABS(XD - DES_X_s(2,Ds))
      WRITE(lUNIT,"(4(3X,F12.8))")S_TIME, XD, DES_X_s(2,Ds), ABS_ERR
      CLOSE(lUNIT)



! Solids Species I:
!.....................................................................//
      FNAME = 'POST_Xs_I.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)
      IF (.NOT.F_EXISTS) THEN
         OPEN(UNIT=lUNIT,FILE=FNAME,STATUS='NEW')

         write(lUNIT,"(7X,A,10X,A,11X,A,8X,A)")&
            'S_TIME','XI','XI-DEM','ABS ERR'
      ELSE
         OPEN(UNIT=lUNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      XI = X_I(T_COMP)
      ABS_ERR = ABS(XI - DES_X_s(2,3))
      WRITE(lUNIT,"(4(3X,F12.8))")S_TIME, XI, DES_X_s(2,Is), ABS_ERR
      CLOSE(lUNIT)

      RETURN
      END SUBROUTINE WRITE_RXN_DATA
