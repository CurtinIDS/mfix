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

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! None

! Local variables
!---------------------------------------------------------------------//
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
      DOUBLE PRECISION, SAVE :: OUT_TIME
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 0.1d0

      IF(FIRST_PASS) THEN
         CALL WRITE_DES_Tp(S_TIME)
         OUT_TIME = OUT_dT
         FIRST_PASS = .FALSE.
      ELSE
         IF(S_TIME >= OUT_TIME) THEN
            CALL WRITE_DES_Tp(S_TIME)
            OUT_TIME = OUT_TIME + OUT_dT
         ENDIF
      ENDIF

      RETURN

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_URS2                                               !
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
      SUBROUTINE WRITE_DES_Tp(lTime)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use fldvar
      Use run
      Use usr
      Use param1

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION lTime

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME1, FNAME2
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS1, F_EXISTS2
! file unit for heat transfer data
      INTEGER, PARAMETER :: TP_UNIT1 = 2030
      INTEGER, PARAMETER :: TP_UNIT2 = 2031

      DOUBLE PRECISION K11, K21, K31, K41
      DOUBLE PRECISION K12, K22, K32, K42
      DOUBLE PRECISION, PARAMETER :: RK4_DT_DEFAULT = 0.00001
      DOUBLE PRECISION RK4_DT
      DOUBLE PRECISION RK4_DT_LAST
      DOUBLE PRECISION TIME_INTERVAL

      DOUBLE PRECISION REL_ERR1, REL_ERR2

      DOUBLE PRECISION, SAVE :: RK4_TIME = 0.0d0
      DOUBLE PRECISION, SAVE :: Tp1_RK4 = 298.15000792d0
      DOUBLE PRECISION, SAVE :: Tp2_RK4 = 453.14999544d0

      DOUBLE PRECISION, PARAMETER :: lC1 = 2.70024635866355*(10.0**(-10))
      DOUBLE PRECISION, PARAMETER :: lC2 = 1.55333291786026*(10.0**(-10))

      INTEGER IJK
      DOUBLE PRECISION EPxTg, hEs, Tenv

      INTEGER I, RK4_STEPS

      LOGICAL NOISY

! Open the files.
      FNAME1 = 'POST_TP1.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=TP_UNIT1,FILE=FNAME1,STATUS='NEW')
         WRITE(TP_UNIT1,1000)
      ELSE
         OPEN(UNIT=TP_UNIT1,FILE=FNAME1,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
 1000 FORMAT(7X,'Time',11X,'Tp1',10X,'Tp1_MFIX',6X,'REL ERR 1')

      FNAME2 = 'POST_TP2.dat'
      INQUIRE(FILE=FNAME2,EXIST=F_EXISTS2)
      IF (.NOT.F_EXISTS2) THEN
         OPEN(UNIT=TP_UNIT2,FILE=FNAME2,STATUS='NEW')
         WRITE(TP_UNIT2,1001)
      ELSE
         OPEN(UNIT=TP_UNIT2,FILE=FNAME2,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
 1001 FORMAT(7X,'Time',11X,'Tp2',10X,'Tp2_MFIX',6X,'REL ERR 2')


! Calculate the value for the RK4 solutions.
      TIME_INTERVAL = lTime - RK4_TIME
      IF(TIME_INTERVAL .LE. RK4_DT_DEFAULT) THEN
         RK4_STEPS = 1
         RK4_DT = TIME_INTERVAL
         RK4_DT_LAST = UNDEFINED
      ELSE
         RK4_STEPS = floor(real(TIME_INTERVAL/RK4_DT_DEFAULT))
         RK4_DT = RK4_DT_DEFAULT
         RK4_DT_LAST = lTime - (RK4_TIME + RK4_STEPS*RK4_DT)
      ENDIF


      NOISY = .FALSE.
      IF(NOISY) THEN
         write(*,*) ''
         write(*,*) ''
         write(*,*) ' Entering RK4 calculations:'
         write(*,*) '    lTime: ', lTime
         write(*,*) '    RK4_TIME: ', RK4_TIME
         write(*,*) '    TIME_INTERVAL: ', TIME_INTERVAL
         write(*,*) '    RK4_STEPS: ', RK4_STEPS
         write(*,*) '    RK4_DT :', RK4_DT
         write(*,*) '    RK4_DT_LAST: ',RK4_DT_LAST
      ENDIF

      IJK = PIJK(1,4)
      EPxTg = EP_g(IJK) * T_g(IJK)
      hEs = 0.5d0*(ONE - EP_g(IJK))

      DO I=1, RK4_STEPS

         Tenv = EPxTg + hEs*(TP1_RK4 + TP2_RK4)
         K11  = RK4_DT * lC1 * (Tenv**4 - TP1_RK4**4)
         K12  = RK4_DT * lC2 * (Tenv**4 - TP2_RK4**4)

         Tenv = EPxTg + hEs*((TP1_RK4 + K11/2.0d0) + (TP2_RK4 + K12/2.0d0))
         K21  = RK4_DT * lC1 * (Tenv**4 - (TP1_RK4 + K11/2.0d0)**4 )
         K22  = RK4_DT * lC2 * (Tenv**4 - (TP2_RK4 + K12/2.0d0)**4)

         Tenv = EPxTg + hEs*((TP1_RK4 + K21/2.0d0) + (TP2_RK4 + K22/2.0d0))
         K31  = RK4_DT*lC1*(Tenv**4 - (TP1_RK4 + K21/2.0d0)**4)
         K32  = RK4_DT*lC2*(Tenv**4 - (TP2_RK4 + K22/2.0d0)**4)

         Tenv = EPxTg + hEs*((TP1_RK4 + K31) + (TP2_RK4 + K32))
         K41  = RK4_DT*lC1*(Tenv**4 - (TP1_RK4 + K31)**4)
         K42  = RK4_DT*lC2*(Tenv**4 - (TP2_RK4 + K32)**4)

         TP1_RK4 = TP1_RK4 + (K11 + 2.0d0*K21 + 2.0d0*K31 + K41)/6.0d0
         TP2_RK4 = TP2_RK4 + (K12 + 2.0d0*K22 + 2.0d0*K32 + K42)/6.0d0

         RK4_TIME = RK4_TIME + RK4_DT

      ENDDO


      IF(RK4_DT_LAST .NE. UNDEFINED) THEN

         Tenv = EPxTg + hEs*(TP1_RK4 + TP2_RK4)
         K11 = RK4_DT_LAST*lC1*(Tenv**4 - TP1_RK4**4)
         K12 = RK4_DT_LAST*lC2*(Tenv**4 - TP2_RK4**4)

         Tenv = EPxTg + hEs*((TP1_RK4 + K11/2.0d0) + (TP2_RK4 + K12/2.0d0))
         K21 = RK4_DT_LAST*lC1*(Tenv**4 - (TP1_RK4+K11/2.0d0)**4)
         K22 = RK4_DT_LAST*lC2*(Tenv**4 - (TP2_RK4+K12/2.0d0)**4)

         Tenv = EPxTg + hEs*((TP1_RK4 + K21/2.0d0) + (TP2_RK4 + K22/2.0d0))
         K31 = RK4_DT_LAST*lC1*(Tenv**4 - (TP1_RK4+K21/2.0d0)**4)
         K32 = RK4_DT_LAST*lC2*(Tenv**4 - (TP2_RK4+K22/2.0d0)**4)

         Tenv = EPxTg + hEs*((TP1_RK4 + K31) + (TP2_RK4 + K32))
         K41 = RK4_DT_LAST*lC1*(Tenv**4 - (TP1_RK4+K31)**4)
         K42 = RK4_DT_LAST*lC2*(Tenv**4 - (TP2_RK4+K32)**4)

         TP1_RK4 = TP1_RK4 + (K11 + 2.0d0*K21 + 2.0d0*K31 + K41)/6.0d0
         TP2_RK4 = TP2_RK4 + (K12 + 2.0d0*K22 + 2.0d0*K32 + K42)/6.0d0

         RK4_TIME = RK4_TIME + RK4_DT_LAST

      ENDIF

      IF(NOISY) THEN
         write(*,*) ' Leaving RK4 calculations:'
         write(*,*) '    RK4_TIME: ', RK4_TIME
         write(*,*) ' Last calculation values:'
         write(*,*)'    K11: ',K11
         write(*,*)'    K21: ',K21
         write(*,*)'    K31: ',K31
         write(*,*)'    K41: ',K41
         write(*,*)'    K12: ',K12
         write(*,*)'    K22: ',K22
         write(*,*)'    K32: ',K32
         write(*,*)'    K42: ',K42
         write(*,*)'    Tp1_RK4: ',Tp1_RK4
         write(*,*)'    Tp2_RK4: ',Tp2_RK4
      ENDIF

      if(Tp1_RK4 .NE. Tp1_RK4 .OR. Tp2_RK4 .NE. TP2_RK4) then
         write(*,*)' NAN found... exiting'
         stop
      endif


! Calculate the relative error.
       REL_ERR1 = (ABS(TP1_RK4 - DES_T_s(1))/ABS(TP1_RK4))*100
       REL_ERR2 = (ABS(TP2_RK4 - DES_T_s(2))/ABS(TP2_RK4))*100

! Write the data to a file.
      WRITE(TP_UNIT1,"(4(3X,F12.8))")lTime,TP1_RK4,DES_T_s(1),REL_ERR1
      CLOSE(TP_UNIT1)

      WRITE(TP_UNIT2,"(4(3X,F12.8))")lTime,TP2_RK4,DES_T_s(2),REL_ERR2
      CLOSE(TP_UNIT2)




      END SUBROUTINE WRITE_DES_Tp

      END SUBROUTINE USR2_DES
