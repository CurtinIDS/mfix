!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: STIFF_CHEM_DEBUG                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STIFF_CHEM_DBG

      PRIVATE

! Variable Access:
!---------------------------------------------------------------------//
      PUBLIC :: ODE_DEBUG_LEVEL

! Subroutine Access:
!---------------------------------------------------------------------//
      PUBLIC :: ALLOCATE_STIFF_CHEM_DBG

      PUBLIC :: UPDATE_ODE_OLD
      PUBLIC :: RESET_ODE

      PUBLIC :: CHECK_ODE_DATA
      PUBLIC :: WRITE_ODE_LOG

! Routine used to compare to values.
      LOGICAL, external :: COMPARE

! Static variables/parameters.
!---------------------------------------------------------------------//
! Debug Level: (Messages)
! 0 - None
! 1 - Limited
! 2 - Aggressive
      INTEGER :: ODE_DEBUG_LEVEL = 1

! File unit for ODE Error Log.
      INTEGER, parameter :: OEL_Unit = 6589

! Original field variables. Used to reset after failed integration.
      DOUBLE PRECISION, allocatable :: ODE_VARS_o(:)

      LOGICAL :: PRINT_ERR_HEADER

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_STIFF_CHEM_STATS                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_STIFF_CHEM_DBG(lODE_DIMN)

      implicit none

! The number ODEs (maximum).
      INTEGER, intent(in) :: lODE_DIMN

      allocate( ODE_VARS_o(lODE_DIMN) )

      RETURN
      END SUBROUTINE ALLOCATE_STIFF_CHEM_DBG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: ODE_UPDATE_OLD                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE UPDATE_ODE_OLD(lODE_DIMN, lODE_VARS)

      implicit none

! The number ODEs (maximum).
      INTEGER, intent(in) :: lODE_DIMN

! MFIX field variables mapped into the ODE format.
      DOUBLE PRECISION, intent(in) :: lODE_VARS(lODE_DIMN)

      ODE_VARS_o = lODE_VARS

      RETURN
      END SUBROUTINE UPDATE_ODE_OLD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: ODE_RESET                                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE RESET_ODE(lDIMN, lVARS, lAtps)

      use compar, only: myPE

      use stiff_chem_stats, only: failedCount

      implicit none

! The number ODEs (maximum).
      INTEGER, intent(in) :: lDIMN

! MFIX field variables mapped into the ODE format.
      DOUBLE PRECISION, intent(inout) :: lVARS(lDIMN)

      INTEGER, intent(in) :: lAtps

      lVARS = ODE_VARS_o

      IF(lAtps == 3) failedCount(myPE) = failedCount(myPE) + 1

      RETURN
      END SUBROUTINE RESET_ODE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_ODE_DATA                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ODE_DATA(lnD, lNEQ, loD, VARS, lUNLIMITED,     &
         lState, lErr)

      use param1,   only : ZERO
      use physprop, only : NMAX, MMAX
      use toleranc, only : TMax, TMin

      implicit none

      INTEGER, intent(in) :: lnD

      INTEGER, intent(in) :: lNEQ(lnD)

! The number ODEs (maximum).
      INTEGER, intent(in) :: loD

! MFIX field variables mapped into the ODE format.
      DOUBLE PRECISION, intent(in) :: VARS(loD)
! Flag for unlimited steps.
      LOGICAL, intent(in) :: lUNLIMITED
! Sate value returned from ODEPACK
      INTEGER, intent(in)  :: lState
! Error Flag
      INTEGER, intent(out) :: lErr


! Loop indicies:
      INTEGER :: M  ! phase
      INTEGER :: N  ! species
      INTEGER :: Node  ! other


      lErr = 0

! Check ODE_VAR for NaNs.
!---------------------------------------------------------------------//
      DO Node=1,loD
         IF(VARS(node).NE.VARS(node)) THEN
            lErr = -8
            return
         ENDIF
      ENDDO


! Check ODEPACK State
!---------------------------------------------------------------------//
      IF(lState == -1) THEN
         lErr = lState
         if(lUNLIMITED) return
      ELSEIF(lState == 2) THEN
         lErr = 0
      ELSE
         lErr = lState
      ENDIF


! Initialize.
      Node = 1

! Check variables for physical values.
!---------------------------------------------------------------------//

! Gas phase density.
      IF(VARS(Node) <= ZERO) THEN
         lErr = -(100 + Node)
         return
      ENDIF
      Node = Node + 1


! Gas phase temperature.
      IF(VARS(Node) > Tmax .OR. VARS(Node) < Tmin) THEN
         lErr = -(100 + Node)
         return
      ENDIF
      Node = Node + 1

! Gas phase species mass fractions.
      DO N=1,NMAX(0)
         IF(VARS(Node) > 1.009d0) THEN
            lErr = -(100 + Node)
            return
         ENDIF
         IF(VARS(Node) < 0.0d0) THEN
            IF(abs(VARS(Node)) > 1.0d-5) THEN
               lErr = -(100 + Node)
               return
            ENDIF
         ENDIF
         Node = Node + 1
      ENDDO

! Solids temperature.
      DO M = 1, MMAX
         IF(VARS(Node) > Tmax .OR. VARS(Node) < Tmin) THEN
            lErr = -(100 + Node)
            return
         ENDIF
         Node = Node + 1
      ENDDO


      DO M = 1, MMAX
         IF(lNEQ(2+M) == 1) THEN

! Solids volume fraction.
            IF(VARS(Node) <= ZERO) THEN
               lErr = -(100 + Node)
               return
            ENDIF
            Node = Node + 1

! Solids phase species mass fractions.
            DO N=1,NMAX(M)

               IF(VARS(Node) > 1.009d0) THEN
                  lErr = -(100 + Node)
                  return
               ENDIF
               IF(VARS(Node) < 0.0d0) THEN
                  IF(abs(VARS(Node)) > 1.0d-5) THEN
                     lErr = -(100 + Node)
                     return
                  ENDIF
               ENDIF
               Node = Node + 1

            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CHECK_ODE_DATA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ODE_ErrorLog                                           !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_ODE_LOG(lErr, lnD, lNEQ, loD, lVARS)

      use compar,   only : myPE
      use run,      only : TIME
      use physprop, only : MMAX

      use indices

      implicit none

      INTEGER, intent(in) :: lErr

      INTEGER, intent(in) :: lnD
      INTEGER, intent(in) :: lNEQ(lnD)

! The number ODEs (maximum).
      INTEGER, intent(in) :: loD
      DOUBLE PRECISION, intent(in) :: lVARS(loD)

      INTEGER :: lc1
      INTEGER :: IJK

      DOUBLE PRECISION :: ddt_lVARS(loD)

      CHARACTER(LEN=255) :: lFile

      LOGICAL :: lExist

      IJK = lNEQ(2)

      lFile = ''; write(lFile,"('StiffChem_',I2.2,'.log')") myPE

      inquire(file=lFile,exist=lExist)
      if(lExist) then
         open(OEL_Unit,file=trim(adjustl(lFile)), position='append',convert='big_endian')
      else
         open(OEL_Unit,file=trim(adjustl(lFile)), status='new',convert='big_endian')
      endif


      if(PRINT_ERR_HEADER) then
         WRITE(OEL_Unit,9000) Time
         PRINT_ERR_HEADER = .FALSE.
      endif


      WRITE(OEL_Unit,9001) IJK, I_OF(IJK), J_OF(IJK), myPE
      WRITE(OEL_Unit,9002) lErr
      WRITE(OEL_Unit,9003) lNEQ(1)

      DO lc1=1,MMAX
         IF(lNEQ(lc1+2) == 1) THEN
            WRITE(OEL_Unit,9004) lc1
         ELSEIF(lNEQ(lc1+2) == 0) THEN
            WRITE(OEL_Unit,9005) lc1
         ELSE
            WRITE(OEL_Unit,9006) lc1, lNEQ(lc1+2)
         ENDIF
      ENDDO

      WRITE(OEL_Unit,"(/5x,'Field Variables:')")
      WRITE(OEL_Unit,9007)
      DO lc1=1, loD
         WRITE(OEL_Unit,9008) lc1, ODE_VARS_o(lc1), lVARS(lc1), &
            (lVARS(lc1) - ODE_VARS_o(lc1))
      ENDDO

      CALL STIFF_CHEM_RRATES(lNEQ, Time, lVARS, ddt_lVARS)

      WRITE(OEL_Unit,"(/5x,'Derivatives:')")
      WRITE(OEL_Unit,9009)
      DO lc1=1, loD
         WRITE(OEL_Unit,9010) lc1, ddt_lVARS(lc1)
      ENDDO


      RETURN

 9000 FORMAT(//3x,'Time: ',g15.8)
 9001 FORMAT(//5x,'Fluid Cell: ',I4,4x,'(',I3,',',I3,')',10x,&
         'Process: ',I4)
 9002 FORMAT(  5x,'Error Flag/ODEPACK Status: ',I4,/)

 9003 FORMAT(  5x,'Number of ODEs Solved: ', I4,/)
 9004 FORMAT(  5x,'Solving solids phase ',I2,&
         ' bulk density and species equations.')

 9005 FORMAT(  5x,'NOT solving solids phase ',I2,&
         ' bulk density and species equations.')

 9006 FORMAT(  5x,'Unknown flag for solids phase ',I2,'  Flag = ',I4)


 9007 FORMAT(  5x,'Var',6x,'Incoming',13x,'Returned',13x,'Difference')
 9008 FORMAT(  5x,I3,3(3x,G18.8))


 9009 FORMAT(  5x,'Var',6x,'  YDOT  ')
 9010 FORMAT(  5x,I3,3x,G18.8)

      END SUBROUTINE WRITE_ODE_LOG


      END MODULE STIFF_CHEM_DBG
