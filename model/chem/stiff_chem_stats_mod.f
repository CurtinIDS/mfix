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
      MODULE STIFF_CHEM_STATS

      PRIVATE

! Variable Access:
!---------------------------------------------------------------------//
      PUBLIC :: failedCount


! Subroutine Access:
!---------------------------------------------------------------------//
      PUBLIC :: ALLOCATE_STIFF_CHEM_STATS

      PUBLIC :: INIT_STIFF_CHEM_STATS
      PUBLIC :: UPDATE_STIFF_CHEM_STATS
      PUBLIC :: WRITE_STIFF_CHEM_STATS


! Routine used to compare to values.
      LOGICAL, external :: COMPARE


! Static variables/parameters.
!---------------------------------------------------------------------//
! Frequency to report the number of steps distribution.
      INTEGER, parameter :: reportNST_Freq = 10

! Variables updated once each call to the stiff solver.
!---------------------------------------------------------------------//
! Frequency to report the number of steps distribution.
      INTEGER :: reportNST

      INTEGER :: failedCount_total
      INTEGER :: countINCPT_total


! Variables updated every IJK loop cycle.
!---------------------------------------------------------------------//

! The minimum number of integrations needed (over all IJK)
      INTEGER, allocatable :: minNST(:)                     ! local
      INTEGER, allocatable :: minNST_all(:)                 ! global

! The maximum number of integrations needed (over all IJK)
      INTEGER, allocatable :: maxNST(:)                     ! local
      INTEGER, allocatable :: maxNST_all(:)                 ! global

! An array that stores the distrubtion of the number of steps needed
! to integrate ODES.
      INTEGER, allocatable :: countNST(:)                   ! local
      INTEGER, allocatable :: countNST_all(:)               ! global

! Number of cells that only have homogeneous chemical reactions.
      INTEGER, allocatable :: Homogns(:)                    ! local
      INTEGER, allocatable :: Homogns_all(:)                ! global

! Number of cells that only have homogeneous and/or heterogeneous
! chemical reactions.
      INTEGER, allocatable :: Hetrgns(:)                    ! local
      INTEGER, allocatable :: Hetrgns_all(:)                ! global

! Number of cells that failed to successfully integration ODEs.
      INTEGER, allocatable :: failedCount(:)                ! local
      INTEGER, allocatable :: failedCount_all(:)            ! global

! Maximum number of attempts to integrate.
      INTEGER, allocatable :: maxAttempts(:)                ! local
      INTEGER, allocatable :: maxAttempts_all(:)            ! global

! Maximum number of incomplete integrations.
      INTEGER, allocatable :: countINCPT(:)                 ! local
      INTEGER, allocatable :: countINCPT_all(:)             ! global

      DOUBLE PRECISION :: ODE_StartTime

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
      SUBROUTINE ALLOCATE_STIFF_CHEM_STATS

      use compar, only: myPE
      use compar, only: numPEs
      use output, only: FULL_LOG

      implicit none

      reportNST = 1
      failedCount_total = 0
      countINCPT_total = 0

! Number of cells that failed to successfully integration ODEs.
      allocate( failedCount(0:numPEs-1) ); failedCount = 0    ! local
      allocate( failedCount_all(0:numPEs-1) )                 ! global

      if(.NOT.FULL_LOG) return

! The minimum number of integrations needed (over all IJK)
      allocate( minNST(0:numPEs-1)); minNST = 0               ! local
      allocate( minNST_all(0:numPEs-1) )                      ! global
      minNST(myPE) = 5000

! The maximum number of integrations needed (over all IJK)
      allocate( maxNST(0:numPEs-1) ); maxNST = 0              ! local
      allocate( maxNST_all(0:numPEs-1) )                      ! global

! An array that stores the distrubtion of the number of steps needed
! to integrate ODES.
      allocate( countNST(5) ); IF(reportNST==1) countNST = 0  ! local
      allocate( countNST_all(5) )                             ! global

! Number of cells that only have homogeneous chemical reactions.
      allocate( Homogns(0:numPEs-1) ); Homogns = 0            ! local
      allocate( Homogns_all(0:numPEs-1) )                     ! global

! Number of cells that only have homogeneous and/or heterogeneous
! chemical reactions.
      allocate( Hetrgns(0:numPEs-1) ); Hetrgns = 0;           ! local
      allocate( Hetrgns_all(0:numPEs-1) )                     ! global

! Maximum number of attempts to integrate.
      allocate( maxAttempts(0:numPEs-1) ); maxAttempts = 0    ! local
      allocate( maxAttempts_all(0:numPEs-1) )                 ! global

! Number of cells that fail to completely integrate the time step
! given the maximum number of steps.
      allocate( countINCPT(0:numPEs-1) ); countINCPT = 0      ! local
      allocate( countINCPT_all(0:numPEs-1) )                  ! global


      RETURN
      END SUBROUTINE ALLOCATE_STIFF_CHEM_STATS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: INIT_ODE_STATS0                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_STIFF_CHEM_STATS

      use compar, only: myPE
      use compar, only: PE_IO
      use output, only: FULL_LOG

      implicit none

      if(.NOT.FULL_LOG) return

      CALL CPU_TIME(ODE_StartTime)

      Hetrgns = 0
      Homogns = 0
      failedCount = 0
      countINCPT = 0

      if(myPE == PE_IO) &
         write(*,"(/3x,'Integrating stiff chemistry...')",ADVANCE="NO")

      RETURN
      END SUBROUTINE INIT_STIFF_CHEM_STATS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: UPDATE_STIFF_CHEM_STATS                                !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE UPDATE_STIFF_CHEM_STATS(lNEQ, lNEQ_DIMN, lNST, &
         lODE_DIMN, lAtps, lIncpt)

      use compar, only: myPE

      implicit none

! The number of steps needed to integrate.
      INTEGER, intent(in) :: lNEQ_DIMN

! (1) :: Number of ODEs
! (2) :: Fluid cell index (IJK) passed into ODEPACK
      INTEGER, dimension(lNEQ_DIMN), intent(in) :: lNEQ
! The number of steps needed to integrate.
      INTEGER, intent(in) :: lNST
! The number ODEs (maximum).
      INTEGER, intent(in) :: lODE_DIMN

! The number of attempts.
      INTEGER, intent(in) :: lAtps

! Flag that the integration is incomplete
      LOGICAL, intent(in) :: lIncpt


      IF(lNEQ(1) == lODE_DIMN) THEN
         Hetrgns(myPE) = Hetrgns(myPE) + 1
      ELSE
         Homogns(myPE) = Homogns(myPE) + 1
      ENDIF

      maxAttempts(myPE) = max(lAtps, maxAttempts(myPE))

      minNST(myPE) = min(minNST(myPE), lNST)
      maxNST(myPE) = max(maxNST(myPE), lNST)

      IF (lNST <           10) THEN
         countNST(1) = countNST(1) + 1
      ELSE IF (lNST <     100) THEN
         countNST(2) = countNST(2) + 1
      ELSE IF (lNST <    1000) THEN
         countNST(3) = countNST(3) + 1
      ELSE IF (lNST <   10000) THEN
         countNST(4) = countNST(4) + 1
      ELSE
         countNST(5) = countNST(5) + 1
      ENDIF

      IF(lIncpt) countINCPT(myPE) = countINCPT(myPE) + 1

      RETURN
      END SUBROUTINE UPDATE_STIFF_CHEM_STATS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: WRITE_ODE_STATS                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_STIFF_CHEM_STATS

      use compar, only: myPE
      use compar, only: PE_IO
      use output, only: FULL_LOG

      use mpi_utility

      implicit none

! Message buffer.
      CHARACTER(LEN=64) :: lMsg0, lMsg1

      DOUBLE PRECISION :: lODE_EndTime, lODE_RunTime

      IF(.NOT.FULL_LOG) return


! Update screen message.
      IF(myPE == PE_IO) WRITE(*,"(2x,'DONE.',/)")


      CALL CPU_TIME(lODE_EndTime)
      lODE_RunTime = lODE_EndTime - ODE_StartTime


! Collect stats on min/max number of steps.
      minNST_all = 0; CALL global_sum(minNST, minNST_all)
      maxNST_all = 0; CALL global_sum(maxNST, maxNST_all)

! Collect stats on the number of cells with pure homogeneous reactions.
      Homogns_all = 0;
      CALL global_sum(Homogns, Homogns_all)

! Collect stats on the number of cells with heterogeneous and
! homogeneous reactions.
      Hetrgns_all = 0
      CALL global_sum(Hetrgns, Hetrgns_all)

! Collect stats on the maximum number of integration attempts.
      maxAttempts_all = 0
      CALL global_sum(maxAttempts, maxAttempts_all)

! Collect stats on the maximum number of incomplete integrations.
      countINCPT_all = 0
      CALL global_sum(countINCPT, countINCPT_all)

! Collect stats on the number of failed integrations.
      failedCount_all = 0
      CALL global_sum(failedCount, failedCount_all)


! Display stiff solver summary.
      IF(myPE == PE_IO) THEN

! Report Min/Max steps:
         lMsg0=''; write(lMsg0,*) minval(minNST_all)
         lMsg1=''; write(lMsg1,*) maxval(maxNST_all)
         write(*,1000)  trim(adjustl(lMsg0)), trim(adjustl(lMsg1))

! Report Homogeneous/Heterogeneous reactions:
         lMsg0=''; write(lMsg0,*) sum(Homogns_all)
         lMsg1=''; write(lMsg1,*) sum(Hetrgns_all)
         write(*,1001) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))

! Report Max attempts:
         lMsg0=''; write(lMsg0,*) maxval(maxAttempts_all)
         write(*,1004)  trim(adjustl(lMsg0))

! Report incomplete integrations:
         countINCPT_total = countINCPT_total + sum(countINCPT_all)

         IF(countINCPT_total > 0) THEN
            lMsg0=''; write(lMsg0,*) sum(countINCPT_all)
            lMsg1=''; write(lMsg1,*) countINCPT_total
            write(*,1002) 'incomplete', trim(adjustl(lMsg0)), trim(adjustl(lMsg1))
         ENDIF

! Report failed integrations:
         failedCount_total = failedCount_total + sum(failedCount_all)

         IF(failedCount_total > 0) THEN
            lMsg0=''; write(lMsg0,*) sum(failedCount_all)
            lMsg1=''; write(lMsg1,*) failedCount_total
            write(*,1002) 'failed', trim(adjustl(lMsg0)), trim(adjustl(lMsg1))
         ENDIF

         IF(lODE_RunTime > 3.6d3) THEN
            lMsg0=''; write(lMsg0,"(f8.4)") lODE_RunTime/3.6d3
            lMsg1='hrs'
         ELSEIF(lODE_RunTime > 6.0d1) THEN
            lMsg0=''; write(lMsg0,"(f8.4)") lODE_RunTime/6.0d1
            lMsg1='min'
         ELSE
            lMsg0=''; write(lMsg0,"(f8.4)") lODE_RunTime
            lMsg1='sec'
         ENDIF
         write(*,1003) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))

      ENDIF


      if(reportNST == reportNST_Freq) then
! Collect the number of steps distributions.
         countNST_all = 0;
         CALL global_sum(countNST, countNST_all)

         countNST_all = int(countNST_all/reportNST_Freq)

         if(myPE == PE_IO) then
            write(*,"(/5x,'Average Integration Distribution:')")
            write(*,"(7x,'NST < 10^1: ', I6)")countNST_all(1)
            write(*,"(7x,'NST < 10^2: ', I6)")countNST_all(2)
            write(*,"(7x,'NST < 10^3: ', I6)")countNST_all(3)
            write(*,"(7x,'NST < 10^4: ', I6)")countNST_all(4)
            write(*,"(7x,'NST > 10^5: ', I6)")countNST_all(5)
         endif
! Reset the reporting counter.
         reportNST = 1
! Clear out old data.
         countNST = 0
         countNST_all = 0
      else
! Increment the reporting counter.
         reportNST = reportNST + 1
      endif

      if(myPE == PE_IO)write(*,"(/' ')")

      RETURN

 1000 Format(5x,'Minimum/Maximum number of steps over all cells: ',A,'/',A)
 1001 Format(5x,'Number of cells with Homogeneous/Heterogeneous reactions: ',A,'/',A)
 1002 Format(5x,'Number of Current/Cumulative ',A,' integrations: ',A,'/',A)
 1003 Format(5x,'CPU Time Used: ',A,' ',A)
 1004 Format(5x,'Maximum number of integration attempts: ',A)

      END SUBROUTINE WRITE_STIFF_CHEM_STATS

      END MODULE STIFF_CHEM_STATS
