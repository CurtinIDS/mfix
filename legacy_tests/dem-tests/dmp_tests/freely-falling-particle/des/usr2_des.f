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
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 1.0d-3

      IF(FIRST_PASS) THEN
         CALL WRITE_DES_OUT(S_TIME)
         OUT_TIME = OUT_dT
         FIRST_PASS = .FALSE.
      ELSE
         IF(S_TIME >= OUT_TIME) THEN
            CALL WRITE_DES_Out(S_TIME)
            OUT_TIME = OUT_TIME + OUT_dT
         ENDIF
      ENDIF

      RETURN

      contains

!......................................................................!
!  Subroutine: WRITE_DES_Out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_DES_Out(lTime)

      use compar
      use desmpi
      use discretelement
      use mpi_utility
      use parallel
      use run
      use sendrecv
      use usr
      use mpi_comm_des, only: des_gather

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION lTime


! Local variables
!---------------------------------------------------------------------//
      character*30 :: filename
! logical used for testing is the data file already exists
      LOGICAL :: exists
! file unit for heat transfer data
      INTEGER, PARAMETER :: uPos = 2030
      INTEGER, PARAMETER :: uVel = 2031

! Absolute relative error between MFIX solution and analytic solution.

      integer :: aStage

      double precision :: lGrav

      integer, allocatable :: lgID(:)

      double precision, allocatable :: lRad(:)

      double precision, allocatable :: lPos_Y(:),  lVel_Y(:)
      double precision              :: aPos_Y   ,  aVel_Y
      double precision              :: Pos_rErr ,  Vel_rErr

! Variables related to gather
      integer llocalcnt
      integer lglocnt
      integer lgathercnts(0:numpes-1)
      integer lproc

      integer lc1, lc2

! Initialize the global count.
      lglocnt = 10

! Calculate the number of real (non-ghost) particles on local process.
!   pip and ighost_cnt come from the discretelement module
      llocalcnt = pip - ighost_cnt
      call global_sum(llocalcnt, lglocnt)

! Allocate the send buffers for each process
      allocate( dprocbuf(llocalcnt) )
      allocate( iprocbuf(llocalcnt) )

! Allocate the receive buffer for root process.
      allocate( drootbuf(lglocnt) )
      allocate( irootbuf(lglocnt) )

! Set the send count variable.
      igath_sendcnt = llocalcnt ! (desmpi_mod)

! Set the gather counts array.
      lgathercnts = 0
      lgathercnts(mype) = llocalcnt
      call global_sum(lgathercnts,igathercnts)

! Set the gather displacement array.
      idispls(0) = 0
      do lproc = 1,numpes-1
         idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
      end do

! Allocate local storage variables.
      allocate( lRad   (lglocnt) )
      allocate( lPos_Y (lglocnt) )
      allocate( lVel_Y (lglocnt) )
      allocate( lgID   (lglocnt) )

! Gather global particle IDs for map.
      call des_gather(iglobal_id)
      if (mype.eq.pe_io) lgID = irootbuf

! Gather particle diameters
      call des_gather(des_radius)
      if (mype.eq.pe_io) lRad = drootbuf

! Gather particle position (Y-axis only)
      call des_gather(des_POS_new(:,2))
      if (mype.eq.pe_io) lPos_Y = drootbuf

! Gather particle position (Y-axis only)
      call des_gather(des_VEL_new(:,2))
      if (mype.eq.pe_io) lVel_Y = drootbuf

! Set local variables.
      lGrav = -grav(2)

      if(mype.eq.pe_io) then

         do lc1=1, lglocnt

            lc2 = outMap(lgID(lc1))
            if(lc1 /= lc2) write(*,"(3x,'Order switched: ',(3x,I4))")  &
               lc1, lc2

! Open the files.
            filename = ''
            write(filename,"('POST_POS_',I1,'.dat')") lc2
            inquire(file=filename, exist=exists)
            if(exists) then
               open(unit=uPos, file=filename,&
                  position="APPEND", status='OLD')
            else
               open(unit=uPos, file=filename, status='NEW')
               write(uPos,"(3X,'Time, Stage, Pos, Pos_MFIX, aErr')")
            endif

            filename = ''
            write(filename,"('POST_VEL_',I1,'.dat')") lc2
            inquire(file=filename, exist=exists)
            if(exists) then
               open(unit=uVel, file=filename,&
                  position="APPEND",status='OLD')
            else
               open(unit=uVel, file=filename, status='NEW')
               write(uVel,"(3X,'Time, Stage, Vel, Vel_MFIX, aErr')")
            endif

! Initialize the stage.
            aStage = 0
! Calculate the position and velocity of the particle
! Stage 1: Free fall
            if(lTime < time_c(lc2)) then
               aStage = 1
               aPos_Y = y_s1(lc2, h0(lc2), lGrav, lTime)
               aVel_Y = dydt_s1(lc2, lGrav, lTime)
! Stage 2: Contact
            elseif( lTime < time_r(lc2)) then
               aStage = 2
               aPos_Y = y_s2(lc2, h0(lc2), lRad(lc2), b_r(lc2),        &
                  w0_r(lc2), lGrav, lTime)
               aVel_Y = dydt_s2(lc2, h0(lc2), lRad(lc2), b_r(lc2),     &
                  w0_r(lc2), lGrav, lTime)
! Stage 3: Rebound
            else
               aStage = 3
               aPos_Y = y_s3(lc2, lRad(lc2), lGrav, lTime)
               aVel_Y = dydt_s3(lc2, lGrav, lTime)
            endif

! Calculate the absolute relative error.
             Pos_rErr = (ABS(aPos_Y - lPos_Y(lc1))/ABS(aPos_Y))*100
             Vel_rErr = (ABS(aVel_Y - lVel_Y(lc1))/ABS(aVel_Y))*100

! Write the results to a file.
            WRITE(uPos,9000) lTime, aStage, aPos_Y, lPos_Y(lc1),Pos_rErr
            CLOSE(uPos)

            WRITE(uVel,9000) lTime, aStage, aVel_Y, lVel_Y(lc1),Vel_rErr
            CLOSE(uVel)

 9000 FORMAT(3x,F15.8,',',1X,I1,3(',',1X,F15.8))

         enddo
      endif

! Dellocate the send buffers for each process
      if(allocated(dprocbuf)) deallocate( dprocbuf )
      if(allocated(iprocbuf)) deallocate( iprocbuf )
      if(allocated(drootbuf)) deallocate( drootbuf )
      if(allocated(irootbuf)) deallocate( irootbuf )
      if(allocated(lRad    )) deallocate( lRad     )
      if(allocated(lPos_Y  )) deallocate( lPos_Y   )
      if(allocated(lVel_Y  )) deallocate( lVel_Y   )
      if(allocated(lgID    )) deallocate( lgID     )


      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
