!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use compar
      use desmpi
      use discretelement
      use exit, only: mfix_exit
      use mpi_comm_des, only: des_gather
      use mpi_utility
      use parallel
      use run
      use sendrecv
      use usr

      implicit none

! Particle properties:
      integer, allocatable :: lPhase(:)         ! solids phase index
      double precision, allocatable :: lRad(:)  ! radius
      double precision, allocatable :: lMass(:) ! mass

! Simulation variables:
      double precision :: lGrav                 ! gravity: y-component

! Iterative solver parameters.
      integer, parameter :: lc_max = 1000
      logical converged

! Variables related to gather
      integer llocalcnt
      integer lglocnt
      integer lgathercnts(0:numpes-1)
      integer lproc

      integer m, lc1, lc2                       ! loop indices
      double precision :: t_r, ly, ldydt        ! temp variables.

      character*30 :: filename

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
      do lproc = 1, numpes-1
         idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
      end do

! Record local data
      filename = ''
      write(filename,"(A,'_PE',I2.2,'.txt')") trim(RUN_NAME), myPE
      open(unit=201, file=filename, status='replace')
      if(llocalcnt > 0) then
         write(201,"(3x,'Local particle count: ',I5)") llocalcnt
         write(201,"(/3x,'Particle Data:')")
         write(201,"(/3x,'Global',3x,'Local',11x,'Position')")
         write(201,"(5x,'ID',6x,'ID',7x,'X',9x,'Y',9x,'Z')")
         do lc1 = 1, max_pip
            if(is_nonexistent(lc1)) cycle
            if(is_ghost(lc1)) cycle
            write(201,"(3x,I6,2x,I6,3(2x,F8.4))") iglobal_id(lc1),     &
               lc1, des_pos_new(lc1,1:3)
         end do
      else
         write(201,"(3x,'No particles on this process.')")
      endif
      close(201)

! This case is limited in its use.
      if(myPE == PE_IO) then
         if(particles /= 4) then
            write(*,"(3x, 'invalid setup for test case')")
            call mfix_exit(0)
         endif
         call init_usr_var(lglocnt, imax_global_id)
      endif

      allocate( lRad   (lglocnt) )
      allocate( lMass  (lglocnt) )
      allocate( lPhase (lglocnt) )

! Gather global particle IDs for map.
      call des_gather(iglobal_id)
      if (mype.eq.pe_io) then
            write(*,"(//3x,'outMap:')")
         do lc1=1, lglocnt
            outMap(irootbuf(lc1)) = lc1
            write(*,"(5x,'gID: ',I4,3x,'lc: ',I2)") irootbuf(lc1), lc1
         enddo
      endif

! Gather particle phase
      call des_gather(pijk(:,5))
      if (mype.eq.pe_io) lPhase = irootbuf

! Gather particle diameters
      call des_gather(des_radius)
      if (mype.eq.pe_io) lRad = drootbuf

! Gather particle mass
      call des_gather(pmass(:))
      if (mype.eq.pe_io) lMass = drootbuf

! Gather particle position (Y-axis only)
      call des_gather(des_pos_new(:,2))
      if (mype.eq.pe_io) h0 = drootbuf

      lGrav = -grav(2)

      if(mype.eq.pe_io) then

         do lc1=1, lglocnt

            m = lPhase(lc1)

! Calculate the start time of particle/wall collision.
            time_c(lc1) =  dsqrt(2.0d0*(h0(lc1) - lRad(lc1))/lGrav)
! Calculate the particle's velocity at the start of the collision.
            vel_c(lc1)  = -dsqrt(2.0d0*lGrav*h0(lc1))

            b_r(lc1)  = DES_ETAN_WALL(m)/(2.0d0*dsqrt(KN_W*lMass(lc1)))
            w0_r(lc1) = dsqrt(KN_W / lMass(lc1))

! Initial loop parameters.
            time_r(lc1) = time_c(lc1) + 0.05d0
            lc2 = 0
            converged = .false.


            write(*,"(//3x,'Particle:',I2)") lc1
            write(*,"(5x,'Kn: ',E9.3)") KN_W
            write(*,"(5x,'En: ',E9.3)") DES_EN_WALL_INPUT(m)

!            write(*,"(/3x,'Iter',7X,'Old',11x,'New',11x,'|Err|')")

! Calculate the start time for the rebound stage.
            do while (.not.converged)
               lc2 = lc2 + 1

               ly = y_s2(lc1, h0(lc1), lRad(lc1), b_r(lc1), w0_r(lc1), &
                  lGrav, time_r(lc1)) - lRad(lc1)

               ldydt = dydt_s2(lc1, h0(lc1), lRad(lc1), b_r(lc1),      &
                  w0_r(lc1), lGrav, time_r(lc1))

               if(ldydt /= 0.0d0) then
                  t_r = time_r(lc1) - ly/ldydt
               else
                  write(*,"(3x, 'Fatal Error: ldydt == 0')")
                  call mfix_exit(0)
               endif

               if(abs(t_r - time_r(lc1)) < 10.0d-16) then
                  converged = .true.
               else
                  converged = .false.
!                  write(*,"(3x,I4,3(3x,F11.8))") lc2, &
!                     t_r, time_r(lc1), abs(t_r - time_r(lc1))
               endif

               time_r(lc1) = t_r

               if(lc2 > lc_max) then
                  write(*,"(3x, 'Fatal Error: lc2 > lc_max')")
                  call mfix_exit(0)
               endif

            enddo
! Calculate the velocity at the end of the rebound stage.
            vel_r(lc1) = dydt_s2(lc1, h0(lc1), lRad(lc1), b_r(lc1),    &
               w0_r(lc1), lGrav, time_r(lc1))

            write(*,"(5x,'Collision:')")
            write(*,"(7x,'Time:     ',F14.8)") time_c(lc1)
            write(*,"(7x,'Velocity: ',F14.8)") vel_c(lc1)
            write(*,"(5x,'Rebound:')")
            write(*,"(7x,'Time: ',F14.8)") time_r(lc1)
            write(*,"(7x,'Velocity: ',F14.8)") vel_r(lc1)
         enddo
         write(*,"(//' ')")
      endif

! Dellocate the send buffers for each process
      if(allocated(dprocbuf)) deallocate( dprocbuf )
      if(allocated(iprocbuf)) deallocate( iprocbuf )
      if(allocated(drootbuf)) deallocate( drootbuf )
      if(allocated(irootbuf)) deallocate( irootbuf )
      if(allocated(lRad    )) deallocate( lRad     )
      if(allocated(lMass   )) deallocate( lMass    )
      if(allocated(lPhase  )) deallocate( lPhase   )


      return

      END SUBROUTINE USR0
