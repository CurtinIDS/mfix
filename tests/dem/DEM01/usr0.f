!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR0                                                   !
!  Author: J.Musser                                   Date: dd-mmm-yy  !
!  Purpose: This routine is called before the time loop starts and is  !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting constants and checking errors in   !
!           data.  This routine is not called from an IJK loop, hence  !
!           all indices are undefined.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use discretelement
      use exit, only: mfix_exit
      use usr

      IMPLICIT NONE

      double precision :: t_r

      double precision :: lRad
      double precision :: lGrav

      double precision :: ly, ldydt


      if(particles /= 1) then
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      endif

      lRad  = des_radius(1)
      h0    = des_pos_new(1,2)
      lGrav = -grav(2)

! Calculate the start time of particle/wall collision.
      time_c =  dsqrt(2.0d0*(h0 - lRad)/lGrav)
! Calculate the particle's velocity at the start of the collision.
      vel_c  = -dsqrt(2.0d0*lGrav*h0)

      b_r  = DES_ETAN_WALL(1)/ (2.0d0 * dsqrt(KN_W * PMASS(1)))
      w0_r = dsqrt(KN_W / PMASS(1))

      write(*,"(//3x,'Parameters:')")
      write(*,"(5x,'Kn ',F14.0)") KN_W
      write(*,"(5x,'En ',F14.8)") DES_EN_WALL_INPUT(1)
      write(*,"(5x,'Br ',F14.8)") b_r
      write(*,"(5x,'Dt ',F14.8)") DTSOLID


! Calculate the reboud time using Brents method.
      if(KN_W >= 5.0d4) then
         time_r = BrentsMethod(0.300d0, 0.350d0)
      elseif(KN_W >= 2.5d4) then
         time_r = BrentsMethod(0.325d0, 0.375d0)
      else
         time_r = BrentsMethod(0.375d0, 0.425d0)
      endif

! Calculate the velocity at the end of the rebound stage.
      vel_r = dydt_s2(h0, lRad, b_r, w0_r, lGrav, time_r)

      write(*,"(//3x,'Collision:')")
      write(*,"(5x,'Time:     ',F14.8)") time_c
      write(*,"(5x,'Velocity: ',F14.8)") vel_c

      write(*,"(/3x,'Rebound:')")
      write(*,"(5x,'Time: ',F14.8)") time_r
      write(*,"(5x,'Velocity: ',F14.8)") vel_r

      write(*,"(//' ')")

      return

      contains

!......................................................................!
!  Function: BrentsMethod                                              !
!  Purpose: Locate a root between two values.                          !
!......................................................................!
      double precision function BrentsMethod(ap, bp)

      integer :: lc
      double precision :: ap, bp
      double precision :: ak, bk, ck, dk, s
      double precision :: fa, fb, fc, fs

      logical :: mflag

      double precision, parameter :: root_tol = 1.0d-8
      double precision, parameter :: delta = 1.0d-8
      logical :: converged

      double precision :: tmp0

      ak = ap
      bk = bp

      fa = y_s2(h0, lRad, b_r, w0_r, lGrav, ak) - lRad
      fb = y_s2(h0, lRad, b_r, w0_r, lGrav, bk) - lRad


      if(abs(fa) < abs(fb)) then
         tmp0 = ak
         ak = bk
         bk = tmp0
         tmp0 = fa
         fa = fb
         fb = tmp0
      endif

! Verify that ak and bk bound the root
      if(fa*fb >= 0) then
         write(*,"(3x,'Fatal Error: fa*fb >= 0')")
         write(*,"(3x,'fa: ',g12.4,3x,'a: ',g12.4)") fa, ak
         write(*,"(3x,'fb: ',g12.4,3x,'b: ',g12.4)") fb, bk
         call mfix_exit(0)
      endif

      ck = ak
      converged = .false.
      mflag = .true.
      lc = 1
      do while (.not.converged .and. lc < 100)
         fc = y_s2(h0, lRad, b_r, w0_r, lGrav, ck) - lRad
         if(fa /= fc .and. fb /= fc) then
            s = ak*fb*fc/((fa-fb)*(fa-fc)) + &
                bk*fa*fc/((fb-fa)*(fb-fc)) + &
                ck*fa*fb/((fc-fa)*(fc-fb))
         else
            s = bk - fb*(bk-ak)/(fb-fa)
         endif

         if((s < ((3.0*ak + bk)/4.0) .OR. s > bk) .OR. &
            (mflag .and. abs(s-bk) >= abs(bk-ck)/2.0d0) .OR. &
            (.not.mflag .and. abs(s-bk) >= (ck-dk)/2.0d0) .OR. &
            (mflag .and. abs(bk-ck) < delta) .OR. &
            (.not.mflag .and. abs(ck-dk) < delta)) then
            s = (ak + bk)/2.0d0
            mflag = .true.
         else
            mflag = .false.
         endif
         fs = y_s2(h0, lRad, b_r, w0_r, lGrav, s) - lRad

         dk = ck
         ck = bk
         if(fa*fs < 0) then
            bk = s
            fb = fs 
         else
            ak = s
            fa = fs
         endif

         if(abs(fa) < abs(fb)) then
            tmp0 = ak
            ak = bk
            bk = tmp0
            tmp0 = fa
            fa = fb
            fb = tmp0
         endif

         converged=(abs(bk-ak) < root_tol .OR. &
                    abs(fb) < root_tol .OR. abs(fs) < root_tol)

         !write(*,"(I3,4(3x,g12.4),3x,L1)") lc, ak, bk, fa, fb, converged
         lc = lc + 1
      enddo

      if(abs(fs) < abs(fb)) then
         BrentsMethod = s
      else
         BrentsMethod = bk
      endif

      end function BrentsMethod


      END SUBROUTINE USR0
