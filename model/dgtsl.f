!
      subroutine dgtsl(n, c, d, e, b, info)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:38:13   1/21/99
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, info
      double precision, dimension(n) :: c, d, e, b
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, kb, kp1, nm1, nm2
      double precision :: t,tc,td,te,tb
!-----------------------------------------------
!
!     dgtsl given a general tridiagonal matrix and a right hand
!     side will find the solution.
!
!     on entry
!
!        n       integer
!                is the order of the tridiagonal matrix.
!
!        c       double precision(n)
!                is the subdiagonal of the tridiagonal matrix.
!                c(2) through c(n) should contain the subdiagonal.
!                on output c is destroyed.
!
!        d       double precision(n)
!                is the diagonal of the tridiagonal matrix.
!                on output d is destroyed.
!
!        e       double precision(n)
!                is the superdiagonal of the tridiagonal matrix.
!                e(1) through e(n-1) should contain the superdiagonal.
!                on output e is destroyed.
!
!        b       double precision(n)
!                is the right hand side vector.
!
!     on return
!
!        b       is the solution vector.
!
!        info    integer
!                = 0 normal value.
!                = k if the k-th element of the diagonal becomes
!                    exactly zero.  the subroutine returns when
!                    this is detected.
!
!     linpack. this version dated 08/14/78 .
!     jack dongarra, argonne national laboratory.
!
!     no externals
!     fortran dabs
!
!     internal variables
!
!     begin block permitting ...exits to 100
!
      info = 0
      c(1) = d(1)
      nm1 = n - 1
      if (nm1 >= 1) then
         d(1) = e(1)
         e(1) = 0.0D0
         e(n) = 0.0D0
!
         do k = 1, nm1
            kp1 = k + 1
!
!              find the largest of the two rows
!
            if (dabs(c(kp1)) >= dabs(c(k))) then
!
!                 interchange row
!
               tc = c(kp1)
               c(kp1) = c(k)
               c(k) = tc

               td = d(kp1)
               d(kp1) = d(k)
               d(k) = td

               te = e(kp1)
               e(kp1) = e(k)
               e(k) = te

               tb = b(kp1)
               b(kp1) = b(k)
               b(k) = tb
            endif
            if (c(k) == 0.0D0) then
               info = k
!     ............exit
               go to 100
            endif
            t = -c(kp1)/c(k)
            c(kp1) = d(kp1) + t*d(k)
            d(kp1) = e(kp1) + t*e(k)
            e(kp1) = 0.0D0
            b(kp1) = b(kp1) + t*b(k)
         end do
      endif
      if (c(n) == 0.0D0) then
         info = n
      else
         nm2 = n - 2
         b(n) = b(n)/c(n)
         if (n /= 1) then
            b(nm1) = (b(nm1)-d(nm1)*b(n))/c(nm1)
            if (nm2 >= 1) then
               do kb = 1, nm2
                  k = nm2 - kb + 1
                  b(k) = (b(k)-d(k)*b(k+1)-e(k)*b(k+2))/c(k)
               end do
            endif
         endif
      endif
  100 continue
      return
      end subroutine dgtsl
