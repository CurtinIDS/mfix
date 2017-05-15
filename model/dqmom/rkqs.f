      SUBROUTINE rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext)

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: N

      DOUBLE PRECISION :: y(N), dydx(N), yscal(N)
      DOUBLE PRECISION x, hTry, EPs, hDid, hNext

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER i

      DOUBLE PRECISION :: errmax, h, htemp, xnew
      DOUBLE PRECISION :: yerr(N), ytemp(N)

      DOUBLE PRECISION, parameter :: SAFETY =  0.9
      DOUBLE PRECISION, parameter :: PGROW  = -0.2
      DOUBLE PRECISION, parameter :: PSHRNK = -0.25
      DOUBLE PRECISION, parameter :: ERRCON =  1.89e-4

      h=htry

      do
         call rkck(y,dydx,n,x,h,ytemp,yerr)
         errmax=0.
         do i=1,n
          errmax=max(errmax,abs(yerr(i)/yscal(i)))
         end do
         errmax=errmax/eps
         if(errmax<=1.0) exit
         htemp=SAFETY*h*(errmax**PSHRNK)
         h=sign(max(abs(htemp),0.1*abs(h)),h)
         xnew=x+h
         if(xnew==x) write(*,*) 'WARNING: stepsize underflow in rkqs'
      enddo

      IF(errmax>ERRCON) THEN
         HNEXT = SAFETY*h*(errmax**PGROW)
      ELSE
         HNEXT = 5.0*h
      ENDIF

      hdid=h
      x=x+h
      y(:)=ytemp(:)

      RETURN
      END SUBROUTINE rkqs
