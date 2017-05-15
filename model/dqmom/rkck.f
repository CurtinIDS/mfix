       SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr)
       IMPLICIT NONE

       INTEGER n,NMAX
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      PARAMETER (NMAX=50)
      INTEGER i
      double precision ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),&
      ak6(NMAX),ytemp(NMAX),A2,A3,A4,A5,A6,B21,&
      B31,B32,B41,B42,B43,B51,B52,B53,&
      B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,&
      B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,&
      B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,&
      B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,&
      C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,&
      DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,&
      DC6=C6-.25)
      do i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
      end do

      call source_population_eq(x+A2*h,ytemp,ak2)
      do i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      end do
      call source_population_eq(x+A3*h,ytemp,ak3)
      do i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      end do
      call source_population_eq(x+A4*h,ytemp,ak4)
      do i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      end do
      call source_population_eq(x+A5*h,ytemp,ak5)
      do i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+&
      B65*ak5(i))
      end do
      call source_population_eq(x+A6*h,ytemp,ak6)
      do i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      end do
      do i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*&
      ak6(i))
      end do
      return
      END
