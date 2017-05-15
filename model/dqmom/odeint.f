
       SUBROUTINE odeint(ystart1,nvar,x1,x2,eps,h1,hmin,nok,nbad)
       IMPLICIT NONE

       INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX

       DOUBLE PRECISION  eps,h1,hmin,x1,x2,ystart1(nvar),TINY
       PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
       INTEGER i,kmax1,kount,nstp
       DOUBLE PRECISION dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),&
       xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
       x=x1
       h=sign(h1,x2-x1)
       nok=0
       nbad=0
       kount=0
!#######################################
!      Don't store the intermediate results.
!#######################################
       kmax1=0.0
       dxsav=0.0
!##############################################
! store the intermediate results. uncooment them.
!      kmax1=KMAXX
!       dxsav=0.0
!##################################################
       do i=1,nvar
        y(i)=ystart1(i)
       end do


       if (kmax1>0) xsav=x-2.*dxsav
       do  nstp=1,MAXSTP
          call source_population_eq(x,y,dydx)

         do i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
         end do
         if(kmax1>0)then
          if(abs(x-xsav)>abs(dxsav)) then
            if(kount<kmax1-1)then
              kount=kount+1
              xp(kount)=x
              do i=1,nvar
                yp(i,kount)=y(i)
              end do
              xsav=x
             endif
           endif
          endif
         if((x+h-x2)*(x+h-x1)>0.) h=x2-x
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext)
         if(hdid==h)then
           nok=nok+1
         else
          nbad=nbad+1
         endif
         if((x-x2)*(x2-x1)>=0.)then
          do i=1,nvar
            ystart1(i)=y(i)
          end do
          if(kmax1/=0)then
            kount=kount+1
            xp(kount)=x
            do i=1,nvar
              yp(i,kount)=y(i)
            end do
          endif
          return
         endif
        if(abs(hnext)<hmin) write(*,*) 'WARNING: stepsize smaller than minimum in odeint'
        h=hnext
       end do
      write(*,*) 'WARNING: too many steps in odeint'
      return
      END
