    SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big,dum,pivinv

        ipiv(:)=0

      do  i=1,n
        big=0.
        do j=1,n
          if(ipiv(j)/=1)then
            do k=1,n
              if (ipiv(k)==0) then
                if (abs(a(j,k))>=big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k)>1) then
                write(*,*) 'WARNING: singular matrix in gaussj'
              endif
           enddo
          endif
        enddo
        ipiv(icol)=ipiv(icol)+1
        if (irow/=icol) then
          do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          enddo
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          enddo
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol)==0.) write(*,*) 'WARNING: singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        enddo
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        enddo
        do ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            enddo
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            enddo
          endif
        enddo
      enddo
      do l=n,1,-1
        if(indxr(l)/=indxc(l))then
          do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
        enddo
        endif
       enddo
      return
      END
